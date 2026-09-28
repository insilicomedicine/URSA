from __future__ import annotations

import logging
import re
import zipfile
from collections.abc import Mapping
from pathlib import Path

from ursa.basic.result import DatasetResult

logger = logging.getLogger("ursa.ursa")

_UNSAFE_FILENAME_CHARS = re.compile(r'[\\/:*?"<>|\s]+')


def _safe_filename(name: str) -> str:
    """Make a target id safe to use as a filename.

    Characters that are not valid in filenames (``/``, ``\\``, ``:``, ``*``,
    ``?``, ``"``, ``<``, ``>``, ``|``) and whitespace are collapsed to
    single underscores.

    :param name: The raw target id.
    :type name: str

    :return: A filesystem-safe version of ``name``.
    :rtype: str
    """
    safe = _UNSAFE_FILENAME_CHARS.sub("_", name).strip("_")
    return safe or "path"


def _prepare_reactions_for_drawing(reactions: list[str]) -> list[str]:
    """Run ChemCensor processing and return canonical reactions without orphans.

    ChemCensor removes auxiliary reactants after atom mapping. Its canonical
    reaction SMILES therefore contains only molecules that participate in the
    reaction and is the representation that should be rendered.

    :param reactions: Original forward reaction SMILES.
    :type reactions: list[str]

    :return: Canonical, orphan-free reaction SMILES in the same order.
    :rtype: list[str]
    """
    from chemcensor.basic import Reaction
    from chemcensor.processing import ReactionProcessor

    processor = ReactionProcessor()
    prepared = []
    for reaction_smiles in reactions:
        reaction = processor.process(Reaction(reaction_smiles=reaction_smiles))
        prepared.append(reaction.canonical_smiles)
    return prepared


def write_schemes_zip(
    result: DatasetResult,
    output_dir: Path | str,
    stem: str = "",
    target_id_by_smiles: Mapping[str, str] | None = None,
) -> Path | None:
    """Render each best path as a CDXML scheme and bundle them into a zip.

    For every :class:`~ursa.PathResult` in ``result`` that has at least one
    reaction step, the steps of its best variant are rendered into a single
    CDXML reaction scheme, with each reaction labelled by its chemcensor
    score. The schemes are written as individual ``.cdxml`` files inside a
    zip archive next to ``best_paths.json``.

    Paths with no steps are skipped. Schemes that cannot be drawn (e.g.
    recursive routes where a molecule is both produced and consumed in
    different steps) are logged and skipped, so a single bad path never
    aborts the whole archive.

    :param result: The dataset evaluation result whose best paths are drawn.
    :type result: DatasetResult
    :param output_dir: Directory to write the archive into (created if needed).
    :type output_dir: Path | str
    :param stem: Optional filename prefix, matching
        :meth:`~ursa.DatasetResult.save` (``{stem}_schemes.zip`` when set,
        otherwise ``schemes.zip``).
    :type stem: str
    :param target_id_by_smiles: Optional mapping from canonical target SMILES
        to benchmark molecule IDs. When supplied, molecule IDs are used as
        CDXML filenames instead of path IDs.
    :type target_id_by_smiles: Mapping[str, str] | None

    :return: Path to the written archive, or ``None`` if no scheme was drawn.
    :rtype: Path | None
    """
    # Imported lazily so that ``import ursa`` does not pull in lxml / the
    # rdkit drawing machinery unless schemes are actually requested.
    from ursa.drawing.cdxml.route_from_rxn_smiles import reaction_smiles_to_cdxml

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    prefix = f"{stem}_" if stem else ""
    zip_path = out / f"{prefix}schemes.zip"

    used_names: set[str] = set()
    written = 0

    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as archive:
        for path_result in result.path_results:
            step_results = path_result.best_variant_solv_2.step_results
            reaction_smiles = [sr.node.reaction_smiles for sr in step_results]
            # Skip paths without any reaction steps (nothing to draw).
            if not reaction_smiles:
                continue

            labels = [f"CC Score: {sr.score_with_fg:.2f}" for sr in step_results]
            path_id = path_result.original_path.path_id
            target_smiles = path_result.original_path.root.canonical_smiles
            target_id = (
                target_id_by_smiles.get(target_smiles, path_id)
                if target_id_by_smiles is not None
                else path_id
            )

            try:
                drawing_reactions = _prepare_reactions_for_drawing(reaction_smiles)
                cdxml = reaction_smiles_to_cdxml(drawing_reactions, labels=labels)
            except Exception as exc:  # noqa: BLE001 - one bad path must not abort
                logger.warning("Could not render CDXML for %s: %s", path_id, exc)
                continue

            base = _safe_filename(target_id)
            filename = f"{base}.cdxml"
            dedup = 1
            while filename in used_names:
                dedup += 1
                filename = f"{base}_{dedup}.cdxml"
            used_names.add(filename)

            archive.writestr(filename, cdxml)
            written += 1

    if not written:
        zip_path.unlink(missing_ok=True)
        logger.warning("No CDXML schemes were rendered; %s not written", zip_path.name)
        return None

    logger.info("Wrote %d CDXML schemes to %s", written, zip_path)
    return zip_path
