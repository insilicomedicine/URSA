from pathlib import Path

import pytest
from retrocast.exceptions import ArtifactFormatError
from retrocast.io import save_collected_candidates
from retrocast.io import save_collected_routes
from retrocast.models.candidates import Candidate
from retrocast.models.candidates import FailureRecord
from retrocast.models.route import Molecule
from retrocast.models.route import Route
from retrocast.typing import ErrorCode
from retrocast.typing import InChIKeyStr
from retrocast.typing import SmilesStr

from ursa.cli import _apply_top_k
from ursa.cli import _build_parser
from ursa.cli import _load_candidates
from ursa.cli import _load_routes


def _route(smiles: str) -> Route:
    """Build a minimal single-node retrocast ``Route`` identified by ``smiles``."""
    return Route(
        target=Molecule(smiles=SmilesStr(smiles), inchikey=InChIKeyStr(f"KEY-{smiles}"))
    )


def test_load_routes_roundtrip_and_top_k(tmp_path: Path):
    archive = tmp_path / "routes.json.gz"
    save_collected_routes({"t1": [_route("C"), _route("CC"), _route("CCC")]}, archive)

    full = _load_routes(str(archive))
    assert [type(r).__name__ for r in full["t1"]] == ["Route", "Route", "Route"]
    assert [r.target.smiles for r in full["t1"]] == ["C", "CC", "CCC"]

    # top_k keeps the first K (list order already encodes rank).
    top2 = _load_routes(str(archive), top_k=2)
    assert [r.target.smiles for r in top2["t1"]] == ["C", "CC"]


def test_load_candidates_sorts_by_rank_and_drops_failures(tmp_path: Path):
    archive = tmp_path / "candidates.json.gz"
    save_collected_candidates(
        {
            "t1": [
                Candidate(rank=3, route=_route("CCC")),
                Candidate(rank=1, route=_route("C")),
                Candidate(
                    rank=4,
                    failure=FailureRecord(code=ErrorCode("io.decode_failed")),
                ),
                Candidate(rank=2, route=_route("CC")),
            ]
        },
        archive,
    )

    routes = _load_candidates(str(archive))
    # Failure candidate dropped; survivors unwrapped to Routes ordered by rank.
    assert all(isinstance(r, Route) for r in routes["t1"])
    assert [r.target.smiles for r in routes["t1"]] == ["C", "CC", "CCC"]

    top2 = _load_candidates(str(archive), top_k=2)
    assert [r.target.smiles for r in top2["t1"]] == ["C", "CC"]


def test_load_candidates_rejects_routes_archive(tmp_path: Path):
    archive = tmp_path / "routes.json.gz"
    save_collected_routes({"t1": [_route("C")]}, archive)
    # The two archive shapes are mutually exclusive; loaders are strict.
    with pytest.raises(ArtifactFormatError):
        _load_candidates(str(archive))


def test_apply_top_k_zero_keeps_all():
    routes = {"t1": [_route("C"), _route("CC")]}
    assert _apply_top_k(routes, 0) == routes
    assert [r.target.smiles for r in _apply_top_k(routes, 1)["t1"]] == ["C"]


_BASE_ARGS = ["--routes", "x.json.gz", "--benchmark", "EXPERT_2026", "--output", "out"]


def test_workers_arg_defaults_to_none():
    args = _build_parser().parse_args(_BASE_ARGS)
    assert args.workers is None


def test_workers_arg_parsed():
    args = _build_parser().parse_args(_BASE_ARGS + ["-j", "4"])
    assert args.workers == 4
