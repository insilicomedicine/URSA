# URSA

[![URSA Paper](https://img.shields.io/badge/URSA%20Paper-arXiv%3A2607.04688-B31B1B?style=flat-square&logo=arxiv&logoColor=white)](https://arxiv.org/abs/2607.04688)
[![ChemCensor Paper](https://img.shields.io/badge/ChemCensor%20Paper-arXiv%3A2602.03554-B31B1B?style=flat-square&logo=arxiv&logoColor=white)](https://arxiv.org/abs/2602.03554)
[![ChemCensor on GitHub](https://img.shields.io/badge/ChemCensor-GitHub-181717?style=flat-square&logo=github&logoColor=white)](https://github.com/insilicomedicine/ChemCensor)
[![RetroCast on GitHub](https://img.shields.io/badge/RetroCast-GitHub-181717?style=flat-square&logo=github&logoColor=white)](https://github.com/ischemist/project-procrustes)


URSA is a framework for evaluating retrosynthetic routes: it checks the structural consistency of the tree, verifies that starting materials are present in the building-block catalog, generates collapsed variants, scores every step with `ChemCensor`, and aggregates dataset-level metrics under the **Solv-N** hierarchy.

## Solv-N metrics

URSA reports a hierarchy of increasingly strict route-validity rates. Every level requires the route to terminate in commercially available building blocks (stock termination). Each rate is `routes_passing / total_molecules`, so targets without a route lower the score.

A *variant* is the original (uncollapsed) route together with every valid collapsed form produced by `PathCollapser` — the original path is always part of the candidate set.

| Level | Field | Requirement |
| --- | --- | --- |
| **Solv-0** | `solv_0` | Stock termination: the tree is consistent (valid SMILES, no breaks) **and** every starting material is in the catalog. |
| **Solv-1** | `solv_1` | Solv-0 **and** some variant (original or collapsed) where every step scores `> 0` with ChemCensor **without** functional-group matching (legal reaction center). |
| **Solv-2** | `solv_2` | Solv-0 **and** some variant (original or collapsed) where every step scores `> 0` with ChemCensor **with** functional-group matching (chemical plausibility). |

Because collapsing a route changes reaction identities (and therefore scores), the best variant is selected **independently for each level** over the full candidate set (original + collapsed): a route passes Solv-N iff *any* of those variants clears that level. `PathResult` exposes `passes_solv_0/1/2` plus `best_variant_solv_1` and `best_variant_solv_2` (the latter is used for display). A target that is itself a building block has no reaction steps and is flagged `is_no_synthesis` (it passes no Solv level).

Each `StepResult` carries both `score_without_fg` (Solv-1) and `score_with_fg` (Solv-2); the dataset metrics also report `mean_score_without_fg` / `mean_score_with_fg` as diagnostics over the per-level best variants, and `routes_solv_0/1/2` / `routes_no_synthesis` as raw counts.

## Installation

URSA depends on [ChemCensor](https://github.com/insilicomedicine/ChemCensor) for
per-step reaction scoring. `uv sync` installs it automatically from GitHub.

### 1. Install uv

We recommend using uv for fast, reliable dependency management.

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

### 2. Install URSA dependencies

```bash
uv sync --extra dev
```

## Data prerequisites

On the first scoring run, URSA downloads required data assets from HuggingFace:

- [ChemCensor database](https://huggingface.co/datasets/insilicomedicine/chemcensor) → `data/chemcensor_db/`
- [Building-block catalog](https://huggingface.co/datasets/insilicomedicine/URSA-BBs) → `data/building_blocks/`
- [Benchmark target sets](https://huggingface.co/datasets/insilicomedicine/URSA-benchmarking-sets) → `data/URSA_benchmarking_sets/` (when a built-in preset is used)

A legacy ChemCensor filename (`ChemCensor_DB_v1_0_0.sqlite`) in the same directory is also accepted.

## Example run

`data/example_data/` ships a ready-to-use bundle of two raw RetroChimera predictions with high ChemCensor scores from the `EXPERT_2026` benchmark:

- `X404-1768-3704.json` — best route score = 3.40, 5 steps
- `X404-1760-0042.json` — best route score = 2.33, 6 steps
- `bundled.json.gz` — gzipped JSON keyed by target SMILES, consumed directly by `ursa-bench`.

Score the bundle:

```bash
uv run ursa-bench \
    --input     data/example_data/bundled.json.gz \
    --adapter   retrochimera \
    --benchmark EXPERT_2026 \
    --output    data/results \
    --stem      example
```

Expected output:

```
Results:
  total_molecules:       100
  molecules_with_route:  2
  Solv-0 (STR):          0.0200  (2 routes)
  Solv-1:                0.0200  (2 routes)
  Solv-2:                0.0200  (2 routes)
  mean_score_without_fg: 3.0000
  mean_score_with_fg:    2.8182
```

Artifacts:

- `data/results/example_metrics.json` — aggregated dataset metrics.
- `data/results/example_best_paths.json` — per target: `passes_solv_0/1/2` and the Solv-1 / Solv-2 best variants with per-step scores.

## RetroCast integration

URSA consumes routes in the [RetroCast](https://pypi.org/project/retrocast/) format. `ursa-bench` exposes three input modes:

- `--routes ROUTES_JSON_GZ` — already-adapted RetroCast collected routes (`dict[target_id, list[Route]]`) serialized with `retrocast.io.save_collected_routes`.
- `--candidates CANDIDATES_JSON_GZ` — already-adapted RetroCast collected candidates; failure records are dropped and survivors are ordered by rank.
- `--input RAW_JSON_GZ --adapter ADAPTER` — raw model predictions keyed by target; URSA dispatches through the matching RetroCast adapter, writes a temporary routes archive, and loads it.

Supported adapter names (passed to `--adapter`):

```
aizynth, askcos, dms, dreamretro, multistepttl, paroutes,
retrochimera, retrostar, synplanner, syntheseus, synllama
```

Skip `--adapter` and use `--routes` when you already have RetroCast-formatted routes on disk — useful for caching the adaptation step across multiple benchmark runs.

## Additional options

- `--top-k K` — keep only the `K` best routes per target (ranked best-first); `0` keeps all. Default: `10`.
- `-j/--workers N` — score reactions in parallel using `N` ChemCensor worker processes (pass `0` to auto-scale to the CPU count). Omit for sequential scoring.

## Built-in benchmark sets

Built-in presets (`EXPERT_2026`, `DRUGS_CLINICALS_2026`) download their CSV files
from [URSA-benchmarking-sets](https://huggingface.co/datasets/insilicomedicine/URSA-benchmarking-sets)
on first use into `data/URSA_benchmarking_sets/`.

A custom set of target molecules can be supplied by passing a CSV path instead of a preset name (default columns: `Structure ID`, `SMILES`; overridable via `--id-col` / `--smiles-col`).

## Custom building blocks and targets

By default, `ursa-bench` downloads the URSA building-block catalog and uses a built-in benchmark preset. You can override either asset with local files.

### Custom target list

Prepare a CSV of product molecules to evaluate (one row per target):

```csv
Structure ID,SMILES
mol-1,CCO
mol-2,c1ccccc1
```

Run against your predictions:

```bash
uv run ursa-bench \
    --input     predictions.json.gz \
    --adapter   retrochimera \
    --benchmark my_targets.csv \
    --output    data/results \
    --stem      custom_targets
```

If your CSV uses different column names:

```bash
uv run ursa-bench \
    --input     predictions.json.gz \
    --adapter   retrochimera \
    --benchmark my_targets.csv \
    --id-col    mol_id \
    --smiles-col smiles \
    --output    data/results
```

### Custom building-block catalog

`BuildingBlockChecker` accepts:

- **CSV** with a `smiles` column (e.g. `bb_id,smiles`)
- **Plain text** — one SMILES per line (`.smi`)

Pass your catalog with `--bb-catalog`:

```bash
uv run ursa-bench \
    --input      predictions.json.gz \
    --adapter    retrochimera \
    --benchmark  my_targets.csv \
    --bb-catalog my_building_blocks.smi \
    --output     data/results
```

### Both custom assets

```bash
uv run ursa-bench \
    --input      predictions.json.gz \
    --adapter    retrochimera \
    --benchmark  my_targets.csv \
    --bb-catalog my_building_blocks.csv \
    --output     data/results \
    --stem       custom
```

The same options are available from Python:

```python
from ursa import Ursa
from ursa.datasets import BenchmarkDataset

ursa = Ursa(bb_catalog_path="my_building_blocks.csv")
dataset = BenchmarkDataset.from_csv(
    "my_targets.csv",
    id_col="mol_id",
    smiles_col="smiles",
)
result = ursa.score_dataset(paths, target_smiles=dataset.target_smiles)
```

---

## License

URSA is released under a license for **independent benchmarking and evaluation purposes only**. Use in products, pipelines, automated workflows, or redistribution requires prior written permission from Insilico. See [LICENSE](LICENSE) for full terms.

---

## Citation

If you use URSA in your work, please cite the URSA paper:

```bibtex
@misc{zagribelnyy2026ursachemistryawarebenchmarkutilitarian,
      title={URSA: Chemistry-Aware Benchmark for Utilitarian Retrosynthesis Assessment},
      author={Bogdan Zagribelnyy and Ivan Ilin and Nikita Bondarev and Anton Morgunov and Arkadii Lin and Maksim Kuznetsov and Rim Shayakhmetov and Vladimir Aladinskiy and Alex Aliper and Alex Zhavoronkov},
      year={2026},
      eprint={2607.04688},
      archivePrefix={arXiv},
      primaryClass={cs.LG},
      url={https://arxiv.org/abs/2607.04688},
}
```

and the ChemCensor paper:

```bibtex
@misc{zagribelnyy2026singleanswerenoughrethinking,
      title={When Single Answer Is Not Enough: Rethinking Single-Step Retrosynthesis Benchmarks for LLMs},
      author={Bogdan Zagribelnyy and Ivan Ilin and Maksim Kuznetsov and Nikita Bondarev and Mathieu Reymond and Roman Schutski and Thomas MacDougall and Rim Shayakhmetov and Zulfat Miftakhutdinov and Mikolaj Mizera and Vladimir Aladinskiy and Alex Aliper and Alex Zhavoronkov},
      year={2026},
      eprint={2602.03554},
      archivePrefix={arXiv},
      primaryClass={cs.LG},
      url={https://arxiv.org/abs/2602.03554},
}
```
