# URSA

URSA is a framework for evaluating retrosynthetic routes: it checks the structural consistency of the tree, verifies that starting materials are present in the building-block catalog, generates collapsed variants, scores every step with `ChemCensor`, and aggregates dataset-level metrics.

## Installation

ChemCensor is not published. Clone the ChemCensor repository into `./chemcensor` at the repo root (this directory is gitignored), then install dependencies:

```bash
git clone https://github.com/<ORG_OR_USER>/chemcensor.git chemcensor
uv sync --extra dev
```

## Data prerequisites

Before working with benchmarks or the full scoring pipeline locally, download the following assets and place them under the paths below (relative to the repo root). Default paths match `ursa.DataConfig` in `src/ursa/configs/data_config.py`:

- `data/building_blocks/URSA_BBs_v1_0_0.csv` — building-block catalog for `BuildingBlockChecker`
- `data/chemcensor_db/ChemCensor_DB_v1_0_0.sqlite` — ChemCensor SQLite database for scoring

Download links: *to be added.*

## Example run

`data/example_data/` ships a ready-to-use bundle of two raw retrochimera predictions with high ChemCensor scores from the `EXPERT_2026` benchmark:

- `X404-1768-3704.json` — best route score = 3.40, 5 steps
- `X404-1760-0042.json` — best route score = 2.33, 6 steps
- `bundled.json.gz` — gzipped JSON keyed by target SMILES, consumed directly by `ursa-bench`.

Score the bundle:

```bash
ursa-bench \
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
  solved_routes:         2
  Solv-2:                0.0200
  mean_chemcensor_score: 2.8182
  passed_steps:          11 / 11
```

Artifacts:

- `data/results/example_metrics.json` — aggregated dataset metrics.
- `data/results/example_best_paths.json` — best variant for each target with per-step scores.

## RetroCast integration

URSA consumes routes in the [RetroCast](https://pypi.org/project/retrocast/) format. `ursa-bench` exposes two input modes:

- `--routes ROUTES_JSON_GZ` — already-adapted RetroCast routes (`dict[target_id, list[Route]]`) serialized with `retrocast.io.save_routes`.
- `--input RAW_JSON_GZ --adapter ADAPTER` — raw model predictions keyed by target; URSA dispatches through the matching RetroCast adapter, writes a temporary routes archive, and loads it.

Supported adapter names (passed to `--adapter`):

```
aizynth, askcos, dms, dreamretro, multistepttl, paroutes,
retrochimera, retrostar, synplanner, syntheseus, synllama
```

Skip `--adapter` and use `--routes` when you already have RetroCast-formatted routes on disk — useful for caching the adaptation step across multiple benchmark runs.

## Programmatic API

```python
from ursa import Ursa, RetrosyntheticPath

ursa = Ursa()

# single path
result = ursa.score(path)

# dataset
dataset_result = ursa.score_dataset(paths, target_smiles=target_smiles)
print(dataset_result.metrics.solv_2)
dataset_result.save("data/results", stem="example")
```

`Ursa` orchestrates the pipeline: `PathConsistencyChecker` → `BuildingBlockChecker` → `PathCollapser` → `PathScorer` → `BestPathSelector` → `DatasetMetricsCalculator`.

## Built-in benchmark sets

Located in `data/URSA_benchmarking_sets/`:

- `EXPERT_2026`
- `DRUGS_CLINICALS_2026`
- `USPTO_190`

A custom set of target molecules can be supplied by passing a CSV path instead of a preset name (default columns: `Structure ID`, `SMILES`; overridable via `--id-col` / `--smiles-col`).

## Tests

```bash
pytest
```
