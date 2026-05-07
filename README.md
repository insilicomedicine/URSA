# URSA

URSA is a framework for evaluating retrosynthetic routes: it checks the structural consistency of the tree, verifies that starting materials are present in the building-block catalog, generates collapsed variants, scores every step with `ChemCensor`, and aggregates dataset-level metrics.

## Installation

### 1. Install ChemCensor dependency
ChemCensor is not published on PyPI. Download the repository archive from
https://anonymous.4open.science/r/ChemCensor-81B0/ (use the **Download ZIP**
button), then unpack it into `./chemcensor` at the URSA repo root
(this directory is gitignored):

```bash
unzip /path/to/downloads/ChemCensor-81B0.zip
mv ChemCensor-81B0 chemcensor
```

### 2. Install uv

We recommend using uv for fast, reliable dependency management.

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

### 3. Install URSA dependencies 

```bash
uv sync --extra dev
```

## Data prerequisites

Before working with benchmarks or the full scoring pipeline locally, download the following assets and place them under the paths below (relative to the repo root). Default paths match `ursa.DataConfig` in `src/ursa/configs/data_config.py`:

- `data/building_blocks/URSA_BBs_v1_0_0.csv` — building-block catalog for `BuildingBlockChecker`
- `data/chemcensor_db/ChemCensor_DB_v1_0_0.sqlite` — ChemCensor SQLite database for scoring

Download links: https://osf.io/wms6r/overview?view_only=073f80629f674f9084d95b7efc9e01ba

```bash
mkdir -p data/building_blocks data/chemcensor_db
unzip /path/to/downloads/'D3. URSA_BBs.csv.zip' -d data/building_blocks
unzip /path/to/downloads/'D6. URSA-minor-0.5.2-U2_database.zip' -d data/chemcensor_db
```

Once the zip contents are extracted, build the ChemCensor SQLite database from
the Parquet export:

```bash
python scripts/import_parquet_to_sqlite.py \
    --parquet-dir data/chemcensor_db/uspto_full_parquet \
    --out-sqlite data/chemcensor_db/ChemCensor_DB_v1_0_0.sqlite
```

## Example run

`data/example_data/` ships a ready-to-use bundle of two raw retrochimera predictions with high ChemCensor scores from the `EXPERT_2026` benchmark:

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

## Built-in benchmark sets

Located in `data/URSA_benchmarking_sets/`:

- `EXPERT_2026`
- `DRUGS_CLINICALS_2026`
- `USPTO_190`

A custom set of target molecules can be supplied by passing a CSV path instead of a preset name (default columns: `Structure ID`, `SMILES`; overridable via `--id-col` / `--smiles-col`).
