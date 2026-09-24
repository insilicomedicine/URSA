# LLM retrosynthesis benchmark

This workflow evaluates an LLM on multi-step retrosynthesis with URSA. For each
target it builds a randomized prompt, samples one or more completions, adapts
the generated `<synthesis_step>` blocks into route trees, and reports Solv-0/1/2.

Use [`benchmark.ipynb`](benchmark.ipynb) for an interactive run or
`scripts/benchmark_llm.py` from the command line.

## Setup

From the repository root:

```bash
uv sync --extra llm-benchmark
cp .env.example .env
```

Put credentials in `.env`; never put secrets in the notebook. Generic
LiteLLM-compatible providers can use:

```dotenv
LLM_API_KEY=...
# LLM_BASE_URL=https://custom-provider.example/v1
```

Azure OpenAI uses:

```dotenv
AZURE_OPENAI_ENDPOINT=https://<resource>.openai.azure.com/
AZURE_OPENAI_API_KEY=...
AZURE_OPENAI_API_VERSION=2024-10-21
```

The model is selected per run, not in `.env`. Examples include
`openai/gpt-5`, `gemini/gemini-2.5-pro`,
`anthropic/claude-sonnet-4-5`, `xai/grok-4`, and
`azure/<deployment-name>`.

## Benchmark sets

Two built-in sets are available:

- `EXPERT_2026` — 100 expert targets.
- `DRUGS_CLINICALS_2026` — drugs and clinical candidates.

A custom CSV can be passed instead. Its default SMILES column is
`product_smiles`; use `--smiles-col` to select another column. Target SMILES
are canonicalized with RDKit before prompting and scoring.

## Run from the CLI

Small smoke test:

```bash
uv run python scripts/benchmark_llm.py \
    --model azure/my-deployment \
    --benchmark EXPERT_2026 \
    --limit 10 \
    --samples 1 \
    --request-workers 8 \
    --score-workers 0 \
    --output data/results/azure-smoke
```

Full best-of-10 protocol:

```bash
uv run python scripts/benchmark_llm.py \
    --model azure/my-deployment \
    --benchmark EXPERT_2026 \
    --samples 10 \
    --request-workers 100 \
    --score-workers 0 \
    --output data/results/azure-best-of-10
```

`--request-workers` controls concurrent inference requests.
`--score-workers` controls ChemCensor worker processes; `0` auto-selects the
process count.

## Resume and rerun

Completions are appended to `OUTPUT/completions.jsonl` immediately. Repeating
the same command requests only missing samples. Use:

- `--sample-only` to stop after inference.
- `--score-only` to rescore an existing completion file.
- `--fresh` to delete existing completions before inference.
- `--dry-run` to print one prompt without calling the model.

Use a separate output directory for each model and protocol. Resume counting
is based on target SMILES, so reusing a directory for another model would mix
their completions.

## Outputs

The output directory contains:

- `completions.jsonl` — raw model responses and target metadata.
- `llm_metrics.json` — Solv-N metrics, completion/adaptation counts,
  `adapt_rate`, and retained route count.
- `llm_best_paths.json` — best route and per-step scores for each target.

Every Solv-N rate uses the full target list as its denominator. Targets with no
parsable route therefore lower the score. `adapt_rate` is the fraction of
completions successfully converted by the RetroCast `ursa` adapter.
