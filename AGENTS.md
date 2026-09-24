# Repository agent guide

This file is tool-agnostic and applies to any coding agent working in URSA.

## Repository orientation (required)

Before running commands or editing files:

1. Read the root `README.md` for project scope, installation, data assets,
   supported inputs, Solv-N semantics, and license restrictions.
2. Check `git status` and preserve unrelated or user-authored changes.
3. Read `pyproject.toml` before changing dependencies or tooling.
4. For changes to the scoring pipeline or domain model, inspect
   `doc/uml/architecture.uml` and the relevant tests first.
5. Read the task-specific documentation listed below before launching an
   expensive or externally billed operation.

Use `uv` for environment and command execution. Do not commit generated data,
credentials, benchmark results, or unrelated untracked files unless the user
explicitly requests it. Do not commit, push, or publish without authorization.

Before handing off code changes, run:

```bash
uv sync --extra dev
uv run pre-commit run --all-files
uv run pytest
```

## LLM benchmark workflow

Use the following sections when a user asks to configure, run, resume, inspect,
or compare an LLM benchmark.

## Read first

- Project overview and Solv-N semantics: `README.md`
- User guide and notebook: `doc/llm_benchmark/`
- Exact prompt templates: `scripts/llm_benchmark/templates.py`
- Scoring implementation: `scripts/llm_benchmark/scoring.py`
- Entry point: `scripts/benchmark_llm.py`

## Gather the run configuration

Before making API calls, determine the following. Ask only for missing choices:

1. Model identifier and provider. Examples: `openai/gpt-5`,
   `gemini/gemini-2.5-pro`, `anthropic/claude-sonnet-4-5`, `xai/grok-4`, or
   `azure/<deployment-name>`.
2. Benchmark: `EXPERT_2026`, `DRUGS_CLINICALS_2026`, or a custom CSV path.
3. Protocol size: target limit and samples per target. The published protocol
   is all targets with 10 samples; a smoke test is 10 targets with 1 sample.
4. Concurrency: inference request workers and ChemCensor scoring workers.
5. Output directory and whether to resume or start fresh.

State the total planned API calls (`targets × samples`) before a large run.
If the user already supplied all choices and explicitly requested the run, do
not ask them again.

## Credentials

Credentials are loaded automatically from the repository `.env`.

- Generic: `LLM_API_KEY`, optionally `LLM_BASE_URL`.
- Azure: `AZURE_OPENAI_ENDPOINT`, `AZURE_OPENAI_API_KEY`, and
  `AZURE_OPENAI_API_VERSION`.
- LiteLLM also understands provider-specific variables such as
  `OPENAI_API_KEY`, `GEMINI_API_KEY`, `ANTHROPIC_API_KEY`, and `XAI_API_KEY`.

Never ask the user to paste a secret into chat, source code, a notebook, or a
command. Ask whether the required variable is configured. If diagnosis is
needed, report only whether variables are set and sanitize endpoint output.
Never commit `.env`.

For Azure, `--model` must be `azure/<deployment-name>`. The deployment name is
the configured Azure deployment, which may differ from the base model name.

## Setup and navigation

```bash
uv sync --extra llm-benchmark
uv run python scripts/benchmark_llm.py --help
```

Use a dry run to verify target loading, canonicalization, and prompt format
without calling a model:

```bash
uv run python scripts/benchmark_llm.py \
  --model MODEL \
  --benchmark EXPERT_2026 \
  --output data/results/RUN_NAME \
  --dry-run
```

Recommended smoke test:

```bash
uv run python scripts/benchmark_llm.py \
  --model MODEL \
  --benchmark EXPERT_2026 \
  --limit 3 \
  --samples 1 \
  --request-workers 8 \
  --score-workers 0 \
  --output data/results/RUN_NAME
```

Full best-of-10:

```bash
uv run python scripts/benchmark_llm.py \
  --model MODEL \
  --benchmark EXPERT_2026 \
  --samples 10 \
  --request-workers 100 \
  --score-workers 0 \
  --output data/results/RUN_NAME
```

## Run behavior

- Target SMILES are canonicalized with RDKit before prompting.
- Completions are appended immediately to `OUTPUT/completions.jsonl`.
- Repeating a command resumes missing samples; `--fresh` discards existing
  completions.
- Use a separate output directory for every model and protocol. Resume counting
  is target-based and does not distinguish models.
- `--sample-only` performs inference only; `--score-only` reuses completions.
- `--request-workers` controls concurrent API calls.
- `--score-workers 0` lets ChemCensor choose the process count automatically.

## Verify and report

After completion, verify the process exit code and report:

- adapted/total completions and `adapt_rate`;
- Solv-0, Solv-1, and Solv-2;
- `molecules_with_route / total_molecules`;
- output paths and elapsed time.

Expected artifacts are `completions.jsonl`, `llm_metrics.json`, and
`llm_best_paths.json`. Do not commit generated results unless the user
explicitly requests it.

CUDA-driver and `pkg_resources` warnings are non-fatal when scoring completes
successfully on CPU. A `KeyboardInterrupt` means a user stopped the run; resume
without `--fresh` to continue from saved completions.
