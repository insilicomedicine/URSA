from __future__ import annotations

import json
import sys
import time
from collections import defaultdict
from concurrent.futures import as_completed
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from .templates import build_prompt


def chat_completion(
    *,
    base_url: str | None,
    api_key: str | None,
    api_version: str | None,
    model: str,
    prompt: str,
    timeout: float,
    temperature: float | None,
    attempts: int = 4,
) -> str:
    """Request one completion from any provider supported by LiteLLM."""
    try:
        import litellm
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "LLM sampling requires LiteLLM. Install it with "
            "`uv sync --extra llm-benchmark`."
        ) from exc
    # Custom Azure deployment names are absent from LiteLLM's public model
    # catalog. Routing still works, but LiteLLM otherwise prints a misleading
    # provider-list hint for every request.
    litellm.suppress_debug_info = True

    request: dict = {
        "model": model,
        "messages": [{"role": "user", "content": prompt}],
        "timeout": timeout,
    }
    if api_key:
        request["api_key"] = api_key
    if base_url:
        request["api_base"] = base_url
    if api_version:
        request["api_version"] = api_version
    if temperature is not None:
        request["temperature"] = temperature

    last_error: Exception | None = None
    for attempt in range(attempts):
        try:
            response = litellm.completion(**request)
            content = response.choices[0].message.content
            if not isinstance(content, str) or not content.strip():
                raise RuntimeError("model returned an empty completion")
            return content
        except Exception as exc:
            last_error = exc
            if attempt == attempts - 1:
                raise
        time.sleep(2**attempt)
    raise RuntimeError(f"chat completion failed: {last_error}")


def sample_completions(
    targets: list[str],
    existing_records: list[dict],
    *,
    n_samples: int,
    completions_path: Path,
    base_url: str | None,
    api_key: str | None,
    api_version: str | None,
    model: str,
    timeout: float,
    temperature: float | None,
    request_workers: int,
    fresh: bool,
) -> None:
    """Append missing best-of-N samples to the completion JSONL."""
    completions_path.parent.mkdir(parents=True, exist_ok=True)
    if fresh and completions_path.exists():
        completions_path.unlink()
        existing_records = []

    have: dict[str, int] = defaultdict(int)
    for record in existing_records:
        have[record["meta"]["product_smiles"]] += 1

    unique_targets = list(dict.fromkeys(targets))
    jobs: list[tuple[str, str]] = []
    for target in unique_targets:
        for _ in range(max(n_samples - have[target], 0)):
            jobs.append((target, build_prompt(target)))

    done = sum(min(have[target], n_samples) for target in unique_targets)
    total = len(unique_targets) * n_samples
    print(f"Sampling {model}: {len(jobs)} requests ({done}/{total} already on disk)")
    if not jobs:
        return

    with completions_path.open("a") as handle:
        # Validate provider configuration before starting the worker pool.
        target, prompt = jobs.pop(0)
        completion = chat_completion(
            base_url=base_url,
            api_key=api_key,
            api_version=api_version,
            model=model,
            prompt=prompt,
            timeout=timeout,
            temperature=temperature,
        )
        record = {
            "completion": completion,
            "meta": {"product_smiles": target},
        }
        handle.write(json.dumps(record, ensure_ascii=False) + "\n")
        handle.flush()
        finished = done + 1
        print(f"  {finished}/{total} {target}", file=sys.stderr)

        with ThreadPoolExecutor(max_workers=request_workers) as pool:
            futures = {
                pool.submit(
                    chat_completion,
                    base_url=base_url,
                    api_key=api_key,
                    api_version=api_version,
                    model=model,
                    prompt=prompt,
                    timeout=timeout,
                    temperature=temperature,
                ): target
                for target, prompt in jobs
            }
            for future in as_completed(futures):
                target = futures[future]
                record = {
                    "completion": future.result(),
                    "meta": {"product_smiles": target},
                }
                handle.write(json.dumps(record, ensure_ascii=False) + "\n")
                handle.flush()
                finished += 1
                print(f"  {finished}/{total} {target}", file=sys.stderr)
