import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from scripts.llm_benchmark import sampling
from scripts.llm_benchmark import templates
from scripts.llm_benchmark.pipeline import canonicalize_targets
from scripts.llm_benchmark.pipeline import load_records
from scripts.llm_benchmark.scoring import completions_to_paths
from scripts.llm_benchmark.scoring import score_records


def test_build_prompt_uses_target_and_format_example(monkeypatch):
    monkeypatch.setattr(
        templates.random,
        "choice",
        lambda choices: choices[0],
    )

    prompt = templates.build_prompt("CCO")

    assert "<smiles>CCO</smiles>" in prompt
    assert "only as an output-format example" in prompt
    assert prompt.endswith("Add no other text.")


def test_canonicalize_targets_preserves_stereochemistry():
    assert canonicalize_targets(["OC(C)", "F/C=C/F"]) == ["CCO", "F/C=C/F"]


def test_canonicalize_targets_rejects_invalid_smiles():
    with pytest.raises(SystemExit, match="invalid target SMILES"):
        canonicalize_targets(["not-a-smiles"])


def test_sampling_resumes_existing_jsonl(tmp_path: Path, monkeypatch):
    output = tmp_path / "completions.jsonl"
    existing = {
        "completion": "existing",
        "meta": {"product_smiles": "CCO"},
    }
    output.write_text(json.dumps(existing) + "\n")

    calls = []

    def fake_completion(**kwargs):
        calls.append(kwargs)
        return "new route"

    monkeypatch.setattr(sampling, "chat_completion", fake_completion)
    sampling.sample_completions(
        ["CCO", "CCN"],
        [existing],
        n_samples=2,
        completions_path=output,
        base_url=None,
        api_key=None,
        api_version=None,
        model="test/model",
        timeout=1,
        temperature=None,
        request_workers=2,
        fresh=False,
    )

    records = load_records(output)
    assert len(calls) == 3
    assert len(records) == 4
    assert [r["meta"]["product_smiles"] for r in records].count("CCO") == 2
    assert [r["meta"]["product_smiles"] for r in records].count("CCN") == 2


def test_format_example_adapts_to_a_route():
    target = "CC1CCN(Cc2ccccc2COC2(C(F)F)CCCC2)CC1"
    records = [
        {
            "completion": templates.FORMAT_INSTRUCTION,
            "meta": {"product_smiles": target},
        }
    ]

    paths, adapted, attempted = completions_to_paths(records, [target], top_k=10)

    assert (adapted, attempted, len(paths)) == (1, 1, 1)
    assert paths[0].root.canonical_smiles
    assert len(paths[0].root.children) == 2


def test_score_records_persists_adaptation_metrics(tmp_path: Path, monkeypatch):
    metrics = SimpleNamespace(
        total_molecules=1,
        molecules_with_route=1,
        routes_solv_0=1,
        routes_solv_1=1,
        routes_solv_2=1,
        solv_0=1.0,
        solv_1=1.0,
        solv_2=1.0,
        mean_score_without_fg=2.0,
        mean_score_with_fg=1.0,
    )

    class FakeResult:
        def __init__(self):
            self.metrics = metrics

        def save(self, output, stem=""):
            output.mkdir(parents=True, exist_ok=True)
            metrics_path = output / f"{stem}_metrics.json"
            paths_path = output / f"{stem}_best_paths.json"
            metrics_path.write_text("{}")
            paths_path.write_text("[]")
            return metrics_path, paths_path

    class FakeUrsa:
        def __init__(self, **kwargs):
            pass

        def score_dataset(self, paths, target_smiles):
            return FakeResult()

    monkeypatch.setitem(sys.modules, "ursa", SimpleNamespace(Ursa=FakeUrsa))
    monkeypatch.setattr(
        "scripts.llm_benchmark.scoring.completions_to_paths",
        lambda records, targets, top_k: ([object()], 8, 10),
    )

    score_records(
        [],
        ["CCO"],
        top_k=10,
        score_workers=None,
        bb_catalog=None,
        chemcensor_db=None,
        output=tmp_path,
        stem="llm",
    )

    saved = json.loads((tmp_path / "llm_metrics.json").read_text())
    assert saved["completions_total"] == 10
    assert saved["completions_adapted"] == 8
    assert saved["adapt_rate"] == 0.8
    assert saved["routes_kept"] == 1
