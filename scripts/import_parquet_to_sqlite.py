#!/usr/bin/env python3
"""
Import the ChemCensor reaction-center Parquet exports back into a SQLite DB.

Expected Parquet files (produced by scripts/export_uspto_full_sqlite_to_csv.py):
- reaction_centers.parquet
- centers_to_reactions.parquet
- reactions.parquet

This recreates the original schema:
- reaction_centers(
    reaction_center_smiles PK,
    fg_signature BLOB,
    is_multi_component INT,
    components TEXT
  )
- centers_to_reactions(
    PK(reaction_center_smiles,reaction_smiles),
    fg_signature BLOB,
    sear_signature BLOB
  )
- reactions(PK(reaction_smiles,document_id))
and the index:
- idx_ctr_reaction_smiles ON centers_to_reactions(reaction_smiles)
"""
from __future__ import annotations

import argparse
import sqlite3
from pathlib import Path
from typing import Iterable


def _require_pyarrow():
    try:
        import pyarrow.parquet as pq  # type: ignore
    except ModuleNotFoundError as e:  # pragma: no cover
        raise ModuleNotFoundError(
            "pyarrow is required for Parquet import. Install it with "
            "`pip install pyarrow` (or add it to your environment)."
        ) from e
    return pq


SCHEMA_SQL = """
CREATE TABLE reaction_centers (
    reaction_center_smiles TEXT PRIMARY KEY,
    fg_signature           BLOB NOT NULL DEFAULT (x''),
    is_multi_component     INTEGER NOT NULL DEFAULT 0,
    components             TEXT NOT NULL DEFAULT '[]'
);
CREATE TABLE centers_to_reactions (
    reaction_center_smiles TEXT NOT NULL
        REFERENCES reaction_centers(reaction_center_smiles),
    reaction_smiles TEXT NOT NULL,
    fg_signature    BLOB NOT NULL DEFAULT (x''),
    sear_signature  BLOB NOT NULL DEFAULT (x''),
    PRIMARY KEY (reaction_center_smiles, reaction_smiles)
);
CREATE TABLE reactions (
    reaction_smiles TEXT NOT NULL,
    document_id     TEXT NOT NULL,
    PRIMARY KEY (reaction_smiles, document_id)
);
CREATE INDEX idx_ctr_reaction_smiles
    ON centers_to_reactions(reaction_smiles);
""".strip()


def _connect(out_sqlite: Path) -> sqlite3.Connection:
    out_sqlite.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(out_sqlite))
    conn.execute("PRAGMA journal_mode=WAL;")
    conn.execute("PRAGMA synchronous=NORMAL;")
    conn.execute("PRAGMA temp_store=MEMORY;")
    return conn


def _init_schema(conn: sqlite3.Connection) -> None:
    conn.executescript(
        """
        DROP TABLE IF EXISTS centers_to_reactions;
        DROP TABLE IF EXISTS reaction_centers;
        DROP TABLE IF EXISTS reactions;
        """
    )
    conn.executescript(SCHEMA_SQL)


def _iter_batches(parquet_path: Path, batch_size: int):
    pq = _require_pyarrow()
    pf = pq.ParquetFile(str(parquet_path))
    # iter_batches yields RecordBatch with columnar arrays (often zero-copy).
    yield from pf.iter_batches(batch_size=batch_size)


def _executemany(conn: sqlite3.Connection, sql: str, rows: Iterable[tuple]) -> None:
    conn.executemany(sql, rows)


def import_reaction_centers(
    conn: sqlite3.Connection,
    parquet_path: Path,
    batch_size: int,
) -> None:
    insert_sql = """
    INSERT OR REPLACE INTO reaction_centers
      (reaction_center_smiles, fg_signature, is_multi_component, components)
    VALUES (?, ?, ?, ?)
    """
    for batch in _iter_batches(parquet_path, batch_size=batch_size):
        cols = batch.to_pydict()
        # Export uses `components_json` column name; map back to `components`.
        rows = zip(
            cols["reaction_center_smiles"],
            cols["fg_signature"],
            cols["is_multi_component"],
            cols["components_json"],
            strict=True,
        )
        _executemany(conn, insert_sql, rows)


def import_reactions(
    conn: sqlite3.Connection, parquet_path: Path, batch_size: int
) -> None:
    insert_sql = """
    INSERT OR REPLACE INTO reactions
      (reaction_smiles, document_id)
    VALUES (?, ?)
    """
    for batch in _iter_batches(parquet_path, batch_size=batch_size):
        cols = batch.to_pydict()
        rows = zip(cols["reaction_smiles"], cols["document_id"], strict=True)
        _executemany(conn, insert_sql, rows)


def import_centers_to_reactions(
    conn: sqlite3.Connection, parquet_path: Path, batch_size: int
) -> None:
    insert_sql = """
    INSERT OR REPLACE INTO centers_to_reactions
      (reaction_center_smiles, reaction_smiles, fg_signature, sear_signature)
    VALUES (?, ?, ?, ?)
    """
    for batch in _iter_batches(parquet_path, batch_size=batch_size):
        cols = batch.to_pydict()
        rows = zip(
            cols["reaction_center_smiles"],
            cols["reaction_smiles"],
            cols["fg_signature"],
            cols["sear_signature"],
            strict=True,
        )
        _executemany(conn, insert_sql, rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--parquet-dir",
        default="data/chemcensor_db/uspto_full_parquet",
        help="Directory containing the Parquet exports.",
    )
    parser.add_argument(
        "--out-sqlite",
        default="uspto_full_roundtrip.sqlite",
        help="Output SQLite path to create/overwrite.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=50_000,
        help="Row batch size for Parquet -> SQLite inserts.",
    )
    args = parser.parse_args()

    parquet_dir = Path(args.parquet_dir)
    out_sqlite = Path(args.out_sqlite)

    reaction_centers_pq = parquet_dir / "reaction_centers.parquet"
    centers_to_reactions_pq = parquet_dir / "centers_to_reactions.parquet"
    reactions_pq = parquet_dir / "reactions.parquet"

    for p in [reaction_centers_pq, centers_to_reactions_pq, reactions_pq]:
        if not p.exists():
            raise FileNotFoundError(p)

    conn = _connect(out_sqlite)
    try:
        _init_schema(conn)
        with conn:
            # Order matters because of FK: reaction_centers first.
            import_reaction_centers(
                conn,
                reaction_centers_pq,
                batch_size=args.batch_size,
            )
            import_reactions(conn, reactions_pq, batch_size=args.batch_size)
            import_centers_to_reactions(
                conn, centers_to_reactions_pq, batch_size=args.batch_size
            )
    finally:
        conn.close()

    print(str(out_sqlite))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
