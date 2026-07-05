#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys
from io import StringIO
from pathlib import Path

from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from aessp.api.cache import FileCache
from aessp.api.ncbi import NCBIClient, sequence_sha256
from aessp.io import ensure_parent, read_csv_or_empty
from aessp.validation.catpred import as_bool, catpred_manifest_from_validation


def _matching_fasta_record(fasta_text: str, accession: str):
    records = list(SeqIO.parse(StringIO(fasta_text), "fasta"))
    if not records:
        return None
    accession_keys = {accession, accession.split(".", 1)[0]}
    for record in records:
        record_keys = {record.id, record.id.split("|")[-1], record.id.split("|")[-1].split(".", 1)[0]}
        if accession_keys & record_keys:
            return record
    return records[0]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Export private FASTA inputs for CatPred.")
    parser.add_argument("--validation-csv", default="data/processed/enzyme_validation_table.csv")
    parser.add_argument("--email", default="", help="Contact email for NCBI E-utilities.")
    parser.add_argument("--out-fasta", default="data/private/catpred/catpred_input.fasta")
    parser.add_argument("--out-manifest", default="data/private/catpred/catpred_input_manifest.csv")
    parser.add_argument("--cache-dir", default="data/cache")
    parser.add_argument("--only-catpred-ready", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    validation = read_csv_or_empty(args.validation_csv)
    selected = validation.fillna("")
    if args.only_catpred_ready and "catpred_ready" in selected.columns:
        selected = selected.loc[selected["catpred_ready"].apply(as_bool)].copy()
    selected = selected.reset_index(drop=True)

    client = NCBIClient(email=args.email, cache=FileCache(args.cache_dir))
    fasta_lines: list[str] = []
    manifest = catpred_manifest_from_validation(selected, only_ready=False)

    for index, row in selected.iterrows():
        accession = str(row.get("protein_accession", "")).strip()
        if not accession:
            continue
        fasta_text = client.efetch_fasta("protein", [accession])
        record = _matching_fasta_record(fasta_text, accession)
        if record is None:
            continue
        sequence = str(record.seq)
        header = f"{row.get('enzyme_validation_id', '')}|{accession}"
        fasta_lines.append(f">{header}")
        fasta_lines.append(sequence)
        manifest.loc[index, "sequence_length"] = len(sequence)
        manifest.loc[index, "sequence_hash"] = sequence_sha256(sequence)

    ensure_parent(args.out_fasta)
    Path(args.out_fasta).write_text("\n".join(fasta_lines) + ("\n" if fasta_lines else ""), encoding="utf-8")
    ensure_parent(args.out_manifest)
    manifest.to_csv(args.out_manifest, index=False)


if __name__ == "__main__":
    main()
