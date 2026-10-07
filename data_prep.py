"""
First draft of data prepoaration techniques using swiss-prot and GO ontology 

Inputs:
  - Swiss-Prot flat file: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz
  - GO ontology: https://purl.obolibrary.org/obo/go/go-basic.obo

Output:
  - data/processed/proteins.parquet — protein_id, sequence, length, go_terms (list)
"""

import argparse
from pathlib import Path

import pandas as pd
from Bio import SwissProt
from goatools.obo_parser import GODag
from tqdm import tqdm


def parse_swissprot(dat_path: str, max_records=None):
    """Yield (protein_id, sequence, raw_go_ids) tuples from a Swiss-Prot .dat file."""
    with open(dat_path) as handle:
        for i, record in enumerate(SwissProt.parse(handle)):
            if max_records and i >= max_records:
                break
            protein_id = record.accessions[0]
            sequence = record.sequence
            go_ids = [xref[1] for xref in record.cross_references if xref[0] == "GO"]
            yield protein_id, sequence, go_ids


def propagate_go_terms(go_ids, godag) -> set:
    """True-path rule: every ancestor of an annotated term is also annotated."""
    propagated = set()
    for go_id in go_ids:
        if go_id not in godag:
            continue
        propagated.add(go_id)
        propagated.update(anc.id for anc in godag[go_id].get_all_parents())
    return propagated


def build_dataset(dat_path, obo_path, out_path, max_records=None):
    print("Loading GO DAG...")
    godag = GODag(obo_path)

    rows = []
    print("Parsing Swiss-Prot records...")
    for protein_id, sequence, go_ids in tqdm(parse_swissprot(dat_path, max_records)):
        if not sequence or not go_ids:
            continue
        propagated = propagate_go_terms(go_ids, godag)
        rows.append({
            "protein_id": protein_id,
            "sequence": sequence,
            "length": len(sequence),
            "go_terms": sorted(propagated),
        })

    df = pd.DataFrame(rows)
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out_path)
    print(f"Saved {len(df)} proteins to {out_path}")
    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--dat", required=True, help="Path to uniprot_sprot.dat")
    parser.add_argument("--obo", required=True, help="Path to go-basic.obo")
    parser.add_argument("--out", default="data/processed/proteins.parquet")
    parser.add_argument("--max-records", type=int, default=None, help="Limit records for quick testing")
    args = parser.parse_args()

    build_dataset(args.dat, args.obo, args.out, args.max_records)
