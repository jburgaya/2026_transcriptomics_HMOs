#!/usr/bin/env python3
import csv
import re
import argparse

parser = argparse.ArgumentParser(description="Extract GO terms from gff file")
parser.add_argument("--gff",
                    required=True,
                    help="Input gff file")
parser.add_argument("--out",
                    required=True,
                    help="Output")

args = parser.parse_args()

gff_file = args.gff
out_file = args.out

locus_go = {}

with open(gff_file) as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 9:
            continue

        feature_type = parts[2]
        attrs = parts[8]

        # Only CDS features contain GO terms
        if feature_type != "CDS":
            continue

        # Extract locus_tag
        locus_match = re.search(r"locus_tag=([^;]+)", attrs)
        if not locus_match:
            continue
        locus_tag = locus_match.group(1)

        # Extract GO terms
        go_match = re.search(r"Ontology_term=([^;]+)", attrs)
        go_terms = go_match.group(1) if go_match else ""

        # Merge GO terms if multiple CDS per locus_tag
        if locus_tag in locus_go:
            existing_terms = set(locus_go[locus_tag].split(",")) if locus_go[locus_tag] else set()
            new_terms = set(go_terms.split(",")) if go_terms else set()
            locus_go[locus_tag] = ",".join(sorted(existing_terms | new_terms))
        else:
            locus_go[locus_tag] = go_terms

with open(out_file, "w", newline="") as outcsv:
    writer = csv.writer(outcsv, quoting=csv.QUOTE_ALL)
    writer.writerow(["Gene", "GO_terms"])
    for locus_tag, go_terms in sorted(locus_go.items()):
        writer.writerow([locus_tag, go_terms])
