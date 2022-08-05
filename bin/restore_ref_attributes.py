#! /usr/bin/env python3
# Retrieve biotypes in new isoforms from input annotation
import argparse
from GTF import GTF, Attributes

parser = argparse.ArgumentParser(description="Restore ref attributes in novel gtf.")
parser.add_argument("-gtf", type=argparse.FileType("r"), help="gtf to restore", required=True)
parser.add_argument(
    "-ref",
    type=argparse.FileType("r"),
    help="ref annotation",
    required=True,
)
args = parser.parse_args()

g_attributes = {}

# Extract gene attributes
for record in GTF.parse_by_line(args.ref):
    if record.feature == "gene" or record.feature == "transcript":
        continue

    g_id = record["gene_id"]
    if g_id in g_attributes:
        continue

    attributes = {}
    attributes["transcript_biotype"] = record["gene_biotype"]
    for k, v in record.attributes.items():
        if "gene" in k:
            attributes[k] = v
    g_attributes[g_id] = attributes

# Merge new isoforms attributes with ref
for record in GTF.parse_by_line(args.gtf):
    record.attributes = Attributes({**record.attributes, **g_attributes[record["gene_id"]]})
    print(record)
