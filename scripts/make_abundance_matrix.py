import csv
import sys
from pathlib import Path


if len(sys.argv) < 3:
    print("Usage: python3 make_abundance_matrix.py <bracken1> <bracken2> ...")
    sys.exit(1)


input_files = [Path(x) for x in sys.argv[1:]]
output_file = Path("abundance_matrix.tsv")


# ------------------------------------------------------------
# Read Bracken files
# ------------------------------------------------------------

sample_tables = {}

for file in input_files:

    sample = file.name.replace(".bracken.tsv", "")

    taxa = {}

    with open(file, newline="") as f:

        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:

            taxon = row["name"]
            abundance = float(row["fraction_total_reads"])

            taxa[taxon] = abundance

    sample_tables[sample] = taxa


# ------------------------------------------------------------
# Collect all taxa
# ------------------------------------------------------------

all_taxa = set()

for taxa in sample_tables.values():
    all_taxa.update(taxa.keys())

all_taxa = sorted(all_taxa)
samples = sorted(sample_tables.keys())


# ------------------------------------------------------------
# Write abundance matrix
# ------------------------------------------------------------

with open(output_file, "w", newline="") as f:

    writer = csv.writer(f, delimiter="\t")

    writer.writerow(["Taxon"] + samples)

    for taxon in all_taxa:

        row = [taxon]

        for sample in samples:

            value = sample_tables[sample].get(taxon, 0.0)

            row.append(f"{value:.6f}")

        writer.writerow(row)


print()
print("==============================================")
print("ABUNDANCE MATRIX")
print("==============================================")
print()
print(f"Samples : {len(samples)}")
print(f"Taxa    : {len(all_taxa)}")
print()
print(f"Output  : {output_file}")
print()
print("Abundance matrix construction complete.")
print("==============================================")
