import cooler
import pandas as pd
import numpy as np
import re
import csv
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--vcf", required=True, help="Input VCF file")
parser.add_argument("--cool", required=True, help="Cooler file path")
parser.add_argument("--promoter", required=True, help="Promoter BED file")
parser.add_argument("--out", required=True, help="Output prefix")

parser.add_argument("--percentile", type=float, default=99.99)
parser.add_argument("--dist", type=int, default=100000)

args = parser.parse_args()

VCF_FILE = args.vcf
COOL_PATH = args.cool
PROMOTER_BED = args.promoter
OUT_PREFIX = args.out

THRESHOLD_PERCENTILE = args.percentile
DIST_FILTER = args.dist

detail_file = OUT_PREFIX + "_detail.tsv"
summary_file = OUT_PREFIX + "_summary.tsv"


# -----------------------------
# read cool
# -----------------------------
clr = cooler.Cooler(COOL_PATH)
bins = clr.bins()[:]
bins["bin_id"] = np.arange(len(bins))

# -----------------------------
# read vcf
# -----------------------------

sv_list = []
with open(VCF_FILE) as f:
    for line in f:
        if line.startswith("#"):
            continue
        row = line.strip().split("\t")
        chrom = row[0]
        pos = int(row[1])
        sv_id = row[2]
        info = row[7]
        m_type = re.search("SVTYPE=([^;]+)", info)
        sv_type = m_type.group(1) if m_type else "NA"
        m_end = re.search("END=([0-9]+)", info)
        end = int(m_end.group(1)) if m_end else pos + 1
        sv_list.append([chrom, pos, end, sv_type, sv_id])

sv_df = pd.DataFrame(sv_list, columns=["chr","start","end","sv_type","sv_id"])

# -----------------------------
# read promoter BED
# -----------------------------

prom_df = pd.read_csv(PROMOTER_BED, sep="\t", header=None,
                      names=["chr","start","end","gene"])

promoter_bins = set()
bin_gene_map = {}

for _, p in prom_df.iterrows():

    sub = bins[bins["chrom"] == p["chr"]]

    overlap = sub[
        (sub["start"] < p["end"]) &
        (sub["end"] > p["start"])
    ]

    for b in overlap["bin_id"]:
        promoter_bins.add(b)

        if b not in bin_gene_map:
            bin_gene_map[b] = []

        bin_gene_map[b].append(p["gene"])


def find_sv_bins(chrom, start, end):
    sub = bins[bins["chrom"] == chrom]
    overlap = sub[(sub["start"] < end) & (sub["end"] > start)]
    return overlap["bin_id"].values


# -----------------------------
# Main function
# -----------------------------

results = []
summary = []
mat = clr.matrix(balance=True,sparse=True)
count = 0
sv_count = 0
for _, sv in sv_df.iterrows():
    sv_count += 1
    if sv["sv_type"] == "TRA":
        continue
    sv_bins = find_sv_bins(sv["chr"], sv["start"], sv["end"])
    promoter_hits = 0
    gene_list = []
    for b in sv_bins:
        
        vec = mat[b,:].toarray()[0]
        vec[b] = 0

        threshold = np.percentile(vec[~np.isnan(vec)], THRESHOLD_PERCENTILE)
        
        if threshold == 0:
            nonzero = vec[vec > 0]
            if len(nonzero) > 0:
                threshold = nonzero.min()
            else:
                continue
        
        sig_bins = np.where(vec >= threshold)[0]
        for sb in sig_bins:
            genes = []
            row = bins.iloc[sb]
            dist = abs(sb - b) * clr.binsize
            if dist < DIST_FILTER:
                continue

            if sb not in promoter_bins:
                continue

            genes = bin_gene_map.get(sb, [])
            row = bins.iloc[sb]
            for g in genes:
                if g not in gene_list:
                    gene_list.append(g)
                promoter_hits += 1
                results.append({
                    "sv_id": sv["sv_id"],
                    "sv_type":sv["sv_type"],
                    "sv_chr": sv["chr"],
                    "contact_bin": sb,
                    "contact_chr": row["chrom"],
                    "contact_start": row["start"],
                    "contact_end": row["end"],
                    "gene": g,
                    "contact_value": vec[sb]
                })
    if promoter_hits > 0:
        count += 1
        
    summary.append({
    "sv_id": sv["sv_id"],
    "has_promoter_contact": int(promoter_hits > 0),
    "gene_list": gene_list,
    "promoter_contact_count": promoter_hits
    })
        
print("Number of SV:",sv_count)
print("Number of SV with high contact promoter:",count)

# -----------------------------
# output
# -----------------------------
pd.DataFrame(results).to_csv(detail_file, sep="\t", index=False)
pd.DataFrame(summary).to_csv(summary_file, sep="\t", index=False)
