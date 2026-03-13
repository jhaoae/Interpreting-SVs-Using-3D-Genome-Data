import pandas as pd
import numpy as np
import random
from scipy.stats import wilcoxon
import argparse

parser = argparse.ArgumentParser(description='Load and merge tumor and control datasets.')
parser.add_argument('tumor', type=str, help='Path to the poreC tumor data insulation score (e.g., HCC1937_score.bedgraph)')
parser.add_argument('normal', type=str, help='Path to the poreC control data insulation score (e.g., HCC1937_BL_score.bedgraph)')
parser.add_argument('sv', type=str, help='Path to the SV Bed file ()')
parser.add_argument('sv', type=str, help='Path to the SV BED file (6 columns: chr, start, end, sv_type, sv_start, sv_end. The first three columns describe a ±50 kb genomic window around the SV breakpoints, and the last three columns describe the SV itself.)')

# load insulation profiles
tumor_ins = pd.read_csv(
    args.tumor,
    sep="\t",
    names=["chr","start","end","score"]
)

bl_ins = pd.read_csv(
    args.normal,
    sep="\t",
    names=["chr","start","end","score"]
)

bins = tumor_ins.merge(
    bl_ins,
    on=["chr", "start", "end"],
    suffixes=("_tumor", "_bl")
)

sv = pd.read_csv(
    args.sv_file,
    sep="\t",
    names=["chr","start","end","sv_type","sv_start","sv_end"]
)

def bins_overlapping_sv(bins_df, chrom, sv_start, sv_end):
    return bins_df[
        (bins_df.chr == chrom) &
        (bins_df.end > sv_start) &
        (bins_df.start < sv_end)
    ]


results = []

for _, row in sv.iterrows():
    chrom, start, end, sv_type, sv_start, sv_end = row

    #sub_sv = bins_overlapping_sv(bins, chrom, sv_start, sv_end)
    sub_win = bins_overlapping_sv(bins, chrom, start, end)

    #if len(sub_sv) == 0:
    #    continue
    if len(sub_win) == 0:
        continue

    win_score_tumor = sub_win.score_tumor
    win_score_bl = sub_win.score_bl

    q = 0.05

    delta_win =  (
        sub_win.score_tumor.quantile(q) -
        sub_win.score_bl.quantile(q)
    )

    results.append({
        "sv_id": f"{sv_type}_{chrom}_{sv_start}_{sv_end}",
        "chr": chrom,
        "start": start,
        "end": end,
        "sv_type": sv_type,
        "delta": delta_win,
        "n_bins": len(sub_win),
        "scores_tumor": sub_win.score_tumor.values,
        "scores_bl": sub_win.score_bl.values
    })


result_df = pd.DataFrame(results)

result_df.to_csv(
    "SV_insulation_change_detail_score.tsv",
    sep="\t",
    index=False
)


stats = []

for r in results:
    stat, p = wilcoxon(
        r["scores_tumor"],
        r["scores_bl"]
    )

    stats.append({
        "sv_id": r["sv_id"],
        "chr": r["chr"],
        "start": r["start"],
        "end": r["end"],
        "sv_type": r["sv_type"],
        "delta": r["delta"],
        "abs_delta": abs(r["delta"]),
        "pvalue": p,
        "n_bins": r["n_bins"]
    })

stat_df = pd.DataFrame(stats)

sig_df = (
    stat_df
    .query("pvalue < 0.05")
    .sort_values("abs_delta", ascending=False)
)

stat_df.to_csv("SV_insulation_change.tsv", sep="\t", index=False)
sig_df.to_csv("SV_insulation_significant_change.tsv",sep="\t",index=False)

