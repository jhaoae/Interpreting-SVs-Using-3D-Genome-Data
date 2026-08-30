import pandas as pd
from scipy.stats import wilcoxon
import argparse

parser = argparse.ArgumentParser(description='Load and merge tumor and control datasets.')
parser.add_argument('tumor', type=str, help='Path to the poreC tumor data insulation score (e.g., HCC1937_score.bedgraph)')
parser.add_argument('normal', type=str, help='Path to the poreC control data insulation score (e.g., HCC1937_BL_score.bedgraph)')
parser.add_argument('sv_file', type=str, help='Path to the SV BED file (6 columns: chr, start, end, sv_type, sv_start, sv_end. The first three columns describe a window around the SV breakpoints, and the last three columns describe the SV itself.)')
parser.add_argument('--out-prefix', default='SV_insulation', help='Output prefix (default: SV_insulation)')
parser.add_argument('--quantile', type=float, default=0.05, help='Insulation-score quantile used to calculate the change (default: 0.05)')

args = parser.parse_args()

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

    delta_win =  (
        sub_win.score_tumor.quantile(args.quantile) -
        sub_win.score_bl.quantile(args.quantile)
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
    f"{args.out_prefix}_change_detail_score.tsv",
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

stat_columns = [
    "sv_id", "chr", "start", "end", "sv_type", "delta",
    "abs_delta", "pvalue", "n_bins"
]
stat_df = pd.DataFrame(stats, columns=stat_columns)

sig_df = (
    stat_df
    .query("pvalue < 0.05")
    .sort_values("abs_delta", ascending=False)
)

stat_df.to_csv(f"{args.out_prefix}_change.tsv", sep="\t", index=False)
sig_df.to_csv(f"{args.out_prefix}_significant_change.tsv",sep="\t",index=False)
