import pandas as pd
import numpy as np
import random
from scipy.stats import wilcoxon
import argparse

parser = argparse.ArgumentParser(description='Load and process genomic data.')
parser.add_argument('sv_file', type=str, help='Path to the SV insulation significant file (SV_insulation_significant_change.tsv)')
parser.add_argument('nosv_file', type=str, help='Path to the no SV region file (3 columns: chr, start, end)')
parser.add_argument('tumor_file', type=str, help='Path to the poreC tumor data insulation score')
parser.add_argument('bl_file', type=str, help='Path to the poreC control data insulation score')

args = parser.parse_args()

sv_df = pd.read_csv(
    args.sv_file,
    sep="\t"
)

nosv = pd.read_csv(
    args.nosv_file,
    sep="\t",
    names=["chr", "start", "end"]
)

tumor_ins = pd.read_csv(
    args.tumor_file,
    sep="\t",
    names=["chr", "start", "end", "score_tumor"]
)

bl_ins = pd.read_csv(
    args.bl_file,
    sep="\t",
    names=["chr", "start", "end", "score_bl"]
)


bins = tumor_ins.merge(
    bl_ins,
    on=["chr", "start", "end"]
)

BIN_SIZE = 10000
N_BG = 500   # number of background windows per SV


def bins_in_window(bins_df, chrom, start, end):
    """Return bins overlapping a genomic window"""
    return bins_df[
        (bins_df.chr == chrom) &
        (bins_df.end > start) &
        (bins_df.start < end)
    ]

def sample_background_window(chrom, win_len, merged_df):
    """Sample one background window from nosv regions"""
    
    candidates = merged_df[merged_df.len >= win_len]
    if len(candidates) == 0:
        return None

    region = candidates.sample(1).iloc[0]
    
    bg_start = random.randrange(
        region.start,
        region.end - win_len+1,
        BIN_SIZE
    )
    bg_end = bg_start + win_len
    return bg_start, bg_end

def delta_quantile(tumor_scores, bl_scores, q=0.05):
    return tumor_scores.quantile(q) - bl_scores.quantile(q)


results = []

for _, sv in sv_df.iterrows():

    chrom = sv["chr"]
    sv_start = sv["start"]
    sv_end = sv["end"]
    sv_id = sv["sv_id"]
    sv_delta = sv["abs_delta"]
    
    win_len = sv_end - sv_start
    # ---- background ----
    nosv_chr = nosv[nosv.chr == chrom]
    bg_deltas = []

    attempts = 0
    
    nosv_chr = nosv_chr.sort_values("start")

    merged = []
    cur_start = None
    cur_end = None

    for _, row in nosv_chr.iterrows():
        s, e = row.start, row.end

        if cur_start is None:
            cur_start, cur_end = s, e
        elif s <= cur_end:  
            cur_end = max(cur_end, e)
        else:
            merged.append((cur_start, cur_end))
            cur_start, cur_end = s, e

    if cur_start is not None:
        merged.append((cur_start, cur_end))

    merged_df = pd.DataFrame(merged, columns=["start", "end"])
    merged_df["len"] = merged_df.end - merged_df.start
    
    
    while len(bg_deltas) < N_BG and attempts < 1000:
        attempts += 1

        bg_win = sample_background_window(
            chrom, win_len, merged_df
        )
        if bg_win is None:
            continue

        bg_start, bg_end = bg_win

        bg_bins = bins_in_window(
            bins, chrom, bg_start, bg_end
        )

        if len(bg_bins) < 3:
            continue

        bg_delta = delta_quantile(
            bg_bins.score_tumor,
            bg_bins.score_bl
        )
        bg_deltas.append(bg_delta)
     
    if len(bg_deltas) < 3:
        continue

    # ---- permutation-style p-value ----
    bg_deltas = np.array(bg_deltas)

    p_empirical = (
        np.sum(np.abs(bg_deltas) >= sv_delta) + 1
    ) / (len(bg_deltas) + 1)

    results.append({
        "sv_id": sv_id,
        "chr": chrom,
        "sv_start": sv_start,
        "sv_end": sv_end,
        "sv_delta": sv_delta,
        "bg_delta": bg_deltas,
        "bg_mean_delta": bg_deltas.mean(),
        "bg_std_delta": bg_deltas.std(),
        "empirical_p": p_empirical,
        "n_bg": len(bg_deltas)
    })

res_df = pd.DataFrame(results)

res_df.to_csv(
    "SV_insulation_permutation_test.tsv",
    sep="\t",
    index=False
)

print("Finished background testing.")
print(res_df.head())

