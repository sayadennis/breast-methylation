# pylint: disable=redefined-outer-name
# pylint: disable=duplicate-code

import json
import os
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter
from scipy.stats import false_discovery_control

######################################
#### Load SeSAMe modeling results ####
######################################

din_dm = "/projects/p30791/methylation/differential_methylation"
din_dv = "/projects/p30791/methylation/differential_variability/missMethylDiffVar"
dout_dm = din_dm
plot_dir = "/projects/p30791/methylation/plots/differential_methylation"

if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)

if not os.path.exists(f"{dout_dm}/all_probes_DM"):
    os.makedirs(f"{dout_dm}/all_probes_DM")

ref_comp_pairs_normals = [  # pairs to compare the normals
    ("CFN", "TU"),
    ("CUB", "TU"),
    ("OQ", "TU"),
    ("AN", "TU"),
]
ref_comp_pairs_tpx = [  # pairs along the tumor proximity axis
    ("CFN", "CUB"),
    ("CUB", "OQ"),
    ("OQ", "AN"),
    ("AN", "TU"),
]

dml_results, dmr_results = {}, {}
for ref, comp in set(ref_comp_pairs_tpx + ref_comp_pairs_normals):
    # DML results
    dml_results[f"{ref} vs {comp}"] = pd.read_csv(
        f"{din_dm}/DML_results_ref{ref}_comp{comp}.csv"
    )
    # DMR results
    dmr_results[f"{ref} vs {comp}"] = pd.read_csv(
        f"{din_dm}/DMR_results_Sample.Region{comp}_ref{ref}.csv"
    )
    # While we're here, write the lists of all probe sets
    dml_results[f"{ref} vs {comp}"].Probe_ID.to_csv(
        f"{dout_dm}/all_probes_DM/all_probes_{ref}_vs_{comp}.txt",
        header=False,
        index=False,
    )

tissue_types = ["CFN", "CUB", "OQ", "AN", "TU"]

with open(
    f"{str(Path.home())}/breast-methylation/pipeline/config.json", "r", encoding="utf-8"
) as f:
    config = json.load(f)

p_thres = config["p_thres"]
effect_thres = config["effect_thres"]

############################################################
#### Count pairwise difference and prep for Sankey Plot ####
############################################################

## Count hyper- and hypo-methylated probes between adjacent tissue categories along TPX
num_diff_probes = pd.DataFrame(
    0.0,
    index=[
        f"{ref} vs {comp}"
        for ref, comp in set(ref_comp_pairs_tpx + ref_comp_pairs_normals)
    ],
    columns=["Higher", "Lower"],
)
trends = pd.DataFrame(
    index=list(
        set().union(*(df["Probe_ID"] for df in dml_results.values()))
    ),  # all probes
    columns=["All probes"]
    + [
        f"{ref} vs {comp}"
        for ref, comp in set(ref_comp_pairs_tpx + ref_comp_pairs_normals)
    ],
)

for ref, comp in set(ref_comp_pairs_tpx + ref_comp_pairs_normals):
    p_vals = dml_results[f"{ref} vs {comp}"][f"Pval_Sample.Region{comp}"]
    p_vals.index = dml_results[f"{ref} vs {comp}"].Probe_ID
    # fix underflow (p=0) by replacing zero with non-zero minimum
    p_nonzero_min = np.min(p_vals.iloc[p_vals.values != 0])
    p_vals.replace(0, p_nonzero_min, inplace=True)
    # replace p=NaN with non-significant p
    p_vals.fillna(0.999999, inplace=True)
    # obtain slope to decipher trend
    slope = dml_results[f"{ref} vs {comp}"][f"Est_Sample.Region{comp}"]
    slope.index = dml_results[f"{ref} vs {comp}"].Probe_ID
    # find hypermethylated probes
    pos_sig = (false_discovery_control(p_vals) < p_thres) & (slope >= effect_thres)
    # find hypomethylated probes
    neg_sig = (false_discovery_control(p_vals) < p_thres) & (
        slope <= (-1) * effect_thres
    )
    # record information
    num_diff_probes.loc[f"{ref} vs {comp}", "Higher"] = pos_sig.sum()
    num_diff_probes.loc[f"{ref} vs {comp}", "Lower"] = neg_sig.sum()
    trends.loc[pos_sig.iloc[pos_sig.values].index, f"{ref} vs {comp}"] = "hyper"
    trends.loc[neg_sig.iloc[neg_sig.values].index, f"{ref} vs {comp}"] = "hypo"

num_diff_probes.to_csv(f"{din_dm}/number_dml_probes_btwn_categories.csv")
print(num_diff_probes.to_markdown())  # printing for Wiki

trends.dropna(axis=0, how="all", inplace=True)
trends = trends.loc[[x.startswith("cg") for x in trends.index], :]
trends.fillna("n.d.", inplace=True)

trends[
    ["All probes"] + [f"{ref} vs {comp}" for ref, comp in ref_comp_pairs_tpx]
].to_csv(f"{din_dm}/dml_hyper_hypo_pairwise_trends.csv")

##############################################
#### Bar graphs of number of DM/DV probes ####
##############################################

probesets = {}
probesets["DM"] = {}
probesets["DV"] = {}

for ref, comp in set(ref_comp_pairs_tpx + ref_comp_pairs_normals):
    for trend in ["hyper", "hypo"]:
        # First get DM probes
        with open(
            f"{din_dm}/probe_set_{trend}_ref{ref}_comp{comp}.txt", "r", encoding="utf-8"
        ) as f:
            probesets["DM"][f"{trend}_ref{ref}_comp{comp}"] = [
                x.strip() for x in f.readlines()
            ]
        # Next get DV probes
        with open(
            f"{din_dv}/{trend}DV_missMethyl_{ref}_vs_{comp}.txt", "r", encoding="utf-8"
        ) as f:
            dv_probes = [x.strip() for x in f.readlines()]
            dv_but_not_dm_probes = list(
                set(dv_probes) - set(probesets["DM"][f"{trend}_ref{ref}_comp{comp}"])
            )
            probesets["DV"][f"{trend}_ref{ref}_comp{comp}"] = dv_but_not_dm_probes


def format_ticks(x, pos):  # pylint: disable=unused-argument
    """
    Function to facilitate changing format of y tick labels.
    """
    return "{:,.0f}".format(x)  # pylint: disable=consider-using-f-string


colors = ["slategray", "cornflowerblue", "forestgreen", "goldenrod"]

for identifier, pairs in zip(
    ["normals", "tpx"], [ref_comp_pairs_normals, ref_comp_pairs_tpx]
):
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(5.5, 7))
    # Plot bar graphs
    for i, analysis in enumerate(["DM", "DV"]):
        for j, trend in enumerate(["Hyper", "Hypo"]):
            axs[i, j].bar(
                np.arange(4),
                [
                    len(probesets[analysis][f"{trend.lower()}_ref{ref}_comp{comp}"])
                    for ref, comp in pairs
                ],
                width=1.0,
                edgecolor="black",
                linewidth=1,
                color=colors,
            )
            axs[i, j].spines["top"].set_visible(False)
            axs[i, j].spines["right"].set_visible(False)
            axs[i, j].set_xticks(np.arange(4))
            axs[i, j].set_xticklabels(
                [f"{ref} vs {comp}" for ref, comp in pairs], ha="right", rotation=30
            )
            axs[i, j].yaxis.set_major_formatter(FuncFormatter(format_ticks))
            if j == 0:
                axs[i, j].set_ylabel("Number of Probes")
    # Add figure labels
    fig.text(0.03, 0.75, "DM", rotation=90, va="center", ha="center", fontsize=16)
    fig.text(0.03, 0.28, "DV", rotation=90, va="center", ha="center", fontsize=16)
    fig.text(0.38, 0.97, "Hyper", rotation=0, va="center", ha="center", fontsize=16)
    fig.text(0.82, 0.97, "Hypo", rotation=0, va="center", ha="center", fontsize=16)
    # Set margins for text
    plt.tight_layout()
    plt.subplots_adjust(left=0.22, top=0.95)
    # Save
    fig.savefig(f"{plot_dir}/bar_num_dm_dv_{identifier}.png")
    plt.close()

#########################################################
#### Group DML patterns across all tissue categories ####
#########################################################

# From here, we only focus on comparisons along the TPX
trends = trends[
    ["All probes"] + [f"{ref} vs {comp}" for ref, comp in ref_comp_pairs_tpx]
]

#### Define set names and their specific patterns ####
set_patterns = {
    # Monotonic increase
    "Monotonic increase A": ["n.d.", "n.d.", "n.d.", "hyper"],
    "Monotonic increase B": ["hyper", "n.d.", "n.d.", "n.d."],
    "Monotonic increase C": ["hyper", "n.d.", "n.d.", "hyper"],
    "Monotonic increase D": ["n.d.", "n.d.", "hyper", "n.d."],
    "Monotonic increase E": ["hyper", "n.d.", "hyper", "n.d."],
    # Monotonic decrease
    "Monotonic decrease A": ["n.d.", "n.d.", "n.d.", "hypo"],
    "Monotonic decrease B": ["hypo", "n.d.", "n.d.", "n.d."],
    "Monotonic decrease C": ["hypo", "n.d.", "n.d.", "hypo"],
    "Monotonic decrease D": ["n.d.", "n.d.", "hypo", "n.d."],
    "Monotonic decrease E": ["hypo", "n.d.", "hypo", "n.d."],
    # Non-monotonic down & up, or "valley"
    "Non-monotonic valley A": ["hypo", "n.d.", "n.d.", "hyper"],
    "Non-monotonic valley B": ["n.d.", "n.d.", "hypo", "hyper"],
    "Non-monotonic valley C": ["hypo", "n.d.", "hypo", "hyper"],
    "Non-monotonic valley D": ["hypo", "n.d.", "hyper", "n.d."],
    "Non-monotonic valley E": ["hypo", "n.d.", "hyper", "hyper"],
    # Non-monotonic up & down, or "hill"
    "Non-monotonic hill A": ["hyper", "n.d.", "n.d.", "hypo"],
    "Non-monotonic hill B": ["n.d.", "n.d.", "hyper", "hypo"],
    "Non-monotonic hill C": ["hyper", "n.d.", "hyper", "hypo"],
    "Non-monotonic hill D": ["hyper", "n.d.", "hypo", "n.d."],
    # Non-monotonic mixed
    "Non-monotonic mixed A": ["hypo", "n.d.", "hyper", "hypo"],
    "Non-monotonic mixed B": ["hyper", "n.d.", "hypo", "hyper"],
}

#### Get probe names with each patter ####
probesets = {}

for setname, pattern in set_patterns.items():
    subdata = trends.iloc[
        [
            np.all(trends.loc[i, :].values.ravel()[1:] == np.array(pattern))
            for i in trends.index
        ],
        :,
    ]
    probesets[setname] = list(subdata.index)

#### Write probe sets to TXT files ####
for setname, probes in probesets.items():
    setname_file = re.sub(r"[-\s]", "_", setname)
    with open(f"{dout_dm}/probe_set_{setname_file}.txt", "w", encoding="utf-8") as f:
        for probe in probes:
            f.write(f"{probe}\n")


#### Plot arrow graph to show global patterns across TPX ####

pattern_names = [
    "Monotonic increase",
    "Monotonic decrease",
    "Non-monotonic \ndown & up",
    "Non-monotonic \nup & down",
    "Non-monotonic mixed",
]
n_patterns = {
    "Monotonic increase": 5,
    "Monotonic decrease": 5,
    "Non-monotonic \ndown & up": 5,
    "Non-monotonic \nup & down": 4,
    "Non-monotonic mixed": 2,
}
max_patterns = np.max(list(n_patterns.values()))

coords = {
    "nd": {"x": -0.25, "y": 0.0, "dx": 0.5, "dy": 0.0},
    "hyper": {"x": -0.25, "y": -0.15, "dx": 0.5, "dy": 0.3},
    "hypo": {"x": -0.25, "y": 0.15, "dx": 0.5, "dy": -0.3},
}

colors = {
    "nd": "gray",
    "hyper": "magenta",
    "hypo": "cyan",
}

fig, axs = plt.subplots(
    nrows=max_patterns,
    ncols=len(pattern_names),
    figsize=(3 * len(pattern_names), 1.1 * max_patterns),
)


def add_num(i, j, text):
    axs[i, j].text(
        0.05,
        0.85,
        f"{text} probes",
        transform=axs[i, j].transAxes,
        fontsize=10,
        verticalalignment="bottom",
        multialignment="center",
        bbox={
            "boxstyle": "round",
            "facecolor": "white",
            "edgecolor": "red",
            "alpha": 0.7,
        },
    )


def plot_arrows_and_num(i: int, j: int, pattern: list, num: int):
    for k, trend in enumerate(pattern):
        axs[i, j].arrow(
            coords[trend]["x"] + 0.5 + k,
            coords[trend]["y"],
            coords[trend]["dx"],
            coords[trend]["dy"],
            color=colors[trend],
            width=0.05,
            head_width=0.2,
            length_includes_head=True,
        )
        add_num(i, j, f"{num:,}")


## Monotonic increase
plot_arrows_and_num(0, 0, ["nd", "nd", "nd", "hyper"], 33698)
plot_arrows_and_num(1, 0, ["hyper", "nd", "nd", "nd"], 739)
plot_arrows_and_num(2, 0, ["hyper", "nd", "nd", "hyper"], 142)
plot_arrows_and_num(3, 0, ["nd", "nd", "hyper", "hyper"], 38)
plot_arrows_and_num(4, 0, ["hyper", "nd", "hyper", "nd"], 6)

## Monotonic decrease
plot_arrows_and_num(0, 1, ["nd", "nd", "nd", "hypo"], 53800)
plot_arrows_and_num(1, 1, ["hypo", "nd", "nd", "nd"], 34882)
plot_arrows_and_num(2, 1, ["hypo", "nd", "nd", "hypo"], 7672)
plot_arrows_and_num(3, 1, ["nd", "nd", "hypo", "hypo"], 5)
plot_arrows_and_num(4, 1, ["hypo", "nd", "hypo", "nd"], 6)

## Non-monotonic down and up
plot_arrows_and_num(0, 2, ["hypo", "nd", "nd", "hyper"], 6379)
plot_arrows_and_num(1, 2, ["nd", "nd", "hypo", "hyper"], 27)
plot_arrows_and_num(2, 2, ["hypo", "nd", "hypo", "hyper"], 45)
plot_arrows_and_num(3, 2, ["hypo", "nd", "hyper", "nd"], 25)
plot_arrows_and_num(4, 2, ["hypo", "nd", "hyper", "hyper"], 3)

## Non-monotonic up and down
plot_arrows_and_num(0, 3, ["hyper", "nd", "nd", "hypo"], 56)
plot_arrows_and_num(1, 3, ["nd", "nd", "hyper", "hypo"], 76)
plot_arrows_and_num(2, 3, ["hyper", "nd", "hyper", "hypo"], 8)
plot_arrows_and_num(3, 3, ["hyper", "nd", "hypo", "nd"], 2)

## Non-monotonic mixed
plot_arrows_and_num(0, 4, ["hypo", "nd", "hyper", "hypo"], 29)
plot_arrows_and_num(1, 4, ["hyper", "nd", "hypo", "hyper"], 4)

for i in range(max_patterns):
    for j, pattern_name in enumerate(pattern_names):
        axs[i, j].spines["right"].set_visible(False)
        axs[i, j].spines["top"].set_visible(False)
        axs[i, j].spines["bottom"].set_visible(False)
        axs[i, j].spines["left"].set_visible(False)
        axs[i, j].set_ylim(-0.5, 0.5)
        axs[i, j].set_yticks([])
        axs[i, j].set_yticklabels([])
        axs[i, j].set_xlim(-0.5, 4.5)
        axs[i, j].set_xticks(np.arange(5))
        # add horizontal lines
        axs[i, j].axhline(
            0.0, linestyle="--", linewidth=1.0, c="grey", zorder=1, alpha=0.5
        )
        if i == (max_patterns - 1):
            axs[i, j].set_xticklabels(
                ["CFN", "CUB", "OQ", "AN", "TU"], fontsize=12, ha="center", rotation=30
            )
        else:
            axs[i, j].set_xticklabels([])
        # add vertical lines
        for x in [0, 1, 2, 3, 4]:
            axs[i, j].axvline(x, linestyle="--", linewidth=0.5, c="grey")
        # set aspect ratio
        axs[i, j].set_aspect("equal")
        # if top row, set pattern titles
        if i == 0:
            axs[i, j].set_title(pattern_name + "\n", fontsize=14)
        # if left-most column, add labels
        if j == 0:
            axs[i, j].set_ylabel(
                ["A", "B", "C", "D", "E"][i] + "  ",
                fontsize=17,
                rotation=0,
                va="center",
                ha="right",
            )

plt.tight_layout()
plt.subplots_adjust(left=0.05, top=0.90)
fig.savefig(f"{plot_dir}/patterns_across_tissues_arrows.png")
plt.close()

###################################################################
#### For smaller probe sets from global patterns, map to genes ####
###################################################################

probemeta = pd.read_csv(
    "/projects/p30791/methylation/copied_from_b1122/data/meta/EPIC-8v2-0_A1.csv",
    skiprows=7,
)

for pattern_name in [
    "Monotonic increase D",
    "Monotonic increase E",
    "Monotonic decrease D",
    "Monotonic decrease E",
]:
    probenames = probesets[pattern_name]
    genes = probemeta.iloc[
        [probe_id in probenames for probe_id in probemeta["Name"]], :
    ]["UCSC_RefGene_Name"].unique()
    print(f"\n\n{pattern_name}")
    print(genes)

##########################################
#### Analyze DML results by biomarker ####
##########################################

for tissue_category in ["CUB", "OQ", "AN", "TU"]:
    for biomarker in ["ER", "HER2"]:
        dml = pd.read_csv(
            f"{din_dm}/DML_results_of_{tissue_category}_by_{biomarker}.csv"
        )
        dml = dml.iloc[[x.startswith("cg") for x in dml.Probe_ID], :]
        # obtain slope and p-values
        slope = dml[f"Est_{biomarker}Positive"]
        slope.index = dml.Probe_ID
        p_vals = dml[f"Pval_{biomarker}Positive"]
        p_vals.index = dml.Probe_ID
        # fix underflow (p=0) by replacing zero with non-zero minimum
        p_nonzero_min = np.min(p_vals.iloc[p_vals.values != 0])
        p_vals.replace(0, p_nonzero_min, inplace=True)
        # replace p=NaN with non-significant p
        p_vals.fillna(0.999999, inplace=True)
        # find hypermethylated probes
        pos_sig = (false_discovery_control(p_vals) < p_thres) & (slope >= effect_thres)
        pos_sig_probes = list(pos_sig.iloc[pos_sig.values].index)
        # find hypomethylated probes
        neg_sig = (false_discovery_control(p_vals) < p_thres) & (
            slope <= (-1) * effect_thres
        )
        neg_sig_probes = list(neg_sig.iloc[neg_sig.values].index)
        if len(pos_sig_probes) > 0:
            with open(
                f"{dout_dm}/probe_set_hyper_in_{tissue_category}_{biomarker}.txt",
                "w",
                encoding="utf-8",
            ) as f:
                for item in pos_sig_probes:
                    f.write(f"{item}\n")
        if len(neg_sig_probes) > 0:
            with open(
                f"{dout_dm}/probe_set_hypo_in_{tissue_category}_{biomarker}.txt",
                "w",
                encoding="utf-8",
            ) as f:
                for item in neg_sig_probes:
                    f.write(f"{item}\n")

#############################################
#### Get more information about segments ####
#############################################

seg_cols = [
    "Seg_ID",
    "Seg_Chrm",
    "Seg_Start",
    "Seg_End",
    "Seg_Est",
    "Seg_Pval",
    "Seg_Pval_adj",
]

for ref, comp in ref_comp_pairs_tpx:
    print(f"#### Comparison - {ref} vs {comp} ####")
    segs = dmr_results[f"{ref} vs {comp}"][seg_cols].drop_duplicates(ignore_index=True)
    segs.set_index("Seg_ID", drop=True, inplace=True)
    # Number of probes per segment
    num_probes = dmr_results[f"{ref} vs {comp}"].groupby("Seg_ID").size()
    # Segment length (kb) distribution
    seglen = pd.Series(segs.Seg_End - segs.Seg_Start, index=segs.index)
    for n in [1, 10, 50, 100, 250, 500, 750]:
        print(f"Number of segments larger than {n}kb: {sum(seglen>n*1000)}")
    # Which chromosomes do larger segments lie in?
    # size_thres = 50e3  # 50kb
    probe_thres = 5
    n_segs_per_chrom = (
        segs.loc[num_probes.iloc[num_probes.values > probe_thres].index, :]
        .groupby("Seg_Chrm")
        .size()
    )
    chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
    # Plot bar graph of number of large segments
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.bar(np.arange(len(chroms)), n_segs_per_chrom.loc[chroms])
    ax.set_xticks(np.arange(len(chroms)))
    ax.set_xticklabels(chroms, rotation=45, ha="right")
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Number of segments")
    ax.set_title(
        # f"Distribution of large (>{int(size_thres/1e+3)}kb) segments across chromosomes"
        f"Distribution of large (>{probe_thres} probes) segments across chromosomes"
    )
    plt.tight_layout()
    fig.savefig(f"{plot_dir}/num_large_segs_by_chrom_ref{ref}_comp{comp}.png")
    plt.close()
