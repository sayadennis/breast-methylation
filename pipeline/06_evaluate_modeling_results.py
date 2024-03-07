import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

######################################
#### Load SeSAMe modeling results ####
######################################

din = "/projects/p30791/methylation/sesame_out/differential_methylation"
dout = "/projects/p30791/methylation/sesame_out/differential_methylation"
plot_dir = "/projects/p30791/methylation/plots"

ref_comp_pairs = [
    ("Normal", "CUB"),
    ("CUB", "OQ"),
    ("OQ", "AN"),
    ("AN", "TU"),
]

dml_results, dmr_results = {}, {}
for ref, comp in ref_comp_pairs:
    # DML results
    dml_results[f"{ref} vs {comp}"] = pd.read_csv(
        f"{din}/DML_results_ref{ref}_comp{comp}.csv"
    )
    # DMR results
    dmr_results[f"{ref} vs {comp}"] = pd.read_csv(
        f"{din}/DMR_results_Sample.Region{comp}_ref{ref}.csv"
    )

tissue_types = ["Normal", "CUB", "OQ", "AN", "TU"]
p_thres = 0.01
effect_thres = 0.1

############################################################
#### Count pairwise difference and prep for Sankey Plot ####
############################################################

## Count up- and down-methylated probes between adjacent tissue categories
num_diff_probes = pd.DataFrame(
    0.0,
    index=[f"{ref} vs {comp}" for ref, comp in ref_comp_pairs],
    columns=["Higher", "Lower"],
)
trends = pd.DataFrame(
    index=list(
        set().union(*(df["Probe_ID"] for df in dml_results.values()))
    ),  # all probes
    columns=["All probes"] + [f"{ref} vs {comp}" for ref, comp in ref_comp_pairs],
)

for ref, comp in ref_comp_pairs:
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
    trends.loc[pos_sig.iloc[pos_sig.values].index, f"{ref} vs {comp}"] = "up"
    trends.loc[neg_sig.iloc[neg_sig.values].index, f"{ref} vs {comp}"] = "down"

num_diff_probes.to_csv(f"{din}/number_dml_probes_btwn_categories.csv")
print(num_diff_probes.to_markdown())  # printing for Wiki

trends.dropna(axis=0, how="all", inplace=True)
trends = trends.loc[[x.startswith("cg") for x in trends.index], :]
trends.fillna("n.d.", inplace=True)

trends.to_csv(f"{din}/dml_up_down_pairwise_trends.csv")

## Print statements for Sankey plots
## https://developers.google.com/chart/interactive/docs/gallery/sankey
for current_node, next_node in [
    ("All probes", "Normal vs CUB"),
    ("Normal vs CUB", "CUB vs OQ"),
    ("CUB vs OQ", "OQ vs AN"),
    ("OQ vs AN", "AN vs TU"),
]:
    if current_node == "All probes":
        num_up = (trends[next_node] == "up").sum()
        num_down = (trends[next_node] == "down").sum()
        num_nd = (trends[next_node] == "n.d.").sum()
        print(f"       [ '{current_node}', '{next_node} - up', {num_up} ],")
        print(f"       [ '{current_node}', '{next_node} - down', {num_down} ],")
        print(f"       [ '{current_node}', '{next_node} - n.d.', {num_nd} ],")
    else:
        for reg1 in ["up", "down", "n.d."]:
            for reg2 in ["up", "down", "n.d."]:
                num = (
                    (trends[current_node] == reg1) & (trends[next_node] == reg2)
                ).sum()
                print(
                    f"       [ '{current_node} - {reg1}', '{next_node} - {reg2}', {num} ],"
                )

# Save interesting probe sets to TXT file
probe_sets = {}
for ref, comp in ref_comp_pairs:
    for trend in ["up", "down"]:
        probe_sets[f"{comp}_{trend}_from_{ref}"] = list(
            trends.iloc[trends[f"{ref} vs {comp}"].values == trend].index
        )

probe_sets["AN_down_and_TU_up"] = list(
    trends.iloc[
        (trends["OQ vs AN"].values == "down") & (trends["AN vs TU"].values == "up")
    ].index
)
probe_sets["AN_up_and_TU_down"] = list(
    trends.iloc[
        (trends["OQ vs AN"].values == "up") & (trends["AN vs TU"].values == "down")
    ].index
)

for key, probelist in probe_sets.items():
    with open(f"{dout}/probe_set_{key}.txt", "w", encoding="utf-8") as f:
        for item in probelist:
            f.write(f"{item}\n")

###########################################################
#### Analyze DML patterns across all tissue categories ####
###########################################################

#### Define set names and their specific patterns ####

set_patterns = {
    # Monotonic increase
    "Monotonic increase A": ["n.d.", "n.d.", "n.d.", "up"],
    "Monotonic increase B": ["up", "n.d.", "n.d.", "n.d."],
    "Monotonic increase C": ["up", "n.d.", "n.d.", "up"],
    "Monotonic increase D": ["n.d.", "n.d.", "up", "n.d."],
    "Monotonic increase E": ["up", "n.d.", "up", "n.d."],
    # Monotonic decrease
    "Monotonic decrease A": ["n.d.", "n.d.", "n.d.", "down"],
    "Monotonic decrease B": ["down", "n.d.", "n.d.", "n.d."],
    "Monotonic decrease C": ["down", "n.d.", "n.d.", "down"],
    "Monotonic decrease D": ["n.d.", "n.d.", "down", "n.d."],
    "Monotonic decrease E": ["down", "n.d.", "down", "n.d."],
    # Non-monotonic down & up, or "valley"
    "Non-monotonic valley A": ["down", "n.d.", "n.d.", "up"],
    "Non-monotonic valley B": ["n.d.", "n.d.", "down", "up"],
    "Non-monotonic valley C": ["down", "n.d.", "down", "up"],
    "Non-monotonic valley D": ["down", "n.d.", "up", "n.d."],
    "Non-monotonic valley E": ["down", "n.d.", "up", "up"],
    # Non-monotonic up & down, or "hill"
    "Non-monotonic hill A": ["up", "n.d.", "n.d.", "down"],
    "Non-monotonic hill B": ["n.d.", "n.d.", "up", "down"],
    "Non-monotonic hill C": ["up", "n.d.", "up", "down"],
    "Non-monotonic hill D": ["up", "n.d.", "down", "n.d."],
    # Non-monotonic mixed
    "Non-monotomic mixed A": ["down", "n.d.", "up", "down"],
    "Non-monotomic mixed B": ["up", "n.d.", "down", "up"],
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
    with open(f"{dout}/probe_set_{setname_file}.txt", "w", encoding="utf-8") as f:
        for probe in probes:
            f.write(f"{probe}\n")

##########################################
#### Analyze DML results by biomarker ####
##########################################

for tissue_category in ["CUB", "OQ", "AN", "TU"]:
    for biomarker in ["ER", "HER2"]:
        dml = pd.read_csv(f"{din}/DML_results_of_{tissue_category}_by_{biomarker}.csv")
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
                f"{dout}/probe_set_hyper_in_{tissue_category}_{biomarker}.txt",
                "w",
                encoding="utf-8",
            ) as f:
                for item in pos_sig_probes:
                    f.write(f"{item}\n")
        if len(neg_sig_probes) > 0:
            with open(
                f"{dout}/probe_set_hypo_in_{tissue_category}_{biomarker}.txt",
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

for ref, comp in ref_comp_pairs:
    segs = dmr_results[f"{ref} vs {comp}"][seg_cols].drop_duplicates(ignore_index=True)
    # Segment length (kb) distribution
    seglen = (segs.Seg_End - segs.Seg_Start).values.ravel()
    for n in [1, 10, 50, 100, 250, 500, 750]:
        print(f"Number of segments larger than {n}kb: {sum(seglen>n*1000)}")
    # Which chromosomes do larger segments lie in?
    size_thres = 50e3  # 50kb
    n_segs_per_chrom = segs.iloc[seglen > size_thres, :].groupby("Seg_Chrm").size()
    chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
    # Plot bar graph of number of large segments
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.bar(np.arange(len(chroms)), n_segs_per_chrom.loc[chroms])
    ax.set_xticks(np.arange(len(chroms)))
    ax.set_xticklabels(chroms, rotation=45, ha="right")
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Number of segments")
    ax.set_title(
        f"Distribution of large (>{int(size_thres/1e+3)}kb) segments across chromosomes"
    )
    plt.tight_layout()
    fig.savefig(f"{plot_dir}/num_large_segs_by_chrom_refCUB.png")
    plt.close()
