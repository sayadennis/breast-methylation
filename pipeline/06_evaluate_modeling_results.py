import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import upsetplot
from scipy.stats import false_discovery_control

# from upsetplot import UpSet

######################################
#### Load SeSAMe modeling results ####
######################################

din = "/projects/p30791/methylation/sesame_out/differential_methylation"
dout = "/projects/p30791/methylation/sesame_out/differential_methylation"
plot_dir = "/projects/p30791/methylation/plots"

dml_results = {
    "Normal": pd.read_csv(f"{din}/DML_results_refNormal.csv"),
    "CUB": pd.read_csv(f"{din}/DML_results_refCUB.csv"),
    "OQ": pd.read_csv(f"{din}/DML_results_refOQ.csv"),
    "AN": pd.read_csv(f"{din}/DML_results_refAN.csv"),
    "TU": pd.read_csv(f"{din}/DML_results_refTU.csv"),
}

tissue_types = ["TU", "AN", "OQ", "CUB", "Normal"]
p_thres = 0.01

dmr_results = {}
for ref, _ in dml_results.items():
    dmr_results[ref] = {}
    # Add results for age
    dmr_results[ref]["Age"] = pd.read_csv(f"{din}/DMR_results_Age_ref{ref}.csv")
    # Add results for tissue regions
    for comparison in tissue_types:
        if os.path.exists(f"{din}/DMR_results_Sample.Region{comparison}_ref{ref}.csv"):
            dmr_results[ref][comparison] = pd.read_csv(
                f"{din}/DMR_results_Sample.Region{comparison}_ref{ref}.csv"
            )

############################################
#### Upset Plot of number of DML probes ####
############################################

## Prepare data for Upset Plot
dml_significance = {}
upset_data = {}

for ref, _ in dml_results.items():
    comparison_columns = (
        dml_results[ref]
        .columns[
            [
                colname.startswith("Est_Sample.Region")
                for colname in dml_results[ref].columns
            ]
        ]
        .tolist()
    )  # elements of this list are like 'Est_Sample.RegionAN'
    comparisons = [
        x[len("Est_Sample.Region") :] for x in comparison_columns
    ]  # elements of this list are like 'AN'
    dml_significance[ref] = pd.DataFrame(
        index=dml_results[ref].Probe_ID, columns=comparisons
    )
    for tissue_type in comparisons:
        # Record significance WITH FDR CORRECTION
        p_vals = dml_results[ref][f"Pval_Sample.Region{tissue_type}"]
        p_vals.replace(0, 1e-320, inplace=True)  # fix underflow
        p_vals.fillna(0.999999, inplace=True)
        dml_significance[ref][tissue_type] = false_discovery_control(p_vals) < p_thres
    dml_significance[ref].rename(
        {
            "TU": "Tumor",
            "AN": "Adjacent normal",
            "OQ": "Opposite quadrant",
            "CUB": "Contralateral",
        },
        axis=1,
        inplace=True,
    )
    upset_data[ref] = (
        dml_significance[ref].groupby(list(dml_significance[ref].columns)).size()
    )

# Plot Upset
for ref, _ in dml_results.items():
    upsetplot.plot(upset_data[ref], sort_categories_by="input")
    plt.title(f"Number of probes differentially methylated compared to {ref}")
    plt.savefig(f"{plot_dir}/upset_ref{ref}.png")
    plt.close()

############################################################
#### Count pairwise difference and prep for Sankey Plot ####
############################################################

effect_thres = 0.2

## Count up- and down-methylated probes between adjacent tissue categories
tissue_types.reverse()
num_diff_probes = pd.DataFrame(0.0, index=tissue_types, columns=tissue_types)
trends = pd.DataFrame(
    index=dml_results["Normal"].Probe_ID,
    columns=["All probes", "Normal vs CUB", "CUB vs OQ", "OQ vs AN", "AN vs TU"],
)

for ref_category, comp_category in [
    ("Normal", "CUB"),
    ("CUB", "OQ"),
    ("OQ", "AN"),
    ("AN", "TU"),
]:
    p_vals = dml_results[ref_category][f"Pval_Sample.Region{comp_category}"]
    p_vals.index = dml_results[ref_category].Probe_ID
    # fix underflow (p=0) by replacing zero with non-zero minimum
    p_nonzero_min = np.min(p_vals.iloc[p_vals.values != 0])
    p_vals.replace(0, p_nonzero_min, inplace=True)
    # replace p=NaN with non-significant p
    p_vals.fillna(0.999999, inplace=True)
    # obtain slope to decipher trend
    slope = dml_results[ref_category][f"Est_Sample.Region{comp_category}"]
    slope.index = dml_results[ref_category].Probe_ID
    # find hypermethylated probes
    pos_sig = (false_discovery_control(p_vals) < p_thres) & (slope >= effect_thres)
    # find hypomethylated probes
    neg_sig = (false_discovery_control(p_vals) < p_thres) & (
        slope <= (-1) * effect_thres
    )
    # record information
    num_diff_probes.loc[ref_category, comp_category] = pos_sig.sum()
    num_diff_probes.loc[comp_category, ref_category] = neg_sig.sum()
    trends.loc[
        pos_sig.iloc[pos_sig.values].index, f"{ref_category} vs {comp_category}"
    ] = "up"
    trends.loc[
        neg_sig.iloc[neg_sig.values].index, f"{ref_category} vs {comp_category}"
    ] = "down"

num_diff_probes.to_csv(f"{din}/number_dml_probes_btwn_categories.csv")

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
probe_sets = {
    "CUB_down_from_Normal": list(trends.iloc[trends.CUB.values == "down"].index),
    "CUB_up_from_Normal": list(trends.iloc[trends.CUB.values == "up"].index),
    "AN_up_from_OQ": list(trends.iloc[trends.AN.values == "up"].index),
    "AN_down_from_OQ": list(trends.iloc[trends.AN.values == "down"].index),
    "AN_up_and_TU_down": list(
        trends.iloc[(trends.AN.values == "up") & (trends.TU.values == "down")].index
    ),
    "AN_up_and_TU_nd_or_up": list(
        trends.iloc[(trends.AN.values == "up") & (trends.TU.values != "down")].index
    ),
    "AN_down_and_TU_nd_or_down": list(
        trends.iloc[(trends.AN.values == "down") & (trends.TU.values != "up")].index
    ),
    "AN_down_and_TU_up": list(
        trends.iloc[(trends.AN.values == "down") & (trends.TU.values == "up")].index
    ),
    "TU_down_from_AN": list(trends.iloc[trends.TU.values == "down"].index),
    "TU_up_from_AN": list(trends.iloc[trends.TU.values == "up"].index),
}

for key, probelist in probe_sets.items():
    with open(f"{dout}/probe_set_{key}.txt", "w", encoding="utf-8") as f:
        for item in probelist:
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
segs = dmr_results["CUB"]["AN"][seg_cols].drop_duplicates(ignore_index=True)

# Segment length (kb) distribution
seglen = (segs.Seg_End - segs.Seg_Start).values.ravel()
for n in [1, 10, 50, 100, 250, 500, 750]:
    print(f"Number of segments larger than {n}kb: {sum(seglen>n*1000)}")

# Number of probes in a segment
n_probes_per_seg = dmr_results["CUB"]["AN"].groupby("Seg_ID").size().values.ravel()
for n in [2, 5, 10, 25, 50, 75]:
    print(f"Segment with >{n} probes: {sum(n_probes_per_seg>n)}")

# Which chromosomes do larger segments lie in?
size_thres = 50e3  # 50kb
n_segs_per_chrom = segs.iloc[seglen > size_thres, :].groupby("Seg_Chrm").size()
chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

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
