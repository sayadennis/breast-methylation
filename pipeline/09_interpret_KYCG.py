import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib_venn import venn2

################################################################
#### Create gene-only TXT files for gene enrichment results ####
################################################################

din = "/projects/p30791/methylation/differential_methylation/KYCG"
dout = din

probe_set_names = [
    "hyper_ER-_refAN_compTU",
    "hyper_ER+_refAN_compTU",
    "hypo_ER-_refAN_compTU",
    "hypo_ER+_refAN_compTU",
    "hyper_ER-_refHDB_compCUB",
    "hyper_ER+_refHDB_compCUB",
    "hypo_ER-_refHDB_compCUB",
    "hypo_ER+_refHDB_compCUB",
    "hyper_refHDB_compCUB",
    "hypo_refHDB_compCUB",
    "hyper_refOQ_compAN",
    "hypo_refOQ_compAN",
    "hyper_refAN_compTU",
    "hypo_refAN_compTU",
]

for setname in probe_set_names:
    genes_df = pd.read_csv(f"{din}/testEnrichment_{setname}_genes.csv")
    genes = genes_df.gene_name.values
    with open(f"{dout}/genes_{setname}.txt", "w") as f:
        for gene in genes:
            f.write(f"{gene}\n")

###########################################################
#### Plot horizontal bar graphs showing gene/probe overlaps
#### driving TFBS hits in HDB/CUB and AN/TU comparison ####
###########################################################

probe_or_gene = "gene"
din = "/projects/p30791/methylation/differential_methylation"
dout = "/projects/p30791/methylation/plots"

data = pd.read_csv(f"{din}/TFBS_hits_{probe_or_gene}_overlaps.csv")

chunks = data[
    ["left_trend", "right_trend", "left_comparison", "right_comparison"]
].drop_duplicates(ignore_index=True)
chunks = [tuple(x) for x in chunks.itertuples(index=False)]
sizes = [sum(all(data.iloc[i, :4].values == x) for i in data.index) for x in chunks]

height_ratios = tuple(2 * size for size in sizes)
figsize = (sum(height_ratios) / 4, sum(height_ratios) / 5)

viridis_cmap = matplotlib.colormaps["viridis"]
colors = [viridis_cmap(index) for index in np.linspace(0, 1, 3)]

left_color = colors[0]
overlap_color = colors[1]
right_color = colors[2]

fig, axs = plt.subplots(
    nrows=len(height_ratios),
    ncols=1,
    height_ratios=height_ratios,
    figsize=figsize,
)

for i, (chunk, size) in enumerate(zip(chunks, sizes)):
    ## Define the beginning and ending rows of `data`
    if i == 0:
        start = 0
        end = size
    elif i == len(sizes) - 1:
        start = sum(sizes[:i])
        end = sum(sizes)
    else:
        start = sum(sizes[:i])
        end = sum(sizes[: i + 1])
    ## Set left and right only labels
    left_only_label = (
        chunk[2].split("_")[0][3:] + " vs. " + chunk[2].split("_")[1][4:] + " only"
    )
    right_only_label = (
        chunk[3].split("_")[0][3:] + " vs. " + chunk[3].split("_")[1][4:] + " only"
    )
    ## Plot
    axs[i].barh(
        data.TF.iloc[start:end].values,
        data.left_only_size.iloc[start:end].values,
        label=left_only_label,
        edgecolor="black",
        color=left_color,
    )
    axs[i].barh(
        data.TF.iloc[start:end].values,
        data.overlap_size.iloc[start:end].values,
        label="Overlap",
        edgecolor="black",
        left=data.left_only_size.iloc[start:end].values,
        color=overlap_color,
    )
    axs[i].barh(
        data.TF.iloc[start:end].values,
        data.right_only_size.iloc[start:end].values,
        label=right_only_label,
        edgecolor="black",
        left=data.overlap_size.iloc[start:end].values
        + data.left_only_size.iloc[start:end].values,
        color=right_color,
    )
    axs[i].legend(
        bbox_to_anchor=(1, 1.2 if i != 0 else 1),
    )
    axs[i].set_xlabel(f"Number of {probe_or_gene}s")

for i in range(len(sizes)):
    axs[i].spines["top"].set_visible(False)
    axs[i].spines["right"].set_visible(False)

fig.suptitle(
    f"Number of unique/overlapping {probe_or_gene}s for each TFBS hit", fontsize=12
)

plt.tight_layout()
plt.subplots_adjust(wspace=0.8, hspace=1.2)
fig.savefig(f"{dout}/TFBS_{probe_or_gene}_overlaps.png")
plt.close()

######################################################################
#### Plot Venn diagrams for ER+ vs ER- hyper/hypo-methylated CpGs ####
######################################################################

din = "/projects/p30791/methylation/differential_methylation"
dout = "/projects/p30791/methylation/plots"

probes = {}
for trend in ["hyper", "hypo"]:
    probes[trend] = {}
    for er_status in ["ER+", "ER-"]:
        with open(f"{din}/probe_set_{trend}_{er_status}_refAN_compTU.txt", "r") as f:
            probes[trend][er_status] = [line.strip() for line in f.readlines()]

for trend in ["hyper", "hypo"]:
    er_pos = set(probes[trend]["ER+"])
    er_neg = set(probes[trend]["ER-"])
    total_size = len(er_pos.union(er_neg))
    pos_only = len(er_pos - er_neg)
    overlap = len(er_pos.intersection(er_neg))
    neg_only = len(er_neg - er_pos)
    venn2(
        subsets=(pos_only, neg_only, overlap),
        set_labels=("Samples from\nER+ patients", "Samples from\nER- patients"),
    )
    plt.savefig(f"{dout}/CpG_overlap_{trend}_ER_pos_vs_neg.png")
    plt.close()
