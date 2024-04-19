import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

probe_or_gene = "probe"

data = pd.read_csv(
    f"/projects/p30791/methylation/plots/TFBS_hits_{probe_or_gene}_overlaps.csv"
)

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
fig.savefig(f"/projects/p30791/methylation/plots/TFBS_{probe_or_gene}_overlaps.png")
plt.close()
