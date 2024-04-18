import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

data = pd.read_csv("/projects/p30791/methylation/plots/TFBS_hits_probe_overlaps.csv")

height_ratios = (12, 3, 2)
# width_ratios = (12, 6)

viridis_cmap = matplotlib.colormaps["viridis"]
colors = [viridis_cmap(index) for index in np.linspace(0, 1, 3)]

left_color = colors[0]
overlap_color = colors[1]
right_color = colors[2]

fig, axs = plt.subplots(
    nrows=3,
    ncols=1,
    height_ratios=height_ratios,
    # width_ratios = width_ratios,
    figsize=(8, 5),
)

## Hypomethylated in both CFN v. CUB and in AN v. TU
axs[0].barh(
    data.TF.iloc[:12].values,
    data.left_only_size.iloc[:12].values,
    label="CFN vs. CUB only",
    edgecolor="black",
    color=left_color,
)
axs[0].barh(
    data.TF.iloc[:12].values,
    data.overlap_size.iloc[:12].values,
    label="Overlap",
    edgecolor="black",
    left=data.left_only_size.iloc[:12].values,
    color=overlap_color,
)
axs[0].barh(
    data.TF.iloc[:12].values,
    data.right_only_size.iloc[:12].values,
    label="AN vs. TU only",
    edgecolor="black",
    left=data.overlap_size.iloc[:12].values + data.left_only_size.iloc[:12].values,
    color=right_color,
)
axs[0].legend(
    bbox_to_anchor=(1, 1),
)

## Hypermethylated in both OQ v. AN and in AN v. TU
axs[1].barh(
    data.TF.iloc[12:15].values,
    data.left_only_size.iloc[12:15].values,
    label="OQ vs. AN only",
    edgecolor="black",
    color=left_color,
)
axs[1].barh(
    data.TF.iloc[12:15].values,
    data.overlap_size.iloc[12:15].values,
    label="Overlap",
    edgecolor="black",
    color=overlap_color,
    left=data.left_only_size.iloc[12:15].values,
)
axs[1].barh(
    data.TF.iloc[12:15].values,
    data.right_only_size.iloc[12:15].values,
    label="AN vs. TU only",
    edgecolor="black",
    color=right_color,
    left=data.overlap_size.iloc[12:15].values + data.left_only_size.iloc[12:15].values,
)
axs[1].legend(
    bbox_to_anchor=(1, 1),
)

## Mixed ones
axs[2].barh(
    data.TF.iloc[15:].values,
    data.left_only_size.iloc[15:].values,
    label="Left only",
    edgecolor="black",
    color=left_color,
)
axs[2].barh(
    data.TF.iloc[15:].values,
    data.overlap_size.iloc[15:].values,
    label="Overlap",
    edgecolor="black",
    color=overlap_color,
    left=data.left_only_size.iloc[15:].values,
)
axs[2].barh(
    data.TF.iloc[15:].values,
    data.right_only_size.iloc[15:].values,
    label="Right only",
    edgecolor="black",
    color=right_color,
    left=data.overlap_size.iloc[15:].values + data.left_only_size.iloc[15:].values,
)
axs[2].legend(
    bbox_to_anchor=(1, 0.5),
)

for i in range(3):
    axs[i].spines["top"].set_visible(False)
    axs[i].spines["right"].set_visible(False)

plt.tight_layout()
fig.savefig("/projects/p30791/methylation/plots/TFBS_probe_overlaps.png")
plt.close()
