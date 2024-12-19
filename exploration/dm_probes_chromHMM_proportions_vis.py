import matplotlib.pyplot as plt
import pandas as pd

input_fn = "/projects/p30791/methylation/differential_methylation/num_dm_probes_by_chromHMM.csv"
dout = "/projects/p30791/methylation/plots/differential_methylation"

cts = pd.read_csv(input_fn, index_col=0)

totals = cts.sum(axis=0)
pcts = 100 * cts.divide(totals)
pcts_cum = pcts.cumsum(axis=0)

color_mapping = {
    "1_TssA": "red",
    "2_TssAFlnk": "tomato",
    "3_TxFlnk": "limegreen",
    "4_Tx": "forestgreen",
    "5_TxWk": "darkgreen",
    "6_EnhG": "lightgreen",
    "7_Enh": "gold",
    "8_ZNF/Rpts": "mediumaquamarine",
    "9_Het": "steelblue",
    "10_TssBiv": "indianred",
    "11_BivFlnk": "lightcoral",
    "12_EnhBiv": "olivedrab",
    "13_ReprPC": "darkgrey",
    "14_ReprPCWk": "lightgrey",
    "15_Quies": "whitesmoke",
}

fig, ax = plt.subplots(figsize=(10, 3.5))

ax.invert_yaxis()

for chromname in pcts.index:
    ax.barh(
        y=[0, 1, 3, 4],
        width=pcts.loc[chromname].values,
        left=pcts_cum.loc[chromname].values - pcts.loc[chromname].values,
        height=0.5,
        label=chromname,
        color=color_mapping[chromname],
        edgecolor="black",
        linewidth=0.75,
    )

ax.spines["top"].set_visible(False)
ax.spines["bottom"].set_visible(False)

ax.set_xlim([0, 100])
ax.set_yticks([0, 1, 3, 4])
ax.set_yticklabels(
    [
        f"Hyper\n({cts.sum(axis=0).values[0]:,})",
        f"Hypo\n({cts.sum(axis=0).values[1]:,})",
        f"Hyper\n({cts.sum(axis=0).values[2]:,})",
        f"Hypo\n({cts.sum(axis=0).values[3]:,})",
    ]
)
ax.set_xlabel("Percentages (%)")
ax.legend(ncols=1, bbox_to_anchor=(1.01, 0.5), loc="center left", fontsize="small")

# Add row labels with genes
fig.text(0.02, 0.78, "CUB vs HDB", rotation=90, va="center", ha="center", fontsize=12)
fig.text(0.02, 0.30, "TU vs AN", rotation=90, va="center", ha="center", fontsize=12)

plt.tight_layout()
plt.subplots_adjust(left=0.10)
fig.savefig(f"{dout}/dm_chromHMM_proportions_bar.png")
plt.close()
