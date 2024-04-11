# pylint: disable=redefined-outer-name

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import false_discovery_control

######################################
#### Load DM/DV results and betas ####
######################################

din = "/projects/p30791/methylation"
plot_dir = "/projects/p30791/methylation/plots/differential_methylation"

dm = pd.read_csv(f"{din}/differential_methylation/DML_results_refAN_compTU.csv")
dv = pd.read_csv(
    f"{din}/differential_variability/missMethylDiffVar/topDV_AN_vs_TU.csv", index_col=0
)

betas = pd.read_csv(
    f"{din}/sesame_data/betas_processed.csv", index_col=0
)  # nrows=2000 for testing purposes
meta = pd.read_csv(f"{din}/raw_data/meta.csv")
meta.set_index("IDAT", drop=False, inplace=True)

###############################
#### Select probes to plot ####
###############################

dm = dm.sort_values("Pval_Sample.RegionTU")
pvals = dm["Pval_Sample.RegionTU"].fillna(0.99)
nonzero_min = min(pvals.iloc[pvals.values != 0])
pvals.replace(to_replace=0, value=nonzero_min, inplace=True)
dm = dm.iloc[false_discovery_control(pvals) < 0.01, :]
hyper_dm = list(
    dm.iloc[dm["Est_Sample.RegionTU"].values > 0.2, :].iloc[:5, :]["Probe_ID"].values
)
hypo_dm = list(
    dm.iloc[dm["Est_Sample.RegionTU"].values < -0.2, :].iloc[:5, :]["Probe_ID"].values
)

# hyper_dm = ['cg00241081', 'cg00034409', 'cg00028595', 'cg00006600']
# hypo_dm = ['cg00169897', 'cg00214628', 'cg00173014', 'cg00046575']

dv = dv.iloc[dv["Adj.P.Value"].values < 0.01, :]
hyper_dv = list(dv.iloc[dv["LogVarRatio"].values > 1, :].iloc[:5, :].index)
hypo_dv = list(dv.iloc[dv["LogVarRatio"].values < -1, :].iloc[:5, :].index)

# hyper_dv = ['cg00146085', 'cg00121634', 'cg00139242', 'cg00062873']
# hypo_dv = ['cg00077264', 'cg00113195', 'cg00242299', 'cg00088464']

###########################
#### Plot swarm plots ####
###########################

num_each = 4

seaborn_data = betas.loc[
    hyper_dm[:num_each] + hypo_dm[:num_each] + hyper_dv[:num_each] + hypo_dv[:num_each],
    meta.iloc[[x in ["AN", "TU"] for x in meta["Sample Region"].values], :].IDAT.values,
]
seaborn_data = pd.concat(
    (
        seaborn_data,
        pd.DataFrame(
            meta.loc[seaborn_data.columns, "Sample Region"].values.reshape(1, -1),
            index=["Sample Region"],
            columns=seaborn_data.columns,
        ),
    ),
    axis=0,
).T

dm_plot_data = seaborn_data[
    hyper_dm[:num_each] + hypo_dm[:num_each] + ["Sample Region"]
]
dm_melted = dm_plot_data.melt(
    id_vars=["Sample Region"], var_name="cg_column", value_name="Beta values"
)

dv_plot_data = seaborn_data[
    hyper_dv[:num_each] + hypo_dv[:num_each] + ["Sample Region"]
]
dv_melted = dv_plot_data.melt(
    id_vars=["Sample Region"], var_name="cg_column", value_name="Beta values"
)

fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(num_each * 3, 8))

sns.swarmplot(
    data=dm_melted,
    x="cg_column",
    y="Beta values",
    hue="Sample Region",
    dodge=True,
    ax=axs[0],
)
axs[0].set_title("DM examples", fontsize=14)
# axs[0].spines["right"].set_visible(False)
# axs[0].spines["top"].set_visible(False)
axs[0].set_ylim(0, 1.05)
axs[0].set_xticklabels(axs[0].get_xticklabels(), rotation=30)
axs[0].set_xlabel("")

sns.swarmplot(
    data=dv_melted,
    x="cg_column",
    y="Beta values",
    hue="Sample Region",
    dodge=True,
    ax=axs[1],
)
axs[1].set_title("DV examples", fontsize=14)
# axs[1].spines["right"].set_visible(False)
# axs[1].spines["top"].set_visible(False)
axs[1].set_ylim(0, 1.05)
axs[1].set_xticklabels(axs[1].get_xticklabels(), rotation=30)
axs[1].set_xlabel("")

plt.tight_layout()
fig.savefig(f"{plot_dir}/dm_dv_example_swarm.png")
plt.close()

###########################
#### Plot violin plots ####
###########################

num_each = 5

seaborn_data = betas.loc[
    hyper_dm[:num_each] + hypo_dm[:num_each] + hyper_dv[:num_each] + hypo_dv[:num_each],
    meta.iloc[[x in ["AN", "TU"] for x in meta["Sample Region"].values], :].IDAT.values,
]
seaborn_data = pd.concat(
    (
        seaborn_data,
        pd.DataFrame(
            meta.loc[seaborn_data.columns, "Sample Region"].values.reshape(1, -1),
            index=["Sample Region"],
            columns=seaborn_data.columns,
        ),
    ),
    axis=0,
).T

dm_plot_data = seaborn_data[
    hyper_dm[:num_each] + hypo_dm[:num_each] + ["Sample Region"]
]
dm_melted = dm_plot_data.melt(
    id_vars=["Sample Region"], var_name="cg_column", value_name="Beta values"
)

dv_plot_data = seaborn_data[
    hyper_dv[:num_each] + hypo_dv[:num_each] + ["Sample Region"]
]
dv_melted = dv_plot_data.melt(
    id_vars=["Sample Region"], var_name="cg_column", value_name="Beta values"
)

fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(num_each * 1.5, 8))

sns.violinplot(
    data=dm_melted,
    x="cg_column",
    y="Beta values",
    hue="Sample Region",
    inner="quart",
    gap=0.1,
    split=True,
    ax=axs[0],
)
axs[0].set_title("DM examples")
axs[0].spines["right"].set_visible(False)
axs[0].spines["top"].set_visible(False)
axs[0].set_ylim(0, 1.05)
axs[0].set_xticklabels(axs[0].get_xticklabels(), rotation=30)
axs[0].set_xlabel("")

sns.violinplot(
    data=dv_melted,
    x="cg_column",
    y="Beta values",
    hue="Sample Region",
    inner="quart",
    gap=0.1,
    split=True,
    ax=axs[1],
)
axs[1].set_title("DV examples")
axs[1].spines["right"].set_visible(False)
axs[1].spines["top"].set_visible(False)
axs[1].set_ylim(0, 1.05)
axs[1].set_xticklabels(axs[1].get_xticklabels(), rotation=30)
axs[1].set_xlabel("")

plt.tight_layout()
fig.savefig(f"{plot_dir}/dm_dv_example_violin.png")
plt.close()
