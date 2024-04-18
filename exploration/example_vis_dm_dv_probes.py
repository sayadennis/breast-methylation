# pylint: disable=redefined-outer-name

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import false_discovery_control

######################################
#### Load DM/DV results and betas ####
######################################

din = "/projects/p30791/methylation"
plot_dir = "/projects/p30791/methylation/plots/differential_methylation"

ref_comp_pairs_tpx = [  # pairs along the tumor proximity axis
    ("CFN", "CUB"),
    # ("CUB", "OQ"),
    ("OQ", "AN"),
    ("AN", "TU"),
]

color_mapping = {
    "All": "darkslategray",
    "CFN": "slategray",
    "CUB": "deepskyblue",
    "OQ": "forestgreen",
    "AN": "goldenrod",
    "TU": "orangered",
}

num_each = 4  # number of probes to plot per trend per category

dm_results = {}
dv_results = {}

for ref, comp in ref_comp_pairs_tpx:
    dm_results[f"{ref}_vs_{comp}"] = pd.read_csv(
        f"{din}/differential_methylation/DML_results_ref{ref}_comp{comp}.csv"
    )
    dv_results[f"{ref}_vs_{comp}"] = pd.read_csv(
        f"{din}/differential_variability/missMethylDiffVar/topDV_{ref}_vs_{comp}.csv",
        index_col=0,
    )

betas = pd.read_csv(
    f"{din}/sesame_data/betas_processed.csv", index_col=0
)  # nrows=2000 for testing purposes
meta = pd.read_csv(f"{din}/raw_data/meta.csv")

# Align the betas and metadata (remove excluded samples etc.)
betas = betas.iloc[:, [x in meta.IDAT.values for x in betas.columns]]
meta = meta.iloc[[x in betas.columns for x in meta.IDAT.values], :]

meta.set_index("IDAT", drop=False, inplace=True)

for ref, comp in ref_comp_pairs_tpx:
    ###############################
    #### Select probes to plot ####
    ###############################

    dm = dm_results[f"{ref}_vs_{comp}"].sort_values(f"Pval_Sample.Region{comp}")
    pvals = dm[f"Pval_Sample.Region{comp}"].fillna(0.99)
    nonzero_min = min(pvals.iloc[pvals.values != 0])
    pvals.replace(to_replace=0, value=nonzero_min, inplace=True)
    dm = dm.iloc[false_discovery_control(pvals) < 0.01, :]
    hyper_dm = list(
        dm.iloc[dm[f"Est_Sample.Region{comp}"].values > 0.2, :]
        .iloc[:5, :]["Probe_ID"]
        .values
    )
    hypo_dm = list(
        dm.iloc[dm[f"Est_Sample.Region{comp}"].values < -0.2, :]
        .iloc[:5, :]["Probe_ID"]
        .values
    )

    # hyper_dm = ['cg00241081', 'cg00034409', 'cg00028595', 'cg00006600']
    # hypo_dm = ['cg00169897', 'cg00214628', 'cg00173014', 'cg00046575']

    dv = dv_results[f"{ref}_vs_{comp}"]
    dv = dv.iloc[dv["Adj.P.Value"].values < 0.01, :]
    hyper_dv = list(dv.iloc[dv["LogVarRatio"].values > 1, :].iloc[:5, :].index)
    hypo_dv = list(dv.iloc[dv["LogVarRatio"].values < -1, :].iloc[:5, :].index)

    # hyper_dv = ['cg00146085', 'cg00121634', 'cg00139242', 'cg00062873']
    # hypo_dv = ['cg00077264', 'cg00113195', 'cg00242299', 'cg00088464']

    #########################
    #### Plot histograms ####
    #########################

    fig, axs = plt.subplots(nrows=2, ncols=2 * num_each, figsize=(num_each * 4, 8))

    probe_list_dm = hyper_dm[:num_each] + hypo_dm[:num_each]
    probe_list_dv = hyper_dv[:num_each] + hypo_dv[:num_each]

    ref_idats = meta.iloc[meta["Sample Region"].values == ref, :].IDAT.values
    comp_idats = meta.iloc[meta["Sample Region"].values == comp, :].IDAT.values

    for i, probe_list in enumerate([probe_list_dm, probe_list_dv]):
        for j, probe_id in enumerate(probe_list):
            # create a list of counts that fall into each bin here
            bins = np.arange(0, 1.01, step=0.05)
            probe_betas_ref = betas.loc[probe_id, ref_idats].values
            probe_betas_comp = betas.loc[probe_id, comp_idats].values
            ref_cts = np.array(
                [
                    sum((probe_betas_ref >= bins[k]) & (probe_betas_ref < bins[k + 1]))
                    for k in range(len(bins) - 1)
                ]
            )
            comp_cts = np.array(
                [
                    sum(
                        (probe_betas_comp >= bins[k]) & (probe_betas_comp < bins[k + 1])
                    )
                    for k in range(len(bins) - 1)
                ]
            )
            ref_cts = ref_cts / sum(ref_cts)  # normalize
            comp_cts = comp_cts / sum(comp_cts)  # normalize
            # plot
            axs[i, j].barh(
                bins[:-1],
                comp_cts,
                color=color_mapping[comp],
                edgecolor="black",
                height=0.05,
            )
            axs[i, j].barh(
                bins[:-1],
                (-1) * ref_cts,
                color=color_mapping[ref],
                edgecolor="black",
                height=0.05,
            )
            axs[i, j].set_ylim(-0.05, 1.05)
            axs[i, j].set_title(probe_id, fontsize=13)
            axs[i, j].spines["right"].set_visible(False)
            axs[i, j].spines["top"].set_visible(False)
            axs[i, j].set_xlabel("Frequencies")
            if j == 0:
                axs[i, j].set_ylabel("Beta values")

    plt.tight_layout()
    fig.savefig(f"{plot_dir}/dm_dv_example_barh_{ref}_vs_{comp}.png")
    plt.close()

    ##########################
    #### Plot swarm plots ####
    ##########################

    seaborn_data = betas.loc[
        hyper_dm[:num_each]
        + hypo_dm[:num_each]
        + hyper_dv[:num_each]
        + hypo_dv[:num_each],
        meta.iloc[
            [x in [ref, comp] for x in meta["Sample Region"].values], :
        ].IDAT.values,
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

    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(num_each * 4, 8))

    sns.swarmplot(
        data=dm_melted,
        x="cg_column",
        y="Beta values",
        hue="Sample Region",
        dodge=True,
        ax=axs[0],
        palette=color_mapping,
    )
    axs[0].set_title(f"DM Example Probes for {ref} vs {comp}", fontsize=14)
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
        palette=color_mapping,
    )
    axs[1].set_title(f"DV Example Probes for {ref} vs {comp}", fontsize=14)
    axs[1].set_ylim(0, 1.05)
    axs[1].set_xticklabels(axs[1].get_xticklabels(), rotation=30)
    axs[1].set_xlabel("")

    plt.tight_layout()
    fig.savefig(f"{plot_dir}/dm_dv_example_swarm_{ref}_vs_{comp}.png")
    plt.close()

    ###########################
    #### Plot violin plots ####
    ###########################

    seaborn_data = betas.loc[
        hyper_dm[:num_each]
        + hypo_dm[:num_each]
        + hyper_dv[:num_each]
        + hypo_dv[:num_each],
        meta.iloc[
            [x in [ref, comp] for x in meta["Sample Region"].values], :
        ].IDAT.values,
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

    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(num_each * 2.5, 8))

    sns.violinplot(
        data=dm_melted,
        x="cg_column",
        y="Beta values",
        hue="Sample Region",
        palette=color_mapping,
        inner="quart",
        gap=0.1,
        split=False,
        bw_adjust=0.5,  # less smoothing
        ax=axs[0],
    )
    axs[0].set_title(f"DM Examples Probes for {ref} vs {comp}")
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
        palette=color_mapping,
        inner="quart",
        gap=0.1,
        split=False,
        bw_adjust=0.5,  # less smoothing
        ax=axs[1],
    )
    axs[1].set_title(f"DV Example Probes for {ref} vs {comp}")
    axs[1].spines["right"].set_visible(False)
    axs[1].spines["top"].set_visible(False)
    axs[1].set_ylim(0, 1.05)
    axs[1].set_xticklabels(axs[1].get_xticklabels(), rotation=30)
    axs[1].set_xlabel("")

    plt.tight_layout()
    fig.savefig(f"{plot_dir}/dm_dv_example_violin_{ref}_vs_{comp}.png")
    plt.close()
