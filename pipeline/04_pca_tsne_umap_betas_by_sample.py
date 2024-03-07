import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import umap
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

## Load data
din = "/projects/p30791/methylation/sesame_out"
dout = "/projects/p30791/methylation/plots"
meta_fn = "/projects/p30791/methylation/data/meta.csv"

betas = pd.read_csv(
    f"{din}/betas_processed.csv", index_col=0
)  # nrows=2000 for testing purposes
meta = pd.read_csv(meta_fn)

with open(f"{din}/exclude_IDATs.txt", "r", encoding="utf-8") as f:
    exclude_idats = [x.strip() for x in f.readlines()]

meta = meta.iloc[[x not in exclude_idats for x in meta.IDAT]]

# Re-code HER2 status so it's Positive/Negative
meta.loc[:, "HER2"] = meta["HER2"].map(
    {
        0: "Negative",
        1: "Negative",
        2: "Negative",
        3: "Positive",
    }
)

## Create mappings between IDAT IDs, tissue category, and integers
idat_to_region = dict(zip(meta["IDAT"], meta["Sample Region"]))
label_int_mapping = {"Normal": 0, "CUB": 1, "OQ": 2, "AN": 3, "TU": 4}
int_label_mapping = {0: "Normal", 1: "CUB", 2: "OQ", 3: "AN", 4: "TU"}
int_color_mapping = {
    0: "slategray",
    1: "cornflowerblue",
    2: "greenyellow",
    3: "goldenrod",
    4: "orangered",
}

## Clean beta
betas = betas[[col for col in betas.columns if col in meta["IDAT"].unique()]]
betas.dropna(axis=0, inplace=True)  # drop columns that contain NA for PCA purposes
betas = betas.iloc[[x.startswith("cg") for x in betas.index], :]

int_labels = np.array([label_int_mapping[idat_to_region[x]] for x in betas.columns])

# ## Take only the top 5000 most variable probes
# betas = betas.iloc[betas.var(axis=1).argsort().tail(5000),:]

print("Input data shape:", betas.T.shape)

fig, axs = plt.subplots(1, 3, figsize=(12, 4))

##########################
#### Run and plot PCA ####
##########################

pca = PCA()
data_tf = pca.fit_transform(betas.T)
data_tf = pd.DataFrame(data_tf, index=betas.columns)

for region, intlabel in label_int_mapping.items():
    subset_data = data_tf.loc[int_labels == intlabel, :]
    axs[0].scatter(
        subset_data[0],
        subset_data[1],
        s=10,
        label=region,
        color=int_color_mapping[intlabel],
    )

axs[0].set_xlabel(f"PC1 ({100 * pca.explained_variance_ratio_[0]:.1f}%)")
axs[0].set_ylabel(f"PC2 ({100 * pca.explained_variance_ratio_[1]:.1f}%)")
axs[0].legend()
axs[0].set_title("PCA", fontsize=14)

# ## Plot cumulative sum of explained variance
#
# num_pcs = 100  # len(pca.explained_variance_ratio_)
# cumvar = np.cumsum(pca.explained_variance_ratio_)
#
# fig, ax = plt.subplots()
#
# ax.plot(np.arange(1, num_pcs + 1), cumvar[:num_pcs])
# ax.set_xlabel("Number of Principal Components")
# ax.set_ylabel("Cumulative Explained Variance")
#
# fig.suptitle("Cumulative Explained variance of PCA")
# plt.tight_layout()
# fig.savefig(f"{dout}/pca_by_sample_cumsum_variance.png")
# plt.close()

##############
#### tSNE ####
##############

tsne = TSNE()
data_tf = tsne.fit_transform(betas.T)
data_tf = pd.DataFrame(data_tf, index=betas.columns)

for region, intlabel in label_int_mapping.items():
    subset_data = data_tf.loc[int_labels == intlabel, :]
    axs[1].scatter(
        subset_data[0],
        subset_data[1],
        s=10,
        label=region,
        color=int_color_mapping[intlabel],
    )

axs[1].set_xlabel("t-SNE Component 1")
axs[1].set_ylabel("t-SNE Component 2")
# axs[1].legend()

axs[1].set_title("tSNE", fontsize=14)

##############
#### UMAP ####
##############

reducer = umap.UMAP()

data_tf = reducer.fit_transform(betas.T)
data_tf = pd.DataFrame(data_tf, index=betas.columns)

for region, intlabel in label_int_mapping.items():
    subset_data = data_tf.loc[int_labels == intlabel, :]
    axs[2].scatter(
        subset_data[0],
        subset_data[1],
        s=10,
        label=region,
        color=int_color_mapping[intlabel],
    )

axs[2].set_xlabel("UMAP-1")
axs[2].set_ylabel("UMAP-2")
# axs[2].legend()

axs[2].set_title("UMAP", fontsize=14)

fig.suptitle("Sample Clustering via Beta Values", fontsize=16)

plt.tight_layout()
fig.savefig(f"{dout}/clustering_by_sample_betas.png")
plt.close()

################################
#### Plot by tumor metadata ####
################################

# Only select cancer patients
idat_cases = list(meta.iloc[meta["Sample Region"].values != "Normal", :].IDAT)
betas_cases = betas.iloc[:, [x in idat_cases for x in betas.columns]]

## Take only the top 5000 most variable probes
betas_cases = betas_cases.iloc[betas_cases.var(axis=1).argsort().tail(5000), :]

int_labels = np.array(
    [label_int_mapping[idat_to_region[x]] for x in betas_cases.columns]
)

print("Input data shape:", betas_cases.T.shape)

tumormeta_to_plot = {
    "Cancer type": ["IDC", "ILC", "DCIS"],
    "Grade": ["I", "II", "III"],
    "ER": ["+", "-"],
    "PgR": ["+", "-"],
    "HER2": ["Positive", "Negative"],
}

# Run PCA only on cases
pca = PCA(n_components=2)
data_pca = pca.fit_transform(betas_cases.T)
data_pca = pd.DataFrame(data_pca, index=betas_cases.columns)

# Run tSNE only on cases
tsne = TSNE(n_components=2)
data_tsne = tsne.fit_transform(betas_cases.T)
data_tsne = pd.DataFrame(data_tsne, index=betas_cases.columns)

# Run UMAP only on cases
reducer = umap.UMAP(n_components=2)
data_umap = reducer.fit_transform(betas_cases.T)
data_umap = pd.DataFrame(data_umap, index=betas_cases.columns)

for feature_name, categories in tumormeta_to_plot.items():
    fig, axs = plt.subplots(
        nrows=len(categories), ncols=3, figsize=(7.5, 2.5 * len(categories))
    )
    for i, category in enumerate(categories):
        idat_focus = meta.iloc[
            meta[feature_name].values == category, :
        ].IDAT.values.ravel()
        background_samples = {
            "PCA": data_pca.loc[
                [idat_id not in idat_focus for idat_id in data_pca.index], :
            ],
            "tSNE": data_tsne.loc[
                [idat_id not in idat_focus for idat_id in data_tsne.index], :
            ],
            "UMAP": data_umap.loc[
                [idat_id not in idat_focus for idat_id in data_umap.index], :
            ],
        }
        foreground_samples = {
            "PCA": data_pca.loc[
                [idat_id in idat_focus for idat_id in data_pca.index], :
            ],
            "tSNE": data_tsne.loc[
                [idat_id in idat_focus for idat_id in data_tsne.index], :
            ],
            "UMAP": data_umap.loc[
                [idat_id in idat_focus for idat_id in data_umap.index], :
            ],
        }
        for j, method in enumerate(["PCA", "tSNE", "UMAP"]):
            # first plot samples that do not belong to focus caetgory
            axs[i, j].scatter(
                background_samples[method][0],
                background_samples[method][1],
                c="k",
                alpha=0.1,
                s=10,
                label=None,
            )
            # next plot samples that do belong to focus category
            for region, intlabel in label_int_mapping.items():
                if region != "Normal":
                    subset_data = foreground_samples[method].loc[
                        [
                            idat_to_region[idat_id] == region
                            for idat_id in foreground_samples[method].index
                        ],
                        :,
                    ]
                    axs[i, j].scatter(
                        subset_data[0],
                        subset_data[1],
                        s=10,
                        label=region,
                        color=int_color_mapping[intlabel],
                    )
            # label axes etc.
            axs[i, j].set_xticklabels([])
            axs[i, j].set_yticklabels([])
            axs[i, j].spines["top"].set_visible(False)
            axs[i, j].spines["right"].set_visible(False)
            # axs[0].set_xlabel(f"PC1 ({100 * pca.explained_variance_ratio_[0]:.1f}%)")
            # axs[0].set_ylabel(f"PC2 ({100 * pca.explained_variance_ratio_[1]:.1f}%)")
            if j == 0:
                axs[i, j].legend()
                axs[i, j].set_ylabel(f"{feature_name}={category}", fontsize=14)
            if i == 0:
                axs[i, j].set_title(method, fontsize=14)
    fig.suptitle(f"Clustering by beta values separated by {feature_name}", fontsize=16)
    plt.tight_layout()
    feature_name = feature_name.replace(
        " ", "_"
    ).lower()  # "Cancer type" -> "cancer_type"
    fig.savefig(f"{dout}/clustering_by_sample_betas_{feature_name}.png")
    plt.close()
