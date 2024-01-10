import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import umap
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

## Load data
din = "/projects/p30791/methylation/sesame_out"
dout = "/projects/p30791/methylation/sesame_out"
meta_fn = "/projects/p30791/methylation/data/meta.csv"

betas = pd.read_csv(
    f"{din}/betas_processed.csv", index_col=0
)  # nrows=2000 for testing purposes
meta = pd.read_csv(meta_fn)

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

betas = betas[[col for col in betas.columns if col in meta["IDAT"].unique()]]
betas.dropna(axis=0, inplace=True)  # drop columns that contain NA for PCA purposes
betas = betas.iloc[[x.startswith("cg") for x in betas.index], :]

int_labels = np.array([label_int_mapping[idat_to_region[x]] for x in betas.columns])

print("Input data shape:", betas.T.shape)
print()

##########################
#### Run and plot PCA ####
##########################

pca = PCA()
data_tf = pca.fit_transform(betas.T)
data_tf = pd.DataFrame(data_tf, index=betas.columns)

fig, ax = plt.subplots()

for region, intlabel in label_int_mapping.items():
    subset_data = data_tf.loc[int_labels == intlabel, :]
    ax.scatter(
        subset_data[0],
        subset_data[1],
        s=14,
        label=region,
        color=int_color_mapping[intlabel],
    )

ax.set_xlabel(f"PC1 ({100 * pca.explained_variance_ratio_[0]:.1f}%)")
ax.set_ylabel(f"PC2 ({100 * pca.explained_variance_ratio_[1]:.1f}%)")
ax.legend()

fig.suptitle("PCA of betas per sample")

plt.tight_layout()
fig.savefig("/projects/p30791/methylation/plots/pca_by_sample.png")
plt.close()

## Plot cumulative sum of explained variance

num_pcs = 100  # len(pca.explained_variance_ratio_)
cumvar = np.cumsum(pca.explained_variance_ratio_)

fig, ax = plt.subplots()

ax.plot(np.arange(1, num_pcs + 1), cumvar[:num_pcs])
ax.set_xlabel("Number of Principal Components")
ax.set_ylabel("Cumulative Explained Variance")

fig.suptitle("Cumulative Explained variance of PCA")
plt.tight_layout()
fig.savefig("/projects/p30791/methylation/plots/pca_by_sample_cumsum_variance.png")
plt.close()

##############
#### tSNE ####
##############

tsne = TSNE()
data_tf = tsne.fit_transform(betas.T)
data_tf = pd.DataFrame(data_tf, index=betas.columns)

fig, ax = plt.subplots()

for region, intlabel in label_int_mapping.items():
    subset_data = data_tf.loc[int_labels == intlabel, :]
    ax.scatter(
        subset_data[0],
        subset_data[1],
        s=14,
        label=region,
        color=int_color_mapping[intlabel],
    )

ax.set_xlabel("t-SNE Component 1")
ax.set_ylabel("t-SNE Component 2")
ax.legend()

fig.suptitle("tSNE of betas per sample")

plt.tight_layout()
fig.savefig("/projects/p30791/methylation/plots/tsne_by_sample.png")
plt.close()

##############
#### UMAP ####
##############

reducer = umap.UMAP()

data_tf = reducer.fit_transform(betas.T)
data_tf = pd.DataFrame(data_tf, index=betas.columns)

fig, ax = plt.subplots()

for region, intlabel in label_int_mapping.items():
    subset_data = data_tf.loc[int_labels == intlabel, :]
    ax.scatter(
        subset_data[0],
        subset_data[1],
        s=14,
        label=region,
        color=int_color_mapping[intlabel],
    )

ax.set_xlabel("UMAP-1")
ax.set_ylabel("UMAP-2")
ax.legend()

fig.suptitle("UMAP of betas per sample")

plt.tight_layout()
fig.savefig("/projects/p30791/methylation/plots/umap_by_sample.png")
plt.close()
