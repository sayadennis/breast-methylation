import matplotlib.pyplot as plt
import pandas as pd
from matplotlib_venn import venn2

dm_din = "/projects/p30791/methylation/differential_methylation"
dv_din = "/projects/p30791/methylation/differential_variability/missMethylDiffVar"
dout = "/projects/p30791/methylation/plots"

comparisons = [
    ("Normal", "CUB"),
    # ("CUB", "OQ"),
    ("OQ", "AN"),
    ("AN", "TU"),
]

cts = {}

for ref, comp in comparisons:
    dm_probes = []
    dv_probes = []
    for trend in ["hyper", "hypo"]:
        # get DM probe set
        dm_probes = dm_probes + list(
            pd.read_csv(
                f"{dm_din}/probe_set_{trend}_{ref}_vs_{comp}.txt",
                header=None,
                index_col=None,
            ).values.ravel()
        )
        # get DV probe set
        dv_probes = dv_probes + list(
            pd.read_csv(
                f"{dv_din}/{trend}DV_missMethyl_{ref}_vs_{comp}.txt",
                header=None,
                index_col=None,
            ).values.ravel()
        )
    # Count overlapping and unique probes
    overlap = set(dm_probes).intersection(set(dv_probes))
    cts[f"{ref} vs {comp}"] = {
        "DV-only": len(dv_probes) - len(overlap),
        "DM-only": len(dm_probes) - len(overlap),
        "overlap": len(overlap),
    }

############################
#### Plot Venn diagrams ####
############################

fig, axs = plt.subplots(nrows=1, ncols=len(comparisons), figsize=(12, 4))

for i, (ref, comp) in enumerate(comparisons):
    venn2(
        subsets=(
            cts[f"{ref} vs {comp}"]["DV-only"],  # left-hand unique
            cts[f"{ref} vs {comp}"]["DM-only"],  # right-hand unique
            cts[f"{ref} vs {comp}"]["overlap"],  # overlap
        ),
        set_labels=("DV", "DM"),
        ax=axs[i],
    )
    axs[i].set_title(f"{ref} vs {comp}", fontsize=14)

plt.tight_layout()
fig.savefig(f"{dout}/venn_DV_vs_DM.png")
plt.close()
