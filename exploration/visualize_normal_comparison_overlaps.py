from itertools import chain

import matplotlib.pyplot as plt
import pandas as pd
import upsetplot

dm_din = "/projects/p30791/methylation/differential_methylation"
dv_din = "/projects/p30791/methylation/differential_variability/missMethylDiffVar"
dout = "/projects/p30791/methylation/plots"

setnames = [
    "Normal_vs_TU",
    "CUB_vs_TU",
    "OQ_vs_TU",
    "AN_vs_TU",
]

#########################
#### Load probe sets ####
#########################

dm_probesets = {setname: [] for setname in setnames}
dv_probesets = {setname: [] for setname in setnames}

# Read DM probe sets
for setname in setnames:
    for direction in ["hyper", "hypo"]:
        with open(
            f"{dm_din}/probe_set_{direction}_{setname}.txt", "r", encoding="utf-8"
        ) as f:
            dm_probesets[setname] = dm_probesets[setname] + [
                x.strip() for x in f.readlines()
            ]

# Read DV probe sets
for setname in setnames:
    for direction in ["hyper", "hypo"]:
        with open(
            f"{dv_din}/{direction}DV_missMethyl_{setname}.txt", "r", encoding="utf-8"
        ) as f:
            dv_probesets[setname] = dv_probesets[setname] + [
                x.strip() for x in f.readlines()
            ]

#########################################
#### Create dataframe for Upset plot ####
#########################################

all_probes_dm = set(chain.from_iterable(dm_probesets.values()))
all_probes_dv = set(chain.from_iterable(dv_probesets.values()))

plot_data_dm = pd.DataFrame(False, index=list(all_probes_dm), columns=setnames)
plot_data_dv = pd.DataFrame(False, index=list(all_probes_dv), columns=setnames)

for setname in setnames:
    plot_data_dm.loc[dm_probesets[setname], setname] = True
    plot_data_dv.loc[dv_probesets[setname], setname] = True

####################
#### Plot Upset ####
####################

# DM probes
upsetplot.plot(plot_data_dm.value_counts(), sort_categories_by="input")
plt.title("Overlap of Differentially Methylated Probes")
plt.savefig(f"{dout}/upset_dm_overlap_normals.png")
plt.close()

# DV probes
upsetplot.plot(plot_data_dv.value_counts(), sort_categories_by="input")
plt.title("Overlap of Differentially Variable Probes")
plt.savefig(f"{dout}/upset_dv_overlap_normals.png")
plt.close()
