from itertools import chain

import matplotlib.pyplot as plt
import pandas as pd
import upsetplot

dm_din = "/projects/p30791/methylation/differential_methylation"
dout = "/projects/p30791/methylation/plots"

#########################
#### Load probe sets ####
#########################

comparison = "refAN_compTU"

dm_probesets = {
    "hyper": {},
    "hypo": {},
}

# Read DM probe sets
for trend, probe_dict in dm_probesets.items():
    for setname, file_pattern in {
        "Pooled all": trend,
        "Stratified ER+": f"{trend}_ER+",
        "Stratified ER-": f"{trend}_ER-",
    }.items():
        with open(
            f"{dm_din}/probe_set_{file_pattern}_{comparison}.txt", "r", encoding="utf-8"
        ) as f:
            probe_dict[setname] = [x.strip() for x in f.readlines()]

#########################################
#### Create dataframe for Upset plot ####
#########################################

all_probes_dm = {
    trend: set(chain.from_iterable(probe_dict.values()))
    for trend, probe_dict in dm_probesets.items()
}

plot_data_dm = {
    trend: pd.DataFrame(
        False, index=list(all_probes_dm[trend]), columns=probe_dict.keys()
    )
    for trend, probe_dict in dm_probesets.items()
}

for trend, probe_dict in dm_probesets.items():
    for setname, probes in probe_dict.items():
        plot_data_dm[trend].loc[probes, setname] = True

####################
#### Plot Upset ####
####################

upsetplot.plot(plot_data_dm["hyper"].value_counts(), sort_categories_by="input")
plt.title("Hyper")
plt.savefig(f"{dout}/upset_effects_of_separating_er_{comparison}_hyper.png")
plt.close()

upsetplot.plot(plot_data_dm["hypo"].value_counts(), sort_categories_by="input")
plt.title("Hypo")
plt.savefig(f"{dout}/upset_effects_of_separating_er_{comparison}_hypo.png")
plt.close()


## Exploratory ##

# Get the sets that are uniquely discovered in stratified comparisons
for trend in ["hyper", "hypo"]:
    df = plot_data_dm[trend]
    pos_only = df.iloc[
        (df["Pooled all"].values is False)
        & (df["Stratified ER+"].values is True)
        & (df["Stratified ER-"].values is False),
        :,
    ].index
    neg_only = df.iloc[
        (df["Pooled all"].values is False)
        & (df["Stratified ER+"].values is False)
        & (df["Stratified ER-"].values is True),
        :,
    ].index
    # write positive only
    with open(f"/home/srd6051/OQ_AN_{trend}_ERpos_only.txt", "w") as f:
        for item in pos_only:
            f.write(f"{item}\n")
    # write negative only
    with open(f"/home/srd6051/OQ_AN_{trend}_ERneg_only.txt", "w") as f:
        for item in neg_only:
            f.write(f"{item}\n")
