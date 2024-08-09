import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from tabulate import tabulate

meta_fn = sys.argv[1]  # e.g. "/projects/p30791/methylation/raw_data/meta.csv"
summary_dir = sys.argv[2]  # e.g. "/projects/p30791/methylation/data_summary"
sesame_dir = sys.argv[3]  # e.g. "/projects/p30791/methylation/sesame_data"
meta = pd.read_csv(meta_fn, index_col=0)

if not os.path.exists(summary_dir):
    os.makedirs(summary_dir)
    print(f"Created output directory: {summary_dir}")

# Summarize mapping of how many samples per patient
n_cases = len(meta.iloc[meta["Case/Control"].values == "case", :]["ID"].unique())
n_controls = len(meta.iloc[meta["Case/Control"].values == "normal", :]["ID"].unique())
print(f"Initial number of cases: {n_cases} // controls: {n_controls}")

###################################################
#### Alleviate age differences in case/control ####
###################################################

# Identify the 10 youngest women among cancer-free normals
cfn = meta.iloc[meta["Case/Control"].values == "normal", :]
youngest10 = cfn.iloc[np.argsort(cfn["Age"]).values[:10], :].ID.values.ravel()

exclude_idat = meta.iloc[[patient_id in youngest10 for patient_id in meta.ID], :].index

with open(f"{sesame_dir}/exclude_IDATs.txt", "w", encoding="utf-8") as f:
    for item in exclude_idat:
        f.write(f"{item}\n")

# Remove these samples from meta dataframe
meta = meta.iloc[[idat_id not in exclude_idat for idat_id in meta.index], :]

###########################################################
#### Distribution of race across case/control category ####
###########################################################

cc_race = meta[["ID", "Case/Control", "Race"]].drop_duplicates(ignore_index=True)
cts = pd.crosstab(
    cc_race["Case/Control"], cc_race["Race"]
)  # pylint: disable=redefined-outer-name
cts["Total"] = cts.sum(axis=1)
cts.index.name = None
cts.columns.name = None
cts = cts.astype("object")

cts_pcts = cts.copy()  # counts and percentages
for i in cts.index:  # pylint: disable=redefined-outer-name
    for j in cts.columns[:-1]:
        pct = 100 * cts.loc[i, j] / cts.loc[i, "Total"]
        cts_pcts.loc[i, j] = f"{cts.loc[i,j]} ({pct:.1f}%)"

tablestring = tabulate(
    cts_pcts[["Asian", "Black or African American", "White", "Other", "Total"]],
    headers="keys",
    tablefmt="pipe",
)
print(tablestring)

cts_pcts.to_csv(f"{summary_dir}/crosstab_casecontrol_race.csv")

########################################################
#### Tumor characteristics of cases (Plot & Tables) ####
########################################################

fig, axs = plt.subplots(nrows=1, ncols=4, figsize=(6.5, 4))

# custom_colors = list(plt.cm.Set2.colors)  # pylint: disable=no-member
histology_colors = [plt.get_cmap("viridis")(i / 4) for i in range(4)]
grade_colors = [plt.get_cmap("plasma")(i / 4) for i in range(4)]
positive_color = "orchid"
negative_color = "darkturquoise"
missing_color = "slategrey"
width = 0.2

tables = {}

patients = meta[
    [
        "ID",
        "Case/Control",
        "Age",
        "Race",
        "BMI",
        "Cancer type",
        "Grade",
        "ER",
        "PgR",
        "HER2",
    ]
].drop_duplicates(ignore_index=True)
patient_cases = patients.iloc[patients["Case/Control"].values == "case", :]
patient_cases.rename({"Cancer type": "Histology"}, axis=1, inplace=True)

## Metadata obtained by manual review from Priyanka
manual_rewrite = {
    "Histology": {
        "CUB-176": "IDC",
        "CUB-194": "IDC/ILC",
    },
    "HER2": {
        "CUB-051": 1.0,
        "CUB-167": 1.0,
        "CUB-314": 1.0,
        "CUB-128": 3.0,
    },
}

## Re-write ambiguous or missing tumor characteristics
for colname, value_maps in manual_rewrite.items():
    for patient_id, value in value_maps.items():
        ix = patient_cases.iloc[patient_cases.ID.values == patient_id, :].index.item()
        patient_cases.loc[ix, colname] = value


def add_pcts(cts: pd.DataFrame):  # pylint: disable=redefined-outer-name
    """
    Given a count table, return the same
    table but with the entries as
    counts and percentages in parentheses.
    """
    cts = cts.astype("object")
    total = cts.sum().item()
    for i in cts.index:  # pylint: disable=redefined-outer-name
        ct = cts.loc[i, "count"]
        pct = 100 * ct / total  # pylint: disable=redefined-outer-name
        cts.loc[i, "count"] = f"{ct} ({pct:.1f}%)"
    return cts


## Histological type
cts = pd.DataFrame(patient_cases.value_counts("Histology"))
cts = cts.loc[["DCIS", "IDC", "ILC", "IDC/ILC"], :]

axs[0].bar(
    [0],
    cts.loc["DCIS"],
    label="DCIS",
    color=histology_colors[0],
    edgecolor="black",
    linewidth=0.5,
    width=width,
)
axs[0].bar(
    [0],
    cts.loc["IDC"],
    bottom=cts.loc["DCIS"],
    label="IDC",
    color=histology_colors[1],
    edgecolor="black",
    linewidth=0.5,
    width=width,
)
axs[0].bar(
    [0],
    cts.loc["ILC"],
    bottom=cts.loc["DCIS"] + cts.loc["IDC"],
    label="ILC",
    color=histology_colors[2],
    edgecolor="black",
    linewidth=0.5,
    width=width,
)
axs[0].bar(
    [0],
    cts.loc["IDC/ILC"],
    bottom=cts.loc["DCIS"] + cts.loc["IDC"] + cts.loc["ILC"],
    label="IDC/ILC",
    color=histology_colors[3],
    edgecolor="black",
    linewidth=0.5,
    width=width,
)
axs[0].spines["right"].set_visible(False)
axs[0].spines["left"].set_visible(False)
axs[0].set_ylim(0, 72)
axs[0].set_ylabel("Number of patients", fontsize=14)
axs[0].legend(bbox_to_anchor=(0.0, -0.05), loc="upper left")
axs[0].set_xticks([])
axs[0].set_xticklabels([])
axs[0].set_title("Histology", fontsize=14)

tables["Histology"] = add_pcts(cts)

## Histological grade
patient_cases["Grade"].fillna("No histological grade", inplace=True)

cts = pd.DataFrame(patient_cases.value_counts("Grade"))
cts = cts.loc[["No histological grade", "I", "II", "III"], :]

axs[1].bar(
    [0],
    cts.loc["No histological grade"],
    label="No grade",
    color=grade_colors[0],
    edgecolor="black",
    linewidth=0.5,
    width=width,
)
axs[1].bar(
    [0],
    cts.loc["I"],
    bottom=cts.loc["No histological grade"],
    label="Grade I",
    color=grade_colors[1],
    edgecolor="black",
    linewidth=0.5,
    width=width,
)
axs[1].bar(
    [0],
    cts.loc["II"],
    bottom=cts.loc["No histological grade"] + cts.loc["I"],
    label="Grade II",
    color=grade_colors[2],
    edgecolor="black",
    linewidth=0.5,
    width=width,
)
axs[1].bar(
    [0],
    cts.loc["III"],
    bottom=cts.loc["No histological grade"] + cts.loc["I"] + cts.loc["II"],
    label="Grade III",
    color=grade_colors[3],
    edgecolor="black",
    linewidth=0.5,
    width=width,
)
axs[1].spines["right"].set_visible(False)
axs[1].spines["left"].set_visible(False)
axs[1].set_ylim(0, 72)
axs[1].set_yticklabels([])
axs[1].legend(bbox_to_anchor=(0.0, -0.05), loc="upper left")
axs[1].set_xticks([])
axs[1].set_xticklabels([])
axs[1].set_title("Grade", fontsize=14)

tables["Grade"] = add_pcts(cts)

## ER
cts = pd.DataFrame(patient_cases["ER"].value_counts())

axs[2].bar(
    [0],
    cts.loc["+"],
    label="Positive",
    color=positive_color,
    edgecolor="black",
    linewidth=0.5,
    width=width,
)
axs[2].bar(
    [0],
    cts.loc["-"],
    bottom=cts.loc["+"],
    label="Negative",
    color=negative_color,
    edgecolor="black",
    linewidth=0.5,
    width=width,
)
axs[2].spines["right"].set_visible(False)
axs[2].spines["left"].set_visible(False)
axs[2].set_ylim(0, 72)
axs[2].set_yticklabels([])
axs[2].legend(bbox_to_anchor=(0.0, -0.05), loc="upper left")
axs[2].set_xticks([])
axs[2].set_xticklabels([])
axs[2].set_title("ER Status", fontsize=14)

tables["ER"] = add_pcts(cts)

## HER2
patient_cases["HER2"] = patient_cases["HER2"].map(
    {0: "Negative", 1: "Negative", 2: "Equivocal", 3: "Positive"}
)
patient_cases["HER2"].fillna("Missing", inplace=True)
cts = pd.DataFrame(patient_cases["HER2"].value_counts())
cts = cts.loc[["Positive", "Negative", "Missing"]]

axs[3].bar(
    [0],
    cts.loc["Positive"],
    label="Positive",
    color=positive_color,
    edgecolor="black",
    linewidth=0.5,
    width=width,
)
axs[3].bar(
    [0],
    cts.loc["Negative"],
    bottom=cts.loc["Positive"],
    label="Negative",
    color=negative_color,
    edgecolor="black",
    linewidth=0.5,
    width=width,
)
axs[3].bar(
    [0],
    cts.loc["Missing"],
    bottom=cts.loc["Positive"] + cts.loc["Negative"],
    label="Missing",
    color=missing_color,
    edgecolor="black",
    linewidth=0.5,
    width=width,
)
axs[3].spines["right"].set_visible(False)
axs[3].spines["left"].set_visible(False)
axs[3].set_ylim(0, 72)
axs[3].set_yticklabels([])
axs[3].legend(bbox_to_anchor=(0.0, -0.05), loc="upper left")
axs[3].set_xticks([])
axs[3].set_xticklabels([])
axs[3].set_title("HER2 Status", fontsize=14)

tables["HER2"] = add_pcts(cts)

plt.tight_layout()
fig.savefig(f"{summary_dir}/tumor_characteristics.png")
plt.close()

for name, df in tables.items():
    print("\n", df.to_markdown())

###################################
#### Plot continuous variables ####
###################################

age_bmi = meta[["ID", "Case/Control", "Age", "BMI"]].drop_duplicates(ignore_index=False)
age_bmi["Case/Control"] = age_bmi["Case/Control"].map(
    {"case": "Cancer cases", "normal": "Cancer-free normals"}
)

label_colors = [plt.get_cmap("winter")(i / 3) for i in range(3)]

fig, ax = plt.subplots(1, 2, figsize=(4, 4))

sns.violinplot(
    data=age_bmi,
    y="Age",
    hue="Case/Control",
    split=True,
    inner="quartile",
    gap=0.1,
    cut=0,
    ax=ax[0],
    legend=False,
)
ax[0].set_xlabel("")
ax[0].spines["right"].set_visible(False)
ax[0].spines["top"].set_visible(False)
ax[0].spines["bottom"].set_visible(False)
ax[0].set_xticklabels([])
ax[0].yaxis.label.set_size(12)

sns.violinplot(
    data=age_bmi,
    y="BMI",
    hue="Case/Control",
    split=True,
    inner="quartile",
    gap=0.1,
    cut=0,
    ax=ax[1],
)
ax[1].set_xlabel("")
ax[1].spines["right"].set_visible(False)
ax[1].spines["top"].set_visible(False)
ax[1].spines["bottom"].set_visible(False)
ax[1].set_xticklabels([])
ax[1].yaxis.label.set_size(12)

handles, labels = ax[1].get_legend_handles_labels()
ax[1].legend_.remove()  # Remove the existing legend from the second plot

# Create a new legend outside the axes space
legend = fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.57, 0.00))
legend.set_title("")  # Remove the title

plt.tight_layout()
plt.subplots_adjust(bottom=0.15)
fig.savefig(f"{summary_dir}/violin_age_bmi.png")
plt.close()
