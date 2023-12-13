import os
import sys

import numpy as np
import pandas as pd
from tabulate import tabulate
import matplotlib.pyplot as plt

meta_fn = sys.argv[1] # metadata file path e.g. "/projects/p30791/methylation/data/meta.csv"
out_dir = sys.argv[2] # where to output the data summary from this script e.g. "/projects/p30791/methylation/data_summary"
meta = pd.read_csv(meta_fn, index_col=0)

if not os.path.exists(out_dir):
    os.makedirs(out_dir)
    print(f'Created output directory: {out_dir}')

# Categorical: "Sample Region", "Case/Control", "Race"
# Continuous: "Age", "BMI"

# Summarize mapping of how many samples per patient
n_cases = len(meta.iloc[meta['Case/Control'].values=='case',:]['ID'].unique())
n_controls = len(meta.iloc[meta['Case/Control'].values=='normal',:]['ID'].unique())
print(f'Number of cases: {n_cases} // controls: {n_controls}')

# num_regions_per_patient = meta.groupby(['ID', 'Case/Control', 'Sample Region']).size()
# # TODO: Add logic to calculate percentages of patients that have the correct number of samples per patient (3 for cases, 1 for control)
# print(f'{pct_case:/1f}% of cases have 3 samples per patient, and {pct_control:.1f}% of controls have 1 sample per patient.')

# Distribution of race across case/control category
cc_race = meta[['ID', 'Case/Control', 'Race']].drop_duplicates(ignore_index=True)
cts = pd.crosstab(cc_race['Case/Control'], cc_race['Race']) # counts
cts['Total'] = cts.sum(axis=1)
cts.index.name = None
cts.columns.name = None

cts_pcts = cts.copy() # counts and percentages
for i in cts.index:
    for j in cts.columns[:-1]:
        pct = 100 * cts.loc[i,j] / cts.loc[i,'Total']
        cts_pcts.loc[i,j] = f'{cts.loc[i,j]} ({pct:.1f}%)'

tablestring = tabulate(cts_pcts[['Asian', 'Black or African American', 'White', 'Other', 'Total']], headers='keys', tablefmt='pipe')
print(tablestring)

cts_pcts.to_csv(f'{out_dir}/crosstab_casecontrol_race.csv')

#### Plot continuous variables ####

num_meta = meta[['ID', 'Case/Control', 'Age', 'BMI']].drop_duplicates(ignore_index=False) # numeric values

age_max, age_min = max(meta['Age']), min(meta['Age'])
bmi_max, bmi_min = max(meta['BMI']), min(meta['BMI'])

fig, ax = plt.subplots(3, 2, figsize=(8,6))

## Age
ax[0,0].hist(num_meta['Age'], bins=np.arange(age_min, age_max), color='tab:gray')
ax[1,0].hist(num_meta.iloc[num_meta['Case/Control'].values=='case',:]['Age'], bins=np.arange(age_min, age_max), color='tab:pink')
ax[2,0].hist(num_meta.iloc[num_meta['Case/Control'].values=='normal',:]['Age'], bins=np.arange(age_min, age_max), color='tab:cyan')
ax[2,0].set_xlabel('Age')

ax[0,0].set_title('Age', fontsize=14)

for i, label in enumerate(['All', 'Cases', 'Normals']):
    ax[i,0].set_xlim(age_min-5, age_max+5)
    ax[i,0].set_ylabel(label, fontsize=14)
    ax[i,0].spines['top'].set_visible(False)
    ax[i,0].spines['right'].set_visible(False)

fig.text(0.03, 0.5, "Number of patients", rotation=90, va='center', ha='center', fontsize=14)

## BMI
ax[0,1].hist(num_meta['BMI'], bins=np.arange(np.floor(bmi_min), np.ceil(bmi_max)), color='tab:gray')
ax[1,1].hist(num_meta.iloc[num_meta['Case/Control'].values=='case',:]['BMI'], bins=np.arange(np.floor(bmi_min), np.ceil(bmi_max)), color='tab:pink')
ax[2,1].hist(num_meta.iloc[num_meta['Case/Control'].values=='normal',:]['BMI'], bins=np.arange(np.floor(bmi_min), np.ceil(bmi_max)), color='tab:cyan')
ax[2,1].set_xlabel('BMI')

ax[0,1].set_title('BMI', fontsize=14)

for i, label in enumerate(['All', 'Cases', 'Normals']):
    ax[i,1].set_xlim(np.floor(bmi_min)-5, np.ceil(bmi_max)+5)
    ax[i,1].spines['top'].set_visible(False)
    ax[i,1].spines['right'].set_visible(False)

fig.suptitle('Distributions of Age and BMI', fontsize=16)

plt.tight_layout()
plt.subplots_adjust(left=0.12)
fig.savefig(f'{out_dir}/histograms_age_bmi.png')
plt.close()

