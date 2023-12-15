import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control
import upsetplot
from upsetplot import UpSet
import matplotlib.pyplot as plt

######################################
#### Load SeSAMe modeling results ####
######################################

din = '/projects/p30791/methylation/sesame_out'
plot_dir = '/projects/p30791/methylation/plots'

dml_results = {
    'refNormal' : pd.read_csv(f'{din}/DML_results_refNormal.csv'),
    'refTU' : pd.read_csv(f'{din}/DML_results_refTU.csv'),
}

tissue_types = ['TU', 'AN', 'OQ', 'CUB', 'Normal']

dmr_results = {}
for ref in ['refNormal', 'refTU']:
    dmr_results[ref] = {}
    # Add results for age
    dmr_results[ref]['Age'] = pd.read_csv(f'{din}/DMR_results_Age_{ref}.csv')
    # Add results for tissue regions
    for comparison in tissue_types:
        if comparison not in ref:  # cannot compare reference to reference
            dmr_results[ref][comparison] = pd.read_csv(f'{din}/DMR_results_Sample.Region{comparison}_{ref}.csv')

p_thres = 0.01

######################################
#### Analyze results ####
######################################

## Differentially methylated probes

dml_significance = {}
upset_data = {}

for ref in dml_results.keys():
    dml_significance[ref] = pd.DataFrame(index=dml_results[ref].Probe_ID, columns=tissue_types)
    for tissue_type in tissue_types:
        if tissue_type not in ref:
            # Record significance WITH FDR CORRECTION
            dml_significance[ref][tissue_type] = (false_discovery_control(dml_results[ref][f'Pval_Sample.Region{tissue_type}']) < p_thres)
        else:
            # remove reference column here
            dml_significance[ref] = dml_significance[ref].drop(tissue_type, axis=1)
    dml_significance[ref].rename({'TU': 'Tumor', 'AN': 'Adjacent normal', 'OQ': 'Opposite quadrant', 'CUB': 'Contralateral'}, axis=1, inplace=True)
    upset_data[ref] = dml_significance[ref].groupby(list(dml_significance[ref].columns)).size()

# Plot
for ref in dml_results.keys():
    upsetplot.plot(upset_data[ref], sort_categories_by='input')
    comparison = 'Tumor' if ref=='refTU' else 'Normal'
    plt.title(f'Number of probes differentially methylated compared to {comparison}')
    plt.savefig(f'{plot_dir}/upset_{ref}.png')
    plt.close()

## Differentially methylated regions

dmr_significance = {}
upset_data = {}

for ref in dmr_results.keys():
    dmr_significance[ref] = pd.DataFrame(index=np.sort(dmr_results[ref]['AN'].Seg_ID.unique()), columns=tissue_types)
    for tissue_type in tissue_types:
        if tissue_type not in ref:
            segs = dmr_results[key][comparison][['Seg_ID', 'Seg_Chrm', 'Seg_Start', 'Seg_End', 'Seg_Est', 'Seg_Pval', 'Seg_Pval_adj']].drop_duplicates(ignore_index=True).sort_values('Seg_ID')
            dmr_significance[ref][tissue_type] = (segs.Seg_Pval_adj < p_thres).values
        else:
            # remove reference column here
            dmr_significance[ref] = dmr_significance[ref].drop(tissue_type, axis=1)
    dmr_significance[ref].rename({'TU': 'Tumor', 'AN': 'Adjacent normal', 'OQ': 'Opposite quadrant', 'CUB': 'Contralateral'}, axis=1, inplace=True)
    upset_data[ref] = dmr_significance[ref].groupby(list(dmr_significance[ref].columns)).size()

## Observation: regions are either significant in all tissue types or in none of the tissue types


#############################################
#### Get more information about segments ####
#############################################

# Segment length (kb) distribution

# Number of probes in a segment

# Look into trends in larger segments

