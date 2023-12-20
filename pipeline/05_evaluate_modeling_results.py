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

#########################
#### Analyze results ####
#########################

## Differentially methylated probes

dml_significance = {}
upset_data = {}

for ref in dml_results.keys():
    dml_significance[ref] = pd.DataFrame(index=dml_results[ref].Probe_ID, columns=tissue_types)
    for tissue_type in tissue_types:
        if tissue_type not in ref:
            # Record significance WITH FDR CORRECTION
            p_vals = dml_results[ref][f'Pval_Sample.Region{tissue_type}']
            p_vals.replace(0, 1e-320, inplace=True)  # fix underflow
            p_vals.fillna(0.999999, inplace=True)
            dml_significance[ref][tissue_type] = false_discovery_control(p_vals) < p_thres
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

#############################################
#### Get more information about segments ####
#############################################

seg_cols = ['Seg_ID', 'Seg_Chrm', 'Seg_Start', 'Seg_End', 'Seg_Est', 'Seg_Pval', 'Seg_Pval_adj']
segs = dmr_results['refNormal']['AN'][seg_cols].drop_duplicates(ignore_index=True)

# Segment length (kb) distribution
seglen = (segs.Seg_End - segs.Seg_Start).values.ravel()
for n in [1,10,50,100,250,500,750]:
    print(f'Number of segments larger than {n}kb: {sum(seglen>n*1000)}')

# Number of probes in a segment
n_probes_per_seg = dmr_results['refNormal']['AN'].groupby('Seg_ID').size().values.ravel()
for n in [2,5,10,25,50,75]:
    print(f'Segment with >{n} probes: {sum(n_probes_per_seg>n)}')

# Which chromosomes do larger segments lie in?
size_thres = 50e+3  # 50kb
n_segs_per_chrom = segs.iloc[seglen > size_thres,:].groupby('Seg_Chrm').size()
chroms = [f'chr{i}' for i in range(1,23)] + ['chrX', 'chrY']

fig, ax = plt.subplots(figsize=(7,5))

ax.bar(np.arange(len(chroms)), n_segs_per_chrom.loc[chroms])
ax.set_xticks(np.arange(len(chroms)))
ax.set_xticklabels(chroms, rotation=45, ha='right')
ax.set_xlabel('Chromosome')
ax.set_ylabel(f'Number of segments')
ax.set_title(f'Distribution of large (>{int(size_thres/1e+3)}kb) segments across chromosomes')

plt.tight_layout()
fig.savefig(f'{plot_dir}/num_large_segs_by_chrom.png')
plt.close()

