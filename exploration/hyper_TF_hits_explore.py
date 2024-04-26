import pandas as pd

din = "/projects/p30791/methylation/differential_methylation"

tfs = ["JARID2", "TCF21", "EZH2", "SUZ12"]

trend = "hyper"
comparison = "refAN_compTU"
genes = {}

for tf in tfs:
    with open(
        f"{din}/genes_{tf}_hits_{trend}_{comparison}.txt",
        "r",
        encoding="utf-8",
    ) as f:
        genes[tf] = [line.strip() for line in f.readlines()]

print(f"\n\n#### {comparison} ####")
relationships = pd.DataFrame(
    0,
    index=tfs,
    columns=list(set(value for sublist in genes.values() for value in sublist)),
)

for tf, genenames in genes.items():
    for genename in genenames:
        relationships.loc[tf, genename] += 1

relationships = relationships.iloc[:, relationships.sum(axis=0).values > 1]

## Print any connections between TFBS hit-causing genes and gene-level enrichment results
for gene in relationships.columns:
    filename = f"{din}/KYCG/testEnrichment_{trend}_{comparison}_genes.csv"
    with open(
        filename,
        "r",
        encoding="utf-8",
    ) as f:
        lines = f.readlines()
    # print any lines that contains the gene name
    for line in lines:
        if gene in line:
            targeting_tfs = list(
                relationships.iloc[relationships[gene].values == 1, :].index
            )
            grep_match = line.rsplit(",", maxsplit=1)[-1].strip()
            print(f"{gene} (matched {grep_match}): Targeted by TFs {targeting_tfs}")
