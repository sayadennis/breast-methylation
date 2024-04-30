import pandas as pd

din = "/projects/p30791/methylation/differential_methylation"

hits = pd.read_csv(f"{din}/TFBS_hits_gene_overlaps.csv")

hits = hits.iloc[:12]  # only focus on hypo overlap for now

trend = "hypo"
left_comparison = "refCFN_compCUB"
right_comparison = "refAN_compTU"

genes = {
    left_comparison: {},
    right_comparison: {},
    "overlap": {},
}

for i in hits.index:
    tf = hits.loc[i, "TF"]
    with open(
        f"{din}/genes_{tf}_hits_{trend}_{left_comparison}_only.txt",
        "r",
        encoding="utf-8",
    ) as f:
        genes[left_comparison][tf] = [line.strip() for line in f.readlines()]
    with open(
        f"{din}/genes_{tf}_hits_{trend}_{right_comparison}_only.txt",
        "r",
        encoding="utf-8",
    ) as f:
        genes[right_comparison][tf] = [line.strip() for line in f.readlines()]
    with open(
        (
            f"{din}/genes_overlap_{tf}_hits_{trend}_{left_comparison}"
            f"_AND_{trend}_{right_comparison}_only.txt"
        ),
        "r",
        encoding="utf-8",
    ) as f:
        genes["overlap"][tf] = [line.strip() for line in f.readlines()]

for comparison in [left_comparison, right_comparison, "overlap"]:
    print(f"\n\n#### {comparison} ####")
    relationships = pd.DataFrame(
        0,
        index=hits.TF.values,
        columns=list(
            set(value for sublist in genes[comparison].values() for value in sublist)
        ),
    )

    for tf, genenames in genes[comparison].items():
        for genename in genenames:
            relationships.loc[tf, genename] += 1

    relationships = relationships.iloc[:, relationships.sum(axis=0).values > 1]

    ## Print any connections between TFBS hit-causing genes and gene-level enrichment results
    for gene in relationships.columns:
        if comparison == "overlap":
            filename_left = (
                f"{din}/KYCG/testEnrichment_{trend}_{left_comparison}_genes.csv"
            )
            filename_right = (
                f"{din}/KYCG/testEnrichment_{trend}_{right_comparison}_genes.csv"
            )
            with open(
                filename_left,
                "r",
                encoding="utf-8",
            ) as f:
                lines_left = f.readlines()
            with open(
                filename_right,
                "r",
                encoding="utf-8",
            ) as f:
                lines_right = f.readlines()
            lines = lines_left + lines_right
        else:
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
