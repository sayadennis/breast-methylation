"""
Usage: python export_KYCG_manuscript.py /home/srd6051/test_kycg_export.csv $(for x in $(ls /projects/p30791/methylation/differential_methylation/KYCG/testEnrichment_hypo_refCUB_compTU_*.csv); do echo -n "$x "; done) # pylint:disable=line-too-long
"""

import re
import sys
from typing import Dict

import pandas as pd

args = sys.argv[1:]  # "testEnrichment_hyper_refAN_compTU_genes.csv"

output_fn = args[0]
input_fn = args[1:]


def kycg_format(
    df: pd.DataFrame, comparison: str, trend: str, db_category: str
) -> pd.DataFrame:
    """
    Function to re-format the raw KYCG output table for ease of handling
    """
    db_colname = "gene_name" if "gene_name" in df.columns else "dbname"
    df = df[[db_colname, "estimate", "FDR", "nQ", "nD", "overlap"]]  # subset columns
    df = df.rename(
        {
            db_colname: "Hit",
            "estimate": "Fold enrichment",
            "FDR": "Adjusted p-value",
            "overlap": "Overlap",
        },
        axis=1,
    )  # rename columns
    n_hits = df.shape[0]
    if n_hits > 20:
        df = df.iloc[:20, :]
        db_category = f"{db_category} ({n_hits} hits; truncated)"
    df.insert(0, "Comparison", comparison)
    df.insert(1, "Trend", trend)
    df.insert(2, "Database Category", db_category)
    df = df[
        [
            "Comparison",
            "Trend",
            "Database Category",
            "Hit",
            "Fold enrichment",
            "Adjusted p-value",
            "nQ",
            "nD",
            "Overlap",
        ]
    ]
    return df


def get_results_attr_from_filename(kycg_base_fn: str) -> Dict:
    """
    Function to take a KYCG result file name, reformat it, and take out the
    necessary attributes and return it as a dictionary.
    """
    attr_dict = {}
    # First match the most basic pattern
    if re.match(
        r"^testEnrichment_(hyper|hypo)_ref[A-Z]{2,3}_comp[A-Z]{2,3}_.*\.csv$",
        kycg_base_fn,
    ):
        # Get comparison
        ref = kycg_base_fn.split("ref")[1].split("_")[0]
        comp = kycg_base_fn.split("comp")[1].split("_")[0]
        attr_dict["comparison"] = f"{comp} vs {ref}"
        # Get trend
        attr_dict["trend"] = kycg_base_fn.split("_")[1].capitalize()
        # Get db category
        db_name = kycg_base_fn.split("_")[-1].split(".")[0]
        db_name = db_name.capitalize() if db_name == "genes" else db_name
        attr_dict["db_category"] = db_name
    elif re.match(
        r"^testEnrichment_(hyper|hypo)_ER[+-]_ref[A-Z]{2,3}_comp[A-Z]{2,3}_.*\.csv$",
        kycg_base_fn,
    ):
        # Get comparison
        ref = kycg_base_fn.split("ref")[1].split("_")[0]
        comp = kycg_base_fn.split("comp")[1].split("_")[0]
        er_status = "ER+" if "ER+" in kycg_base_fn else "ER-"
        attr_dict["comparison"] = f"{er_status} {comp} vs {ref}"
        # Get trend
        attr_dict["trend"] = kycg_base_fn.split("_")[1].capitalize()
        # Get db category
        db_name = kycg_base_fn.split("_")[-1].split(".")[0]
        db_name = db_name.capitalize() if db_name == "genes" else db_name
        attr_dict["db_category"] = db_name
    elif re.match(
        r"^testEnrichment_(hyper|hypo)_(ER|HER2)_neg_vs_pos_in_[A-Z]{2,3}_.*\.csv$",
        kycg_base_fn,
    ):
        # Get comparison
        tissue = kycg_base_fn.split("_")[7]
        biomarker = "HER2" if "HER2" in kycg_base_fn else "ER"
        attr_dict["comparison"] = f"{biomarker}+ vs {biomarker}- in {tissue}"
        # Get trend
        attr_dict["trend"] = kycg_base_fn.split("_")[1].capitalize()
        # Get db category
        db_name = kycg_base_fn.split("_")[-1].split(".")[0]
        db_name = db_name.capitalize() if db_name == "genes" else db_name
        attr_dict["db_category"] = db_name
    elif re.match(
        r"^relaxed_testEnrichment_(hyper|hypo)_[A-Z]{2,3}_vs_[A-Z]{2,3}_.*\.csv$",
        kycg_base_fn,
    ):
        # Get comparison
        ref = kycg_base_fn.split("_")[3]
        comp = kycg_base_fn.split("_")[5]
        attr_dict["comparison"] = f"Relaxed {comp} vs {ref}"
        # Get trend
        attr_dict["trend"] = kycg_base_fn.split("_")[2].capitalize()
        # Get db category
        db_name = kycg_base_fn.split("_")[-1].split(".")[0]
        db_name = db_name.capitalize() if db_name == "genes" else db_name
        attr_dict["db_category"] = db_name
    elif re.match(
        r"^testEnrichment_(hyper|hypo)DV_missMethyl_[A-Z]{2,3}_vs_[A-Z]{2,3}_.*\.csv$",
        kycg_base_fn,
    ):
        # Get comparison
        ref = kycg_base_fn.split("_")[3]
        comp = kycg_base_fn.split("_")[5]
        attr_dict["comparison"] = f"DV {comp} vs {ref}"
        # Get trend
        attr_dict["trend"] = kycg_base_fn.split("DV")[0].split("_")[-1].capitalize()
        # Get db category
        db_name = kycg_base_fn.split("_")[-1].split(".")[0]
        db_name = db_name.capitalize() if db_name == "genes" else db_name
        attr_dict["db_category"] = db_name
    elif re.match(
        r"^testEnrichment_hypoDV_missMethyl_Non_monotonic_valley_A_.*\.csv$",
        kycg_base_fn,
    ):
        # Get comparison name
        attr_dict["comparison"] = "Non-monotonic valley A"
        # Get trend
        attr_dict["trend"] = "Hypo"
        # Get db category
        db_name = kycg_base_fn.split("_")[-1].split(".")[0]
        db_name = db_name.capitalize() if db_name == "genes" else db_name
        attr_dict["db_category"] = db_name
    # elif re.match(f"^testEnrichment_.*onotonic_.*_[A-F]_.*.csv$", kycg_base_fn):
    #     # Get comparison
    #     setname = re.findall(r"_.*onotonic_.*_[A-F]_", kycg_base_fn)[0]
    #     setname = setname.replace("_", " ").strip()
    #     tissue = kycg_base_fn.split("_")[7]
    #     biomarker = "HER2" if "HER2" in kycg_base_fn else "ER"
    #     attr_dict["comparison"] = f"{biomarker}+ vs {biomarker}- in {tissue}"
    #     # Get trend
    #     attr_dict["trend"] = kycg_base_fn.split("_")[1].capitalize()
    #     # Get db category
    #     db_name = kycg_base_fn.split("_")[-1].split(".")[0]
    #     db_name = db_name.capitalize() if db_name == "genes" else db_name
    #     attr_dict["db_category"] = db_name
    return attr_dict


export_list = []

for fn in input_fn:
    data = pd.read_csv(fn)
    attrs = get_results_attr_from_filename(fn.rsplit("/", maxsplit=1)[-1])
    export = kycg_format(data, **attrs)
    export_list.append(export)

export = pd.concat(export_list, axis=0)
export["Fold enrichment"] = export["Fold enrichment"].round(2)
export["Adjusted p-value"] = export["Adjusted p-value"].apply(lambda x: f"{x:.2e}")


export.to_csv(output_fn, index=False)
