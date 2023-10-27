# Methylation in Breast Cancer 

The goal of this project is to analyze the patterns in methylation across breast tumor, tumor-adjacent, opposite quadrant, and opposite breast tissue.

## Data exploration

Documenting the data file exploration since this data was shared by a collaborator.

* Files copied from `/projects/b1042/BartomLab/SeemaMethylation/` (copied to `copied_from_b1042/`)
  * Directory `Clare_Project_003/`
    * `Clare_Project_003_ControlDashboard.csv` -- columns: `Category,Control,BeadType,Sample_ID,Sentrix_Label,Section 1 X,Section 1 Y` // 39,953 rows 
    * `Clare_Project_003_Group_Meth_Profile.txt` -- columns: `Index,TargetID,ProbeID_A,ProbeID_B,Default Group.AVG_Beta,Default Group.Intensity` // 865,919 rows 
    * `Clare_Project_003_Sample_Meth_Profile.txt` -- columns: `Index,TargetID,ProbeID_A,ProbeID_B,NA10859_1.AVG_Beta` // 865,919 rows
    * `Clare_Project_003_Sample_Sheet.csv` -- columns: `Sample_Name,Sample_Well,BCD_Well,Sample_Plate,Sample_Group,Pool_ID,Sentrix_ID,Sentrix_Position` // 96 rows
    * `Clare_Project_003_SamplesTable.txt` -- columns: `% Loci Detected Index   Sample ID   Sample Group    Sentrix Barcode Sample Section  Detected CpG (0.01) Detected Cp    G (0.05) Signal Average GRN  Signal Average RED  Signal P05 GRN  Signal P05 RED  Signal P25 GRN  Signal P25 RED  Si    gnal P50 GRN  Signal P50 RED  Signal P75 GRN  Signal P75 RED  Signal P95 GRN  Signal P95 RED  Sample_Well BCD_Well        Sample_Plate    Pool_ID` // 89 rows
    * `infinium-methylationepic-v-1-0-b5-manifest-file (1).bpm` No column header available? // 867,197 rows
    * Directory `Data/`: has files of unknown formats (`*.analysis`, `*.analysis.config`, `*.analysis.user`)
    * Directory `idat_files/`: 11 directories named with numerical sequences
  * `ReNormalized.MethylData.txt` -- header section (8 rows) followed by a tab-separated data section with 2533 columns (`TargetID  ProbeID_A  ProbeID_B  206949970109_R01C01.AVG_Beta  ...  GENOME_BUILD  CHR    MAPINFO`) and 866k rows
  * Directory `SusanClareSFC03357200/` has file formats including `*.project`, `*.project.bin`, `*.analysis`, and `*.analysis.config`
  * `TakaSehl_meta_ordered.withAdj.sorted.csv` -- Sample metadata? with columns `IDAT,ID,Sample Region,Case/Control,Age,Race,BMI` and 472 rows
* Files copied from `/projects/b1122/gannon/CUB/` (copied to `copied_from_b1122`)
  * `CUB.Rproj`
  * `data/`
  * `datamanifests`
  * `Elizabeth/`
  * `results/`
  * `scripts/`

