# Notes from initial data exploration and project handover

I inherited this project from collaborators who previously worked on this project. Here, I'm recording the results of my initial data exploration (data includes original data files, metadata files, and analysis results of previous members who worked on this project) as well as some pointers I got from Gannon.

## Data exploration

Documenting the data file exploration since this data was shared by a collaborator.

* Files copied from `/projects/b1122/gannon/CUB/` (copied to `/projects/p30791/methylation/copied_from_b1122/`) - mostly Gannon's work
  * `CUB.Rproj`
  * `data/`
    * `archive/` - Gannon let me know that I can mostly disregard these. Files associated with samples that did not pass filtering etc.
    * `IDAT_Files/` - Important data files!
      * 52 directories named with numerical sequences
      * `IDAT_only/` directory that contains only `*.idat` files - this will be the input to SeSaMe
      * `Normal_CUB_Project.RData`
    * `manifests/` - empty directory 
    * `meta/`
      * `EPIC-8v2-0_A1.csv` - 937,699 rows, seems to be information about the probes (ID, sequence, genomic loci, gene, SNP IDs etc.) There are 7 header rows before data section.
      * `infinium-methylationepic-v-1-0-b5-manifest-file.csv` - similar to `EPIC-8v2-0_A1.csv` but with 866,562 rows.
      * `TakaSehl_BMIFiltered_meta.csv` - 388 rows with columns `IDAT,ID,Sample Region,Case/Control,Age,Race,BMI`
      * `Taka_Sehl_comboFiltered_meta.csv` - 409 rows with columns `IDAT,ID,Sample Region,Case/Control,Age,Race,BMI`
      * `TakaSehl_meta_ordered.csv` - 409 rows with columns `IDAT,ID,Sample Region,Case/Control,Age,Race,BMI` - Gannon created this file. Likely does not include adjacent normals, which were processed later.
      * `TakaSehl_meta_ordered.noSlashR.csv` - 409 rows with columns `IDAT,ID,Sample Region,Case/Control,Age,Race,BMI` - Elizabeth created this file. Not sure what `noSlashR` means at this time.
      * `TakaSehl_meta_ordered.withAdj.sorted.csv` - 472 rows with columns `IDAT,ID,Sample Region,Case/Control,Age,Race,BMI` - Elizabeth created this file. Likely includes adjacent normals but double-check.
    * `results/` - results from SeSaMe. For each run, there is a reference and comparison, and the reference is indicated like `refTU` (tumor is the referemce) and the comparison like `CUB` after `merged_` etc.
    * `Sehl_data/` - Sehl is the collaborator that provided the normals (I think?)
  * `datamanifests` - file with 867,197 rows
  * `Elizabeth/TakaSehl_meta_ordered.withAdj.sorted.csv` - file with columns: `IDAT,ID,Sample Region,Case/Control,Age,Race,BMI` and 472 rows. Same as the one from b1042?
  * `results/Region_Comparisons/` - CSV, Rdata, JPEG files etc. - results from SeSaMe. For each run, there is a reference and comparison, and the reference is indicated like `refTU` (tumor is the referemce) and the comparison like `CUB` after `merged_` etc.
  * `scripts/`
    * `CUB_sesame.R`
    * `openSesame.sh`
    * `Sesame_Normal_CUB_script.R`
    * `Sesame_redo_script.R`
    * `sesame_tutorial.R`
* Files copied from `/projects/b1042/BartomLab/SeemaMethylation/` (copied to `/projects/p30791/methylation/copied_from_b1042/`) - mostly Elizabeth's work
  * Directory `Clare_Project_003/`
    * `Clare_Project_003_ControlDashboard.csv` -- columns: `Category,Control,BeadType,Sample_ID,Sentrix_Label,Section 1 X,Section 1 Y` // 39,953 rows 
    * `Clare_Project_003_Group_Meth_Profile.txt` -- columns: `Index,TargetID,ProbeID_A,ProbeID_B,Default Group.AVG_Beta,Default Group.Intensity` // 865,919 rows 
    * `Clare_Project_003_Sample_Meth_Profile.txt` -- columns: `Index,TargetID,ProbeID_A,ProbeID_B,NA10859_1.AVG_Beta` // 865,919 rows
    * `Clare_Project_003_Sample_Sheet.csv` -- columns: `Sample_Name,Sample_Well,BCD_Well,Sample_Plate,Sample_Group,Pool_ID,Sentrix_ID,Sentrix_Position` // 96 rows
    * `Clare_Project_003_SamplesTable.txt` -- columns: `% Loci Detected Index   Sample ID   Sample Group    Sentrix Barcode Sample Section  Detected CpG (0.01) Detected Cp    G (0.05) Signal Average GRN  Signal Average RED  Signal P05 GRN  Signal P05 RED  Signal P25 GRN  Signal P25 RED  Si    gnal P50 GRN  Signal P50 RED  Signal P75 GRN  Signal P75 RED  Signal P95 GRN  Signal P95 RED  Sample_Well BCD_Well        Sample_Plate    Pool_ID` // 89 rows
    * `infinium-methylationepic-v-1-0-b5-manifest-file (1).bpm` No column header available? // 867,197 rows
    * Directory `Data/`: has files of unknown formats (`*.analysis`, `*.analysis.config`, `*.analysis.user`)
    * Directory `idat_files/`: 11 directories named with numerical sequences. These seem to be the adjacent normal ones!!
  * `ReNormalized.MethylData.txt` -- header section (8 rows) followed by a tab-separated data section with 2533 columns (`TargetID  ProbeID_A  ProbeID_B  206949970109_R01C01.AVG_Beta  ...  GENOME_BUILD  CHR    MAPINFO`) and 866k rows
  * Directory `SusanClareSFC03357200/` has file formats including `*.project`, `*.project.bin`, `*.analysis`, and `*.analysis.config`
  * `TakaSehl_meta_ordered.withAdj.sorted.csv` -- Sample metadata? with columns `IDAT,ID,Sample Region,Case/Control,Age,Race,BMI` and 472 rows

