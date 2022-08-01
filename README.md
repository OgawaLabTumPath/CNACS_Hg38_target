# CNACS_hg38_target
CNACS_hg38_target is a program for detecting copy number changes.

## Information
CNACS_hg38_target was built from [original CNACS](https://github.com/OgawaLabTumPath/CNACS) by adding following two types of changes:
1. Liftover the reference version from Hg19 to Hg38
2. Bug fix

## Manual
For more details such as analysis flows, directory settings, and arguments, please refer to [original CNACS manual](https://github.com/OgawaLabTumPath/CNACS/blob/main/info/CNACS_Manual_072919.docx).
In CNACS_hg38_target, following steps are compatible with analysis in the reference version Hg38:  
### Steps compatible with Hg38 
* Step2: Establish a control pool
* Step3: Analysis of tumor samples
### Steps incompatible with Hg38 
* Step1: Design SNP probe
* Step4: Adjustments of results
* Step5: Summarization of results

## Scripts
The information on the scripts are as follows:
* subscript_target
    * The reference version is Hg38.
    * Bugs were fixed.
    * The scripts for drawing plot figures, `plot.sh` and  `plot_all.R` were modified.
* create_pool.sh
    * Establishing a control pool steps originaly performed by `ref_count.sh` and `ref_install.sh` in [original CNACS](https://github.com/OgawaLabTumPath/CNACS) are combined into this script.
    * The analysis mode is target sequence only.
    * Thresholds `(Threshold.txt)` to adjust regions to be analyzed are automatically calculated in this script. If the thresholds need to be changed, please edit `Threshold.txt` and re-run the step of `ref_install.sh`
    * The sex information of control samples should be listed in `sex_info.txt` and put it in the appropriate directory.
* cn_analysis.sh
    * This scripts includes the step to CNA analysis correspodning to `cnacs.sh` as in [original CNACS](https://github.com/OgawaLabTumPath/CNACS).
    * The analysis mode is target sequence only.

## Dependency
Dependency is as the same as in [original CNACS](https://github.com/OgawaLabTumPath/CNACS). Please intall the appropriate libraries.

## Database
CNACS refers to the following database files. Please create and put in the directories accordingly.
* [db](https://github.com/OgawaLabTumPath/CNACS/tree/main/db) in [original CNACS](https://github.com/OgawaLabTumPath/CNACS).
    * The reference version should be liftover to Hg38.
* Publicly available databases
    * The human genome reference (GRCh38 Genome Reference Consortium Human Reference 38)
    * Pseudoautosomal regions
    * Cytoband 
        * db/cytoBand in [original CNACS](https://github.com/OgawaLabTumPath/CNACS) includes the database files need to be edited. Please edit by referring to the database file as in [original CNACS](https://github.com/OgawaLabTumPath/CNACS).
* Replication timing
    * `hg19.interval.bed` in db/cytoBand directory of [original CNACS](https://github.com/OgawaLabTumPath/CNACS) has following format.
    ```
    Chromosome# location1 location2 Replication_timing_degree_loc#1 Replication_timing_degree_loc#2
    ```
    * Please create a text file in above format by using the replication timing database publicly available and make certain the reference version is Hg38.

