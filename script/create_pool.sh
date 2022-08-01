#! /bin/bash
# set -u : a reference to any variable you haven't previously defined is an error, and causes the program to immediately exit
# set -o pipefail : If any command in a pipeline fails, that return code will be used as the return code of the whole pipeline.
set -uo pipefail

write_usage() {
    echo "Usage: $(basename "$0")"
}

PROBE_DIR=$(dirname "${PROBE_BED}")
PROBE_DIR=$(cd "${PROBE_DIR}" && pwd)
PROBE_FILE="${PROBE_BED##*/}"
PROBE_BED=${PROBE_DIR}/${PROBE_FILE}


# import configuration files
# shellcheck disable=SC1090
source "${UTIL}"


# check input files
check_file_exists "${PROBE_BED}"

# check format of a BED file
${PERL_PATH} "${COMMAND_CNACS}"/subscript_design/check_bed.pl "${PROBE_BED}" "${AUTOSOME}"

# # check a bait size
# BAIT_SIZE=$("${PERL_PATH}" "${COMMAND_CNACS}"/subscript_target/bait_size.pl "${PROBE_BED}" "${MODE_THRES}")
# if [ "${BAIT_SIZE}" -lt "${MODE_THRES}" ]; then
#     MODE=Targeted-seq
#     SUBSCRIPT=subscript_target
# else
#     MODE=Exome-seq
#     SUBSCRIPT=subscript_exome
# fi
# echo "${MODE} mode will be applied."
readonly SUBSCRIPT=subscript_target

# check an input directory
readonly ORGDIR="${CNACSDIR}"/control/"${TAG}"

if [ ! -e  "${ORGDIR}" ]; then
    echo "${ORGDIR} does not exist"
    write_usage
    exit 1
fi

# check a file with information on sex of the samples
readonly SEX_INFO=$(find "${ORGDIR}" -name "sex_info.txt")
check_file_exists "${SEX_INFO}"


# collect information on sample names
find "${ORGDIR}" -maxdepth 2 -name "*.bam" > "${ORGDIR}"/bam_list.txt
FILECOUNT=$(wc -l "${ORGDIR}"/bam_list.txt | cut -d " " -f 1)

for i in $(seq 1 "${FILECOUNT}")
do
    SEQBAM=$(head -n "${i}" "${ORGDIR}"/bam_list.txt | tail -n 1)
    TMP_ID="${SEQBAM##*/}"
    ID="${TMP_ID//\.bam/}"
    ID2="${ID//s_/}"

    SEX_STATUS_ALL=$("${PERL_PATH}" "${COMMAND_CNACS}"/subscript_target/check_sex.pl "${SEX_INFO}" "${ID2}")
    SEX_STATUS=$(echo "${SEX_STATUS_ALL}" | cut -d ":" -f 1)
    if [ ! "${SEX_STATUS}" = "Done" ]; then
        echo "${SEX_STATUS_ALL}"
        exit 1
    fi
done


### define job names
# PREPROCESSING
readonly PREPROCESS=preprocess.ref.${TAG}
readonly PROC_BAM=proc_bam.ref.${TAG}
readonly DIVIDE_BED=divide_bed.ref.${TAG}
readonly BAIT_GC_WGA=bait_gc_wga.ref.${TAG}

# PROCESSING INFORMATION ON BAF
readonly BAM2HETERO=bam2hetero.ref.${TAG}
readonly CAT_BAF_INFO=cat_baf_info.ref.${TAG}

# COMPENSATION FOR LOW CAPTURE EFFICIENCY DUE TO MINOR ALLELE
readonly PROBE2SCALE=probe2scale.ref.${TAG}
readonly CORRECT_BAF=correct_baf.ref.${TAG}

# GC BIAS CORRECTION
readonly CORRECT_GC=correct_gc.ref.${TAG}
readonly PLOT_GC=plot_gc.ref.${TAG}

# CORRECTION OF FRAGMENT LENGTH BIAS
readonly COUNT_DUP=count_dup.ref.${TAG}
readonly CORRECT_LENGTH=correct_length.ref.${TAG}

# CORRECTION OF BIAS FROM WGA
readonly CORRECT_WGA=correct_wga.ref.${TAG}

# CALCULATION OF MEAN DEPTH FOR EACH GENE
readonly GENE_DEPTH=gene_depth.ref.${TAG}

# PROCESSING INFORMATION ON DEPTH
readonly CAT_DEPTH=cat_depth.ref.${TAG}

process_name=""

##### PREPROCESSING #####

# select SNPs in the target regions
# obtain sequences of target regions (including flanking regions)
# obtain information on replication timing

process_name="select SNPs in the target regions"
echo "START: ${process_name}"
"${COMMAND_CNACS}"/subscript_target/preprocess.sh "${ORGDIR}" "${PROBE_BED}"
ret=$?
if [ ${ret} -ne 0 ]; then
    echo "ERROR: ${process_name}"
    exit 1
fi


# process BAM files
process_name="process BAM files"
echo "START: ${process_name}"
for j in $(seq 1 "${FILECOUNT}");
do
    export SGE_TASK_ID=$j
    "${COMMAND_CNACS}"/"${SUBSCRIPT}"/proc_bam.sh "${ORGDIR}" "${PROBE_BED}"
    ret=$?
    if [ ${ret} -ne 0 ]; then
        echo "ERROR: ${process_name}"
        exit 1
    fi
done


# divide BED files according to fragments' length
process_name="divide BED files according to fragments' length"
echo "START: ${process_name}"
for j in $(seq 1 "${FILECOUNT}");
do
    export SGE_TASK_ID=$j
    "${COMMAND_CNACS}"/subscript_target/divide_bed.sh "${ORGDIR}"
    ret=$?
    if [ ${ret} -ne 0 ]; then
        echo "ERROR: ${process_name}"
        exit 1
    fi
done


# # calculate GC content for correction of bias from whole-genome amplification
# if [ "${flg_repliG}" = "TRUE" ]; then
#     "${COMMAND_CNACS}"/subscript_target/bait_gc_wga.sh "${ORGDIR}" "${PROBE_BED}"
# fi

##### END OF PREPROCESSING #####



##### PROCESSING INFORMATION ON BAF #####

# SNP typing
process_name="SNP typing"
echo "START: ${process_name}"
for j in $(seq 1 "${FILECOUNT}");
do
    export SGE_TASK_ID=$j
    "${COMMAND_CNACS}"/subscript_target/bam2hetero_ref.sh "${ORGDIR}"
    ret=$?
    if [ ${ret} -ne 0 ]; then
        echo "ERROR: ${process_name}"
        exit 1
    fi
done


# summarize information on BAF
process_name="summarize information on BAF"
echo "START: ${process_name}"
for j in $(seq 1 "${FILECOUNT}");
do
    export SGE_TASK_ID=$j
    "${COMMAND_CNACS}"/subscript_target/cat_baf_info.sh "${ORGDIR}" "${PROBE_BED}"
    ret=$?
    if [ ${ret} -ne 0 ]; then
        echo "ERROR: ${process_name}"
        exit 1
    fi
done

##### END OF PROCESSING INFORMATION ON BAF #####



##### COMPENSATION FOR LOW CAPTURE EFFICIENCY DUE TO MINOR ALLELES  #####

readonly BAF_FACTOR=${ORGDIR}/stats/baf_factor.bed

# decide scaling factors for SNP-overlapping fragments
process_name="decide scaling factors for SNP-overlapping fragments"
echo "START: ${process_name}"
for j in $(seq 1 "${FILECOUNT}");
do
    export SGE_TASK_ID=$j
    "${COMMAND_CNACS}"/subscript_target/probe2scale.sh "${ORGDIR}" "${BAF_FACTOR}"
    ret=$?
    if [ ${ret} -ne 0 ]; then
        echo "ERROR: ${process_name}"
        exit 1
    fi
done


# compensate for low capture efficiency of fragments containing minor alleles
process_name="compensate for low capture efficiency of fragments containing minor alleles"
echo "START: ${process_name}"
for i in $(seq 1 4)
do
    for j in $(seq 1 "${FILECOUNT}");
    do
        export SGE_TASK_ID=$j
        "${COMMAND_CNACS}"/subscript_target/correct_baf_ref.sh "${ORGDIR}" "${PROBE_BED}" "${i}"
        ret=$?
        if [ ${ret} -ne 0 ]; then
            echo "ERROR: ${process_name}"
            exit 1
        fi
    done
done

##### END OF COMPENSATION FOR LOW CAPTURE EFFICIENCY DUE TO MINOR ALLELES  #####



##### GC BIAS CORRECTION #####

process_name="GC BIAS CORRECTION"
echo "START: ${process_name}"
for i in $(seq 1 4)
do
    for j in $(seq 1 "${FILECOUNT}");
    do
        export SGE_TASK_ID=$j
        "${COMMAND_CNACS}"/"${SUBSCRIPT}"/correct_gc_ref.sh "${ORGDIR}" "${PROBE_BED}" "${SEX_INFO}" "${i}"
        ret=$?
        if [ ${ret} -ne 0 ]; then
            echo "ERROR: ${process_name}"
            exit 1
        fi
    done
done

FILECOUNT
echo "${COMMAND_CNACS}"/subscript_target/plotGC.sh "${ORGDIR}"
"${COMMAND_CNACS}"/subscript_target/plotGC.sh "${ORGDIR}"

##### END OF GC BIAS CORRECTION #####



##### CORRECTION OF FRAGMENT LENGTH BIAS #####

# count %duplicate for fragments with binned length
process_name="count %duplicate for fragments with binned length"
echo "START: ${process_name}"
for j in $(seq 1 "${FILECOUNT}");
do
    export SGE_TASK_ID=$j
    "${COMMAND_CNACS}"/subscript_target/count_dup.sh "${ORGDIR}"
    ret=$?
    if [ ${ret} -ne 0 ]; then
        echo "ERROR: ${process_name}"
        exit 1
    fi
done


# correct fragment length bias
process_name="correct fragment length bias"
echo "START: ${process_name}"
for j in $(seq 1 "${FILECOUNT}");
do
    export SGE_TASK_ID=$j
    "${COMMAND_CNACS}"/subscript_target/correct_length_ref.sh "${ORGDIR}" "${PROBE_BED}"
    ret=$?
    if [ ${ret} -ne 0 ]; then
        echo "ERROR: ${process_name}"
        exit 1
    fi
done

##### END OF FRAGMENT LENGTH BIAS CORRECTION #####



##### CORRECTION OF BIAS FROM WGA #####

# if [ ${flg_repliG} = "TRUE" ]; then
#     process_name="CORRECTION OF BIAS FROM WGA"
#     echo "START: ${process_name}"
#     for j in $(seq 1 "${FILECOUNT}");
#     do
#         export SGE_TASK_ID=$j
#         "${COMMAND_CNACS}"/subscript_target/correct_wga_ref.sh "${ORGDIR}" "${ORGDIR}"/bait_gc.txt
#         ret=$?
#         if [ ${ret} -ne 0 ]; then
#             echo "ERROR: ${process_name}"
#             exit 1
#         fi
#     done
# fi

##### END OF CORRECTION OF BIAS FROM WGA #####



##### CALCULATION OF MEAN DEPTH FOR EACH GENE #####

# if [ ${MODE} == "Exome-seq" ]; then
#     readonly GENE_BED=${PROBE_DIR}/gene2exon.bed
#     check_file_exists ${GENE_BED}
# 
#     process_name=" CALCULATION OF MEAN DEPTH FOR EACH GENE"
#     echo "START: ${process_name}"
#     for j in $(seq 1 "${FILECOUNT}");
#     do
#         export SGE_TASK_ID=$j
#         "${COMMAND_CNACS}"/subscript_exome/gene_depth.sh "${ORGDIR}" "${GENE_BED}"
#         ret=$?
#         if [ ${ret} -ne 0 ]; then
#             echo "ERROR: ${process_name}"
#             exit 1
#         fi
#     done
# fi

##### END OF MEAN DEPTH CALCULATION #####



##### PROCESSING INFORMATION ON DEPTH #####

process_name="PROCESSING INFORMATION ON DEPTH"
echo "START: ${process_name}"
"${COMMAND_CNACS}"/"${SUBSCRIPT}"/cat_depth.sh "${ORGDIR}" "${PROBE_BED}"
ret=$?
if [ ${ret} -ne 0 ]; then
    echo "ERROR: ${process_name}"
    exit 1
fi


##### END OF PROCESSING INFORMATION ON DEPTH #####


##### PROCESSING FINALISE THRESHOLDS FOR POOL OF NORMALS #####

# define a file specifying thresholds
readonly THRESHOLD=${ORGDIR}/stats/threshold.txt

# # override copy a file specifying thresholds 
# echo "cp ${THRESHOLD_FILE_PATH} ${THRESHOLD}"
# cp "${THRESHOLD_FILE_PATH}" "${THRESHOLD}"
# ret=$?
# if [ ${ret} -ne 0 ]; then
#     echo "ERROR: override copy a file specifying thresholds "
#     exit 1
# fi

# check a file specifying thresholds
check_file_exists "${THRESHOLD}"

# a directory for plots of bait distribution
check_mkdir "${ORGDIR}"/stats/bait_dist

#####  Install information on probes
process_name="Install information on probes"
echo "START: ${process_name}"
"${COMMAND_CNACS}"/"${SUBSCRIPT}"/ref_install_main.sh "${ORGDIR}" "${PROBE_BED}" "${THRESHOLD}"
ret=$?
if [ ${ret} -ne 0 ]; then
    echo "ERROR: ${process_name}"
    exit 1
fi

##### END OF PROCESSING FINALISE THRESHOLDS FOR POOL OF NORMALS #####
