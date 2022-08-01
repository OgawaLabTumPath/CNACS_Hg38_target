#! /bin/bash
# set -u : a reference to any variable you haven't previously defined is an error, and causes the program to immediately exit
# set -o pipefail : If any command in a pipeline fails, that return code will be used as the return code of the whole pipeline.
set -uo pipefail

write_usage() {
    echo "Usage: $(basename "$0") [options]"
}

# shellcheck disable=SC1090
source "${UTIL}"

# select CONTROL (default: Intact)
CONTROL=${CONTROL_INTACT}
if [ "${SPECIMEN_CLASSIFICATION:-unknown}" == "FFPE" ]; then
    CONTROL=${CONTROL_FFPE}
    echo "CONTROL is FFPE."
else
    echo "CONTROL is Intact."
fi

# check BED files
PROBE_DIR=$(dirname "${PROBE_BED}")
PROBE_DIR=$(cd "${PROBE_DIR}" && pwd)
PROBE_FILE="${PROBE_BED##*/}"
PROBE_BED=${PROBE_DIR}/${PROBE_FILE}
check_file_exists "${PROBE_BED}"

readonly SNP_BED=${CNACSDIR}/control/${CONTROL}/overlapping_snp.bed
readonly REP_TIME=${CNACSDIR}/control/${CONTROL}/repli_time.txt
check_file_exists "${SNP_BED}"
check_file_exists "${REP_TIME}"

# check format of a BED file
${PERL_PATH} "${COMMAND_CNACS}"/subscript_design/check_bed.pl "${PROBE_BED}" "${AUTOSOME}"

# # check a bait size
# readonly BAIT_SIZE=`${PERL_PATH} ${COMMAND_CNACS}/subscript_target/bait_size.pl ${PROBE_BED} ${MODE_THRES}`
# if [ ${BAIT_SIZE} -lt ${MODE_THRES} ]; then
#     MODE=Targeted-seq
#     SUBSCRIPT=subscript_target
# else
#     MODE=Exome-seq
#     SUBSCRIPT=subscript_exome
# fi
# echo "${MODE} mode will be applied."
readonly SUBSCRIPT=subscript_target

# check required files
readonly BAIT_FA=${CNACSDIR}/control/${CONTROL}/sequence.fa
readonly BAF_INFO=${CNACSDIR}/control/${CONTROL}/stats/baf_stats.txt
readonly BAF_FACTOR_ALL=${CNACSDIR}/control/${CONTROL}/stats/baf_factor.all.bed
readonly BAF_FACTOR=${CNACSDIR}/control/${CONTROL}/stats/baf_factor.bed
readonly ALL_DEPTH=${CNACSDIR}/control/${CONTROL}/stats/all_depth.txt
check_file_exists "${BAIT_FA}"
check_file_exists "${BAF_INFO}"
check_file_exists "${BAF_FACTOR_ALL}"
check_file_exists "${BAF_FACTOR}"
check_file_exists "${ALL_DEPTH}"


# check an input directory
readonly INPUTDIR=${INDIR}/${TAG}

if [ ! -e  "${INPUTDIR}" ]; then
    echo "${INPUTDIR} does not exist"
    write_usage
    exit 1
fi

# make an output directory
readonly OUTPUTDIR=${OUTDIR}/${TAG}
check_mkdir "${OUTPUTDIR}"

# collect information on sample names
find "${INPUTDIR}" -maxdepth 2 -name "*.bam" > "${OUTPUTDIR}"/bam_list.txt
readonly FILECOUNT=$(wc -l "${OUTPUTDIR}"/bam_list.txt | cut -d " " -f 1)
if [ "${FILECOUNT}" = "0" ]; then
    echo "bam file does not exist in ${INPUTDIR}"
    write_usage
    exit 1
fi

### define job names
# PREPROCESSING
readonly PROC_BAM=proc_bam.${TAG}
readonly DIVIDE_BED=divide_bed.${TAG}

# PROCESSING INFORMATION ON BAF
readonly BAM2HETERO=bam2hetero.${TAG}

# COMPENSATION FOR LOW CAPTURE EFFICIENCY DUE TO MINOR ALLELE
readonly PROBE2SCALE=probe2scale.${TAG}
readonly CORRECT_BAF=correct_baf.${TAG}

# GC BIAS CORRECTION
readonly CORRECT_GC=correct_gc.${TAG}
readonly PLOT_GC=plot_gc.${TAG}

# CORRECTION OF FRAGMENT LENGTH BIAS
readonly COUNT_DUP=count_dup.${TAG}
readonly CORRECT_LENGTH=correct_length.${TAG}

# CORRECTION OF BIAS FROM WGA
readonly CORRECT_WGA=correct_wga.${TAG}

# CALCULATION OF MEAN DEPTH FOR EACH GENE
readonly GENE_DEPTH=gene_depth.${TAG}

# MAIN
readonly CNACS_MAIN=cnacs_main.${TAG}
readonly PLOT=plot.${TAG}

# SUMMARIZE
readonly CAT_RESULT=cat_result.${TAG}

process_name=""

##### PREPROCESSING #####

# process BAM files
process_name="process BAM files"
echo "START: ${process_name}"
for j in $(seq 1 "${FILECOUNT}");
do
    export SGE_TASK_ID=$j
    "${COMMAND_CNACS}"/"${SUBSCRIPT}"/proc_bam.sh "${OUTPUTDIR}" "${PROBE_BED}"
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
    "${COMMAND_CNACS}"/subscript_target/divide_bed.sh "${OUTPUTDIR}"
    ret=$?
    if [ ${ret} -ne 0 ]; then
        echo "ERROR: ${process_name}"
        exit 1
    fi
done
##### END OF PREPROCESSING #####



##### PROCESSING INFORMATION ON BAF #####

# SNP typing
process_name="SNP typing"
echo "START: ${process_name}"
for j in $(seq 1 "${FILECOUNT}");
do
    export SGE_TASK_ID=$j
    "${COMMAND_CNACS}"/subscript_target/bam2hetero.sh "${OUTPUTDIR}" "${SNP_BED}" "${BAF_INFO}"
    ret=$?
    if [ ${ret} -ne 0 ]; then
        echo "ERROR: ${process_name}"
        exit 1
    fi
done
##### END OF PROCESSING INFORMATION ON BAF #####



##### COMPENSATION FOR LOW CAPTURE EFFICIENCY DUE TO MINOR ALLELES  #####

# decide scaling factors for SNP-overlapping fragments
process_name="decide scaling factors for SNP-overlapping fragments"
echo "START: ${process_name}"
for j in $(seq 1 "${FILECOUNT}");
do
    export SGE_TASK_ID=$j
    "${COMMAND_CNACS}"/subscript_target/probe2scale.sh "${OUTPUTDIR}" "${BAF_FACTOR}"
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
        "${COMMAND_CNACS}"/subscript_target/correct_baf.sh "${OUTPUTDIR}" "${PROBE_BED}" "${i}" "${BAF_FACTOR_ALL}"
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
        "${COMMAND_CNACS}"/"${SUBSCRIPT}"/correct_gc.sh "${OUTPUTDIR}" "${PROBE_BED}" "${i}" "${BAIT_FA}"
        ret=$?
        if [ ${ret} -ne 0 ]; then
            echo "ERROR: ${process_name}"
            exit 1
        fi
    done
done

process_name="GC PLOT"
echo "START: ${process_name}"
for j in $(seq 1 "${FILECOUNT}");
do
    export SGE_TASK_ID=$j
    "${COMMAND_CNACS}"/subscript_target/plotGC.sh "${OUTPUTDIR}"
    ret=$?
    if [ ${ret} -ne 0 ]; then
        echo "ERROR: ${process_name}"
        exit 1
    fi
done
##### END OF GC BIAS CORRECTION #####



##### CORRECTION OF FRAGMENT LENGTH BIAS #####

# count %duplicate for fragments with binned length
process_name="count %duplicate for fragments with binned length"
echo "START: ${process_name}"
for j in $(seq 1 "${FILECOUNT}");
do
    export SGE_TASK_ID=$j
    "${COMMAND_CNACS}"/subscript_target/count_dup.sh "${OUTPUTDIR}"
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
    "${COMMAND_CNACS}"/subscript_target/correct_length.sh "${OUTPUTDIR}" "${PROBE_BED}"
    ret=$?
    if [ ${ret} -ne 0 ]; then
        echo "ERROR: ${process_name}"
        exit 1
    fi
done
##### END OF FRAGMENT LENGTH BIAS CORRECTION #####



# ##### CALCULATION OF MEAN DEPTH FOR EACH GENE #####
# if [ ${MODE} == "Exome-seq" ]; then
#     process_name="CALCULATION OF MEAN DEPTH FOR EACH GENE"
#     echo "START: ${process_name}"
#     GENE_BED=${PROBE_DIR}/gene2exon.bed
#     check_file_exists ${GENE_BED}
#     ${COMMAND_CNACS}/subscript_exome/gene_depth.sh ${OUTPUTDIR} ${GENE_BED}
# 
#     if [ $? -ne 0 ]; then
#         echo "ERROR: ${process_name}"
#         exit 1
#     fi
# fi
# ##### END OF MEAN DEPTH CALCULATION #####



##### MAIN #####

# filter out low quality probes
# make appropriate signals for control from control samples
# calculate signals for CBS
# perform CBS

process_name="filter out low quality probes"
echo "START: ${process_name}"
for j in $(seq 1 "${FILECOUNT}");
do
    export SGE_TASK_ID=$j
    "${COMMAND_CNACS}"/"${SUBSCRIPT}"/cnacs_main.sh "${OUTPUTDIR}" "${BAF_INFO}" "${BAF_FACTOR}" "${BAF_FACTOR_ALL}" "${ALL_DEPTH}" "${REP_TIME}"
    ret=$?
    if [ ${ret} -ne 0 ]; then
        echo "ERROR: ${process_name}"
        exit 1
    fi
done


# draw a plot
process_name="draw a plot"
echo "START: ${process_name}"
for j in $(seq 1 "${FILECOUNT}");
do
    export SGE_TASK_ID=$j
    # "${COMMAND_CNACS}"/"${SUBSCRIPT}"/plot.sh "${OUTPUTDIR}"
    "${COMMAND_CNACS_CUSTOM}"/"${SUBSCRIPT}"/plot.sh "${OUTPUTDIR}"
    ret=$?
    if [ ${ret} -ne 0 ]; then
        echo "ERROR: ${process_name}"
        exit 1
    fi
done
##### END OF MAIN #####



##### SUMMARIZATION #####
process_name="SUMMARIZATION"
echo "START: ${process_name}"
"${COMMAND_CNACS}"/subscript_target/cat_result.sh "${OUTPUTDIR}"
ret=$?
if [ ${ret} -ne 0 ]; then
    echo "ERROR: ${process_name}"
    exit 1
fi
##### END OF SUMMARIZATION #####


exit 0