#!/bin/bash
#$ -S /bin/bash
#$ -cwd
readonly ORGDIR=$1
readonly PROBE_BED=$2
readonly LENG_NUM=$3

source ${CONFIG}
source ${UTIL}

check_num_args $# 3

readonly SEQBAM=`head -n ${SGE_TASK_ID} ${ORGDIR}/bam_list.txt | tail -n 1`
TMP_ID="${SEQBAM##*/}"
ID=`echo ${TMP_ID} | sed -e "s/\.bam//"`

LINE_NUM=`expr ${LENG_NUM} \\* 2 - 1`
readonly LENG=`head -n ${LINE_NUM} ${ORGDIR}/${ID}/tmp/length_stats.txt | tail -n 1`
echo ${LENG}


# calculate depth from SNP-overlapping fragments
echo "cut -f 4-7 ${ORGDIR}/stats/baf_factor.all.bed | \
${BEDTOOLS_PATH}/intersectBed -a ${ORGDIR}/${ID}/tmp/mapped_leng${LENG_NUM}.bed -b stdin -wa | \
${BEDTOOLS_PATH}/intersectBed -a stdin -b ${PROBE_BED} -wo | \
${PERL_PATH} ${COMMAND_CNACS}/subscript_target/baf_adjusted_depth.pl ${ORGDIR}/${ID}/tmp/scaling_factor.txt \
> ${ORGDIR}/${ID}/tmp/overlapping_depth_leng${LENG_NUM}.txt"
cut -f 4-7 ${ORGDIR}/stats/baf_factor.all.bed | \
${BEDTOOLS_PATH}/intersectBed -a ${ORGDIR}/${ID}/tmp/mapped_leng${LENG_NUM}.bed -b stdin -wa | \
${BEDTOOLS_PATH}/intersectBed -a stdin -b ${PROBE_BED} -wo | \
${PERL_PATH} ${COMMAND_CNACS}/subscript_target/baf_adjusted_depth.pl ${ORGDIR}/${ID}/tmp/scaling_factor.txt \
> ${ORGDIR}/${ID}/tmp/overlapping_depth_leng${LENG_NUM}.txt
check_error $?


# calculate depth from SNP-nonoverlapping fragments
echo "cut -f 4-7 ${ORGDIR}/stats/baf_factor.all.bed | \
${BEDTOOLS_PATH}/intersectBed -v -a ${ORGDIR}/${ID}/tmp/mapped_leng${LENG_NUM}.bed -b stdin -wa | \
${BEDTOOLS_PATH}/intersectBed -a stdin -b ${PROBE_BED} -wo | \
${PERL_PATH} ${COMMAND_CNACS}/subscript_target/actual_depth.pl \
> ${ORGDIR}/${ID}/tmp/nonoverlapping_depth_leng${LENG_NUM}.txt"
cut -f 4-7 ${ORGDIR}/stats/baf_factor.all.bed | \
${BEDTOOLS_PATH}/intersectBed -v -a ${ORGDIR}/${ID}/tmp/mapped_leng${LENG_NUM}.bed -b stdin -wa | \
${BEDTOOLS_PATH}/intersectBed -a stdin -b ${PROBE_BED} -wo | \
${PERL_PATH} ${COMMAND_CNACS}/subscript_target/actual_depth.pl \
> ${ORGDIR}/${ID}/tmp/nonoverlapping_depth_leng${LENG_NUM}.txt
check_error $?


# add depth from SNP-overlapping and SNP-nonoverlapping fragments
echo "${PERL_PATH} ${COMMAND_CNACS}/subscript_target/add_depth.pl ${PROBE_BED} ${ORGDIR}/${ID}/tmp/overlapping_depth_leng${LENG_NUM}.txt ${ORGDIR}/${ID}/tmp/nonoverlapping_depth_leng${LENG_NUM}.txt \
> ${ORGDIR}/${ID}/tmp/actual_depth_leng${LENG_NUM}.txt"
${PERL_PATH} ${COMMAND_CNACS}/subscript_target/add_depth.pl ${PROBE_BED} ${ORGDIR}/${ID}/tmp/overlapping_depth_leng${LENG_NUM}.txt ${ORGDIR}/${ID}/tmp/nonoverlapping_depth_leng${LENG_NUM}.txt \
> ${ORGDIR}/${ID}/tmp/actual_depth_leng${LENG_NUM}.txt
check_error $?

echo "rm ${ORGDIR}/${ID}/tmp/overlapping_depth_leng${LENG_NUM}.txt"
echo "rm ${ORGDIR}/${ID}/tmp/nonoverlapping_depth_leng${LENG_NUM}.txt"
rm ${ORGDIR}/${ID}/tmp/overlapping_depth_leng${LENG_NUM}.txt
rm ${ORGDIR}/${ID}/tmp/nonoverlapping_depth_leng${LENG_NUM}.txt

: <<'#__COMMENT_OUT__'
#__COMMENT_OUT__
