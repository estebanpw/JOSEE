#!/bin/bash


if [ $# != 6 ]; then
   echo " ==== ERROR ... you called this script inappropriately."
   echo ""
   echo "   usage:  $0 <multi_frags> <out.csv> <n_trim_itera> <n_files> <distance> <min_len_to_trim>"
   echo ""
   exit -1
fi

dirNameX=$(readlink -f $1 | xargs dirname)
seqXName=$(basename "$1")
extensionX="${seqXName##*.}"
seqXName="${seqXName%.*}"

#seqXName=`basename $1 .fasta`

BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

outfile=${2}
iteras=${3}
nfiles=${4}
distance=${5}
minlen=${6}


echo "${BINDIR}/JOSEE -multifrags $1 -out $2 -min_trim_itera $3 --debug"
${BINDIR}/JOSEE -multifrags $1 -out $2 -min_trim_itera $iteras -min_len_trimming $minlen --debug

echo "Generating plots"
echo "Rscript --vanilla ${BINDIR}/fragsTrim.R $1 $nfiles $distance"
Rscript --vanilla ${BINDIR}/fragsTrim.R $1_${iteras}_${minlen} $nfiles $distance



