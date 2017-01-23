#!/bin/bash


if [ $# != 6 ]; then
   echo " ==== ERROR ... you called this script inappropriately."
   echo ""
   echo "   usage:  $0 <paths> <out_dna> <out_class> <blocks_breakpoints_data> <min_len> <mode [0/1]>"
   echo ""
   exit -1
fi


BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

FILEPATHS=$1
OUTDNA=$2
OUTCLASS=$3
DATANAME=$4
LEN=$5
MODE=$6


echo "Cutting out"
${BINDIR}/cutter $FILEPATHS $OUTDNA $OUTCLASS $DATANAME $LEN $MODE

echo "Shuffling"
paste -d ':' $OUTDNA $OUTCLASS | shuf | awk -v FS=":" '{ print $1 > "shuffled.dna" ; print $2 > "shuffled.tclass" }'
mv shuffled.dna $OUTDNA.shuffled
mv shuffled.tclass $OUTCLASS.shuffled






