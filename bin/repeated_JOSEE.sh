#!/bin/bash


if [ $# != 2 ]; then
   echo " ==== ERROR ... you called this script inappropriately."
   echo ""
   echo "   usage:  $0 <paths> <configs>" 
   echo ""
   exit -1
fi

GECKODIR=/home/esteban/github/multigecko/gecko/bin

BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

FASTAS=$1
CONFIGS=$2

mkdir TEMP
mkdir RESULTS

MINSIM=( 40 50 60 70 80 )
MINLEN=( 20 100 500 1000 5000 )
MINITERA=5000
KSIZE=16

for ((i=0 ; i < ${CONFIGS} ; i++))
do
	for ((j=0 ; j < ${CONFIGS} ; j++))
	do
		# Run each gecko on it
		${GECKODIR}/gecko_all_vs_all $FASTAS TEMP/frags_${MINLEN[$j]}_${MINSIM[$i]}_${KSIZE}.frags ${MINLEN[$j]} ${MINSIM[$i]} ${KSIZE}
		
	done

done



for ((i=0 ; i < ${CONFIGS} ; i++))
do
	for ((j=0 ; j < ${CONFIGS} ; j++))
	do
		# Run each JOSEE on it
		$BINDIR/JOSEE -multifrags TEMP/frags_${MINLEN[$j]}_${MINSIM[$i]}_${KSIZE}.frags -pathfiles $FASTAS -out RESULTS/ee_${MINLEN[$j]}_${MINSIM[$i]}_${KSIZE}.ee -min_trim_itera $MINITERA
		
	done

done


