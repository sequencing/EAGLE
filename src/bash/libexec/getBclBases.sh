#!/bin/bash

if [ $# -ne 5 ]; then
  echo "Usage:" 
  echo "      ${0##*/} <laneNum(usually 1-8)> <tileNAME(1/2/3/4/21/22/...)> <posInTile(zero-based)> <firstCycle> <lastCycle>"
  exit 1
fi

LANE=$1
TILE=$2
POS_IN_TILE=$3
CYCLE1=$4
CYCLE2=$5
if [ $CYCLE1 -gt $CYCLE2 ] ; then
  REVERSE_COMPLEMENT=1
  STEP=-1
  BASES=TGCATGCA
else
  REVERSE_COMPLEMENT=0
  STEP=1
  BASES=ACGTACGT
fi

let BYTENUM=POS_IN_TILE+4

echo "Extracting read $LANE:$TILE:$POS_IN_TILE from cycle $CYCLE1 to $CYCLE2: reverse-complement=${REVERSE_COMPLEMENT}"

for i in `seq $CYCLE1 $STEP $CYCLE2` ;
do
  (hexdump -s $BYTENUM -n 1 -b RunFolder/Data/Intensities/BaseCalls/L00${LANE}/C$i.1/s_${LANE}_${TILE}.bcl | head -1 | head -c 11 |tail -c 1 | tr "01234567" "${BASES}");
done
echo ""
for i in `seq $CYCLE1 $STEP $CYCLE2` ;
do
  BCL_BASE=`hexdump -s $BYTENUM -n 1 -b RunFolder/Data/Intensities/BaseCalls/L00${LANE}/C$i.1/s_${LANE}_${TILE}.bcl | head -1 | head -c 11 |tail -c 3`

  let QUAL=16*\(BCL_BASE/100%10\)+2*\(BCL_BASE/10%10\)+\(BCL_BASE%10/4\)
  echo -n "$BCL_BASE($QUAL)."
done

echo ""

