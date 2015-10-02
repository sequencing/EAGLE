#!/bin/bash

if [ $# -ne 3 ]; then
  echo "Usage:" 
  echo "      ${0##*/} <laneNum(usually 1-8)> <tileNAME(1/2/3/4/21/22/...)> <posInTile(zero-based)>"
  exit 1
fi


convertGlobalPosToContig() {
  GlobalPos=$1
  echo -n "Global pos $GlobalPos = "
  ChrAndPos=`grep totalBases sample_genome/genome_size.xml | cut -d '"' --output-delimiter=" " -f 4,6 | awk "BEGIN { POS=$GlobalPos } { if (POS<\\$2) { CHR=\\$1; exit } else { POS -= \\$2 } } END { print CHR \":\" POS+1 }"`
  ChrAndPosArr=(${ChrAndPos//:/ })
  Chr=${ChrAndPosArr[0]}
  Pos=${ChrAndPosArr[1]}
  echo "$Chr:$Pos"

  echo -n "101 ref bases at this position: "
  tail -n +2 sample_genome/${Chr}.fa |tr -d "\n"|tail -c +${Pos} |head -c 101
  echo ""
}

abscommand="$(cd "${0%/*}" 2>/dev/null; echo "$PWD"/"${0##*/}")"
abspath=`dirname "$abscommand"`


LANE=$1
TILE=$2
POS_IN_TILE=$3

let POS1=POS_IN_TILE+1

echo "Extracting info for read $LANE:$TILE:$POS_IN_TILE"

if [ $LANE -ne 0 ]; then
  let TILENUM=32*\(LANE-1\)+16*\(TILE/1000%10\)-16+8*\(TILE/100%10\)-8+\(TILE%10\)-1
  HexTileNum=`printf "%02X" $TILENUM`
  echo "Zero-based sequential tile num: $TILENUM, 0x$HexTileNum"

  #tilenum("lane4 tile62")=3*32+6*4+2-1=121=0x79
  let FragmentNum=`hexdump -v -e '1/2 "%02X\n"' fragments.tile |grep -n $HexTileNum -m $POS1 |tail -1|cut -f 1 -d ':'`

else

  let FragmentNum=${POS_IN_TILE}+1

  # Find in which tile and tile position this fragment got allocated
  let BytesToSkip=2*POS_IN_TILE
  FullTileNum=`hexdump -v -e '1/2 "%d\n"' fragments.tile -s $BytesToSkip -n 2`
  FullTileHexNum=`hexdump -v -e '1/2 "%02X\n"' fragments.tile -s $BytesToSkip -n 2`
  let LANE=1+FullTileNum/32
  let TILENUM=1+FullTileNum%32
  TILE=`ls RunFolder/Data/Intensities/BaseCalls/s_1_*.filter | sed 's/.*_0*//' | sed 's/.filter//'|tail -n +$TILENUM |head -1`
echo "LANE=$LANE, TILE=$TILE (TileNum=$TILENUM, FullTileNum=$FullTileNum=0x$FullTileHexNum)"
  POS_IN_TILE=`hexdump -v -e '1/2 "%02X\n"' fragments.tile -n $BytesToSkip | grep -c $FullTileHexNum`
  echo "PosInTile=$POS_IN_TILE"
  echo " => $LANE:$TILE:$POS_IN_TILE"
fi
echo "FragmentNum: $FragmentNum"

let posInFile=2*$FragmentNum-2
let FragmentLength=`hexdump -v -e '1/2 "%d\n"' -s $posInFile fragments.length -n 2`
echo "FragmentLength: $FragmentLength"

let GlobalPos=`hexdump -v -e '1/2 "%d\n"' fragments.pos|head -$FragmentNum| awk '{ if ($1 == 65536) { print "Error: extended 0xffff fragments.pos value not supported in this tool"; exit 1 } ; SUM += $1 } END { print SUM }'`
echo "GlobalPos: $GlobalPos"

echo -n "Left-most read should start at: "
convertGlobalPosToContig $GlobalPos

let Read2Start=GlobalPos+FragmentLength-101
echo -n "Right-most read should start at: "
convertGlobalPosToContig $Read2Start

echo ""
echo "If reads are ForwardReverse+:"
echo -n "Read 1 from BCL: "
${abspath}/getBclBases.sh $LANE $TILE $POS_IN_TILE 1 101
echo -n "Read 2 from BCL: "
${abspath}/getBclBases.sh $LANE $TILE $POS_IN_TILE 209 109

echo ""
echo "If reads are ReverseForward-:"
echo -n "Read 2 from BCL: "
${abspath}/getBclBases.sh $LANE $TILE $POS_IN_TILE 109 209
echo -n "Read 1 from BCL: "
${abspath}/getBclBases.sh $LANE $TILE $POS_IN_TILE 101 1

echo ""
