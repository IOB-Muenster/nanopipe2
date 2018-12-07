#!/bin/sh

#
# ========================================================================
# Script to simplify creation of a target
# ========================================================================
#

if [ $# != 2 ]; then
	cat <<EOF
Usage: $0 targetname fastafile
EOF
	exit 0
fi

PWD=`pwd`

TARGET=$1
FASTA=$2

[ "$NANOPIPE" = "" ] && echo "Please set environment variable NANOPIPE!" && exit 1
[ ! -f "$FASTA" ] && echo "Cannot find fasta file $FASTA!" && exit 1
[ ! -d $NANOPIPE/targets ] && echo "Missing $NANOPIPE/targets directory!" && exit 1
[ -d $NANOPIPE/targets/$TARGET ] && echo "Target $TARGET exists already!" && exit 1

echo "Create target folder..."
mkdir $NANOPIPE/targets/$TARGET

echo "Copy target file..."
cp $FASTA $NANOPIPE/targets/$TARGET/target.fasta

echo "Create target database..."
cd $NANOPIPE/targets/$TARGET
$NANOPIPE/tools/bin/lastdb target target.fasta
cd -

echo "Target $TARGET is created."
