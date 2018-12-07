#!/bin/sh

#
# ========================================================================
# Nanopipe2 installation script
# ========================================================================
#

NANOPIPE=`pwd`
OS=`uname -s`

if [ "$OS" != "FreeBSD" && "$OS" != "Linux" ]; then
	echo "This script runs only on FreeBSD and Linux"
	exit 0
fi

# ------------------------------------------------------------------------
# Check prerequites
# ------------------------------------------------------------------------

[ ! -d calculate && ! -d modules ] && "Please run inside the extracted nanopipe2 package!" && exit 1
[ `which wget` = "" ] && "Please install wget!" && exit 1
[ `which bunzip2` = "" ] && "Please install bunzip2!" && exit 1
[ `which cc` = "" ] && "Please install cc or gcc (C compiler)!" && exit 1

# ------------------------------------------------------------------------
# Correct PROJDIR in paths.pm
# ------------------------------------------------------------------------

PWD=`pwd`
FILE=$NANOPIPE/modules/nanopipe2/paths.pm
sed -i .tmp -e "s|Please enter your nanopipe2 directory here!|$PWD|" $FILE
rm $FILE.tmp

# ------------------------------------------------------------------------
# Set executable
# ------------------------------------------------------------------------

chmod a+x calculate/*

# ========================================================================
# Install tools
# ========================================================================

# ------------------------------------------------------------------------
# Check for directories
# ------------------------------------------------------------------------

[ ! -d tools ] && mkdir tools
cd tools
[ ! -d bin ] && mkdir bin

# ------------------------------------------------------------------------
# Install last
# ------------------------------------------------------------------------

# Remove
[ -d last ] && rm -rf last

# Get index.html to get filename
wget http://last.cbrc.jp
FILE=`cat index.html | awk '{if(match($0,/last-[0-9\.]+\.zip/))print substr($0,RSTART,RLENGTH)}'`
rm index.html*

# Get last*zip and extract
wget http://last.cbrc.jp/$FILE
unzip -q $FILE; rm $FILE; mv last-* last

if [ $OS = "FreeBSD" ]; then
	sed -i.tmp -e '{
s/= g++/= c++/;
s/= gcc/= cc/
}'  last/makefile
	cd last; gmake; cd ..
fi
if [ $OS = "Linux" ]; then
	cd last; make; cd ..
fi

cp -av last/src/last-split bin
cp -av last/src/lastal bin
cp -av last/src/lastdb bin
cp -av last/scripts/maf-convert bin
cp -av last/scripts/last-train bin

rm -rf last*

# ------------------------------------------------------------------------
# Install samtools
# ------------------------------------------------------------------------

# Remove samtools folder
[ -d samtools ] && rm -rf samtools

# Get index.html to get filename
wget http://www.htslib.org/download/
FILE=`cat index.html | awk '{if(match($0,/samtools-[0-9\.]+\.tar\.bz2/))print substr($0,RSTART,RLENGTH)}'`
rm index.html*

# Get samtools*tar.bz2 (with a trick to get version)
VERSION=`echo $FILE | awk '{if(match($0,/[0-9\.]+/))print substr($0,RSTART,RLENGTH-1)}'`
wget https://github.com/samtools/samtools/releases/download/$VERSION/$FILE
bunzip2 -c $FILE | tar -xf -; rm $FILE; mv samtools-* samtools

if [ $OS = "FreeBSD" ]; then
	sed -i.tmp -e '{
s/= gcc/= cc/
}'  samtools/Makefile
	sed -i.tmp -e '{
s/= gcc/= cc/
}'  samtools/htslib*/Makefile
	cd samtools; gmake; cd ..
fi
if [ $OS = "Linux" ]; then
	cd samtools; make; cd ..
fi

cp samtools/samtools bin

rm -rf samtools

# ------------------------------------------------------------------------
# Settins in .bashrc
# ------------------------------------------------------------------------

echo -n "Add settings in $HOME/.bashrc? [Y/N] "
read ANSWER
if [ "$ANSWER" = "Y" ]; then
	cat <<EOF >> $HOME/.bashrc
# Settings for nanopipe2
export NANOPIPE=$NANOPIPE
export PATH=${PATH}:$NANOPIPE/calculate
export PERL5LIB=$NANOPIPE/modules
EOF
fi

# ------------------------------------------------------------------------
# Add modules (as root)
# ------------------------------------------------------------------------

echo -n "Install perl (system) modules (you must have admin priviledges)? [Y/N] "
read ANSWER
if [ "$ANSWER" = "Y" ]; then
	sudo perl -MCPAN -e "install File::Basename"
	sudo perl -MCPAN -e "install File::Copy"
	sudo perl -MCPAN -e "install File::Path"
	sudo perl -MCPAN -e "install Getopt::Long"
	sudo perl -MCPAN -e "install Time::HiRes"
	sudo perl -MCPAN -e "install JSON::Tiny"
	sudo perl -MCPAN -e "install Proc::ProcessTable"
	sudo perl -MCPAN -e "install File::Touch"
fi
