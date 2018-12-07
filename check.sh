#!/bin/sh

# ========================================================================
# Check perl modules
# ========================================================================

PERL=`which perl`
if [ "$PERL" = "" ]; then
	echo "Missing perl interpreter"
else
	for MODULE in File::Basename File::Copy File::Path Getopt::Long Time::HiRes JSON::Tiny Proc::ProcessTable File::Touch; do
		$PERL -M$MODULE -e 'exit' 2>/dev/null
		if [ "$?" != 0 ]; then
			echo "Missing perl modul $MODULE"
		fi
	done
fi

# ========================================================================
# Check python modules
# ========================================================================

PYTHON=`which python`
[ "$PYTHON" = "" ] && PYTHON=`which python2.7`
if [ "$PYTHON" = "" ]; then
	echo "Missing python interpreter!"
else
	for MODULE in re subprocess urllib urllib2; do
		$PYTHON -c "import $MODULE" 2>/dev/null
		if [ "$?" != 0 ]; then
			echo "Missing python modul $MODULE!"
		fi
	done
fi

# ========================================================================
# Check paths to nanopipe, last aligner and samtools
# ========================================================================

[ "`which nanopipe_calc.pl` 2>/dev/null; echo $?" != 0 ] && echo "Cannot find nanopipe_calc.pl - maybe missing PATH setting!"
[ "$PERL5LIB" = "" ] && echo "Missing PERL5LIB setting to modules directory!"
CONFIG_PATH=`which nanopipe_calc.pl 2>/dev/null`/../modules/nanopipe2/config.pm
[ -f CONFIG_PATH -a "`grep 'Please enter' $CONFIG_PATH 2>/dev/null`" = "" ] && echo "Parameter PROJDIR not changed in config.pm!"
