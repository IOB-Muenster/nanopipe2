#!/usr/bin/env perl

#
# ========================================================================
# This is the Nanopipe script for analyzing Nanopore data
#
# This script takes a maf fil from last and generates following files:
#
# For Metagenomics:
#   calc.targetscounts
#   calc.querylens
# Else
#   calc.targetscounts
#   calc.querylens
#   calc.nuccounts
#   calc.nuccounts.help
#   calc.consensus
#
# Changes
# [2018-02-14] Integrate "hot spots" / locations
# [2017-04-25] Set /tmp/nanopipe as extraction folder
# ========================================================================
#

use strict;

use Config;
use Cwd;
use Fcntl;
use File::Basename;
use File::Path qw(mkpath);
use Getopt::Long;

$| = 1;

use nanopipe2::paths;
use nanopipe2::config;
use nanopipe2::messages;
use nanopipe2::calculate::query;
use nanopipe2::calculate::analyze;
use nanopipe2::calculate::last;

# ------------------------------------------------------------------------
# The paths to additional programs
# ------------------------------------------------------------------------

my $MAFCONVERT = "$nanopipe2::paths::TOOLSDIR/bin/maf-convert";
my $SAMTOOLS   = "$nanopipe2::paths::TOOLSDIR/bin/samtools";
$ENV{PYTHONPATH} = "$nanopipe2::paths::TOOLSDIR/lib/python2.7/site-packages";

# ------------------------------------------------------------------------
# Generated files and their names
# ------------------------------------------------------------------------

# Input
my $EMAILFILE      = qq(input.email);
my $LASTPARAMSFILE = qq(input.lastparams);
my $QUERYFILE      = qq(input.query);
my $TARGETFILE     = qq(input.target);

# Calculate
my $LASTFILE       = qq(calc.lastalign.maf);
my $PIDFILE        = qq(calc.pid);
my $STATISTICSFILE = qq(calc.statistics);

# ------------------------------------------------------------------------
# Global parameters
# ------------------------------------------------------------------------

# The data directory where all generated files are found (by default
# the current directory)
my $datadir;

# ------------------------------------------------------------------------
# Command line parameters
# ------------------------------------------------------------------------

# The query filename
my $query;

# The internal target or a filename
my $target;

# The email address
my $email;

# The user's last parameters
my $lastparams;

#
# ------------------------------------------------------------------------
# Help
# ------------------------------------------------------------------------
#
if (!@ARGV) {
	my $pgm = basename $0;
	print <<EOS;
Usage: $pgm -q query -t target [-e emailaddr] [-p lastparams]

Where
   query      - the query file (either a single file or archive like *.tgz)
   target     - the target for alaignment (dengue, plasmodium or a filename)
   emailaddr  - the email address (send if job finishes)
   lastparams - parameters for lastal call, surrounded with quotes (like "-A 13 -b 3")
EOS

	exit 1;
}

#
# ------------------------------------------------------------------------
# Send a mail to the user if request is finished
# ------------------------------------------------------------------------
#
sub sendMail {
	$datadir =~ m/\/.(\d+)$/;
	my $id = $1;

	my $title = join("", nanopipe2::utils::readFile(qq($datadir/input.title)));
	my $fulltitle = $id . ($title ? ": $title" : "");

	my ($mailmsg, $mailtitle);
	if (@nanopipe2::messages::messages) {
		$mailtitle = qq(NanoPipe job '$fulltitle' had problems);
		$mailmsg   = qq(Problems with request id $id);
	}
	else {
		$mailtitle = qq(NanoPipe job '$fulltitle' finished);
		my $addr = qq(www.bioinformatics.uni-muenster.de/tools/nanopipe2/generate/index.pl);
		$mailmsg = qq(Click here http://$addr?id=$id to see results);
	}

	my $mail = qq(To: $email\nFrom: bioinfo\@uni-muenster.de\nSubject: $mailtitle\n\n$mailmsg);

	#if (system(qq(echo "$mailmsg" | mail -s "$mailtitle" $email -fbioinfo\@uni-muenster.de)) != 0) {
	if (system(qq(echo "$mail" | ssmtp $email)) != 0) {
		print qq(Cannot send mail for id $id to $email);
	}
	else {
		print qq(Send mail for id $id to $email);
	}
}

#
# ------------------------------------------------------------------------
# Initialize parameters
# ------------------------------------------------------------------------
#
sub init {
	print "==> Initialize\n";

	# Read the command line options
	return 0 if (
		!GetOptions(
			'q=s' => \$query,
			't=s' => \$target,
			'e=s' => \$email,
			'p=s' => \$lastparams
		)
	);

	$datadir = getcwd();
	if (!-w $datadir) {
		print STDERR qq(Cannot write into data directory $datadir!\n);
		return 0;
	}

	$email      = nanopipe2::utils::readFile($EMAILFILE)      if (!$email);
	$lastparams = nanopipe2::utils::readFile($LASTPARAMSFILE) if (!$lastparams);

	if (!$query) {
		print STDERR qq(Missing query file!\n);
		return 0;
	}
	elsif (!-r $query) {
		print STDERR qq(Cannot read query file $query!\n);
		return 0;
	}

	if (!$target) {
		print STDERR qq(Missing target!\n);
		return 0;
	}
	else {
		$target = nanopipe2::utils::readFile($TARGETFILE);
		$target =~ s/^\s+//;
		$target =~ s/\s+$//;
	}

	nanopipe2::config::init($target);

	return 1;
}

#
# ------------------------------------------------------------------------
# Generate bam and bai files
# ------------------------------------------------------------------------
#
sub calcBam {
	print "==> Calculate bam and bai files\n";

	my $start = time;

	if (!-f $LASTFILE) {
		print qq(Missing file $LASTFILE from a previous step!\n);
		exit 1;
	}

	my $command = qq($MAFCONVERT sam -d $LASTFILE >calc.sam
$SAMTOOLS view -bS -o calc.bam calc.sam >/dev/null 2>&1
$SAMTOOLS sort -o calc.bam-sorted calc.bam >/dev/null 2>&1
mv calc.bam-sorted calc.bam
$SAMTOOLS index calc.bam >/dev/null 2>&1);
	my ($res, $error) = nanopipe2::utils::execute($command);
	if ($res > 0 || $error) {
		nanopipe2::utils::printError($command, $res, $error);
	}

	print "Time: " . (time - $start) . " seconds\n";
}

#
# ------------------------------------------------------------------------
# Actions on finnish
# ------------------------------------------------------------------------
#
sub finnish {
	my ($error) = @_;

	nanopipe2::messages::save();
	sendMail() if ($email);
	unlink($PIDFILE) if (-f $PIDFILE);
	sleep(1);

	exit 1 if ($error);
}

#
# ------------------------------------------------------------------------
# Terminate all subprocesses
# ------------------------------------------------------------------------
#
sub termSubprocs {
	for my $signal (15, 9) {
		for my $pid (nanopipe2::utils::getSubpids($$, Proc::ProcessTable->new())) {
			print "Terminate subprocess $pid...\n";
			kill($signal, $pid) if (kill(0, $pid));
			sleep 1;
		}
		sleep 5;
	}
}
#
# ------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------
#
sub main {
	my $start = time;

	return if (!init());

	nanopipe2::utils::writeFile($PIDFILE, $$);

	local $SIG{ALRM} = sub {
		print "==> Terminate " . basename(__FILE__) . "...\n";
		nanopipe2::messages::add(qq(Timeout is reached!));
		termSubprocs();
		finnish(1);
	};
	local $SIG{INT} = sub {
		print "==> Terminate " . basename(__FILE__) . "...\n";
		nanopipe2::messages::add(qq(Request was stopped!));
		termSubprocs();
		finnish(1);
	};
	local $SIG{TERM} = sub {
		print "==> Terminate " . basename(__FILE__) . "...\n";
		nanopipe2::messages::add(qq(Request was stopped!));
		termSubprocs();
		finnish(1);
	};

	alarm $nanopipe2::config::values{common}->{timeout} * 3600;

	eval {
		my $target = nanopipe2::utils::readFile(qq($TARGETFILE));
		$target =~ s/^\s+//;
		$target =~ s/\s+$//;

		nanopipe2::calculate::query::run(
			{
				query  => $query,
				minlen => $nanopipe2::config::values{input}->{minlen}
			}
		);

		nanopipe2::calculate::last::run(
			{
				metagenomics => $nanopipe2::config::values{common}->{metagenomics},
				query        => $QUERYFILE,
				target       => $target,
				threads      => $nanopipe2::config::values{last}->{threads},
				params       => $nanopipe2::config::values{last}->{params},
				substmatrix  => $nanopipe2::config::values{last}->{substmatrix}
			}
		);

		if (!nanopipe2::calculate::last::check()) {
			nanopipe2::messages::add(qq(There are no results, because the basic alignment failed!));
		}
		else {
			nanopipe2::calculate::analyze::run(
				{
					metagenomics => $nanopipe2::config::values{common}->{metagenomics},
					maxmem       => $nanopipe2::config::values{analyze}->{maxmem},
					mincount     => $nanopipe2::config::values{analyze}->{mincount},
					minlen       => $nanopipe2::config::values{analyze}->{minlen},
					maxgap       => $nanopipe2::config::values{analyze}->{maxgap},
					maxN         => $nanopipe2::config::values{analyze}->{maxN},
					equal        => $nanopipe2::config::values{analyze}->{equal},
					hscore       => $nanopipe2::config::values{analyze}->{hscore},
				}
			);
			if ($nanopipe2::config::values{common}->{metagenomics} ne "Y") {
				calcBam();

				# Check if tidmap has size > 0, meaning no data was generated
				my $size = (stat(qq(calc.tidmap)))[7];
				if ($size == 0) {
					nanopipe2::messages::add(qq(No results had been generated!  Maybe the input data was too weak?));
				}
				else {
					my $param_s = $target =~ m/^(hg|Human)/ ? "-s human" : "";
					$param_s = $target =~ m/^plasmodium/ ? "-s plasf" : "" if (!$param_s);
					my $command = qq($nanopipe2::paths::CALCDIR/nanopipe_calc_polymorphism.py -q $param_s);
					my ($res, $error) = nanopipe2::utils::execute($command);
					if ($res > 0 || $error) {
						nanopipe2::utils::printError($command, $res, $error);
					}
				}
			}
		}

		my $duration = time - $start;
		open(F, ">", $STATISTICSFILE);
		print F "Time: $duration seconds\n";
		close(F);

		print "==> Total time: $duration seconds\n";
		finnish(0);
	};
	if ($@) {
		nanopipe2::messages::add(qq(An error occured in the calculation of results!));
		print STDERR $@;
		finnish(1);
	}
}

main();
