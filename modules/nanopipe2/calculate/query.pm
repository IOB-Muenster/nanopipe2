package nanopipe2::calculate::query;

#
# ========================================================================
# Prepare query data
# ========================================================================
#

use strict;

use Cwd;
use File::Copy;
use File::Touch;
use File::Path 'rmtree';
use Time::HiRes qw(gettimeofday);

use nanopipe2::paths;
use nanopipe2::messages;
use nanopipe2::utils;

# ------------------------------------------------------------------------
# System programs
# ------------------------------------------------------------------------

my $BUNZIP = "bunzip2";
my $GUNZIP = "gunzip";
my $TAR    = "tar";
my $UNZIP  = "unzip";

# ------------------------------------------------------------------------
# Files
# ------------------------------------------------------------------------

my $MINLENFILE    = qq(input.minlen);
my $QUERYFILE     = qq(input.query);
my $QUERYDONEFILE = qq(input.query.done);

# ------------------------------------------------------------------------
# Parameters
# ------------------------------------------------------------------------

my $query;
my $minlen = 10;

#
# ------------------------------------------------------------------------
# Skep query sequences with lower than a minimum length
# ------------------------------------------------------------------------
#
sub skipMinlen {
	$minlen = nanopipe2::utils::readFile($MINLENFILE) if (-f $MINLENFILE);
	return if ($minlen < 1);

	print "==> Skip minimum query length\n";

	my $type = nanopipe2::utils::readPraefix($QUERYFILE, 1) eq ">" ? 1 : 2;

	rename($QUERYFILE, "$QUERYFILE.tmp");

	open(IN,  "<", "$QUERYFILE.tmp");
	open(OUT, ">", $QUERYFILE);
	if ($type == 1) {
		my ($header, $seq);
		while (my $line = <IN>) {
			if ($line =~ m/^>/) {
				if ($seq) {
					my $len = $seq =~ tr/A-Za-z//;
					print OUT $header, $seq if ($len >= $minlen);
				}
				$header = $line;
				$seq    = "";
			}
			else {
				$seq .= $line;
			}
		}
		if ($seq) {
			my $len = $seq =~ tr/A-Za-z//;
			print OUT $header, $seq if ($len >= $minlen);
		}
	}
	else {
		my $index = 0;
		my $firstline;
		my $print = 0;
		while (my $line = <IN>) {
			if ($index == 0) {
				$firstline = $line;
			}
			elsif ($index == 1) {
				my $len = $line =~ tr/A-Za-z//;
				$print = $len >= $minlen;
				print OUT $firstline, $line if ($print);
			}
			else {
				print OUT $line if ($print);
			}
			$index = ($index + 1) % 4;
		}
	}

	close(OUT);
	close(IN);

	unlink("$QUERYFILE.tmp");
}

#
# ------------------------------------------------------------------------
# Convert a file from fastq to fasta format
# ------------------------------------------------------------------------
#
sub convert2fasta {
	my ($file) = @_;

	print "--> convert2fasta $file\n";

	my $tmpfile = qq($file.tmp);
	rename($file, $tmpfile);
	nanopipe2::utils::execute(qq(awk '(NR - 1) % 4 < 2' $tmpfile | sed 's/@/>/' > $file));
}

#
# ------------------------------------------------------------------------
# Prepare the query data
#
# What means: if query is a plain fastq file, just make a link to
# QUERYFILE.  If it is a fast5 file: convert and make a link.  If it
# is an archive file, extract and put all fastq files into one file
# called QUERYFILE.  Extraction is done in a temporary folder and all
# the extracted content is deleted after preparation.
# ------------------------------------------------------------------------
#
sub run {
	my ($params) = @_;

	my $query  = $params->{query};
	my $minlen = $params->{minlen};

	return if (-f $QUERYDONEFILE);

	print "==> Initialize query\n";

	my $start = time;

	die qq(Query $query is not a file!) if (!-f $query);

	# If query has the name of QUERYFILE: rename
	if ($query eq $QUERYFILE) {
		move($query, "$query.orig");
		$query = "$query.orig";
	}

	my $praefix = nanopipe2::utils::readPraefix($query, 16);

	# Query is a fasta file
	if ($praefix =~ m/^\>/) {
		symlink($query, $QUERYFILE);
	}

	# Query is a fastq file
	elsif ($praefix =~ m/^\@/) {
		convert2fasta($query);
		symlink($query, $QUERYFILE);
	}

	# Extract files into tmp folder
	else {
		my $basedir = qq(/tmp/nanopipe);
		mkdir($basedir) if (!-d $basedir);
		my $tmpdir = sprintf("%s/%010d%02d", $basedir, (gettimeofday)[0], rand(100));

		my $datadir = getcwd();
		my $command;
		my $type = join("", qx(file $query));
		if ($type =~ m/:\s*POSIX tar archive/) {
			$command = "cd $tmpdir; $TAR -xf $datadir/$query";
		}
		elsif ($type =~ m/:\s*bzip2 compressed data, was/) {
			$command = "cd $tmpdir; $BUNZIP -c $datadir/$query >$QUERYFILE.tmp";
		}
		elsif ($type =~ m/:\s*gzip compressed data, was/) {
			$command = "cd $tmpdir; $GUNZIP -c $datadir/$query >$QUERYFILE.tmp";
		}
		elsif ($type =~ m/:\s*bzip2 compressed data/) {
			$command = "cd $tmpdir; $BUNZIP -c $datadir/$query | $TAR -xf -";
		}
		elsif ($type =~ m/:\s*gzip compressed data/) {
			$command = "cd $tmpdir; $GUNZIP -c $datadir/$query | $TAR -xf -";
		}
		elsif ($type =~ m/:\s*Zip archive data/) {
			$command = "cd $tmpdir; $UNZIP $datadir/$query";
		}

		# Nothing to do
		die qq(Cannot extract query file!) if (!$command);

		mkdir($tmpdir);

		# Clean tmp folder and start extraction
		my ($res, $error) = nanopipe2::utils::execute($command);
		if ($res > 0 || $error) {
			nanopipe2::utils::printError($command, $res, $error);
			return;
		}

		# Get and convert fastq files
		my @queryfiles;
		my @files = split(/\n/, qx(find $tmpdir -type f));
		for my $file (@files) {
			my $praefix = nanopipe2::utils::readPraefix($file, 16);

			# Fasta file
			if ($praefix =~ m/^\>/) {
				push(@queryfiles, $file);
			}

			# Fastq file
			elsif ($praefix =~ m/^\@/) {
				convert2fasta($file);
				push(@queryfiles, $file);
			}
		}

		# Concat all fastq files into a global one
		open(QUERY, ">", $QUERYFILE);
		for my $queryfile (@queryfiles) {
			open(F, "<", $queryfile);
			print QUERY join("", <F>);
			close(F);
		}
		close(QUERY);

		rmtree($tmpdir);
	}

	skipMinlen();

	# Mark process to be finished
	touch($QUERYDONEFILE);

	print "Time: " . (time - $start) . " seconds\n";
}

1;
