package nanopipe2::calculate::last;

#
# ========================================================================
# Run the first alignment with last
# ========================================================================
#

use strict;

use nanopipe2::paths;
use nanopipe2::utils;
use nanopipe2::messages;

my $LASTAL    = "$nanopipe2::paths::TOOLSDIR/bin/lastal";
my $LASTDB    = "$nanopipe2::paths::TOOLSDIR/bin/lastdb";
my $LASTSPLIT = "$nanopipe2::paths::TOOLSDIR/bin/last-split";

my $LASTPARAMSFILE    = qq(input.lastparams);
my $QUERYFILE         = qq(input.query);
my $SUBSTMATRIXFILE   = qq(input.substmatrix);
my $TARGETFILE_UPLOAD = qq(input.targetfile);

my $LASTFILE        = qq(calc.lastalign.maf);
my $TARGETDB_UPLOAD = qq(calc.targetdb);

my $targetdb;

# ------------------------------------------------------------------------
# Parameters
# ------------------------------------------------------------------------

my $metagenomics;
my $query;
my $target;
my $threads;
my $lastparams;
my $substmatrix;

#
# ------------------------------------------------------------------------
# Convert a string with parameters to a hash
# ------------------------------------------------------------------------
#
sub s2hash {
	my ($params) = @_;

	my %h;
	my @a = split(/\s+/, $params);
	for (my $i = 0 ; $i < @a ; $i++) {
		if ($a[$i] =~ m/^\-/) {
			if (!($a[$i + 1] =~ m/^\-/)) {
				$h{$a[$i]} = $a[$i + 1];
				$i++;
			}
			else {
				$h{$a[$i]} = "";
			}
		}
	}

	return \%h;
}

#
# ------------------------------------------------------------------------
# Add new parameters in hash to another hash
# ------------------------------------------------------------------------
#
sub addParams {
	my ($params, $newparams) = @_;
	map {$params->{$_} = $newparams->{$_} if (!$params->{$_})} (keys %{$newparams});
}

#
# ------------------------------------------------------------------------
# Convert a parameter hash to a string
# ------------------------------------------------------------------------
#
sub params2s {
	my ($params) = @_;
	return join(" ", map {"$_ $params->{$_}"} (sort keys %{$params}));
}

#
# ------------------------------------------------------------------------
# Check if target already exists and if not: create it
# ------------------------------------------------------------------------
#
sub initTargetdb {
	my $path = qq($nanopipe2::paths::TARGETSDIR/$target/target);

	# Check if target is already prepared or not
	if (-f qq($path.prj)) {
		$targetdb = $path;
	}
	else {
		print "==> Initialize target database\n";
		$targetdb = $TARGETDB_UPLOAD;
		my $command = qq($LASTDB $targetdb $TARGETFILE_UPLOAD);
		my ($res, $error) = nanopipe2::utils::execute($command);
		if ($res > 0 || $error) {
			nanopipe2::utils::printError($command, $res, $error);
		}
	}
}

#
# ------------------------------------------------------------------------
# Checks if there is an algnment - otherwise there is an error
# ------------------------------------------------------------------------
#
sub check {
	my $ok = 0;

	open(F, "<", $LASTFILE);
	while (my $line = <F>) {
		if ($line =~ m/^a /) {
			$ok = 1;
			last;
		}
	}
	close(F);

	return $ok;
}

#
# ------------------------------------------------------------------------
# Create a matrixfile - if not there
# ------------------------------------------------------------------------
#
sub createMatrixfile {
	if (!-f $SUBSTMATRIXFILE && $substmatrix) {
		my @a = split(/\]\s*,\s*\[/, $substmatrix);
		map {$_ =~ s/\s*,\s*/\t/g} @a;
		open(F, ">", $SUBSTMATRIXFILE);
		printf F "\tA\tC\tG\tT\n";
		printf F "A\t$a[0]\n";
		printf F "C\t$a[1]\n";
		printf F "G\t$a[2]\n";
		printf F "T\t$a[3]\n";
		close(F);
	}
}

#
# ------------------------------------------------------------------------
# Set config parameter (if available)
# ------------------------------------------------------------------------
#
sub setConfig {
	$_[0] = $_[1] if ($_[1]);
}

#
# ------------------------------------------------------------------------
# Execute lastal and last-split command to generate a maf/tab file.
# ------------------------------------------------------------------------
#
sub run {
	my ($params) = @_;

	print "==> Align with last\n";

	my $start = time;

	setConfig($metagenomics, $params->{metagenomics} eq "Y");
	setConfig($query,        $params->{query});
	setConfig($target,       $params->{target});
	setConfig($threads,      $params->{threads});
	setConfig($lastparams,   -f $LASTPARAMSFILE ? nanopipe2::utils::readFile($LASTPARAMSFILE) : $params->{params});
	setConfig($substmatrix,  $params->{substmatrix});

	$threads = 2 if (!$threads);
	$threads = 4 if ($threads > 4);

	createMatrixfile();
	initTargetdb();

	my $p = {-P => $threads, -k => "2"};
	$p->{-p} = $SUBSTMATRIXFILE if (-f $SUBSTMATRIXFILE);
	$p->{-Q} = nanopipe2::utils::readPraefix($query, 1) eq ">" ? "0" : "1";
	addParams($p, s2hash($lastparams));

	my $splitcommand = $LASTSPLIT;
	$splitcommand .= " -m1" if ($metagenomics);
	$splitcommand .= " -m1 -n" if ($target =~ m/Human_transcriptome/);
	my $command = "$LASTAL " . params2s($p) . " $targetdb $query | $splitcommand >$LASTFILE";
	my ($res, $error) = nanopipe2::utils::execute($command);
	if ($res > 0 || $error) {
		nanopipe2::utils::printError($command, $res, $error);
	}
	print "Time: " . (time - $start) . " seconds\n";
}

1;
