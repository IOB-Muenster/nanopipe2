package nanopipe2::calculate::analyze;

#
# ========================================================================
# This script analyzses a maf file and extracts information about loci
# (locations / hot spots) of target alignments.
#
# Changes
# [2018-01-10] Act with vec
# [2018-01-09] Count gaps
# [2018-01-09] Add highest score
# [2018-01-09] Four equal nucelotid counts are now "X"
# [2018-01-09] Rename extended sids by adding "#" + position
# [2017-12-30] equal nucleotide number by similar numbers
# [2017-12-26] No sort for TWONUCS and THREENUCS
# ========================================================================
#

use strict;

use nanopipe2::paths;
use nanopipe2::messages;
use nanopipe2::utils;

my $QUERYFILE = qq(input.query);

my $CONSENSUSFILE = qq(calc.consensus);
my $DEBUGFILE     = qq(calc.debug);
my $LASTFILE      = qq(calc.lastalign.maf);
my $NUCCOUNTSFILE = qq(calc.nuccounts);
my $TIDMAPFILE    = qq(calc.tidmap);

# ------------------------------------------------------------------------
# Static parameters
# ------------------------------------------------------------------------

# The (second) block size
my $BLOCKSIZE = 1000;

# The bit size for cells/counting
my $BITSIZE = 32;

# The byte size for the cells
my $BYTESIZE = $BITSIZE / 8;

# The number of cells/countings (of ACGT-), the target nucleotid is the 6-th cell!
my $CELLCOUNT = 6;

# The size of the second (inner) block
my $BLOCKBYTES = $BLOCKSIZE * $CELLCOUNT * $BYTESIZE;

# The index for every nuceotid, mapping nuc -> index
my %NUCINDEX = ("A" => 0, "C" => 1, "G" => 2, "T" => 3, "-" => 4);

# If 2 or 3 nucleotids have the same count, write a special symbol in
# IUPAC notation.
my %TWONUCS = (
	"AC" => "M",
	"CA" => "M",
	"AG" => "R",
	"GA" => "R",
	"AT" => "W",
	"TA" => "W",
	"CG" => "S",
	"GC" => "S",
	"CT" => "Y",
	"TC" => "Y",
	"GT" => "K",
	"TG" => "K"
);
my %THREENUCS = (
	"ACG" => "V",
	"CAG" => "V",
	"GAC" => "V",
	"AGC" => "V",
	"CGA" => "V",
	"GCA" => "V",
	"ACT" => "H",
	"CAT" => "H",
	"TAC" => "H",
	"ATC" => "H",
	"CTA" => "H",
	"TCA" => "H",
	"AGT" => "D",
	"GAT" => "D",
	"TAG" => "D",
	"ATG" => "D",
	"GTA" => "D",
	"TGA" => "D",
	"CGT" => "B",
	"GCT" => "B",
	"TCG" => "B",
	"CTG" => "B",
	"GTC" => "B",
	"TGC" => "B"
);

# ------------------------------------------------------------------------
# Parameters
# ------------------------------------------------------------------------

# The maximum memory usable in MB
my $maxmem = 4 * 1024;

# The minimum count of nucleotids
my $mincount = 10;

# The minimum length of the sequence in the results
my $minlen = 40;

# The maximum acceptable gap width between fragments
my $maxgap = 20;

# The maximum acceptable percentage of N
my $maxN = 10;

# Nucleotide counts are considered equal, if the numbers are very
# similar: (n1 * equal <= n2).  It is evaluated after sorting the
# nucelotide counts.
my $equal = 0.8;

# Take highest score?
my $hscore;

# ------------------------------------------------------------------------
# Runtime
# ------------------------------------------------------------------------

my $newtid   = 1;
my $tidcount = 0;

my $metagenomics;

# ------------------------------------------------------------------------
# (Exported) data structures
# ------------------------------------------------------------------------

# The data hash mapping tid -> array -> bitarray/string
our %nuccountsdata = ();

# The mapping tid -> highest score
our %scores;

# The mappend genes: (tid+qid) -> count
our %targetcounts;

# The length of the queries, reagrding the targets
our %querylens;

# The (unique) target count 
our $targetcount = 0;

#
# ------------------------------------------------------------------------
# Reads the maf file and find/fill highest scores
# ------------------------------------------------------------------------
#
sub fillScores {
	print "==> Fill Scores\n";
	my $starttime = time;

	my ($index, $score);

	open(F, "<", $LASTFILE);
	while (my $line = <F>) {
		chomp($line);

		# Skip comments and blank lines
		next if ($line =~ m/^(#|\s*$)/);

		# Take score from "a" line
		if ($line =~ m/^a.*score=(\d+)/) {
			$index = 1;
			$score = $1;
		}

		# Take tid from first "s" line
		elsif ($line =~ m/^s\s+/) {
			if (++$index == 2) {
				my $tid = (split(/\s+/, $line))[1];
				my $x = $scores{$tid};
				$scores{$tid} = $score if (!$x || $x < $score);
			}
		}
	}
	close(F);

	print "Time: " . (time - $starttime) . " seconds\n";
}

#
# ------------------------------------------------------------------------
# Add sequence counts of query nucleotids ""ACGT-" "to nuccountsdata.
# And additionally set the target nucleotid.
# ------------------------------------------------------------------------
#
sub add {
	my ($tid, $pos, $tseq, $qseq) = @_;

	my $len = length($tseq);
	my $a   = $nuccountsdata{$tid};
	$nuccountsdata{$tid} = $a = [] if (!$a);

	for (my $i = 0 ; $i < $len ; $i++) {
		my $tnuc = substr($tseq, $i, 1);

		# If target has a gap: skip
		next if ($tnuc eq "-");

		my $qnuc = uc(substr($qseq, $i, 1));
		if ($qnuc =~ m/[ACGT\-]/) {
			my $index = int($pos / $BLOCKSIZE);

			# Create inner block data (if not done)
			$a->[$index] = "\0" x $BLOCKBYTES if (!defined($a->[$index]));

			my $index2 = ($pos % $BLOCKSIZE) * $CELLCOUNT;

			# Increase the query nucleotid count
			vec($a->[$index], $index2 + $NUCINDEX{$qnuc}, $BITSIZE)++;

			# Check and set target nucleotid
			#my $n = vec($a->[$index], ($index2 + 5) * $BYTESIZE, 8);
			#if ($n) {
			#	my $c = chr($n);
			#	printf("Wrong nuc in %s, pos=%d: %s != %s\n", $tid, $pos, $c, $tnuc) if ($c != $tnuc);
			#}
			vec($a->[$index], ($index2 + 5) * $BYTESIZE, 8) = ord($tnuc);
		}
		$pos++;
	}
}

#
# ------------------------------------------------------------------------
# Find sequences (locations) and count reads and lengths
# ------------------------------------------------------------------------
#
sub analyze {
	print "==> Analyze\n";
	my $starttime = time;

	my ($index, $score);
	my ($tid, $tstart, $tseq);
	my $prevtime;
	my (%queryintarget, %queryintarget2);

	open(F, "<", $LASTFILE);
	while (my $line = <F>) {
		my $time = time;
		if ($time != $prevtime && !nanopipe2::utils::checkMemory($maxmem)) {
			nanopipe2::messages::add("Request is only partly done, because the internal memory limit is reached!");
			last;
		}
		$prevtime = time;

		# Skip comments and blank lines
		next if ($line =~ m/^(#|\s*$)/);

		if ($line =~ m/^a.*score=(\d+)/) {
			$index = 1;
			$score = $1;
		}
		elsif ($line =~ m/^s\s+/) {
			$index++;
			chomp($line);

			# Target
			if ($index == 2) {
				my @a = split(/\s+/, $line);
				($tid, $tstart, $tseq) = ($a[1], $a[2], $a[6]);

				# Skip if highest score is active but the score does not fit
				$index = 4 if ($hscore && $scores{$tid} != $score);
			}

			# Query
			elsif ($index == 3) {
				my @a = split(/\s+/, $line);
				if (!$metagenomics) {
					add($tid, $tstart, $tseq, $a[6]) if ($tid);
				}

				# Set target counts (per target and query id)
				my $tqid = qq($tid-$a[1]);
				if (!$queryintarget{$tqid}) {
					$targetcounts{$tid}++;
					$queryintarget{$tqid} = 1;
				}

				# Increment total target count
				if (!$queryintarget2{$a[1]}) {
					$targetcount++;
					$queryintarget2{$a[1]} = 1;
				}

				# Set query lengths
				my $qlens = $querylens{$tid};
				$querylens{$tid} = {} if (!$qlens);
				my $s = $a[6];
				$s =~ s/-//g;
				$qlens->{length($s)}++;
			}
		}
	}
	close(F);

	print "Time: " . (time - $starttime) . " seconds\n";
}

#
# ------------------------------------------------------------------------
# Get the consensus nucleotide
# ------------------------------------------------------------------------
#
sub getConsensus {
	my ($A, $C, $G, $T) = @_;

	# Sort nucleotid numbers decreasing
	my @n = sort {$b->[0] <=> $a->[0]} ([$A, "A"], [$C, "C"], [$G, "G"], [$T, "T"]);

	# The highest count is the reference
	my $base = $n[0]->[0] * $equal;

	# Single nucleotid
	if ($base > $n[1]->[0]) {
		return $n[0]->[1];
	}

	# Two maximum counts
	if ($base > $n[2]->[0]) {
		return $TWONUCS{join("", ($n[0]->[1], $n[1]->[1]))};
	}

	# Three maximum counts
	if ($base > $n[3]->[0]) {
		return $THREENUCS{join("", ($n[0]->[1], $n[1]->[1], $n[2]->[1]))};
	}

	# Four equal counts
	return "X";
}

#
# ------------------------------------------------------------------------
# Save nucleotid counts helper data
# ------------------------------------------------------------------------
#
sub saveNuccountsHelp {
	print "==> Save Nucleotid Counts Helper\n";
	my $starttime = time;

	my @files = <$NUCCOUNTSFILE.*>;
	for my $file (@files) {
		open(NUCCOUNTS,     "<", $file);
		open(NUCCOUNTSHELP, ">", "$file.help");

		my ($tid, $offset2);
		my ($min, $max, $entries, $count) = (100000000, 0, "", 0);
		my @entries;
		my @a;

		my $line   = <NUCCOUNTS>;
		my $offset = tell NUCCOUNTS;

		while ($line = <NUCCOUNTS>) {
			my $len = length($line);

			if (!($line =~ m/^>([^\s]+)/)) {
				chomp($line);

				@a = split(/\t/, $line);

				# Set minimum and maximum nucleotid count
				my $sum = $a[1] + $a[2] + $a[3] + $a[4];
				$min = $sum if ($sum < $min);
				$max = $sum if ($sum > $max);

				# Every 1000-th necleotid position
				if (($count++ % 1000) == 0) {
					push(@entries, qq($a[0]:$offset));
					@a = ();
				}

				$offset2 = $offset;
			}

			$offset += $len;
		}

		if (@entries) {
			push(@entries, qq($a[0]:$offset2)) if (@a);
			print NUCCOUNTSHELP qq($min $max ), join(" ", @entries), "\n";
		}

		close(NUCCOUNTSHELP);
		close(NUCCOUNTS);
	}

	print "Time: " . (time - $starttime) . " seconds\n";
}

#
# ------------------------------------------------------------------------
# Save data for a specific sequence id and start position
# ------------------------------------------------------------------------
#
sub saveFragment {
	my ($tid, $start, $end) = @_;

	return if ($start == -1);
	my $n = $_[3] =~ tr/N//;
	my $l = length($_[3]);
	return if ($l < $minlen || int($n * 100.0 / $l) > $maxN);

	if ($newtid) {
		if ($tidcount > 0) {
			close(CONSENSUS);
			close(NUCCOUNTS);
		}
		$tidcount++;
		print TIDMAP qq($tid\t$tidcount\n);
		open(NUCCOUNTS, ">", "$NUCCOUNTSFILE.$tidcount");
		open(CONSENSUS, ">", "$CONSENSUSFILE.$tidcount");
		$newtid = 0;
	}

	# Transfer from zero based to one based counting
	my ($start1, $end1) = ($start + 1, $end + 1);
	my $len = $end - $start + 1;

	print CONSENSUS qq(>$tid ($start1:$end1)\n);
	print CONSENSUS $_[3], qq(\n);

	print NUCCOUNTS qq(>$tid\n);
	print NUCCOUNTS $_[4];
}

#
# ------------------------------------------------------------------------
# Save the output
#
# Note: start and stop are zero based - but intializied by -1!
# ------------------------------------------------------------------------
#
sub saveData {
	print "==> Save Data\n";
	my $starttime = time;

	my @totalkeys = keys %nuccountsdata;

	open(TIDMAP, ">", $TIDMAPFILE);

	for my $tid (sort @totalkeys) {
		$newtid = 1;

		my ($start, $stop, $gap, $consseq, $consseqtmp, $nuccount, $nuccounttmp) = (-1, -1, 0, "", "", "", "");

		my @a = @{$nuccountsdata{$tid}};
		for (my $i = 0 ; $i < @a ; $i++) {

			# No info at the block?
			next if (!$a[$i]);

			# Set the position.
			my $pos = $i * $BLOCKSIZE;

			my $b = $a[$i];

			for (my $k = 0 ; $k < $BLOCKSIZE * $CELLCOUNT ; $k += $CELLCOUNT, $pos++) {
				my $A     = vec($b, $k,     $BITSIZE);
				my $C     = vec($b, $k + 1, $BITSIZE);
				my $G     = vec($b, $k + 2, $BITSIZE);
				my $T     = vec($b, $k + 3, $BITSIZE);
				my $qgaps = vec($b, $k + 4, $BITSIZE);

				my $nucsum = $A + $C + $G + $T;

				# The minimum count of nucleotids is reached and the
				# query gaps are less than the number of nucleotids found
				if ($nucsum >= $mincount && $nucsum > $qgaps) {

					# Set start only if not already started (so
					# minimal gaps are considered)
					$start = $pos if ($start == -1);
					$stop  = $pos;
					$gap   = 0;
				}

				# Else - we have a gap
				else {
					$gap++;
				}

				my $cons = $qgaps > $nucsum ? "-" : "";
				$cons = $nucsum >= $mincount ? getConsensus($A, $C, $G, $T) : "N" if (!$cons);

				# Get target nucleotid
				my $tnuc = chr(vec($b, ($k + 5) * $BYTESIZE, 8));

				next if ($start == -1);

				# If the gap is too wide: the sequence is at end
				if ($gap > $maxgap) {
					saveFragment($tid, $start, $stop, $consseq, $nuccount);

					# Initialize to enable start again
					($start, $stop, $consseq, $consseqtmp, $nuccount, $nuccounttmp) = (-1, -1, "", "", "", "");
				}

				# We are still in the sequence
				else {
					# There is no gap
					if ($gap == 0) {
						my $pos1 = $pos + 1;
						$consseq  .= $consseqtmp  if ($consseqtmp);
						$consseq  .= $cons;
						$nuccount .= $nuccounttmp if ($nuccounttmp);
						$nuccount .= qq($pos1\t$A\t$C\t$G\t$T\t$cons\t$tnuc\t$qgaps\n);
						($consseqtmp, $nuccounttmp) = ("", "");
					}

					# There is a gap: store in temporary
					# variables, you don't know if this gap should
					# be included into the consensus sequence.
					else {
						my $pos1 = $pos + 1;
						$consseqtmp  .= $cons;
						$nuccounttmp .= qq($pos1\t$A\t$C\t$G\t$T\t$cons\t$tnuc\t$qgaps\n);
					}
				}
			}
		}

		saveFragment($tid, $start, $stop, $consseq, $nuccount);
	}

	if ($tidcount > 0) {
		close(CONSENSUS);
		close(NUCCOUNTS);
	}
	close(TIDMAP);

	print "Time: " . (time - $starttime) . " seconds\n";
}

#
# ------------------------------------------------------------------------
# Save a debuuging output from nuccountsdata
# ------------------------------------------------------------------------
#
sub debug {
	print "==> debug\n";
	my $starttime = time;

	my @totalkeys = keys %nuccountsdata;

	open(DEBUG, ">", $DEBUGFILE);

	for my $tid (sort @totalkeys) {
		print DEBUG qq(>$tid\n);

		my @a = @{$nuccountsdata{$tid}};
		for (my $i = 0 ; $i < @a ; $i++) {

			# No info at the block?
			next if (!$a[$i]);

			# Set the position.
			my $pos = $i * $BLOCKSIZE;

			my $b = $a[$i];

			for (my $k = 0 ; $k < $BLOCKSIZE * $CELLCOUNT ; $k += $CELLCOUNT, $pos++) {
				my $A = vec($b, $k,     $BITSIZE);
				my $C = vec($b, $k + 1, $BITSIZE);
				my $G = vec($b, $k + 2, $BITSIZE);
				my $T = vec($b, $k + 3, $BITSIZE);

				my $qgaps = vec($b, $k + 4, $BITSIZE);
				my $cons = $A + $C + $G + $T >= $mincount ? getConsensus($A, $C, $G, $T) : "N";
				my $tnuc = chr(vec($b, ($k + 5) * $BYTESIZE, 8));

				print DEBUG qq($pos\t$A\t$C\t$G\t$T\t$cons\t$tnuc\t$qgaps\n);
			}
		}
	}

	close(DEBUG);

	print "Time: " . (time - $starttime) . " seconds\n";
}

#
# ------------------------------------------------------------------------
# Save target counts
# ------------------------------------------------------------------------
#
sub saveTargetCounts {
	print "==> Save Target Counts\n";

	open(TARGETCOUNTS, ">", "calc.targetcounts");

	my $counts = 0;
	open(QUERYFILE, "<", $QUERYFILE);
	while (my $line = <QUERYFILE>) {
		$counts++ if ($line =~ m/^>/);
	}
	close(QUERYFILE);
	print TARGETCOUNTS qq($counts\tNumber of reads in the query\n);

	my @tcounts;
	for my $tid (sort keys %targetcounts) {
		push(@tcounts, sprintf("%d\t%s", $targetcounts{$tid}, $tid));
	}
	print TARGETCOUNTS qq($targetcount\tNumber of reads mapped\n);
	print TARGETCOUNTS join("\n", @tcounts);
	close(TARGETCOUNTS);
}

#
# ------------------------------------------------------------------------
# Save query lengths
# ------------------------------------------------------------------------
#
sub saveQueryLens {
	print "==> Save Query Lengths\n";

	open(QUERYLENS, ">", "calc.querylens");

	my $counts = 0;
	$counts = 0;
	for my $tid (sort keys %querylens) {
		my @lens = keys %{$querylens{$tid}};
		if (@lens) {
			print QUERYLENS qq($tid\t),
			  join("\t", map {$_ = qq($_:) . $querylens{$tid}->{$_}} sort ({$a <=> $b} @lens)),
			  "\n";
		}
	}

	close(QUERYLENS);
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
# The main part
# ------------------------------------------------------------------------
#
sub run {
	my ($params) = @_;

	setConfig($metagenomics, $params->{metagenomics} eq "Y" ? 1 : 0);
	setConfig($maxmem,       $params->{maxmem});
	setConfig($mincount,     $params->{mincount});
	setConfig($minlen,       $params->{minlen});
	setConfig($maxgap,       $params->{maxgap});
	setConfig($maxN,         $params->{maxN});
	setConfig($equal,        $params->{equal});
	setConfig($hscore,       $params->{hscore} eq "Y"       ? 1 : 0);

	if ($hscore) {
		fillScores();
	}

	analyze();
	saveTargetCounts();
	saveQueryLens();
	if (!$metagenomics) {
		saveData();
		saveNuccountsHelp();
	}
}

1;
