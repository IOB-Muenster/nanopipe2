#!/usr/bin/env perl


#===========================================================================================================================================#
# Loads maf file into memory. The key of the global hash is the chromosome/transcriptome ID. The value is a hash of all alignment starts
# on that chromosome. It contains several sublists with end and quality symbols of the respective alignments. Needs ca. 1GB Ram.
#-------------------------------------------------------------------------------------------------------------------------------------------#

sub loadMAF {
	
	my $s_counter = 0;
	my $align_start = 0;
	my $align_end = 0;
	my $chr = undef;
	my $cwd = ${$_[0]};
	my $match = ${$_[1]};
	my %hash = ();
	my $prevM = "";
	
	#Search for SNP positions in maf file
	open(DATA, "<", $cwd."/calc.lastalign.maf") or die "Couldn't open maf file, $!";
	while(<DATA>) {
		# Skip comments
		if($_ !~m"^#") {
			# Choose lines staring with s... 
			if($_ =~m"^s"){ 
				# ...which belong to a chromosome
				if($_ =~m"($match)"){ 
					
					# If target id and immediately following query id match the reg expression,
					# I only want to use the target data.
					if ($1 eq $prevM) {
						$prevM = "";
						next;
					}
					$prevM = $1;
					$s_counter += 1;
					# Split by any number of white spaces
					my @line_data = split ' ', $_; 
					my $align_len = length $line_data[-1];
					$chr = $line_data[1]." ";
					$align_start = $line_data[2];
					$align_end = $align_len + $align_start -1;
				}
				# If target and following query id differ, $prevM needs to be resetted for the nex pair.
				else {
					$prevM = "";
				}
			}
			elsif($_ =~m"^p"){
				if ($s_counter == 1){
					$s_counter = 0;
					# Split by any number of white spaces
					my @quality_data = split ' ', $_; 
					my $quality = $quality_data[-1];
					if (exists ($hash{$chr})) {
						if (exists ($hash{$chr}{$align_start})) {
							push(@{$hash{$chr}{$align_start}}, [$align_end, $quality]);
						}
						else {
							$hash{$chr}{$align_start} = [[$align_end, $quality]];
						}	
					}
					else {
						$hash{$chr} = {$align_start => [[$align_end, $quality]]}
					}
					
				}
				
			}
				
		}
	
	}
	close(DATA);
	return %hash;
}


#============================================================================================================================================#
# Takes an ordered array and performes a search for array elements smaller than the
# query element without looking at every single item.
#--------------------------------------------------------------------------------------------------------------------------------------------#
	

sub searchStarts {
	
	# Input array has to be sorted
	my $listRef = $_[1];
	
	
	# Query
	my $query = ${$_[0]};
	
	# Variables
	my @result = (); #undef won't work
	my $lenList = @{$listRef};
	$lenList = sprintf("%.0f", $lenList);
	
	# Position of putative match
	my $posList = $lenList/2;
	$posList = sprintf("%.0f", $posList);
	
	# Variables for controlling the amount of searches
	my $isRep = 1;
	my $lastValue = undef;
	my $isNeighbor = 0;
	
	
	# Check once if the query could be found in the list (smaller than or equal to first)
	if ($query >= $listRef -> [0]) {
	}
	else {
		$isRep =0;
		die "Query not in list!";
	}
	
	
	# If the item has not been found yet...
	while ($isRep == 1) {
		
		# Check if the array index is still valid.
		if ($posList >= 0 and $posList <= $#$listRef) {
			
			# The query is greater or equal to entry at that position.
			if ($listRef -> [$posList] <= $query) {
							
				if ($posList+1 < $#$listRef) {
					
					# The next item should be bigger than the query.
					if ($listRef -> [$posList+1] > $query) {
						$isRep = 0;
						@result = @{$listRef}[0..$posList];
					}
					# Else: Increase the index by 50 %.
					else {
						$isNeighbor = 1;
						$lastValue = $posList;
						$posList = $posList + 0.5 * $posList;
						$posList = sprintf("%.0f", $posList);
						
					}
				}
				# Avoids infinite loops for query also bigger last element
				else {
					
					if ($listRef -> [$#$listRef] <= $query) {
						$isRep = 0;
						@result = @{$listRef};
					}
					else{
						$isRep = 0;
						@result = @{$listRef}[0..$#$listRef-1];
						
					}
					
				}
				
			}
			# The query is smaller than the item at the position. Decrease the index by 50 %.
			else {
				
					# A match was found, but the index was increased to find matches at further positions. 
					# Now, the element at the list index is too big.
					if ($isNeighbor == 1) {
						for my $a (($lastValue..$posList)) {
							if ($listRef -> [$a] > $query) {
								@result = @{$listRef}[0..$a-1];
								$isRep = 0;
								last;
							}
						}
						if (not @result) {
							$isRep = 0;
							die "Not found";
						}
					}
					$lastValue = $posList;
					$posList = $posList - 0.5 * $posList;
					$posList = sprintf("%.0f", $posList);			
			}
			
		}
		# The array index is invalid...
		else {
			
			
			# ... because it is smaller than 0.
			if ($posList < 0) {
				for my $a ((0..$lastValue)) {
					if ($listRef -> [$a] > $query) {
						if (not $a == 0) {
							@result = @{$listRef}[0..$a-1];
							$isRep = 0;
							last;
						}
					}
				}
			}
			# ... because it is bigger than the last index.
			else {
				if ($listRef -> [$#$listRef] <= $query) {
					@result = @{$listRef};
					$isRep = 0;
				}
				else {
					for my $a (($lastValue..$#$listRef)) {
						# The query is smaller, break and save previous
						if ($listRef -> [$a] > $query) {
							@result = @{$listRef}[0..$a-1];
							$isRep = 0;
							last;
						}
					}
				}
			}
			
			if (not @result) {
				$isRep = 0;
				die "Not found";
			}
		}
	}

	return @result
}



#===========================================================================================================================================#
# Checks if the SNP is smaller than the end of candidate alignments. Extracts quality region around SNP, which is at max 20.
# This depends on the position of the SNP in  the alignment. The SNP quality is not extracted. THe quality symbols are translated to a
# p-error and the average p-error for a SNP over all alignments is returned.
#-------------------------------------------------------------------------------------------------------------------------------------------#


sub getQual {
		
	my $quality_str = "";
	
	my $queryChr = ${$_[0]};
	my $snp = ${$_[1]};
	# $hashRef points to a hash of hashes of hashes of arrays of arrays
	my $hashRef = $_[2];
	my $resultsRef = $_[3];
	
	# For every alignment start smaller than the SNP
	foreach my $alignStart (@{$resultsRef}) {
		
		# For every alignment array (alignment end, quality symbols)
		foreach my $alignData (@{$hashRef -> {$queryChr} -> {$alignStart}}) { #@{$hash{$queryChr}{$alignStart}}
			
			my $end = @{$alignData}[0];
			
			# If the alignment ends after or with the SNP
			if ( $end >= $snp) {
				
				my $back_space = $end - $snp;
				my $qualPos = $snp - $alignStart;
				my $snp_start =0;
				my $snp_stop = 0;
				my $string_len = 0;
				my $quality = @{$alignData}[1];
				
				# Grep coordinates of a region of at max 20 nucleotides around the SNP
				if ($qualPos >= 10) {
					$snp_start = $qualPos - 10;
					$string_len = 10;
				}
				else {
					# Start at beginning of quality string
					$snp_start = 0; 
					$string_len = $qualPos;
				}
				
				if ($back_space > 10) {
					# Substr will not return the very last position --> +11
					$snp_stop = $qualPos + 11; 
				}
				else {
					# Trigger, that string to the end is taken
					$snp_stop = 0; 
				}
				
				# Get quality symbols for the region around the SNP; SNP itself is skipped
				my $quality_symb = undef;
				if ($snp_stop != 0) {
					$quality_symb = substr($quality, $snp_start, $string_len);
					$quality_symb = $quality_symb.substr($quality, $qualPos + 1, 10);
					$quality_str = $quality_str.$quality_symb;
				}
				else {
					$quality_symb = substr($quality, $snp_start, $string_len);
					$quality_symb = $quality_symb.substr($quality, $qualPos + 1);
					$quality_str = $quality_str.$quality_symb;
					
				
				}
				
			}
			
		}
		
	}
	
	# Convert ascII symbols to floats (p-error)
	my $ave_p_err = 0;
	my $qual_len = length $quality_str;
	foreach my $ascII (split('', $quality_str)) {
		my $unicode = ord($ascII);
		my $p_error = 10 ** -(($unicode -33) / 10);
		$ave_p_err = $ave_p_err + $p_error;
	}
	
	# Average p_error over all alignments of a single SNP
	$ave_p_err = $ave_p_err/$qual_len;
	$ave_p_err = sprintf("%.4f", $ave_p_err);
	$quality_str = "";
	
	return $ave_p_err;
}

#============================================================ Main Function ================================================================#
# Read in .poly files from Python. Analyze the alignment quality for each SNP position and append the information to the rows. Then, the file
# is rewritten.
#-------------------------------------------------------------------------------------------------------------------------------------------#

use strict;
use Cwd;
use JSON::Tiny qw(decode_json encode_json);

# Use unbuffered stdout
select(STDOUT);
$| = 1;

my %print_dict = ();
my %tidmap = ();
my %enc = ();
my $match = undef;
my $matchRef = undef;
my $chr = undef;
my $cwd = getcwd;
my $cwdRef = \$cwd;
opendir my $dir,$cwd , or die "Cannot open directory: $!";
my @files = readdir $dir;
closedir $dir;


# Get chromosome encodings from tidmap file
open(CHR_HASH, "<", $cwd."/calc.tidmap") or die "Couldn't open tidmap file, $!";
while(<CHR_HASH>) {
	$_ =~ s/^\s+|\s+$//g;
	my @line_split = split('\t', $_);
	$chr =$line_split[0]." ";
	my $encode =$line_split[1];
	$tidmap{$encode} = $chr;
	$enc{$chr} = $encode
	
}

my $tidLen = %tidmap;

# Tidmap file is empty: Data is not good enough. Raise no exception. 
if ($tidLen == 0) {
	exit();
}
# Check whether chromosome or transcriptome
elsif ($tidmap{"1"} =~ m"^chr"){
	print "Chromosome data!\n";
	$match = 'chr'
}
elsif ($tidmap{"1"} =~ m"[A-Z]{2}_\d+[.]\d+"){
	print "Transcriptome data!\n";
	$match = '[A-Z]{2}_\d+[.]\d+'
}
elsif ($tidmap{"1"} =~ m"Pf\d{1}D\d{1}_\d+_v\d+"){
	print "Plasmodium data!\n";
	$match = 'Pf\d{1}D\d{1}_\d+_v\d+'
}
else {
	die "Unrecognized data for encodings in tidmap."
}

$matchRef =\$match;

# Load MAF life
my %hash = loadMAF($cwdRef, $matchRef);

# Read from stdin of python; pre-selectetd SNPs/db results
my $pyInput = <STDIN>;
my %resHash = %{decode_json($pyInput)};
my @chromosomes = keys(%resHash);

foreach my $queryChr (@chromosomes) {
	
	# Variable outfile headers for PlasmoDB and dbSNP searches
	if ($queryChr eq "**header**") {
		next;
	}
	
	my @snps = keys(%{$resHash{$queryChr}});
	
	# resHash only accepts chromosomes like "1", hashes created in this script need "chr1 "
	my $resChr = $queryChr;
	
	# Use references to keep datatype for subs
	if ($match eq 'chr') {
		$queryChr = "chr".$queryChr
	}
	$queryChr = $queryChr." ";
	my $ChrRef = \$queryChr;
	
	# List of all alignment starts for the query chromosome
	my @alignList = keys(%{$hash{$queryChr}});
	@alignList =  sort {$a <=> $b} @alignList;
	my $hashRef = \%hash;
	
	foreach my $snp (@snps) { 
		
		my @results = ();
		my $snpRef =\$snp;
		
		# Look for alignment starts smaller than the SNP
		@results = searchStarts($snpRef, [@alignList]);
		my $resultsRef = \@results;
		
		# Check if the alignment ends are bigger than the SNP. If so, extract quality symbols of a region around it 
		# and calculate the average p-error. At max, the region can span 20 nucleotides, 10 before and after the SNP.
		# This depends on the position of the SNP in the alignment,
		my $avPerr = getQual ($ChrRef, $snpRef, $hashRef, $resultsRef);
		
		# Save to resHash
		splice (@{$resHash{$resChr}{$snp}}, -4, 0, $avPerr."\t");
		
	}
	my @positions = keys(%{$resHash{$resChr}});
	
	@positions =  sort {$a <=> $b} @positions;
	
	# Get variable outfile header, depending on database used
	my $printstr = $resHash{"**header**"};
	
	foreach my $pos (@positions) {
		$printstr = $printstr.$pos."\t".join("",@{$resHash{$resChr}{$pos}});
	}
	
	my $outpath = $cwd."/"."calc.nuccounts.".$enc{$queryChr}.".poly";
	
	open(OUTFILE, ">", $outpath);
	print OUTFILE $printstr;
	close OUTFILE;
}

print "Quality analysis finished!"