package nanopipe2::config;

#
# ========================================================================
# Manage config data
# ========================================================================
#

use strict;

use File::Basename;

use nanopipe2::paths;

# The config values/parameters
our %values;

#
# ------------------------------------------------------------------------
# Returns a hash mapping name -> id (like human -> hg38)
# ------------------------------------------------------------------------
#
sub getTargets {
	my %targets;

	my @dirs = <$nanopipe2::paths::TARGETSDIR/*>;
	for my $dir (@dirs) {
		my $file = qq($dir/config);
		next if (!-f $file);
		my $topic;
		open(F, "<", $file);
		while (my $line = <F>) {
			$line =~ s/^\s+//;
			$line =~ s/#.*$//;
			$line =~ s/\s+$//;
			next if (!$line);
			if ($line =~ m/^\[(.*?)\]$/) {
				$topic = $1;
			}
			elsif ($topic eq "common") {
				my ($name, $value) = split(/=/, $line);
				if ($name eq "name" && !($value =~ m/^\s*$/)) {
					$targets{$value} = basename($dir);
					last;
				}
			}
		}
		close(F);
	}

	return %targets;
}

#
# ------------------------------------------------------------------------
# Load config parameters into values
# ------------------------------------------------------------------------
#
sub load {
	my ($file) = @_;

	if (-f $file) {
		my $topic;

		open(F, "<", $file);
		while (my $line = <F>) {
			$line =~ s/^\s+//;
			$line =~ s/#.*$//;
			$line =~ s/\s+$//;
			next if (!$line);
			if ($line =~ m/^\[(.*?)\]$/) {
				$topic = $1;
			}
			elsif ($topic) {
				my ($name, $value) = split(/=/, $line);
				if ($name && $value) {
					my $v = $values{$topic};
					$values{$topic} = $v = {} if (!$v);
					$v->{$name} = $value;
				}
			}
		}
		close(F);
	}
}

#
# ------------------------------------------------------------------------
# Read the whole config for a target
# ------------------------------------------------------------------------
#
sub init {
	my ($target) = @_;

	load(qq($nanopipe2::paths::TARGETSDIR/config));
	load(qq($nanopipe2::paths::TARGETSDIR/$target/config));
}

1;
