package nanopipe2::utils;

#
# ========================================================================
# Basic functions for nanopipe2
# ========================================================================
#

use strict;

use Fcntl;
use Proc::ProcessTable;

#
# ------------------------------------------------------------------------
# Reads a file and returns the content
# ------------------------------------------------------------------------
#
sub readFile {
	my ($file) = @_;

	my $content;
	if (-r $file) {
		open(F, "<", $file);
		$content = join("", <F>);
		close(F);
	}

	return $content;
}

#
# ------------------------------------------------------------------------
# Write content to a file
# ------------------------------------------------------------------------
#
sub writeFile {
	my ($file, $content) = @_;

	open(F, ">", $file);
	print F $content;
	close(F);
}

#
# ------------------------------------------------------------------------
# Checks memory used
# ------------------------------------------------------------------------
#

sub checkMemory {
	my ($maxmem) = @_;

	my $p = Proc::ProcessTable->new;
	my %info = map {$_->pid => $_} @{$p->table};

	if ($info{$$}->rss / (1024 * 1024) > $maxmem) {
		return 0;
	}

	return 1;
}

#
# ------------------------------------------------------------------------
# Get the sub process ids from a parent process
# ------------------------------------------------------------------------
#
sub getSubpids {
	my ($ppid, $procs) = @_;

	my @pids;
	foreach my $p (@{$procs->table()}) {
		if ($p->ppid == $ppid) {
			push(@pids, $p->pid);
			push(@pids, getSubpids($p->pid, $procs));
		}
	}

	return @pids;
}

#
# ------------------------------------------------------------------------
# Reads the first size bytes from a file
# ------------------------------------------------------------------------
#
sub readPraefix {
	my ($filename, $size) = @_;

	my $praefix;
	sysopen(F, $filename, O_RDONLY) or die qq(Cannot read praefix for $filename!\n);
	read(F, $praefix, $size);
	close(F);

	return $praefix;
}

#
# ------------------------------------------------------------------------
# Execute command
# ------------------------------------------------------------------------
#
sub execute {
	my ($command) = @_;

	print qq(--> Execute: $command\n);

	my $errorfile = qq(calc.execute.error);
	my $res       = system("$command 2>$errorfile");
	my $error     = readFile($errorfile) if ((stat($errorfile))[7] > 0);

	return ($res, $error);
}

#
# ------------------------------------------------------------------------
# Print the error
# ------------------------------------------------------------------------
#
sub printError {
	my ($command, $res, $error) = @_;
	print STDERR "Command: $command\nfailed with exit status $res\nError: $error\n";
}

1;
