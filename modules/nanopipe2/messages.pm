package nanopipe2::messages;

#
# ========================================================================
# Central messages
# ========================================================================
#

use strict;

my $MESSAGESFILE = qq(calc.messages);

our @messages;

sub add {
	push(@messages, $_[0]);
}

sub save {
	open(MESSAGES, ">", $MESSAGESFILE);
	print MESSAGES @messages;
	close(MESSAGES);
}

1;
