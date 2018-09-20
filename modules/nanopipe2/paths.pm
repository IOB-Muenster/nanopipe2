package nanopipe2::paths;

#
# ========================================================================
# System paths
# ========================================================================
#

use Config;

# The base/project directory for the software
our $PROJDIR = "Please enter your nanopipe2 directory here!";

# The calculate step programs
our $CALCDIR = "$PROJDIR/calculate";

# The systems targets directory
our $TARGETSDIR = "$PROJDIR/targets";

# The tools directory
our $TOOLSDIR = "$PROJDIR/tools";

1;
