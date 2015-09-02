#!/usr/local/bin/perl

use strict;
die ("usage $0 coverage") unless (scalar(@ARGV)==1);

# A "hole" has this much coverage or less.
my ($cvh) = $ARGV[0];
die ("parameter must be non-negative") unless ($cvh>=0);

print STDERR "Reading STDIN...\n";
print STDERR "Writing STDOUT...\n";
&input();
print STDERR "DONE\n";

# example of expected input (must be sorted by columns 1,2)
# SCAFFOLD  WINDOW_BEGING_POSITION  WINDOW_AVERAGE_COVERAGE
#
# 100 0 2  
# 100 1000 5
# 100 2000 4  
# 100 3000 0

my ($leading_hole);
my ($internal_hole);
my ($num_gaps)=0;

sub input () {
    my ($line);
    my ($prev_scaf)=-1;
    my ($prev_window)=-1;
    my ($prev_cov)=-1;
    my ($next_scaf)=-1;
    my ($next_window)=-1;
    my ($next_cov)=-1;
    my ($hole_start)=-1;
    my ($hole_end)=-1;
    $num_gaps=0;
    while (<STDIN>) {
	chomp;
	$line = $_;
	if ($line =~ /^(\w+)\s(\d+)\s(\d+)/) {
	    $next_scaf =        $1;
	    $next_window =       int($2);
	    $next_cov =         int($3);
	    if ($next_scaf != $prev_scaf) {
		# Previous scaffold may have ended in hole. Ignore it.
		# Start processing next scaffold.
		$prev_scaf = $next_scaf;
		$prev_window = -1;
		$prev_cov = -1;
		$num_gaps=0;
		$hole_start = -1;
		$leading_hole = ($next_cov <= $cvh);
		$internal_hole = 0;
	    } elsif ($next_cov != $prev_cov) {
		# We only care about windows where coverage changes.
		if ($next_cov > $cvh) {
		    # High coverage. Possible end of a hole.
		    if ($leading_hole) {
			# End of leading hole.
			$leading_hole = 0;
			$hole_start = -1;
		    } elsif ($internal_hole) {
			# End of internal hole.
			&output($next_scaf,$hole_start,$next_window);	
			$hole_start = -1;
			$internal_hole = 0;
		    }
		} else {
		    # Low coverage. Start of a hole or continuation of a hole.
		    if ($leading_hole) {
			# Continuation. Not interesting.
		    } elsif ($internal_hole) {
			# Continuation. Not interesting.
		    } else {
			$internal_hole = 1;
			$hole_start = $next_window;
		    }
		}
	    }
	    # Get ready for next input
	    $prev_scaf = $next_scaf;
	    $prev_window = $next_window;
	    $prev_cov = $next_cov;
	} else {
	    warn ("Bad line, $line");
	}
    }
}

# process bug: we cannot report the length of the zero coverage at the scaffold end
# because the input file stops at the last window with coverage

sub output () {
    my ($show_scf) = shift;
    my ($show_start) = shift;
    my ($show_end) = shift;
    print STDOUT "$show_scf $show_start $show_end\n";
}

