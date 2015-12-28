#! /usr/bin/env perl -w
use strict;

die ("Usage: $0 <min_N_len>") unless (scalar(@ARGV)==1);
my ($MIN_N_LEN) = $ARGV[0];
die ("Usage: $0 <positive_integer>") unless ($MIN_N_LEN>=1 && $MIN_N_LEN==int($MIN_N_LEN));

# Data problem:
# Our simulated mates can be partially Ns from gaps.
# Those reads won't map as accurately as others.
# We don't want to use those reads for scaffolding.

# Data description:
# The pairs from the simulator have IDs like so:
#   PE_<SCAFFOLD>@<START>-<END>#<READ> where SCAFFOLD, START, END, READ are integers.
# Note the simulation number is not encoded.
# Thus we cannot guarrantee that IDs are unique across multiple simulations.

# Solution:
# Prepare a list of pair IDs that have Ns in either read.
# Operationally, require at least 5 consecutive Ns to believe the read entered a gap.
# Note that ALLPATHS gap sizes are enriched for multiples of 5.
# Must do this once per simulation (10K, 40K, 160K).

# Example usage:
# cat SIMULATION.10K.?.fastq | ./list_mates_with_n.pl > SIMULATION.10K.mates_with_5N

my ($line);
my ($pair_id);
my ($cnt)=0;
my (%IDS);
my ($count_in)=0;
my ($count_out)=0;


print STDERR "$0\n";
print STDERR "Reading STDIN...\n";
print STDERR "Writing STDOUT...\n";
print STDERR "Min length for N-run to flag: $MIN_N_LEN\n";

my ($PATTERN)="";
for (my ($i) = 1; $i<=$MIN_N_LEN; $i++) {
    $PATTERN .= "N";
}
print STDERR "Flagging pairs with: $PATTERN\n";

while (<STDIN>) {
    $line = $_;
    # parsing FASTQ in groups of 4 lines
    if (++$cnt==5) {
	$cnt=1; 
    }
    if ($cnt==1) { # parse FASTQ def line
	my ($position);
	$count_in++;
	die ("Input is not FASTQ. Expected @ but got this: $line")
	    unless (substr($line,0,1) eq "@");
	$position = index ($line,'#');
	die ("Expected # in this ID: $line")
	    unless ($position >= 2);   # shortest possible line would be @X#1 
	$pair_id = substr($line,1,$position-1);
    } elsif ($cnt==2) { # parse FASTQ sequence line
	my ($position);
	$position = index ($line,$PATTERN);
	if ($position >= 0) {
	    # This sequence is one we want to avoid. 
	    # Output its ID if we haven't already.
	    if (!defined($IDS{$pair_id})) {
		print STDOUT "$pair_id\n";
		$count_out++;
	    }
	    $IDS{$pair_id}++;
	}
    }
}
print STDERR "Reads input: $count_in\n";
print STDERR "Reads output: $count_out\n";
print STDERR "Done.\n";
