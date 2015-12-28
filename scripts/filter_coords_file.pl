#! /usr/bin/env perl -w
use strict;

die ("Usage: $0 <list_file>>") unless (scalar(@ARGV)==1);
my ($LIST_FILE) = $ARGV[0];

# Data problem:
# We want to filter out bad mappings.
# Bad type 1: the mapped read has Ns
# Bad type 2: the mate of the mapped read has Ns
# Bad type 3: the mapping does not cover the full length of the read

# Data stream:
# Coords files have mappings for individual reads from invalid pairs, such as pairs with N.
# In coords file, the read id is field #11.
# A read ID has format <pair>#<read>. Example: PE_1_389@86802-96802#2

# Algorithm
# Load the IDs from the exclusion list.
# Filter STDIN to remove lines containing IDs in the exclusion list.
# If the IDs are PAIR IDs (like ABC), then this will remove either read from the pair (like ABC#1).
# Then, require aligned length >= read length

# Example usage:
# cat SIMULATION.10K.?.coords | ./exclude_from_coords.pl  SIMULATION.mates_with_5N > SIMULATION.10K.filter.coords

my ($line);
my ($pair_id,$read_id);
my ($cnt)=0;
my (%IDS);
my ($count_in)=0;
my ($count_out)=0;
my (@FIELDS);

print STDERR "$0\n";

print STDERR "Reading file: $LIST_FILE...\n";
open (LIST, "<$LIST_FILE") or die ("Cannot open $LIST_FILE");
while (<LIST>) {
    $line = $_;
    @FIELDS = split (' ', $line);
    die ("Expect one field per line, got $line") unless (scalar(@FIELDS==1));
    $pair_id = $FIELDS[0];
    $IDS{$pair_id} = 1;
    $count_in++;
}
close (LIST);
print STDERR "Size of exclude list: $count_in\n";
die ("Incomplete hash") unless (scalar(keys(%IDS))==$count_in);
$count_in=0;

print STDERR "Reading STDIN...\n";
print STDERR "Writing STDOUT...\n";

while (<STDIN>) {
    $line = $_;
    ++$count_in;
    my ($print_it) = 1;
    my ($read_length,$aligned_length);
    @FIELDS = split (' ', $line);
    # Print any non-data lines (headers, comments)
    # Print any data lines unless the ID field is on the exclude list
    if (scalar(@FIELDS)==11) {
	$aligned_length = $FIELDS[5];  # coords file reports alignment length in column 6
	$read_length = $FIELDS[7];    # coords file reports read length in column 8
	die ("Invalid aligned length in $line") unless ($aligned_length<=$read_length);
	$read_id = $FIELDS[10];      # coords file reports read ID in column 11
	@FIELDS = split ('#', $read_id);  # expect read id has suffix of #1 or #2
	if (scalar(@FIELDS)==2) {
	    $pair_id = $FIELDS[0];  # expect read id has pair id as its prefix
	    # Filter out bad types
	    if (defined($IDS{$pair_id}) &&
		$aligned_length==$read_length) {
		$print_it = 0;
	    }
	}
    }
    if ($print_it == 1) {
	print STDOUT $line;
	++$count_out;
    }
}
print STDERR "Lines input: $count_in\n";
print STDERR "Lines output: $count_out\n";
print STDERR "Done.\n";
