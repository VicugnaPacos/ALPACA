#! /usr/bin/env perl 
use warnings;

# Input: Integer
# This process breaks old scaffolds into new ones.
# The new ones will get sequential UIDs, starting with this integer.

# Input: holes file
# List of <scaffold_uid> <start> <stop>
# where the {start,stop} regions need to be removed.
# Assume sorted numerically so coordinates are in order.
# Assume UID like 7180000028980.
# Assume holes are distinct and non-overlapping.

# Input: fasta file
# The scaffolds including scaffold_uid on deflines.
# Assume the FASTA has one line per scaffold sequence (Ns are ok).
# Assume holes do not correspond to Ns (by definition).
# Assume holes are not bounded by Ns (by definition). 
# Otherwise, we could output a scaffold starting with Ns.
# Assume UID like scf7180000028980.

# Output: fasta

use strict;
use warnings;

my ($MIN_SCAFF_LEN) = 0;  # do not output smaller than this
my ($PRINT_DOTS) = 0;  # if set, new scaffold id = old scaffold id _dot_ increment
my ($PRINT_HOLES) = 0;   # usually we do not output sequence from the hole

die ("Usage: $0 <base> <fasta> <holes>") unless (scalar(@ARGV)==3);
my ($base_id)  = $ARGV[0];
my ($in_fasta) = $ARGV[1];
my ($in_holes) = $ARGV[2];
my (%COLLECTION);

print STDERR "Reading holes...\n";
&input_holes;
print STDERR "Writing to STDOUT...\n";
my ($num_scafs_with_holes) = scalar(keys(%COLLECTION));
print STDERR "Loaded $num_scafs_with_holes scaffolds with holes.\n";
print STDERR "Processing FASTA...\n";
&process_fasta;
print STDERR "Done\n";

sub process_fasta {
    my ($line,$in_id,$in_sequence,$doprint);
    my ($alternate)=0;
    my ($current_id) = $base_id;
    open (FA, "<$in_fasta") or die ("Cannot open $in_fasta");
    while (<FA>) {
	chomp;
	$line = $_;
	if (++$alternate == 3) {
	    $alternate = 1;
	}
	if ($alternate == 1) {
	    die ("Expected defline, got: $line") unless ($line =~ /^>scf(\w+)/);
	    $in_id = $1;
	    die ("Existing ID exceeds base for new ID: $line") unless ($in_id < $base_id);
	} else {
	    $in_sequence = $line;
	    if (defined($COLLECTION{$in_id})) {
		# Modify the input scaffold.
		my ($tot,$xx,$stored_string,$stored_id,$stored_start,$stored_stop);
		my ($in_seq_len,$out_seq_len,$current_pos,$out_seq,$dot);
		$in_seq_len = length ($in_sequence);
		$current_pos = 1;
		$dot = 1;
		$out_seq_len = 0;
		$tot = scalar( @ { $COLLECTION{$in_id} } );
		for ($xx=0; $xx<$tot; $xx++) {
		    $stored_string = @ { $COLLECTION{$in_id} } [$xx];
		    die ("Stored garbage: $stored_string") unless 
			$stored_string =~ /(\d+)\s(\d+)\s(\d+)/;
		    $stored_id = $1;
		    $stored_start = $2;
		    $stored_stop = $3;
		    die ("Out of sync: $stored_string") unless ($stored_id eq $in_id);
		    #
		    $out_seq = substr($in_sequence,$current_pos-1,$stored_start-$current_pos);
		    $out_seq_len += length($out_seq);
		    if ($stored_start-$current_pos >= $MIN_SCAFF_LEN) {
			print STDOUT ">scf${current_id}\n" unless ($PRINT_DOTS);
			print STDOUT ">scf${in_id}.${dot}\n" if ($PRINT_DOTS);
			print STDOUT "$out_seq\n";
		    }
		    $current_id++;
		    $dot++;
		    #
		    $out_seq = substr($in_sequence,$stored_start-1,$stored_stop-$stored_start+1);
		    $out_seq_len += length($out_seq);
		    if ($stored_stop-$stored_start+1 >= $MIN_SCAFF_LEN && $PRINT_HOLES) {
			print STDOUT ">scf${current_id}\n" unless ($PRINT_DOTS);
			print STDOUT ">scf${in_id}.${dot}\n" if ($PRINT_DOTS);
			print STDOUT "$out_seq\n";
		    }
		    $current_id++;
		    $dot++;
		    $current_pos = $stored_stop + 1;
		}
		$out_seq = substr($in_sequence,$current_pos-1);
		$out_seq_len += length($out_seq);
		if ($current_pos + $MIN_SCAFF_LEN < $in_seq_len) {
		    print STDOUT ">scf${current_id}\n" unless ($PRINT_DOTS);
		    print STDOUT ">scf${in_id}.${dot}\n" if ($PRINT_DOTS);
		    print STDOUT "$out_seq\n";
		}
		$current_id++;
	    } else {
		# This is the usual case.
		# Just mirror the input to the output.
		print STDOUT ">scf$in_id\n";
		print STDOUT "$in_sequence\n";
	    }
	}
    }
    close (FA);
}

sub input_holes {
    my ($line);
    my (@FIELDS);
    my ($scaffold,$start,$stop);
    open (HF, "<$in_holes") or die ("Cannot open $in_holes");
    while (<HF>) {
	chomp;
	$line = $_;
	@FIELDS = split (' ', $line);
	die ("Cannot parse 3 fields: $line") unless (scalar(@FIELDS)==3);
	$scaffold = $FIELDS[0];
	$start = $FIELDS[1];
	$stop = $FIELDS[2];
	die ("Bad coordinates: $line") unless ($stop >= $start && $start >= 0);
	if (!defined($COLLECTION{$scaffold})) {
	    @{$COLLECTION{$scaffold}}=();
	}
	push (@{$COLLECTION{$scaffold}}, $line);
    }
    close (HF);
}

