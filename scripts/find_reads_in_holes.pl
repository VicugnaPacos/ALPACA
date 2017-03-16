#! /usr/bin/env perl 
use warnings;

use strict;
use warnings;

die ("usage $0 <read_holes> <pair_holes> <posmap>") unless (scalar(@ARGV)==3);
my ($read_holes_file) = $ARGV[0];
my ($pair_holes_file) = $ARGV[1];
my ($posmap_file) = $ARGV[2];
my (@READ_HOLES,%READ_INDEX);
my (@PAIR_HOLES,%PAIR_INDEX);
my ($SAFETY)=5;

print STDERR "Report all reads spanning overlapping coverage holes.\n";
print STDERR "For safety, expand pair holes by SAFETY=$SAFETY bases per side.\n";
print STDERR "The report makes most sense applied to LongReads (ID starts with 1).\n";
print STDERR "Any SimuMates reported must be unmated in the mapping.\n";

print STDERR "Reading $read_holes_file...\n";
&load_holes(\$read_holes_file, \@READ_HOLES, \%READ_INDEX);
print STDERR "Loaded " . scalar(@READ_HOLES) . " read holes in " . 
    scalar(keys(%READ_INDEX)) . " scaffolds.\n";

print STDERR "Reading $pair_holes_file...\n";
&load_holes(\$pair_holes_file, \@PAIR_HOLES, \%PAIR_INDEX);
print STDERR "Loaded " . scalar(@PAIR_HOLES) . " pair holes in " . 
    scalar(keys(%PAIR_INDEX)) . " scaffolds.\n";

print STDERR "Reading $posmap_file...\n";
&process_posmaps(\$posmap_file);
print STDERR "Done\n";

sub process_posmaps () {
    my ($file_ref) = shift;
    my ($line);
    my ($readID,$scfID,$startPos,$stopPos);
    my ($pair_index,$pair_max);
    my ($read_index,$read_max);
    my ($intersections)=0;
    my ($posmaps)=0;
    $pair_max = scalar(@PAIR_HOLES);
    $read_max = scalar(@READ_HOLES);
    my ($read_holes)=0;
    my ($pair_holes)=0;
    my (%CHIMER_READS);
    my ($DEBUG)=0;
    open (POSMAP, "<$$file_ref") or die ("Cannot open input");
    while (<POSMAP>) {
	# Read in one read-to-scaffold mapping.
	# The input need not be sorted.
	# The input may contain SimuMates though the algorithm should pick only unmated LongRead chimera.
	chomp;
	$line = $_;
	if ($line =~ /^(\w+)\s(\w+)\s(\d+)\s(\d+)/) {
	    $posmaps++;
	    $readID=$1;
	    $scfID=$2;
	    $startPos=$3;
	    $stopPos=$4;
	    $pair_index = $PAIR_INDEX{$scfID};
	    print STDERR "DEBUG read maps $startPos-$stopPos id=$readID\n" if ($DEBUG>0);
	    if (defined($pair_index)) {
		# We have previously loaded one or more pair holes for this scaffold.
		# A pair hole is a region of 0X coverage, counting pairs as end-to-end intervals.
		# For each pair hole...
		while ($pair_index < $pair_max) {
		    my ($tmp1) = $PAIR_HOLES[$pair_index];
		    my ($pairHole_Scf,$pairHole_Start,$pairHole_Stop);
		    if ($tmp1 =~ /^(\w+)\s(\d+)\s(\d+)$/) {
			$pairHole_Scf=$1;
			$pairHole_Start=$2;
			$pairHole_Stop=$3;
		    } else {
			die ("Cannot parse saved pair string, $tmp1");
		    }
		    last if ($pairHole_Scf ne $scfID);
		    ++$pair_holes;
		    $pair_index++;
		    $pairHole_Start -= $SAFETY;
		    $pairHole_Stop += $SAFETY;
		    if ($pairHole_Start < $stopPos && $pairHole_Stop > $startPos) {
			# Now, this pair hole overlaps this read mapping
			print STDERR "DEBUG pair hole $pairHole_Start-$pairHole_Stop\n" if ($DEBUG>0);
			$read_index = $READ_INDEX{$scfID};
			if (defined($read_index)) {
			    # We have previously loaded one or more read holes for this scaffold.
			    # A read hole is a region of 1X coverage, counting only long reds.
			    # For each read hole...
			    while ($read_index < $read_max) {
				my ($tmp2) = $READ_HOLES[$read_index];
				my ($readHole_Scf,$readHole_Start,$readHole_Stop);
				if ($tmp2 =~ /^(\w+)\s(\d+)\s(\d+)$/) {
				    $readHole_Scf=$1;
				    $readHole_Start=$2;
				    $readHole_Stop=$3;
				} else {
				    die ("Cannot parse saved read string, $tmp2");
				}
				last if ($readHole_Scf ne $scfID);
				++$read_holes;
				$read_index++;
				if ($readHole_Start < $stopPos && $readHole_Stop > $startPos) {
				    # Now, this read hole overlaps this read mapping
				    print STDERR "DEBUG read hole $readHole_Start-$readHole_Stop\n" if ($DEBUG>0);
				    $intersections++;
				    if ($readHole_Start < $pairHole_Stop && $readHole_Stop > $pairHole_Start) {
					# Now, this read hole also overlaps this pair hole.
					# This read crosses both holes and the holes overlaps.
					# This can only be a LongRead since it spans a pair hole.
					# This read is therefore 1X, the only thing supporting this join.
					# Add this read to the chimer pile.
					# Add some information about the chimeric location.
					die unless ($scfID eq $readHole_Scf);
					die unless ($scfID eq $pairHole_Scf);
					$CHIMER_READS{$readID} .= 
					    "$scfID $pairHole_Start-$pairHole_Stop $readHole_Start-$readHole_Stop ";
					print STDERR "DEBUG yes chimer\n" if ($DEBUG>0); 
				    } else {
					print STDERR "DEBUG yes chimer\n" if ($DEBUG>0); 
				    }
				}					    
			    }
			}
		    }
		}
	    }
	} else {
	    die ("Cannot parse $line");
	}
    }
    close (POSMAP);
    my ($tot_reads) = scalar ( keys (%CHIMER_READS) );
    my ($avg);
    print STDERR "Processed $posmaps read mappings from posmap\n";
    print STDERR "... and process $pair_holes pair holes\n";
    print STDERR "... and process $read_holes read holes\n";
    print STDERR "Found $intersections intersections of posmap + pair hole + read hole\n";
    print STDERR "Total reads in intersections = $tot_reads\n";

    # Output
    my ($ii);
    foreach $ii (keys(%CHIMER_READS)) {
	my ($detail_string) = $CHIMER_READS{$ii};
	my ($chimer_per_read) = 0;
	while ($detail_string =~ /^(\d+)\s(\d+-\d+)\s(\d+-\d+)/) {
	    my ($scaf) = $1;
	    my ($read_pos) = $2;
	    my ($chimer_pos) = $3;
	    ++$chimer_per_read;
	    print STDOUT "Chimer $ii $scaf $read_pos $chimer_pos\n";
	    my ($prefix_len) = length($scaf)+1 + length($read_pos)+1 + length($chimer_pos)+1;
	    $detail_string = substr($detail_string,$prefix_len);
	}
	die ("Cannot parse $CHIMER_READS{$ii}") unless ($chimer_per_read>0); 
    }
}

sub load_holes () {
    my ($file_ref) = shift;
    my ($array_ref) = shift;
    my ($hash_ref) = shift;
    die ("Array not empty") unless (0==scalar(@$array_ref));
    die ("Hash not empty") unless (0==scalar(keys(%$hash_ref)));
    open (HOLES, "<$$file_ref") or die ("Cannot open input");
    my ($line);
    my ($scf,$start,$stop);
    my ($prev_scf) = "";
    my ($prev_start) = -1;
    my ($index);
    while (<HOLES>) {
	chomp;
	$line = $_;
	die ("Cannot parse $line") unless ($line =~ /^(\w+) (\d+) (\d+)$/);
	$scf = $1;
	$start = $1;
	$stop = $1;
	die ("Input not sorted by scf, $line") 
	    if ($scf lt $prev_scf);
	die ("Input not sorted by start, $line") 
	    if ($prev_scf eq $scf && $prev_start > $start);
	if ($scf gt $prev_scf) {
	    $index = scalar(@$array_ref);
	    $hash_ref->{$scf} = $index;
	}
	push (@$array_ref, $line);
	$prev_scf = $scf;
	$prev_start = $start;
    }
    close (HOLES);
}

