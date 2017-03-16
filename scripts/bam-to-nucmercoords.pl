#! /usr/bin/env perl 
use warnings;

use strict;
my ($reffile) = $ARGV[0];
my ($bamfile) = $ARGV[1];
die ("Usage: $0 <ref_fasta> <bam>") unless (defined($reffile));
die ("Usage: $0 <ref_fasta> <bam>") unless (defined($bamfile));
print STDERR "$0 $reffile $bamfile\n";


#  Coords are space-based.

#  From simulation, cigar has been verified to mean:
#    I  - inserted base in the read
#    D  - deleted base from the read
#
#  The spec is contradictory and/or not clear on terminology:
#    I  - "insertion to the reference"
#    D  - "deletion from the reference"

#  action = 0 -- compute the alignment length on the reference
#  action = 1 -- compute the alignment length on the query
#  action = 2 -- compute the length of the query
#  action = 3 -- compute the percent identity, needs the tags in the array

my ($READ_NAME);
my ($CIGAR_VALIDITY);

sub cigarToLength ($$@) {
    my $cigar   = shift @_;
    my $action  = shift @_;

    my $len     = 0;
    my $count   = 0;
    my $label   = undef;
    my $nm      = undef;
    my $id      = 0;

    if ($CIGAR_VALIDITY == 0) {
	if ($action == 4) {
	    $len = 0;
	    # If cigar starts with S, then first number is clip left.
	    if ($cigar =~ m/^(\d+)S/) {
		$len=$1;
	    }
	    return ($len);  # return just the clip left
	}	
	if ($action == 3) {
	    $nm = -1;
	    foreach my $i (@_) {
		if ($i =~ m/^NM:i:(\d+)/) {
		    $nm = $1;
		    last;
		}
	    }
	    $action = 0;
	}
	#print "$cigar\t$len\n";
	while ($cigar =~ m/^(\d+)([MIDNSHP=X])(.*)$/) {
	    $count = $1;
	    $label = $2;
	    $cigar = $3;	    
	    #  From the SAM spec, sum of M, I, S, = and X shall equal the length of the query
	    #  We want the length of the reference covered by the mapped query.
	    if      ($label eq "M") {
		#  Alignment, either match of mismatch
		$len += $count;
	    } elsif ($label eq "I") {
		#  Insertion to the reference
		$len += $count if ($action == 1);
		$len += $count if ($action == 2);
		$id  += $count;
	    } elsif ($label eq "D") {
		#  Deletion from the reference
		$len += $count if ($action == 0);
		$len += $count if ($action == 3);
		$id  += $count;
	    } elsif ($label eq "N") {
		#  Skipped in the reference (intron)
		$len += $count if ($action == 0);
	    } elsif ($label eq "S") {
		#  Soft clipped from the query (present in query)
		#  Length of the read = (clip left) + (aligned span) + (clip right)
		$len += $count if ($action == 2);
		# print STDERR "SOFTCLIPPING detected, not supported.\n";
	    } elsif ($label eq "H") {
		#  Hard clipped from the query (not present in query)
		#$len += $count;
		print STDERR "HARDCLIPPING detected in read $READ_NAME, not supported.\n";		
	    } elsif ($label eq "P") {
		#  Padding - silent deletion from padded reference
		#???
		print STDERR "PADDING detected in read $READ_NAME, not supported.\n";		
	    } elsif ($label eq "=") {
		#  Alignment, match
		$len += $count;
	    } elsif ($label eq "X") {
		#  Alignmment, mismatch
		$len += $count;		
	    } else {
		$CIGAR_VALIDITY = 1;  # invalid label
	    }	    
	    #print "$cigar\t$len\n";
	}
	if (length($cigar) > 0) {
	    $CIGAR_VALIDITY = 2;  # invlaid remnant
	}
	if (defined($nm)) {
	    return(0.0)  if ($nm == -1);	    
	    $nm = $id  if ($nm < $id);
	    #print "nm $nm id $id len $len\n";	    
	    return(int(10000 * (1 - $nm / $len)) / 100);
	}	
    }
    return($len);
}


#  LEN 1 - length of alignment region in reference
#  LEN 2 - length of alignment region in query
#  LEN R - length of reference
#  LEN Q - length of query

#  1 - reference
#  2 - read

my $refNam = undef;
my $refLen = 0;
my %refLen;


open(F, "< $reffile") or die "Failed to open reference sequences $reffile.\n";
print STDERR "Reading $reffile...\n";
while (<F>) {
    if (m/^>(\S+)\s*/) {
        $refLen{$refNam} = $refLen  if (defined($refNam));
        $refNam = $1;
        $refLen = 0;
    } else {
        $refLen += length($_) - 1;
    }
}
close(F);
$refLen{$refNam} = $refLen  if (defined($refNam));
print STDERR "Found ", scalar(keys %refLen), " reference sequences in '$ARGV[0]'.\n";

open (PIPE1, "samtools view $bamfile |") or die "Failed to open pipe from $bamfile.\n";
print STDERR "Reading $bamfile...\n";
print "file1 file2\n";
print "NUCMER\n";
print "\n";
print "[S1]    [E1]    [S2]    [E2]    [LEN 1] [LEN 2] [LEN R] [LEN Q] [% IDY] [TAGS]\n";

while (<PIPE1>) {
    next if ($_ =~ m/^@/);

    my @v = split '\s+', $_;

    my $n2 = $v[0];  #  Read name
    my $fl = $v[1];  #  Flags
    my $n1 = $v[2];  #  Reference name

    $READ_NAME = $n2;

    my $s1 = $v[3];  #  Reference begin position, 1-based
    my $s2 = 0;      #  Read begin position, not known, maybe in cigar
    my $mq = $v[4];  #  Mapping quality

    $CIGAR_VALIDITY = 0;  # cigarToLength() can change this flat
    #  The read start position might be encoded in the cigar (as S, soft clipped from query).
    #  BLASR reports this as XE:i:## in the features.  (MAYBE. SEE BELOW)
    my $e1 = $s1 + cigarToLength($v[5], 0) - 1;  #  Position in reference, in base-based coords
    $s2 = cigarToLength($v[5], 4) + 1;           #  Position on read. Function returns # bases in clip left. Convert to base-based coords 
    my $e2 = $s2 + cigarToLength($v[5], 1) - 1;  #  Position in query read, in base-based coords

    my $l1 = $refLen{$n1};               # length of  reference, not used for math, included in output
    my $l2 = cigarToLength($v[5], 2);    # length of  read, not used for math, included in output 
    my $pi = cigarToLength($v[5], 3, @v);   # percent identity (same for reference or read)

    # THIS IS NOT THE CASE, AT LEAST ON CALMD BAM FILES FROM BLASR WITH SOFT CLIP.
    if (0) {
	#  BLASR XS:i:# and XE:i:# seems to be the start and end of the alignment in the query read.
	#  BLASR XL:i:# is the length of this region.    
	foreach my $i (@v) {
	    $s2 = $1      if ($i =~ m/^XS:i:(\d+)/);
	    $e2 = $1 - 1  if ($i =~ m/^XE:i:(\d+)/);
	}
    }

    if ($CIGAR_VALIDITY == 0) {
	#  Reverse the coords, if the read maps reversed.
	if (($fl & 16) == 16) {
	    ($e2, $s2) = ($s2, $e2);
	}
	my $al1 = ($e1 < $s1) ? $s1 - $e1 : $e1 - $s1;
	my $al2 = ($e2 < $s2) ? $s2 - $e2 : $e2 - $s2;
	$al1 += 1;  # aligned length on assembly = end - start + 1
	$al2 += 1;  # aligned length on read = end - start + 1  
	print "$s1\t$e1\t$s2\t$e2\t$al1\t$al2\t$l1\t$l2\t$pi\t$n1\t$n2\n";
    } elsif ($CIGAR_VALIDITY == 1) {
	warn "Invalid CIGAR label, read $READ_NAME, " .$v[5]. "\n";
    } elsif ($CIGAR_VALIDITY == 2) {
	warn "Invalid CIGAR remnant, read $READ_NAME, " .$v[5]. "\n";
    } else {   # unexpected other condition
	warn "Invalid CIGAR string, read $READ_NAME, " .$v[5]. "\n";
    }
}
close (PIPE1);
print STDERR "Done\n";
