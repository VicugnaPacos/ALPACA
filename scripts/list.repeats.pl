#! /usr/bin/perl -w
use strict;

my ($DEBUG) = 0;
# Only credit scaffolds for capturing repeats up to this far apart.
# Set to 0 to disable this check.
my ($MAX_REF_SEPARATION) = 1000000;  

my ($line)="";
my ($msg)="";

die ("Usage: $0 <repeats> <aligns> <min_ident_rep> <min_ident_aln>")
    unless (@ARGV == 4);
my ($REPS) = shift (@ARGV);
my ($ALNS) = shift (@ARGV);
my ($MIN_IDENT_REP) = shift (@ARGV);
my ($MIN_IDENT_ALN) = shift (@ARGV);
print STDERR "repeat_content $REPS $ALNS $MIN_IDENT_REP $MIN_IDENT_ALN\n";

my (@AFLD_chr_start,@AFLD_chr_end,@AFLD_scf_start,@AFLD_scf_end,@AFLD_chr,@AFLD_len,@AFLD_scf);
my ($NUM_ALNS) ;
my (%CHR_OFFSET);

my (@COVERS1,@COVERS2,@COVERS_BOTH);
my ($UNDEFINED,
    $DIFF_SCAFFOLD,
    $SAME_SCAFFOLD_OVERLAP,
    $SAME_SCAFFOLD_DISTINCT,
    $BOTH,
    $NEITHER,
    $ONLY_ONE)
    =(0,1,2,3,4,5,6,7,8);
my (@STATUS)=(
    "UNDEFINED",             # Not expected in output.
    "DIFFERENT_SCAFFOLDS",   # The repeat copies align to different scaffolds.
    "SAME_POSITION",         # The repeat copies align to exactly one position (collapse).
    "SAME_SCAFFOLD",         # The repeat copies have distinct alignments to the same scaffold (almost perfect).
    "BOTH_ALIGNED",          # One alignment spans both copies (perfect).
    "NEITHER_ALIGNED",       # Neither repeat aligned to a scaffold (min 50% coverage required).
    "ONE_ALIGNED");          # Only one repeat copy aligned to a scaffold (min 50% coverage at source). 

my($datestring);

$datestring = scalar localtime();
print STDERR "Start time $datestring\n";
###&load_alignments();

$datestring = scalar localtime();
print STDERR "Ready time $datestring\n";
&process_repeats();

$datestring = scalar localtime();
print STDERR "Finish time $datestring\n";

# Load alignments of scaffold interval to reference interval.
# Only some of these cover reference repeats.
sub load_alignments () { 
    my (@ALNS) ;
    open (ALNS, "<$ALNS") or die ("Cannot open $ALNS");
    # slurp the file
    @ALNS = <ALNS>;
    close (ALNS);
    my($loaded) = scalar(@ALNS);
    print STDERR "$loaded alignments loaded\n";
    my ($prev_chr)="";
    my ($prev_start)=0;
    my ($this_chr);
    my ($this_start);
    my (@FIELDS);
    my ($one_aln);
    my ($ii)=0;
    foreach $one_aln (@ALNS) {
	chomp ($one_aln);
	@FIELDS = split ('\t',$one_aln);
	# validate each line
	$this_chr = $FIELDS[11];
	$this_start = $FIELDS[0];
	if ($this_chr gt $prev_chr) {
	    $CHR_OFFSET{$this_chr} = $ii;
	}
	die ("Alignments not sorted\n$prev_chr $prev_start\n$this_chr $this_start\n") 
	    unless ($this_start >= $prev_start || $this_chr gt $prev_chr);
	$prev_chr = $this_chr;
	$prev_start = $this_start;
	# filter
	next unless ($FIELDS[6]>=$MIN_IDENT_ALN);
	# store each line
	$AFLD_chr_start[$ii] = $FIELDS[0];
	$AFLD_chr_end[$ii]   = $FIELDS[1];
	$AFLD_scf_start[$ii] = $FIELDS[2];
	$AFLD_scf_end[$ii]   = $FIELDS[3];
	$AFLD_len[$ii]   = $FIELDS[4];
	$AFLD_chr[$ii]   = $FIELDS[11];
	$AFLD_scf[$ii]  = $FIELDS[12];
	$ii++;
    }
    $NUM_ALNS = $ii;
    print STDERR "$NUM_ALNS alignments processed.\n";
}

# Stream the alignments of reference to itself.
# These are the reference repeats.
# Process each by comparing to the alignments in memory.
# Known limitation, deemed unimportant: 
#   outputs nothing for a reference sequence that has no alignments whatsoever.
sub process_repeats () {
    my ($line);
    my (@FIELDS);
    my ($a_start,$a_end,$a_chr,$a_len);
    my ($r_start,$r_chr,$r_len);
    my ($midpoint1,$midpoint2);
    my ($half_len,$midpoint_separation);
    my ($good_start)=0;
    my ($display);
    my ($r_count)=0;
    my ($prev_chr)=("");
    my ($prev_start)=0;
    my ($DO_NOT_OPTIMIZE,$PLEASE_OPTIMIZE,$OPTIMIZED)=(1,2,3);
    my ($optimize);
    my ($min1,$max1,$min2,$max2);
    open (REPS, "<$REPS") or die ("Cannot open $REPS");
    while (<REPS>) {
	chomp;
	$line = $_;
	#	print STDERR "\nExamine this repeat\n$line\n" if ($DEBUG);
	$r_count++;
	@COVERS1 = (); # start at zero covering alignments for each repeat
	@COVERS2 = ();
	@COVERS_BOTH = ();
	@FIELDS= split ('\t', $line);
	$r_chr = $FIELDS[11];
## REMOVE THIS TEST WHILE LISTING ALL REPEATS
##	next unless (defined($CHR_OFFSET{$r_chr}));  # limitation: ignore repeats on unknown reference sequences
	next unless ($FIELDS[11] eq $FIELDS[12]);  # only process nearby alignments of sequence to itself
	next unless ($FIELDS[6]>=$MIN_IDENT_REP);
	$r_start = $FIELDS[0];
	$r_len = $FIELDS[4]; # length of the first repeat copy
	$half_len = int ($r_len / 2);  # half_len of the first repeat copy, rounded down
	
	# Only process a repeat once. Keep align(j,k) but discard align(k,j), by insisting j upstream of k.
	($min1,$max1) = ($FIELDS[0] < $FIELDS[1]) ? ($FIELDS[0],$FIELDS[1]) : ($FIELDS[1],$FIELDS[0]);
	($min2,$max2) = ($FIELDS[2] < $FIELDS[3]) ? ($FIELDS[2],$FIELDS[3]) : ($FIELDS[3],$FIELDS[2]);
	$midpoint1 = $min1 + $half_len;
	$midpoint2 = $min2 + $half_len; # assume half_len is about equal in both repeat copies  
	$midpoint_separation = $midpoint2 - $midpoint1;
	$display = "$r_count $r_len $midpoint_separation";
	my ($repeat_status)= $UNDEFINED;
	next if ($MAX_REF_SEPARATION>0 && $midpoint_separation>$MAX_REF_SEPARATION); # only process nearby repeats
	next if ($min2 < $min1);
	# Ignore this repeat if the two copies overalap.
	# Interpretation of overlapping repeats would be too complicated.
	next if ($max1 > $min2);

	# just print the repeats that this program will process
	print STDOUT "$line\n";
	next;
	
	# For optimization, prepare to jump to next chr in array of alignments
	if ($r_chr ne $prev_chr) {
	    die ("Unrecognized chr $r_chr") unless (defined($CHR_OFFSET{$r_chr}));
	    $good_start=$CHR_OFFSET{$r_chr};
	}
	$optimize=$PLEASE_OPTIMIZE;  # update the good_start point on every repeat
	# Assert repeats are sorted by chr, start
	die ("Repeats not sorted\nprev_chr $prev_start\n$r_chr $r_start\n") 
	    unless ($r_start >= $prev_start || $r_chr gt $prev_chr );
	$prev_start = $r_start;
	$prev_chr = $r_chr;
	#
	# FOR ONE REPEAT, LOOP THROUGH ALL ALIGNMENTS
	# (Use good_start optimization to skip over upstream alignments.)
	# (Use midpoint optimization to skip over downstream alignments.)
	# (Exploit the requirement that alignment must cover at least half of a repeat.)
	#
	for (my($ii)=$good_start; $ii<$NUM_ALNS; $ii++) {
	    #	    print STDERR "Examine this alignment ($ii)\n" if ($DEBUG);
	    $a_start = $AFLD_chr_start[$ii];
	    $a_end =   $AFLD_chr_end[$ii];
	    $a_chr =   $AFLD_chr[$ii];
	    $a_len =   $AFLD_len[$ii];
	    #	    print STDERR "Same chr?\n" if ($DEBUG);
	    # good_start optimization: for each repeat, skip over upstream alignments
	    last unless ($a_chr eq $r_chr);
	    #	    print STDERR "Starts upstream?\n" if ($DEBUG);
	    # midpoint optimization: for each repeat, skip over downstream alignments
	    last unless ($a_start <= $midpoint2+1); 
	    #	    print STDERR "Start here next time?\n" if ($DEBUG);
	    # update the good_start optimization
	    if ($optimize==$PLEASE_OPTIMIZE && $a_end >= $r_start) {
		# On next repeat, can skip ahead to first alignment that ends in this repeat.
		# All previous alignments guarranteed not to reach next repeat.
		# That is, next repeat must start downstream of here, and this is first alignment to reach here.
		# This is safe but not optimal; one very long alignment followed by shorties won't get optimized.
		$good_start = $ii;
		$optimize = $OPTIMIZED;   # only update good_start once per repeat
	    }
	    #	    print STDERR "Long enough to consider?\n" if ($DEBUG);
	    # ignore this alignment unless it reaches at least half-way across the first repeat copy
	    next unless ($a_end >= $midpoint1);
	    #	    print STDERR "Is it long enough?\n" if ($DEBUG);
	    # ignore this alignment unless it is long enough to cover half of a repeat copy
	    if ($a_len >= $half_len-1) {
		# Save any alignment that covers either midpoint (and is at least half as long as the repeat).
		# If one alignment spans at least half of both repeat copies then we are done, so
		# stop searching alignments, output this status, continue to next repeat.
		# Else save this alignment if it covers either repeat copy, but continue searching.
		if ($a_start <= $midpoint1 && $a_end >= $midpoint2) {
		    push (@COVERS_BOTH, $ii);
		    last;  
		} elsif ($a_start <= $midpoint1 && $a_end >= $midpoint1) {
		    push (@COVERS1, $ii);
		} elsif ($a_start <= $midpoint2 && $a_end >= $midpoint2) {
		    push (@COVERS2, $ii);
		}
	    }
	}

	if (scalar(@COVERS_BOTH)>0) {
	    $repeat_status = $BOTH;  # both repeat copies covered by one alignment
	} elsif (scalar(@COVERS1)==0 && scalar(@COVERS2)==0) {
	    $repeat_status = $NEITHER;   # neither repeat copy is covered
	} elsif (scalar(@COVERS1)==0 || scalar(@COVERS2)==0) {
	    $repeat_status = $ONLY_ONE;    # exactly one repeat copy is covered
	} else {
	    # Both repeat copies are covered, but not by one alignment.
	    # Unfortunately, if either repeat copy was covered more than once, 
	    # we must run costly all vs all test 
	    # in order to determine if any alignment pair goes to same scaffold.
	    $repeat_status = &many_to_many();
	}
	&output_repeat_status($display,$repeat_status);
    }
    print STDERR "$r_count repeats processed\n";
}

sub output_repeat_status () {
    my ($input) = shift;
    my ($output) = shift;
    my ($verbal) = $STATUS[$output];
    print STDOUT "$input $verbal\n";
}

sub many_to_many() {
    my ($size1) = scalar(@COVERS1);
    my ($size2) = scalar(@COVERS2);
    my ($j,$k);
    my (@JFLD,@KFLD);
    my ($j_start,$j_end,$j_scf,$j_min,$j_max);
    my ($k_start,$k_end,$k_scf,$k_min,$k_max);
    my ($status)=0;
    my ($different_scaffolds)=0;
    my ($same_scaffold_distinct)=0;
    my ($same_scaffold_overlap)=0;
    my ($ii);
    for ($j=0; $j<$size1; $j++) {
	$ii = $COVERS1[$j];
	$j_start = $AFLD_scf_start[$ii];
	$j_end =   $AFLD_scf_end[$ii];
	$j_scf =   $AFLD_scf[$ii];
	($j_min,$j_max) = ($j_start < $j_end) ? ($j_start,$j_end) : ($j_end,$j_start);	
	die ("Bad coords") unless ($j_min<=$j_max);	    
	for ($k=0; $k<$size2; $k++) {
	    $ii = $COVERS2[$k];
	    $k_start = $AFLD_scf_start[$ii];
	    $k_end =   $AFLD_scf_end[$ii];
	    $k_scf =   $AFLD_scf[$ii];
	    ($k_min,$k_max) = ($k_start < $k_end) ? ($k_start,$k_end) : ($k_end,$k_start);
	    die ("Bad coords") unless ($k_min<=$k_max);	 
	    # We are looking at one repeat in the reference.
	    # We filtered such that both repeat copies are nearby in the reference.
	    # We are examining two scaffold-to-reference alignments.
	    # Alignment j is from scaffold j to the upstream repeat copy.
	    # Alignment k is from scaffold k to the downstream repeat copy.
	    # Best outcome is if j and k are same scaffold but non-overlapping positions.
	    if ($j_scf eq $k_scf) {
		if ($j_max <= $k_min || $k_max <=$j_min) {
		    $same_scaffold_distinct++;
		    last;  # Best possible outcome. No need to search further.
		} else {
		    $same_scaffold_overlap++;
		}
	    } else {
		$different_scaffolds++;  
	    }
	}
	# Short circuit the search if best possible outcome is found.
	# In practice, this is the most common outcome.
	# Otherwise, keep searching.
	last if ($same_scaffold_distinct > 0);
    }
    if ($same_scaffold_distinct > 0) {
	# Repeats 1 and 2 mapped without overlap to same scaffold at least once. 
	# Interpretation: tandem duplicate was captured in a scaffold.
	$status = $SAME_SCAFFOLD_DISTINCT;
    } elsif ($different_scaffolds > 0) {
	# Repeats 1 and 2 never mapped to same scaffold twice, but both repeats are present in scaffolds.
	# Interpretation: tandem duplicates got separated but not collapsed.
	$status = $DIFF_SCAFFOLD;
    } else {
	# Repeats 1 and 2 mapped to same location, and that is their only mapping. 
	# Interpretation: tandem duplicate is collapsed in a scaffold.
	$status = $SAME_SCAFFOLD_OVERLAP;
    }
    return ($status);
}
