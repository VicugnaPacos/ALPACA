#! /usr/bin/env perl 
use warnings;

use strict;
die ("usage $0 <WindowSize>") unless (scalar(@ARGV)==1);
my ($WindowSize)=$ARGV[0];

my (@STARTS,@ENDS);
my (@COVERAGE);
my ($VERBOSE)=0;  # turn on for warnings
my ($bad_coords_input)=0;

print STDERR "Reading STDIN...\n";
&input();
print STDERR "DONE\n";
print STDERR "WARN: Ignored $bad_coords_input inputs with bad coords.\n" if ($bad_coords_input>0);

sub compute () {
    my ($IGNORE_SCAFFOLD) = shift;
    my (@SORT_START,@SORT_END);

    @COVERAGE = ();
    my ($covidx)=0;
    my ($current_coverage)=0;

    @SORT_START = sort {$a <=> $b} @STARTS;
    @SORT_END = sort {$a <=> $b} @ENDS;
    my ($binner);
    my ($s_index,$e_index)=(0,0);
    my ($inputs)=scalar(@SORT_START);
    my ($max_pos)=$SORT_END[$inputs-1];
    my ($position);
    my (@BIN)=();
    my ($modulo)=0;
    my ($sum)=0;
    for ($position=0;$position<=$max_pos;$position++) {
	while ($s_index<$inputs && $SORT_START[$s_index] <= $position) {
	    $current_coverage++;
	    $s_index++;
	}
	while ($e_index<$inputs && $SORT_END[$e_index] <= $position) {
	    $current_coverage--;
	    $e_index++;
	}
	$sum += $current_coverage;
	if (++$modulo == $WindowSize) {
	    my ($average) = ($sum / $WindowSize);
	    $average = int (0.5 + $average);
	    $COVERAGE[$covidx++]=$average;
	    $modulo = 0;
	    $sum=0;
	    # Do not re-initialze $current_coverage. 
	    # Instead, let previous window coverage carry into next window.
	}
    }
}

sub input () {
    my ($line);
    my ($ref_start,$ref_end);
    my ($swap);
    my ($prev_scaf)=-1;
    my ($scaffold);
    @STARTS=();
    @ENDS=();
    while (<STDIN>) {
	chomp;
	$line = $_;
# example of expected input (evidence could be pair ID or read ID; it is ignored)
# SCF BEG  END   EVIDENCE
#
# 100 1696 2673  1111
# 130 3456 2341  1111
	if ($line =~ /^(\w+)\s(\d+)\s(\d+)/) {
	    $scaffold =        $1;
	    $ref_start =       int($2);
	    $ref_end =         int($3);
	    if ($ref_start>=0 && $ref_end>=0 && $ref_start != $ref_end) {
		if ($ref_start > $ref_end) {
		    $swap=$ref_start;   # swap beg & end if necessary
		    $ref_start=$ref_end;
		    $ref_end=$swap;
		}
		if ($scaffold != $prev_scaf) {
		    if ($prev_scaf != -1 && scalar(@STARTS)>0) {
			&compute($prev_scaf);
			&output($prev_scaf);
		    }
		    $prev_scaf = $scaffold;
		    @STARTS = ();
		    @ENDS = ();
		}
		push (@STARTS, $ref_start);
		push (@ENDS, $ref_end);
	    } else {
		# We were given a mapping with start==end. It should not happen but it does.
		warn ("Bad coords, $line") if ($VERBOSE==1);
		$bad_coords_input++;
	    }
	} else {
	    die ("Cannot parse, $line");
	}
    } # WHILE STDIN
    if ($prev_scaf != -1 && scalar(@STARTS)>0) {
	&compute($prev_scaf);
	&output($prev_scaf);
    }
}

sub output () {
    my ($SHOW_SCAFFOLD) = shift;
    die ("Expected scaffold, got $SHOW_SCAFFOLD") unless (length($SHOW_SCAFFOLD)>1);
    my ($ii);
    my ($bin);
    my ($lines)=scalar(@COVERAGE);
    for ($ii=0; $ii<$lines; $ii++) {       
	$bin = $ii * $WindowSize;
	print STDOUT "$SHOW_SCAFFOLD $bin $COVERAGE[$ii]\n";
    }
}

