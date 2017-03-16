#! /usr/bin/env perl 
use warnings;
use strict;

my ($line);
my ($state);
my ($count);
my ($big_line)="";

my ($PREFIX_FASTA) = '>';
my ($MAX_CONSECUTIVE_N) = 20;
if (scalar(@ARGV)>0) {
    $MAX_CONSECUTIVE_N = $ARGV[0];
}
die ("Argument must be non-negative") unless ($MAX_CONSECUTIVE_N>=0);
print STDERR "INFO MaxNs=$MAX_CONSECUTIVE_N\n";
print STDERR "INFO Reading STDIN, Writing STDOUT\n";

my ($prev_ID) ;
my ($next_ID) ;

while (<STDIN>) {
    chomp;
    $line = $_;
    if ($line =~ /^>(\w+)/) {
	$next_ID = $1;
	&output($prev_ID) if (defined($prev_ID));
	$prev_ID = $next_ID;
	$big_line = "";
    } else {
	$big_line .= $line;
    }
}
&output($prev_ID);
print STDERR "Done.\n";

sub output () {
    my ($name) = pop (@_);
    my ($position,$prefix,$suffix);
    my ($nrun);
    my ($contig_per_scaffold) = 0;
    my ($contig_len) = 0;
    my ($seq_len);
    my ($run_len);

    # Ignore N prefix.
    if ($big_line =~ /^(N+)/) {
	$nrun = $1;
	$run_len = length($nrun);
	$big_line = substr ($big_line, $run_len);
    }
    # Print in loop removing Ns.	
    if (length($big_line) > 0) {
	print STDOUT ">$name.$contig_per_scaffold\n";
	while ($big_line =~ /^([^N]+)(N+)/) {
	    $prefix = $1;
	    $nrun = $2;
	    $seq_len = length($prefix);
	    $run_len = length($nrun);
	    $big_line = substr ($big_line, $seq_len+$run_len);
	    print STDOUT "$prefix";
	    if ($run_len <= $MAX_CONSECUTIVE_N) {
		# Print short N-runs but not terminal N-runs.
		print STDOUT "$nrun" unless (length($big_line)==0);
	    } elsif (length($big_line)>0) {
		# Long N-runs, if not terminal, force us to start next contig.
		$contig_per_scaffold++;
		print STDOUT "\n>$name.$contig_per_scaffold\n"
	    }
	}
	# Print non-N suffix. OK if big_line is empty string.
	print STDOUT "$big_line\n";  
    }
}
