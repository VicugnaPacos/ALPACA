#! /usr/bin/env perl 
use warnings;

use strict;

die ("Usage: $0 <UIDtoIIDfile>") unless (scalar(@ARGV)==1);
my ($CONVERSION_FILE) = shift @ARGV;

print STDERR "$0\n";
print STDERR "Reading $CONVERSION_FILE...\n";

my %NAMEtoIID;

open(MAP, "< $CONVERSION_FILE") or die ("Failed to open $CONVERSION_FILE");
while (<MAP>) {
    chomp;
    my @v = split '\s+', $_;
    
    if      (scalar(@v) == 3) {
	# read oriented file has IID count UID
        $NAMEtoIID{$v[2]} = $v[1];
    } elsif (scalar(@v) == 6) {
	# pair oriented file has IID1 count UID1 IID2 count UID2
        $NAMEtoIID{$v[2]} = $v[1];
        $NAMEtoIID{$v[5]} = $v[4];
    } else {
        die ("unknown format '$_'");
    }
}
close (MAP);

print STDERR "Reading STDIN...\n";
print STDERR "Writing STDOUT...\n";
print STDOUT "Header.\n";

my ($count_in)=0;
my ($count_out)=0;
my ($line);
while (<STDIN>) {
    chomp;
    $line = $_;
    # Skip first 4 lines of any coords file. Those are header lines.
    if (++$count_in > 4) {
        my ($s1, $e1, $s2, $e2, $al1, $al2, $rl1, $rl2, $id, $utg, $name) = split '\s+', $line;
	die ("Cannot parse: $line") unless (defined($name));
        my $iid;
        my $orient = ($s2 < $e2) ? "f" : "r";
	$iid = $NAMEtoIID{$name};
	die ("Failed to find IID for this UID: $name") unless (defined($iid));
        print STDOUT "$name,$iid\t.\t.\t.\t$s2\t$e2\t$utg\t.\t$s1\t$e1\n";
	++$count_out;
    }
}
print STDERR "Input lines: $count_in\n";
print STDERR "Output lines: $count_out\n";
exit(0);
