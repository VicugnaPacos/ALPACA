#!/usr/bin/perl

use strict;

my $hdr;
my $seq;

while (!eof(STDIN)) {
    $_ = <STDIN>;

    if (m/^>/) {
        if (length($seq) > 0) {
            print "$hdr$seq\n";
        }

        $hdr = $_;
        $seq = "";

    } else {
        chomp;
        $seq .= $_;
    }
}

print "$hdr$seq\n";
