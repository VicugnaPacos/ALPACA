#!/usr/local/bin/perl

use strict;

my ($READ_MIN)=15;
my ($MATE_MIN)=1;

my ($first_mate_id) = $ARGV[0];
my ($posmap_file) = $ARGV[1];
my ($fasta_file) = $ARGV[2];
my ($output_prefix)= $ARGV[3];

print STDERR ("Usage: $0 <first_mate_id> <posmap> <fasta> <prefix>") unless (scalar(@ARGV)==4);
print STDERR "Running $0\n first_mate_id=$first_mate_id\n posmap_file=$posmap_file\n fasta_file=$fasta_file\n output_prefix=$output_prefix\n";
die ("Detected only ".scalar(@ARGV)." parameters") unless (scalar(@ARGV)==4);
print STDERR "Running $0\n first_mate_id=$first_mate_id\n posmap_file=$posmap_file\n fasta_file=$fasta_file\n output_prefix=$output_prefix\n";
print STDERR "Fixed parameters\n READ_MIN=$READ_MIN\n MATE_MIN=$MATE_MIN\n"; 

my (%READ_COUNTS);
my (%MATE_COUNTS);
my ($scaffold_counts)=0;
print STDERR "Process posmap file\n";
print STDERR " Expect line-oriented space-delimted format:^<readID> <scfID> ...\n";
print STDERR " Reading $posmap_file ...\n";
open (POSMAP, "<$posmap_file") or die ("Cannot read $posmap_file");
while (<POSMAP>) {
    chomp;
    my ($line)=$_;
    my ($scaffold_id, $read_id);
    my ($is_read);
    if ($line =~ /^(\d+)\s(\d+)/) {
	$scaffold_id = $2;
	$read_id = $1;
    } else {
	die ("Expect line-oriented space-delimted format:^readID scfID ...");
    }
    $is_read = ($read_id < $first_mate_id);
    if ($is_read) {
	$READ_COUNTS{$scaffold_id}++;
    } else {
	$MATE_COUNTS{$scaffold_id}++;
    }
}
close (POSMAP);
print STDERR " Loaded read counts for " . scalar ( keys ( %READ_COUNTS) ) . " scaffolds\n";
print STDERR " Loaded mate counts for " . scalar ( keys ( %MATE_COUNTS) ) . " scaffolds\n";

my ($scaffold_sequences)=0;
my ($is_trusted) = 0;
my ($trusted_scaffolds)=0;
my ($output_suspect) = $output_prefix . ".suspect.fasta";
my ($output_trusted) = $output_prefix . ".trusted.fasta";
open (FASTA, "<$fasta_file") or die ("Cannot read $fasta_file");
open (SUSPECT, ">$output_suspect") or die ("Cannot write $output_suspect");
open (TRUSTED, ">$output_trusted") or die ("Cannot write $output_trusted");
print STDERR "Process fasta file\n";
print STDERR " Reading $fasta_file ...\n";
print STDERR " Writing $output_suspect ...\n";
print STDERR " Writing $output_trusted ...\n";
while (<FASTA>) {
    chomp;
    my ($line)=$_;
    my ($scfid);
    my ($mcnt,$rcnt);
    if (substr($line,0,1) eq ">") {
	++$scaffold_sequences;
	# Internally, strip the "scf" prefix found in fasta files but not posmap files. 
	$scfid = substr($line,4); 
	$rcnt = $READ_COUNTS{$scfid};
	$mcnt = $MATE_COUNTS{$scfid};
	die ("Zero reads for scaffold $scfid") if (!defined($rcnt));
	$mcnt = 0 if (!defined($mcnt));
	if ($mcnt >= $MATE_MIN || $rcnt >= $READ_MIN) {
	    $is_trusted = 1;
	    ++$trusted_scaffolds;
	} else {
	    $is_trusted = 0;
	}
	# add some information to the fasta defline
	$line .= " reads=$rcnt mates=$mcnt";
    }
    if ($is_trusted) {
	print TRUSTED "$line\n";
    } else {
	print SUSPECT "$line\n";
    }
}
close (FASTA);
close (TRUSTED);
close (SUSPECT);
print STDERR " Loaded sequence for $scaffold_sequences scaffolds\n";
print STDERR " Wrote sequence for $trusted_scaffolds trusted scaffolds\n";


