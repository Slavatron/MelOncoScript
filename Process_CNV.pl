#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;



die "Usage: $0 <Copy Number Table> <Gene List>\n" unless @ARGV == 2;

my $file = shift @ARGV;
my $ref = shift @ARGV;


# POPULATE LIST OF TARGET GENES
#############################################################################
my %target; 	# HASH OF TARGET GENES;  HASH{GENE}{CHROMOSOME} = {START-STOP}
open (REF, $ref);
while ( <REF> ) {
	my $line = $_;
	chomp $line;
	my @a = split("\t", $line);
	my $gene = $a[0];
	my $chrom = $a[1];
	my $slice = $a[2];
	$target{$gene}{$chrom} = $slice;
}
close REF;


# PARSE COPY-NUMBER FILE POPULATE LIST OF ALL CHANGES
#############################################################################
my $counter = 0;
my %all;	# HASH OF CHANGES;  HASH{SAMPLE}{CHROMSOME}{START-STOP} = EVENT
open ( FILE, $file );
while (<FILE>) {
	$counter++;
	my $line = $_;
	chomp $line;
	my @a = split("\t", $line);
	if ($counter == 1 ) { next; }
	my $sample = $a[0];
	my $chrom = $a[1];
	if ( $chrom =~ /^chr/ ) { $chrom = substr($chrom, 3); }
	my $start = $a[2];
	my $stop = $a[3];
	my $seg_mean = $a[4];
	my $event = $a[8];
### SKIP COPY-NUMBER CHANGES THAT DO NOT MEET THRESH-HOLD
	if ( $event =~ /^deletion$/ && $seg_mean > -0.5987 ) { next; }
	if ( $event =~ /^amplification$/ && $seg_mean < 0.4946 ) { next; }
	my $slice = $start . "-" . $stop;
#	my @b = ($slice, $event);
	$all{$sample}{$chrom}{$slice} = $event;
}
close FILE;
#print Dumper(\%target);
#print Dumper(\%all);
#exit;



# LOOP THROUGH %all HASH THEN LOOP THROUGH %target HASH, PRINT OVERLAPS
#############################################################################
for my $s ( keys %all ) {
	for my $c ( keys $all{$s} ) {
		for my $slice ( keys $all{$s}{$c} ) {
			my @a = split("-", $slice);
			my $start = $a[0];
			my $stop = $a[1];
			my $e = $all{$s}{$c}{$slice};
			for my $g (keys %target) {
				if ( exists $target{$g}{$c} ) {
					my @b = split("-", $target{$g}{$c});
					my $left = $b[0];
					my $right = $b[1];
					if ( $start <= $left ) {
						if ( $stop >= $right ) {
							print "$s\t$g\t$e\n";
						}
					}
				}
			}
		}
	}
}
