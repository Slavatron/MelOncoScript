#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
##############################################################################
# INPUT: 
# - RAW MUTATIONS FILE
# - COPY NUMBER DATA
# - INDEL DATA
# - TARGET GENES LIST
# - HOTSPOTS LIST (OPTIONAL)
# 
# OUTPUT:	
# 1. PRINTS TO STANDARD OUTPUT  GENE-SAMPLE MATRIX FOR THE ONCOPRINT
# 2. PRINTS MUTATION LOAD FOR EACH SAMPLE TO A FILE CALLED "Mut_Count.pl"
##############################################################################


die "Usage: $0 <Point Mutations> <Copy Number> <Indels> <Target Genes> <Hotspots (optional)>\n" unless @ARGV >= 4;

# DEFINE FILE NAMES
my $pm_file = shift @ARGV;
my $cp_file = shift @ARGV;
my $in_file = shift @ARGV;
my $tg_file = shift @ARGV;

# COMPILE LIST OF TARGET GENES
my %target_genes;
my @o_genes;
open (G_FILE, $tg_file);
while ( <G_FILE> ) {
	my $line = $_;
	chomp $line;
	$target_genes{$line} = 1;
	push(@o_genes, $line);
}
close G_FILE;
##############################################################################


##############################################################################
# READ POINT MUTATION FILE AN COMPILE SAMPLE{GENE} HASH OF MUTATIONS
my %out;
my %all_samples; # HASH OF ALL SAMPLES
my %all_genes;	 # HASH OF ALL GENES
my %mut_count;	 # HASH OF SAMPLES PAIRED TO THE NUMBER OF MUTATIONS EACH HAS
my $gs;	# GENE_SAMPLE COLUMN
my $cp;	# CHROM_POS COLUMN
my $counter = 0;
open (PM_FILE, $pm_file);
while (<PM_FILE>) {
	$counter++;
	my $line = $_;
	chomp $line;
	my @a = split("\t", $line);
	if ( $counter == 1 ) {
		my $cur = 0;
		foreach ( @a ) {
			if ( $_ =~ /^GENE_SAMPLE$/ ) { $gs = $cur; }
			if ( $_ =~ /^CHROM_POS$/ ) { $cp = $cur; }
			$cur++;
		}
	next;
	}
	my @b = split("_", $a[$gs]);
	my $gene = $b[0];
	my $sample = $b[1];
	my $coord = $a[$cp];
	$mut_count{$sample} += 1;
	$all_samples{$sample} = 1;
	if ( not exists $target_genes{$gene} ) { next; }
	$all_genes{$gene} = 1;
	$out{$sample}{$gene} = "MUT;";	# DEFINE MUTATION TYPES LATER
	
}
close PM_FILE;
# FILL OUT THE STRUCTURE OF %out SO IT INCLUDES ALL SAMPLES AND ALL GENES;
for my $s ( keys %all_samples ) {
	for my $g ( keys %all_genes ) {
		if ( not exists $out{$s}{$g} ) {
			$out{$s}{$g} = "";
		}
	}
}
#print Dumper(\%out);
#exit;
##############################################################################


##############################################################################
# READ INDEL DATA TO MODIFY %out AND %s_g_coords ACCORDING TO ADD LOF MUTATIONS
$counter = 0;
my $imp;	# IMPACT COLUMN
my $g;		# GENE COLUMN
my $s;		# SAMPLE COLUMN
open ( INDEL, $in_file );
while (<INDEL>) {
	$counter++;
	my $line = $_;
	chomp $line;
	my @a = split("\t", $line);
	if ( $counter == 1 ) {
		my $cur = 0;
		foreach ( @a ) {
			if ( $_ =~ /^SNPEFF_GENE_NAME$/ ) { $g = $cur; }
			if ( $_ =~ /^SNPEFF_EFFECT$/ ) { $imp = $cur;  }
			if ( $_ =~ /^Sample$/ ) { $s = $cur; }
			if ( $_ =~ /^CHROM_POS$/ ) { $cp = $cur; }
			$cur++;
		}
	next;
	}
	my $gene = $a[$g];
	my $sample = $a[$s];
	my $impact = $a[$imp];
	my $coord = $a[$cp];
	# SKIP NON-FRAME SHIFT MUTATIONS
	if ( $impact !~ /^FRAME_SHIFT$/ ) { next; }
## CLEAN UP SAMPLE NAMES;
	if ( $sample =~ /_PE_/ ) {
		my @b = split("_PE_", $sample);
		$sample = $b[0] . $b[1];
	}
	if ( $sample =~ /_1$|_2$/ ) {
		$sample = substr($sample, 0,-2);
	}
#	print "$sample\t$gene\t$impact\n";
## MODIFY %out TO INCLUDE LOF VALUES
	if ( exists $out{$sample}{$gene} ) {
		$out{$sample}{$gene} = "LOF;";
	}

}
close INDEL;
#print Dumper(\%out);
#exit;
##############################################################################



##############################################################################
# ADD COPY-NUMBER INFORMATION
my %s_g_copy;	# HASH OF SAMPLES{GENES} = COPY-NUMBER-EVENT
open ( CP_FILE, $cp_file);
while ( <CP_FILE> ) {
	my $line = $_;	
	chomp $line;
	my @a = split("\t", $line);
	my $sample = $a[0];
	my $gene = $a[1];
	my $event = $a[2];
	$s_g_copy{$sample}{$gene} = $event;
}
close CP_FILE;

# LOOP THROUGH %out ADDING COPY NUMBER EVENTS
for my $s ( keys %s_g_copy ) {
	for my $g ( keys $s_g_copy{$s} ) {
		my $e = $s_g_copy{$s}{$g};
		$e = uc substr($e, 0, 3);	# CHANGE TO UPPER CASE ABBREVIATION
		
		$out{$s}{$g} .= $e;
	}
}
#print Dumper(\%out);
#exit;
##############################################################################



##############################################################################
# COMPILE ORDERED LIST OF PATIENTS (why though...?)
my @o_pts;
for my $i (keys %all_samples) {
	push(@o_pts, $i);
}
#print Dumper(\@o_pts);
#print Dumper(\%all_samples);
#print Dumper(\%mut_count);
#print Dumper(\%out);
#exit;

##############################################################################


##############################################################################;
# LOOP THROUGH OUT IN ORDER OF @o_genes); TO POPULATED TABLE
print "Sample ID";
foreach (@o_genes) { print "\t$_"; }
print "\n";

for my $s (keys %all_samples) {
	print "$s";
	foreach ( @o_genes ) {
		if ( exists $out{$s}{$_} ) {
			print "\t$out{$s}{$_}";
		}
		else { print "\t"; }
	}
	print "\n";
}

##############################################################################


##############################################################################
# PRINT MUTATION LIST TO SEPARATE FILE
open ( my $tmp, '>', 'Mut_Count.anno' );
print $tmp "Sample\tMut_Count\n";
for my $s ( keys %mut_count ) {
	print $tmp "$s\t$mut_count{$s}\n";
}
close $tmp;
