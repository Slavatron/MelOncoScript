#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;

# PREPARES LIST OF HOTSPOT ANNOTATIONS FOR MELANOMA DATA
# INPUT:  HOTSPOT LIST & RAW MUTATIONS FILE & INDELS
# OUTPUT: HOTSPOT ANNOTATIONS LIST & MODIFICATIONS INPUT

die "Usage: $0 <Hotspots.lst> <Raw Mutations> <Indels>\n" unless @ARGV == 3;

my $hots_file = shift @ARGV;
my $file = shift @ARGV;
my $indel = shift @ARGV;
################ AMINO ACID HASH FOR TRANSLATING INCONVENIENT AMINO ACID DATA FORMAT #############
my %codons = ('CYS'=> 'C', 'ASP'=> 'D', 'SER'=> 'S', 'GLN'=> 'Q', 'LYS'=> 'K',
     'ILE'=> 'I', 'PRO'=> 'P', 'THR'=> 'T', 'PHE'=> 'F', 'ASN'=> 'N', 
     'GLY'=> 'G', 'HIS'=> 'H', 'LEU'=> 'L', 'ARG'=> 'R', 'TRP'=> 'W', 
     'ALA'=> 'A', 'VAL'=>'V', 'GLU'=> 'E', 'TYR'=> 'Y', 'MET'=> 'M', '*'=>'*');
my %subs;	# HASH OF SAMPLE = SUBTYPE
my %mods;	# HASH OF SAMPLE = GENE \t HOT;

#####################################################################
# POPULATE HOTSPOTS HASH
my %hots;	# GENE{AA CHANGE}
open ( HOT, $hots_file );
while (<HOT>) {
	my $line = $_;
	chomp $line;
	my @a = split("_", $line);
	my $g = $a[0];
	my $c = $a[1];
	$hots{$g}{$c} = 1;
}
close HOT;
#print Dumper(\%hots);
#exit;

#####################################################################
# PARSE MUTATIONS FILE
my %changes;	# ALL SAMPLES{GENE}{AA CHANGE}
my %all_samples;
my $counter = 0;
my $s;		# FOR SAMPLE COLUMN
my $g;		# FOR GENE COLUMN
my $aa;		# FOR AA COLUMN
open (FILE, $file);
while (<FILE>) {
	$counter++;
	my $line = $_;
	chomp $line;
	my @a = split("\t", $line);
	if ( $counter == 1 ) {
		my $cur = 0;
		foreach ( @a ) {
			if ( $_ =~ /^Sample$/ ) { $s = $cur; }
			if ( $_ =~ /^SNPEFF_GENE_NAME$/ ) { $g = $cur; }
			if ( $_ =~ /^SNPEFF_AMINO_ACID_CHANGE$/ ) { $aa = $cur; }
			$cur++;
		}
	next;
	}
	my $sample = $a[$s];
	my $gene = $a[$g];
	my $amino = $a[$aa];
	$all_samples{$sample} = 1;
### DEFINE AMINO ACIDS IF THERE IS A CHANGE 
	if ( $amino =~ /\// ) {
		$amino = substr($amino, 2);
		my @b = split("/", $amino);
		my $num = $b[0]; #print "$num\t";
### 	SKIP WEIRD CHARACTERS IN AMINO ACID CHANGE... THESE ARE PROBABLY "HIGH IMPACT CHANGES
		if ( $num !~ /^[\*|\w]+[\d]+[\*|\w]+$/ ) { next; }
		$num =~ s/\W{^\*]//g; #print "$num\n";
		my @c = split(/\d+/, $num);
		$num =~ s/\D+//g;
		my $start = $c[0];
		my $stop = $c[1];
		if ( $start !~ /\*/ ) { $start = uc $start; }
		if ( $stop !~ /\*/ ) { $stop = uc $stop; }
		$amino = $codons{$start} . $num . $codons{$stop};
		$changes{$sample}{$gene}{$amino} = 1;
	}
}
	
close FILE;
#exit;
#print "\n\n\n";
#print Dumper(\%hots);
#print Dumper(\%changes);


#####################################################################
# PARSE INDEL FILE TO CAPTURE NF1 CHANGES
$counter = 0;
my $eff;	# COLUMN FOR FINDING FRAMESHIFTS
open ( INDEL, $indel);
while (<INDEL>) {
	$counter++;
	my $line = $_;
	chomp $line;
	my @a = split("\t", $line);
	if ( $counter == 1 ) {
		my $cur = 0;
		foreach ( @a ) {
			if ( $_ =~ /^Sample$/ ) { $s = $cur; }
			if ( $_ =~ /^SNPEFF_GENE_NAME$/ ) { $g = $cur; }
			if ( $_ =~ /^SNPEFF_EFFECT$/ ) { $eff = $cur; }
			$cur++;
		}
	next;
	}
	my $sample = $a[$s];
	my $gene = $a[$g];
	my $effect = $a[$eff];
	if ( $gene =~ /^NF1$/ && $effect =~ /^FRAME_SHIFT$/ ) {
		$subs{$sample} = 3;
		$mods{$sample} = "$gene\tLOF;";
	}
}
close INDEL;
#print Dumper(\%subs);


#####################################################################
# LOOP THROUGH %hots AND %changes TO POPULATE OUTPUT HASHES
for my $s ( keys %changes ) {  
	for my $g ( keys $changes{$s} ) {
		if ( exists $hots{$g} ) {
			for my $a ( keys $changes{$s}{$g} ) {
				if ( exists $hots{$g}{$a} ) {
					$mods{$s} = "$g\tHOT;";
					if ( $g =~ /^BRAF$/ ) { $subs{$s} = 1; }
					else { $subs{$s} = 2; }
				}
			}
		}
	}
} 

#print Dumper(\%mods);
#print Dumper(\%subs);



########################################################################
# LOOP THROUGH ALL SAMPLES AND CALL THEM WILD-TYPES IF THEY AREN'T DEFINED YET
for my $s ( keys %all_samples ) {
	if ( not exists $subs{$s} ) {
		$subs{$s} = 4;
	}
}

#print Dumper(\%subs);
#print Dumper(\%mods);


#########################################################################
# PRINT OUTPUT FILES

# SUBTYPES FILE
open ( my $tmp, '>', 'Melanoma_Subtypes.anno' );
print $tmp "Sample\tSubtype\n";
for my $i ( keys %subs ) {
	print $tmp "$i\t$subs{$i}\n";
}

# MODIFY INPUT FILE (should NOT have a header)
open ( my $tmp2, '>', 'Melanoma_Modifications.mods' );
for my $j ( keys %mods ) {
	print $tmp2 "$j\t$mods{$j}\n";
}
