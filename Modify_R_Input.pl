#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;

# SCRIPT READS INPUT FILES AND USES THEM TO MAKE CHANGES TO R INPUT FILE

die "Usage: $0 <R Input File> <Changes File A> ...\n" unless @ARGV > 1;

my $file = shift @ARGV;

my %header;	# HASH OF HEADER; FORMAT = Gene{Column}
my %out;        # HASH TO BE PRINTED; FORMAT = Sample{Gene}{Cell Contents}
my @g_order;
my %all_samples;
my $counter = 0;
open ( FILE, $file );
while (<FILE>) {
	$counter++;
	my $line = $_;
	chomp $line;
	my @a = split("\t", $line);
	my $cur = 0;
	if ( $counter == 1 ) {
		@g_order = @a;
		foreach ( @a ) {
			$header{$cur} = $_;
			$cur++;
		}
	next;
	}
	my $sample = shift @a;
	$all_samples{$sample} = 1;
	foreach ( @a ) {
		$cur++;
		my $gene = $header{$cur};
		$out{$sample}{$gene} = $_;
	}
}
close FILE;

#print Dumper(\%header);
#	print Dumper(\%out);
#for my $i ( keys %out ) { print "$i\n"; }
#exit;
#########################################################################



#########################################################################
# LOOP THROUGH ALL REMAINING FILES TO POPULATE %modify HASH
my %modify;
while ( @ARGV > 0 ) {
	$file = shift @ARGV;
	open ( FILE, $file);
	while ( <FILE> ) {
		my $line = $_;
		chomp $line;
		my @a = split("\t", $line);
		my $sample = $a[0];
		my $gene = $a[1];
		my $change = $a[2];
		$modify{$sample}{$gene} = $change;
	}
}
#print Dumper(\%modify);
#print Dumper(\%out);
#exit;
#########################################################################



########################################################################
# this is fucked up somehow... and it probably has to do with the way %out is defined
for my $s ( keys %all_samples ) {
	foreach ( @g_order) {
		my $g = $_;
		if ( exists $modify{$s}{$g} ) {
			my $new = $modify{$s}{$g};
			my $old = $out{$s}{$g};
	#		print "$s\t$_\t$old\t$new\n";
		### DOUBLY-NEW DATA REPLACES EVERYTHING
			if ( $new =~ /[\w]+\;[\w]+/ ) { $out{$s}{$g} = $new; }
		### SITUATIONS WHERE MUT AND CNV DATA ARE ALREADY PRESENT
			elsif ( $old =~ /[\w]+\;[\w]+/ ) {
				my @a = split(";", $old);
				if ( $new =~ /\;/ ) { $out{$s}{$g} = $new . $a[1]; }	# NEW MUT DATA
				else { $out{$s}{$g} = $a[0] . ";" . $new; }	# NEW CNV DATA
			}
		## SITUATIONS WHERE ONLY MUT DATA IS PRESENT
			elsif ( $new =~ /\;$/ && $old =~ /\;$/ ) { $out{$s}{$g} = $new; } # NEW MUT DATA
			elsif ( $new !~ /\;$/ && $old =~ /\;$/ ) { $out{$s}{$g} = $old . $new;  } # NEW CNV DATA
		## SITUATIONS WHERE ONLY CNV DATA IS PRESENT
			elsif ( $new =~ /\;$/ && $old !~ /\;$/ ) { $out{$s}{$g} = $new . $old; } # NEW MUT DATA
			elsif ( $new !~ /\;$/ && $old !~/\;$/ ) { $out{$s}{$g} = $new; } # NEW CNV DATA
		}
	}
}
#print Dumper(\%out);
#exit;
###################################	MIGHT BE WORTH CHECKING IF THIS CODE WORKS NOW THAT THE FORMATTING OF %out IS FIXED
# MODIFY %out BASED ON INPUT FILES
#for my $s ( keys %out ) {
#	for my $g ( keys $out{$s} ) {
#		if ( exists $modify{$s}{$g} ) {
#			my $old = $out{$s}{$g};
#			my $change = $modify{$s}{$g};
#			# FOR CHANGING MUTATION AND COPY NUMBER DATA
#			if ( $change !~ /^;/ && $change !~ /;$/ ) {
#				$out{$s}{$g} = $change;
#			}
#			# FOR CHANGING MUTATION OR COPY NUMBER DATA
#			else { 
#				my @b = split(";", $change);
#				my @c = split(";", $old);
#				if ( $change =~ /^;/ ) { $out{$s}{$g} = $c[0] . ";" . $b[1];  }
#				if ( $change =~ /;$/ ) { $out{$s}{$g} = $b[0] . ";" . $c[1];  }
#			}
##			my @b = split(";", $change);
##			print "$b[0]\t$b[1]\n";
##			if ( $change = /[A-Z]+\;[A-Z]+/ ) { $out{$s}{$g} = $modify{$s}{$g}; }
#		}
#	}
#}

#print Dumper(\%out);
#########################################


#########################################################################



#########################################################################
# LOOP THROUGH @g_order AND USE %out TO PRINT TABLE

my $topleft = shift @g_order;
print "$topleft";
foreach (@g_order) { print "\t$_"; }
print "\n";
# DEFINE @order USING %header THEN LOOP THROUGH TO RE-POPULATE THE TABLE
#for my $i ( keys %out) { print "$i\n"; }
for my $s ( keys %all_samples ) {
	print "$s";
	foreach ( @g_order ) {
		my $g = $_;
		print "\t$out{$s}{$g}";
	}
	print "\n";
}
