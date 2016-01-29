use warnings;
use strict;
use Data::Dumper;

# READS ANNOTATION FILE AND MUT_COUNT LIST TO PREPARE INITIAL ANNOTATIONS FILE

die "Usage: $0 <Raw Annotations> <Mut_Count.lst>\n" unless @ARGV >= 2;


########################################################################
# DEFINE COLUMNS DESIRED FROM RAW ANNOTATIONS FILE
# NOTE: EACH KEY WILL BE PAIRED WITH THE COLUMN IT WILL HOLDS IN THE INPUT FILE
my %col = (
	'Sample' => '1',
	'Response' => '1',
);
#print Dumper(\%col);


########################################################################
# READ ANNOTATIONS FILE TO COMPILE OUT HASH
my $file = shift @ARGV;
my %out;	# HASH PAIRED WITH THE FULL LINE IN THE OUTPUT FILE
my $counter = 0;
open ( FILE, $file );
while ( <FILE> ) {
	$counter++;
	my $line = $_;
	chomp $line;
	my @a = split("\t", $line);
## DEFINE WHICH COLUMNS TO COLLECT INTO ACTUAL{COLUMN NAME}
	if ( $counter == 1 ) {
		my $cur = 0;
		foreach ( @a ) {
			if ( exists $col{$_} ) { $col{$_} = $cur; }
			$cur++;
		}
	}
## POPULATE %out WITH ALL COLUMNS DEFINED IN %col
	my $sample = $a[$col{'Sample'}];
	my $response = $a[$col{'Response'}];
## ...
	$out{$sample} = "$sample\t$response";
	

}
#print Dumper(\%out);
# HACK TO MODIFY HEADER
$out{'Sample'} = "Sample\tRespType";

########################################################################
# READ OTHER FILES TO COMPILE ADDITIONAL ANNOTATIONS
# NOTE THAT ANNOTATIONS ARE ADDED IN THE ORDER THEY ARE FED TO THIS PERL SCRIPT
while ( @ARGV > 0 ) {
	my $new = shift @ARGV;
	my %add = ();        # HASH OF SAMPLE{ANNOTATION} TO ADD TO %out
	open ( NEW, $new );
	while ( <NEW> ) {
		my $line = $_;
		chomp $line;
		my @a = split("\t", $line);
		my $sample = $a[0];
		my $anno = $a[1];
		$add{$sample} = $anno;
	}
	close NEW;
	for my $s ( keys %out ) {
		if ( exists $add{$s} ) {
			$out{$s} .= "\t$add{$s}";
		}
		else { $out{$s} .= "\t";
		}
	}
}
#print Dumper(\%out);


########################################################################
# PRINT OUTPUT
print "$out{'Sample'}\n";
delete $out{'Sample'};
for my $s ( keys %out ) {
	print "$out{$s}\n";
}

