#!/usr/bin/env perl

# process genbank flat file from stdin and create a
#     mini bedfile for aligator

my $chrom = "";
my $start = 0;
my $end = 0;
my $strand = "+";
my $func = "pseudogene"; #doesn't seem like they really have any ORFs in these sets?
my $gene = "";

foreach my $line (<STDIN>){
	if ($line =~ /^LOCUS/) {
		my @arr = split( /\s+/, $line );
		$chrom = $arr[1];
	} elsif ( $line =~ /^\s+gene\s+(complement\()?<?(\d+)\.\.>?(\d+)/ ){

		#write previous gene and reset
		print( "$chrom\t$start\t$end\t$gene $func\t.\t$strand\n" ) if $gene ne "";
		($start, $end, $strand, $func, $gene) = (0,0,"+","pseudogene","");

		#capture data
		$strand = "-" if $1;
		$start = $2 - 1;
		$end = $3;
	} elsif ($line =~ m{/gene="(.+)"}) {
		$gene = $1;
	} elsif ($line =~ /CDS/) {
		$func = "gene"
	}
}

#write the last line
print( "$chrom\t$start\t$end\t$gene $func\t.\t$strand\n" ) if $gene ne "";

