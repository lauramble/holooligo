#! /usr/bin/perl

#program converts vcf file to fastPHASE input format
#September 13th, 2012

$phase = $ARGV[0];
$vcf = $ARGV[1];
$output = $ARGV[2];

unless ($#ARGV==2) {
    print STDERR "Please provide name of input vcf file, filename for fastPHASE formatted output, filename for positions file, and sample size on command line\n\n";
    die;
} #end unless

open(PHASE, $phase);

@positions = ();
@names = ();
$loop_size = $sample_size + 8;
$position = "P";
$type = "";
$prev_breed = "";
$breedID = -1;
%genotypes = ();
%haplotype1 = ();
%haplotype2 = ();

$ind = 0;
@alleles = ();
@lines = ();



print STDERR "Reading in PHASE file...\n";

while(<PHASE>) {
    chomp;
    if (/BEGIN BESTPAIRS1/../END BESTPAIRS1/) {
    	next if /BEGIN BESTPAIRS1/ || /END BESTPAIRS1/;
		if (/\#/) {
			next;
		} else {
			$hap = $ind % 2; #haplotype (0 or 1)
			$index = int($ind/2); #index of individual
			s/[\(\)\[\]]//g; #remove () and []
			@line = split(/\s/, $_); #get the variants
			for ($var = 0; $var <= $#line; $var++) {
				$alleles[$index][$var][$hap] = $line[$var];
				$lines[$var] .= $hap ? "|$line[$var]\t" : $line[$var]; #add to the line for the variant
			}
			$ind++;
			print "$hap $index\n";
		}
	}
}

open(OUTPUT, ">$output");

$ind = 0;

open(VCF, $vcf);

s/\s+$//g for @lines;

while(<VCF>) {
	print;
	if (/\#/) {
		print STDOUT "hello\n";
		print OUTPUT;
	} else {
		print STDOUT "hi\n";
		m/(.*?\s)[0-9]+\|[0-9]+/;
		print OUTPUT "$1$lines[$ind]\n";
		$ind++;
	}
}


print "$alleles[0][0][0] $alleles[0][0][1] $alleles[0][19][0] $alleles[0][19][1]\n";
print "$lines[0]\n$lines[1]\n";
print "$vcf\n";
