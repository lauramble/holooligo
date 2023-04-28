#! /usr/bin/perl

$vcf = $ARGV[0];
$output = $ARGV[1];

unless ($#ARGV==1) {
    print STDERR "Please provide name of input vcf file and filename for PHASE formatted output on command line\n\n";
    die;
} #end unless

open(VCF, $vcf);

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

print STDERR "Reading in VCF file...";

while(<VCF>) {
    chomp;
    if ($_=~/\#\#/) {
	next;
    } elsif ($_=~/\#/) {
		@input_line = split(/\s+/, $_);
		$loop_size = @input_line-1;
		for ($a=9; $a<=$loop_size; $a++) {
			@breed = split(/\-/, $input_line[$a]);
			if ($breed[0] ne $prev_breed){
			print "Changing breed: $prev_breed to $breed[0]\n";
			$breedID += 1;
			}
			$prev_breed = $breed[0];
			push(@names, "$breedID #$input_line[$a]");
		} #end for
		next;
    } #end elsif
    
    @input_line = split(/\s+/, $_);
    $ref = $input_line[3];
    $alt = $input_line[4];
    $type .= "S";

	print STDERR "Position: $input_line[1]\n";
	$position .= " $input_line[1]";

    for ($i=9; $i<=$loop_size; $i++) {
		$o = $i - 9;
		$hap1 = $names[$o];
		@genotype = split(":", $input_line[$i]);
	#	print STDERR "Hap name: $hap1; Original genotype: $input_line[$i]; Next step: $genotype[0];";
	#	$genotype[0] =~ s/0/$ref/g;
	#	$genotype[0] =~ s/1/$alt/g;
		$genotype[0] =~ s/\./\?/g;
		@geno = split(/\|/, $genotype[0]);
	#	print STDERR " Final: $genotype[0]\n";
		push @{$haplotype1{$hap1}}, $geno[0];
		push @{$haplotype2{$hap1}}, $geno[1];
    } #end for
    $length = @{$haplotype1{$hap1}}
} #end while

print STDERR "done.\nNumber of loci: $length.\nNow printing output...";
$hap_size = @names;

open(OUTPUT, ">$output");
print OUTPUT "$hap_size\n$length\n$position\n$type";

for ($c=0; $c<=$#names; $c++) {
    print OUTPUT "\n$names[$c]\n";

    $hap1 = $names[$c];
    @hap1_geno = @{$haplotype1{$hap1}};
    @hap2_geno = @{$haplotype2{$hap1}};
    $hap1_line = "";
    $hap2_line = "";

    for ($d=0; $d<=$#hap1_geno; $d++) {
	$hap1_line .= "$hap1_geno[$d] ";
	$cnt2 = length($hap1_line);
	if ($cnt2<500000) {
	    print OUTPUT "$hap1_geno[$d] ";
	} elsif ($cnt2>=500000) {
	    print OUTPUT "\n$hap1_geno[$d] ";
	    $hap1_line = "";
	} #end elsif
	} #end for
	
	print OUTPUT "\n";
	
	for ($d=0; $d<=$#hap1_geno; $d++) {
	$hap2_line .= "$hap2_geno[$d] ";
	$cnt2 = length($hap2_line);
	if ($cnt2<500000) {
	    print OUTPUT "$hap2_geno[$d] ";
	} elsif ($cnt2>=500000) {
	    print OUTPUT "\n$hap2_geno[$d] ";
	    $hap2_line = "";
	} #end elsif
	
    } #end for

#    print OUTPUT "\n$hap1_line";

} #end for

print OUTPUT "\n";

print STDERR "done.\n";



