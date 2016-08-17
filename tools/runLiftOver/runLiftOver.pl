#!/usr/bin/perl
use strict;


my $in_vcf = $ARGV[0];
my $out_vcf = $ARGV[1];
my $is_maf = $ARGV[2];


#store the header,
#squish columns together to make bed file
#1 226712475 . G A 4.61 PASS PR=0.654;TR=26;TA=5;NR=23;NA=0;TC=GGC
#becomes
#chr1 226712475 226712476 .|G|A|4.61|PASS|PR=0.654;TR=26;TA=5;NR=23;NA=0;TC=GGC
#then run liftover and convert output back to vcf

#for MAF format: 
#CUBN    ENSG00000107611 BCGSC   hg19    10      16996400        16996401        +       Missense_Mutation       SNP     G       C       G                       RG003   RG003N  G       G                                                       Somatic                                                 NA      NA1
# need to put column 5-8 as first 3 columns and squish the others together
#


my @header_lines;
my $hg18_bed = $in_vcf . ".bed";
my $hg19_bed = $out_vcf . ".bed";
my $unmapped = $out_vcf . ".unmapped";
my $error = $out_vcf . "error";

open IN, $in_vcf or die "$! $in_vcf\n";
open OUT, ">$hg18_bed" or die "$! $hg18_bed\n";
while(<IN>){
    chomp;
    if(/^\#.+/){
	push @header_lines, $_;
	next;
    }
    my @cols = split;
    if($is_maf){
	my $name;
	for(@cols){
	    $name = $name . "|";
	}
	chop($name);
	
	print OUT "chr$cols[4]\t$cols[5]\t$cols[6]\t$name\n";
    }
    else{
	my $name = "$cols[2]|$cols[3]|$cols[4]|$cols[5]|$cols[6]|$cols[7]";
	my $end = $cols[1]+1;
	unless($cols[0] =~ /chr/){
	    $cols[0] = "chr$cols[0]";
	}
	print OUT "$cols[0]\t$cols[1]\t$end\t$name\n";
    }

}
close IN;
close OUT;
#run liftOver
#liftOver test.bed /projects/rmorin/common/genomes/hg18ToHg19.over.chain out.bed unmapped.bed

my $cmd = "liftOver $hg18_bed /depot/data2/galaxy/hg18/liftOver/hg18ToHg19.over.chain $hg19_bed $unmapped 2> $error";
print "running $cmd\n";
system($cmd);
#convert back to original VCF with header

open NEWBED, $hg19_bed or die "$! $hg19_bed\n";
open NEWVCF, ">$out_vcf" or die "$! $out_vcf\n";
for(@header_lines){
    print NEWVCF "$_\n";
}
my $n;
while(<NEWBED>){
    chomp;
    $n++;
    if($is_maf){
	my @cols = split /\t/, $_;
	my @maf_cols = split /|/, $cols[3];
	$maf_cols[4] = $cols[1];
	$maf_cols[5] = $cols[2]; #replace genome start/end with hg19 coordinate in MAF columns
	my $to_print = "";
	for (@maf_cols){
	    $to_print = $to_print . "\t";
	}
	chomp($to_print);
	print NEWVCF "$to_print\n";
    }
    else{
	s/\|/\t/g;
	s/^chr//;
	my @cols = split /\t/, $_;
	print NEWVCF "$cols[0]\t$cols[1]\t$cols[3]\t$cols[4]\t$cols[5]\t$cols[6]\t$cols[7]\t$cols[8]\n";
    }
}
close NEWBED;
close NEWVCF;
print "wrote $n hg19 lines to $out_vcf\n";
