#!/usr/bin/env perl

use strict;
use warnings;

# Modify the variables to be compatible with dsub routine
my $refGenome = $ENV{REFGEN};
my $vcfFile = $ENV{VCF};
my $outputDir = $ENV{OUTPUT_DIR};

my %cpgs;

print STDERR "Reading variants...";
open(VCF,$vcfFile);
my %vars1;
my %vars2;
my %hetSnps;

while(my $line = <VCF>){
    if($line =~ /^#/){
	next;
    }
    chomp($line);
    my @fields = split("\t",$line);
    my $chr = $fields[0];
    my $pos = $fields[1];
    my $refAl = $fields[3];
    my @altAl;
    if($fields[4] =~ /,/){
	@altAl = split(",",$fields[4]);
    }
    else{
	push(@altAl,$fields[4]);
    }

    my $otherInfo = $fields[9];
    my @otherInfo = split(":",$otherInfo);
    my $gt = $otherInfo[0];
    $gt =~ s/0/$refAl/g;
    for(my $i = 1; $i <= @altAl; $i++){
	$gt =~ s/$i/$altAl[$i-1]/g;
    }

    my ($al1,$al2) = split("/",$gt);
    push(@{$vars1{$chr}},[$pos,$al1]);
    push(@{$vars2{$chr}},[$pos,$al2]);
    if($al1 ne $al2){
	my $key = $chr."\t".($pos-1)."\t".$pos;
	$hetSnps{$key} = [0,$al1,$al2];
    }
}
print STDERR "Done!\n";

open(GENOME,$refGenome);
$/ = ">";
<GENOME>;
while(my $input = <GENOME>){
    chomp($input);
    my @lines = split("\n",$input);
    my $chr = shift(@lines);
    my $seq1 = join("",@lines);
    my $length = length($seq1);
    my $seq2 = $seq1;
    print STDERR "Building diploid ${chr} ($length)...";
    
    foreach my $var (@{$vars1{$chr}}){
	substr($seq1,$var->[0]-1,1,$var->[1]);
    }
    foreach my $var (@{$vars2{$chr}}){
	substr($seq2,$var->[0]-1,1,$var->[1]);
    }
    print STDERR "Done!\n";
    print STDERR "Finding CpG positions in ${chr}...";
    findCpGs($chr,$seq1);
    findCpGs($chr,$seq2);
    print STDERR "Done!\n";
}
close(GENOME);

open(HOMOCPG,">","${outputDir}/homoCpGs.txt");
foreach my $cpg (keys(%cpgs)){
  if($cpgs{$cpg} == 2){
    print HOMOCPG $cpg,"\n";
  }
}
close(HOMOCPG);

open(HETSNPS,">","${outputDir}/hetSnps.txt");
foreach my $hetSnp (keys(%hetSnps)){
  print HETSNPS $hetSnp,"\t",$hetSnps{$hetSnp}->[1],"\t",$hetSnps{$hetSnp}->[2],"\t",$hetSnps{$hetSnp}->[0],"\n";
}
close(HETSNPS);

sub findCpGs {
    my ($chr, $seq) = @_;
    while ($seq =~ /cg/gi) {
	my $startF = $-[0];
	my $endF = $-[0] + 1;
	my $startR = $+[0] - 1;
	my $endR = $+[0];
	$cpgs{"$chr\t$startF\t$endF\t+"}++;
	$cpgs{"$chr\t$startR\t$endR\t-"}++;
	if(exists($hetSnps{"$chr\t$startF\t$endF"})){
	  $hetSnps{"$chr\t$startF\t$endF"}->[0]++
	}
	if(exists($hetSnps{"$chr\t$startR\t$endR"})){
	  $hetSnps{"$chr\t$startR\t$endR"}->[0]++
	}
    }
}
