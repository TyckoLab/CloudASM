#!/usr/bin/env perl

use strict;
use warnings;
use Bio::DB::Sam;

# Modify how the variables are defined in the script
my $bamFile = $ENV{BAM};
my $refGenome = $ENV{REF_GENOME};
my $cpgPosFile = $ENV{CPGPOS};
my $hetSnpPosFile = $ENV{SNPS};
my $outputDir = $ENV{OUTPUT_DIR};


my $sam = Bio::DB::Sam->new(-bam=>$bamFile,
			    -fasta=>$refGenome,
			    -autoindex=>1);

print STDERR "Creating hash of CpG positions...";
open(CPGS,$cpgPosFile);
my %cpgs;
while(my $line = <CPGS>){
  chomp($line);
  if($line =~ /^\s*$/){
    next;
  }
  my ($chr,$start,$pos,$strand) = split("\t",$line);
  $cpgs{"$chr:$pos"} = $strand;
}
close(CPGS);
print STDERR " Done!\n";

my $processedSnps = 0;
open(SNPS,$hetSnpPosFile);
print "Chromosome\tStart\tEnd\tAllele1/Methylated\tAllele1/Unmethylated\tAllele2/Methylated\tAllele2/Unmethylated\tNumber.of.good.reads\tis.on.heterozygous.CpG\tAllele1\tAllele2\n";
while(my $snpLine = <SNPS>){
  chomp($snpLine);
  if($snpLine =~ /^\s*$/){
    next;
  }
  my ($snpChr,$start,$snpPos,$al1,$al2,$isHetCpg) = split("\t",$snpLine);
  my $noNegAs = 0;
  my $noPosTs = 0;
  my $orAl1 = $al1;
  my $orAl2 = $al2;
  if(($al1 =~ /A/i and $al2 =~ /C/i) or ($al2 =~ /A/i and $al1 =~ /C/i)){
      $al1 = "A";
      $al2 = "[CT]";
      $orAl1 = "A";
      $orAl2 = "C";
  }
  elsif(($al1 =~ /T/i and $al2 =~ /G/i) or ($al2 =~ /T/i and $al1 =~ /G/i)){
      $al1 = "T";
      $al2 = "[GA]";
      $orAl1 = "T";
      $orAl2 = "G";
  }
  elsif(($al1 =~ /C/i and $al2 =~ /G/i) or ($al2 =~ /C/i and $al1 =~ /G/i)){
      $al1 = "[CT]";
      $al2 = "[GA]";
      $orAl1 = "C";
      $orAl2 = "G";
  }
  elsif(($al1 =~ /A/i and $al2 =~ /G/i) or ($al2 =~ /A/i and $al1 =~ /G/i)){
      $noNegAs = 1;
  }
  elsif(($al1 =~ /C/i and $al2 =~ /T/i) or ($al2 =~ /C/i and $al1 =~ /T/i)){
      $noPosTs = 1;
  }
      
  my %reads;
  my $getHetSnpReads = sub {
    my ($chr,$pos,$pileups) = @_[0..2];
    if($pos != $snpPos && !exists($cpgs{"$chr:$pos"})){
      return;
    }
    elsif($pos == $snpPos){
      for my $pileup (@$pileups) {
	if($pileup->indel or $pileup->is_refskip){
	  next; # Ignore reads with indels in this position
	}
	my $alignment = $pileup->alignment;
	my $qscore = $alignment->qscore->[$pileup->qpos];
	if($qscore < 20){
	  next; # Ignore reads with base quality less than 20 in this position
	}
	my $qbase  = substr($alignment->qseq,$pileup->qpos,1);
	my $strand = $alignment->strand;
	
	if($noPosTs and (($qbase =~ /T/i) && $strand == 1)){
	  next; # Ignore reads with T basecall from forward when genotype is C/T
	}
	elsif($noNegAs and  (($qbase =~ /A/i) && $strand == -1)){
	  next; # Ignore reads with basecall A from reverse strand when genotype is A/G
	}
	my $readName = $alignment->qname;
	unless(exists($reads{$readName})){
	  $reads{$readName} = [0,0,0,0];
	} 
	if($qbase =~ /$al1/i){
	  $reads{$readName}->[0]++;
	}
	elsif($qbase =~ /$al2/i){
	  $reads{$readName}->[1]++;
	}
      }
    }
    else{
      my $cStrand = $cpgs{"$chr:$pos"};
      for my $pileup (@$pileups){
	if($pileup->indel or $pileup->is_refskip){
	  next;
	}
	my $alignment = $pileup->alignment;
	my $qscore = $alignment->qscore->[$pileup->qpos];
	if($qscore < 20){
	  next;
	}
	my $qbase  = substr($alignment->qseq,$pileup->qpos,1);
	my $strand = $alignment->strand;
	my $readName = $alignment->qname;
	unless(exists($reads{$readName})){
	  $reads{$readName} = [0,0,0,0];
	} 
	if(($qbase=~ /C/i && $strand == 1 && $cStrand eq "+") or ($qbase=~ /G/i && $strand == -1 && $cStrand eq "-")){
	  $reads{$readName}->[2]++;
	}
	elsif(($qbase=~ /T/i && $strand == 1 && $cStrand eq "+") or ($qbase=~ /A/i && $strand == -1 && $cStrand eq "-")){
	  $reads{$readName}->[3]++;
	}
      }
    }
  };
  $sam->fast_pileup("$snpChr:$snpPos-$snpPos",$getHetSnpReads);
  my @allelicMethTab = (0,0,0,0,0);
  for my $snpRead (values(%reads)){
    $allelicMethTab[0] += $snpRead->[0] * $snpRead->[2];
    $allelicMethTab[1] += $snpRead->[0] * $snpRead->[3];
    $allelicMethTab[2] += $snpRead->[1] * $snpRead->[2];
    $allelicMethTab[3] += $snpRead->[1] * $snpRead->[3];
    if(($snpRead->[0] > 0 or $snpRead->[1] > 0) and ($snpRead->[2] > 0 or $snpRead->[3] > 0)){
      $allelicMethTab[4]++;
    }
  }
  
  #if($allelicMethTab[4] > 0){
  print $snpChr,"\t",$snpPos-1,"\t",$snpPos,"\t",$allelicMethTab[0],"\t",$allelicMethTab[1],"\t",$allelicMethTab[2],"\t",$allelicMethTab[3],"\t",$allelicMethTab[4],"\t",$isHetCpg,"\t",$orAl1,"\t",$orAl2,"\n";
  #}
  $processedSnps++;
  if(($processedSnps % 500 == 0) && $processedSnps != 0){
    print STDERR "$processedSnps SNPs processed.\n";
  }
}
