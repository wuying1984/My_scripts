#!/usr/bin/perl
use strict;
use warnings;

#intersect.augustus_masked.txt
#intersect.est2genome.txt
#intersect.protein2genome.txt
#intersect.repeatmasker.txt
#CDS.maker.gff
#gene.maker.gff
# part.intron.bed

my @intersects=glob "intersect*.txt";
my $cds_gff = "CDS.maker.gff";
my $gene_gff = "gene.maker.gff";
my $intron_bed = "part.intron.bed";
my $genome_gff = "part.genome.gff";

my %gene2scafold;
my %scafold2len;
my %gene2len;
my %cds2len;
my %cds2min;
my %intron2len;
my %intron2min;
my %gene2end;
my %gene2intersect;

open R, "$genome_gff";
my @lines = <R>;
close R;
chomp @lines;
foreach my $line (@lines){
	my @eles = split /\t/,$line;
	my $scafold = $eles[0];
	my $len = $eles[4];
#	print "$line\n$scafold\t$len\n";
	$scafold2len{$scafold}=$len;
}

open R, "$gene_gff";
my @lines_gene = <R>;
close R;
chomp @lines_gene;
foreach my $line (@lines_gene){
        my @eles = split /\t/,$line;
        my $scafold = $eles[0];
	my $gs = $eles[3];
	my $ge = $eles[4];
	my $desc =$eles[8];
	my @eles2=split /\;/,$desc;
	my $geneID=$eles2[0];
	$geneID=~ s/ID\=//;
	my $len = $ge - $gs +1;
	my $scafold_len = $scafold2len{$scafold};
	my $distance1 = $gs;
	my $distance2= $scafold_len - $ge;
	my $distance = $distance1;
	$distance = $distance2 if $distance2 < $distance1;
	$gene2end{$geneID}=$distance;
	$gene2len{$geneID}=$len;
	$gene2scafold{$geneID}=$scafold;
#	print "$line\n$scafold\t$geneID\t$len\t$distance1\t$distance2\t$distance\n";
#	sleep 2;
}



open R, "$cds_gff";
my @lines_cds = <R>;
close R;
chomp @lines_cds;
foreach my $line (@lines_cds){
        my @eles = split /\t/,$line;
        my $scafold = $eles[0];
        my $gs = $eles[3];
        my $ge = $eles[4];
        my $desc =$eles[8];
        my @eles2=split /\;/,$desc;
        my $geneID=$eles2[0];
        $geneID=~ s/ID\=//;
	$geneID=~ s/\-RA\:cds//;
        my $len = $ge - $gs +1;
	if (!$cds2len{$geneID}){
	        $cds2len{$geneID}=$len;
		$cds2min{$geneID}=$len;
	} else {
		$cds2min{$geneID}=$len if $len < $cds2min{$geneID};
		$cds2len{$geneID} += $len;
	}
#        print "$line\n$scafold\t$geneID\t$len----$cds2len{$geneID}\t$cds2min{$geneID}\n";
#        sleep 2;
}



open R, "$intron_bed";
my @lines_intron = <R>;
close R;
chomp @lines_intron;
foreach my $line (@lines_intron){
        my @eles = split /\t/,$line;
        my $scafold = $eles[0];
        my $gs = $eles[1];
        my $ge = $eles[2];
        my $desc =$eles[3];
        my @eles2=split /\_\_/,$desc;
        my $geneID=$eles2[0];
        $geneID=~ s/ID\=//;
        $geneID=~ s/\-RA//;
        my $len = $ge - $gs ;
        if (!$intron2len{$geneID}){
                $intron2len{$geneID}=$len;
                $intron2min{$geneID}=$len;
        } else {
                $intron2min{$geneID}=$len if $len < $intron2min{$geneID};
                $intron2len{$geneID} += $len;
        }
#       print "$line\n$scafold\t$geneID\t$len----$intron2len{$geneID}\t$intron2min{$geneID}\n";
#        sleep 2;
}



foreach my $file (@intersects){
	open R, "$file";
	my @lines_in = <R>;
	close R;
	chomp @lines_in;
	foreach my $line (@lines_in){
	        my @eles = split /\t/,$line;
	        my $intersect  = $eles[19];
        	my $desc =$eles[9];
        	my @eles2=split /\;/,$desc;
       		my $geneID=$eles2[0];
        	$geneID=~ s/ID\=//;
        	$geneID=~ s/\-RA\:cds//;
        	my $symbol = $file;
		$symbol =~ s/intersect\.//;
		$symbol =~ s/\.txt//;
		if (!$gene2intersect{$geneID}{$symbol}){
			$gene2intersect{$geneID}{$symbol} = $intersect;
		} else {
			$gene2intersect{$geneID}{$symbol} += $intersect;
		}
	#	print "$line\n$geneID\t$symbol----$intersect---$gene2intersect{$geneID}{$symbol}\n";
       #		 sleep 2;
	}
}

my @symbols = qw /augustus_masked est2genome  protein2genome repeatmasker/;
open W, ">gene_list_for_checkup.tab";
print W "Gene\tlen\tscafold\tdistance\tcds_len\tintron_len\tcds_min\tintron_min\taugustus_masked\test2genome\tprotein2genome\trepeatmasker\n";
foreach my $geneID (sort {$a cmp $b} keys %gene2len){
#	print "$geneID\n";
	my $gene_len = $gene2len{$geneID};
	my $scafold = $gene2scafold{$geneID};
	my $distance = $gene2end{$geneID};
	my $cds_len = $cds2len{$geneID};
	my $cds_min = $cds2min{$geneID};
	my $intron_len = 0;
	$intron_len = $intron2len{$geneID} if $intron2len{$geneID};
	my $intron_min = 0;
	$intron_min =  $intron2min{$geneID} if $intron2min{$geneID};
	print W "$geneID\t$gene_len\t$scafold\t$distance\t$cds_len\t$intron_len\t$cds_min\t$intron_min";
#	print "$geneID\t$scafold\t$distance\t$cds_len\t$intron_len\t$cds_min\t$intron_min";
	foreach my $symbol (@symbols){
		my $intersect = 0;
		$intersect =  $gene2intersect{$geneID}{$symbol} if $gene2intersect{$geneID}{$symbol};
		my $percentage = $intersect / $cds_len;
		print W "\t$percentage";
#		print "\t$symbol---$percentage";
	}
	print W "\n";
#	print "\n"; sleep 2;	
}
close W;

