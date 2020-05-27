#!/usr/bin/perl
use strict;
use warnings;

my $CDS_file = "Bfra_individual_CDS.fa";
my $out_CDS_fasta="Bfra_CDS.fa";
open R, "$CDS_file";
my @lines = <R>;
close R;
chomp @lines;

my %name2seq;
my $name;
foreach my $line (@lines){
	if ($line =~ />(.*)-RA.*/){
		$name= $1;
	}
	else {
		if (!$name2seq{$name}){
			$name2seq{$name}=$line;
			#print $name,"\n";
		} else {
			$name2seq{$name} .= $line;
			#print $name2seq{$name}, "\n";
			#sleep 1;
		}
	}
}


open W, ">$out_CDS_fasta";
foreach my $gene (sort {$a cmp $b} keys %name2seq){
	print W ">$gene\n";
	print W $name2seq{$gene}, "\n";
}
close W;
