#!/usr/bin/perl
use strict;
use warnings;

my @proteins = glob "fasta/*_maker_final.proteins_rename_rmhost.fa";
my @cdss=glob "fasta/*cds_rmhost.fa";

my %name2protein;
my %name2cds;

foreach my $file (@proteins){
	(my $isolate = $file) =~ s/_maker_final.proteins_rename_rmhost.fa//;	
	$isolate =~ s/fasta\///g;
	open R, $file;
	my @lines = <R>;
	close R;
	chomp @lines;
	my $name;
	my $seq;
	foreach my $line (@lines){
		if ($line =~ />(.*) protein.*/){
			$name = $1;
		}
		else {
			$seq =$line;
			$name2protein{$isolate}{$name}=$seq;
		}
	}
}

foreach my $file (@cdss){
        (my $isolate = $file) =~ s/_maker_only.cds_rmhost.fa//;
        $isolate =~ s/fasta\///g;
	open R, $file;
        my @lines = <R>;
        close R;
        chomp @lines;
        my $name;
        my $seq;
        foreach my $line (@lines){
                if ($line =~ />(.*)\:cds/){
                        $name = $1;
                }
                else {
                        $seq =$line;
                        $name2cds{$isolate}{$name}=$seq;
                }
        }
}


open R, "seq_need.txt";
my @names = <R>;
chomp @names;
close R;

open PRO, ">similar_protein.fasta";
open CDS, ">similar_cds.fasta";
foreach my $name (@names){
	if ($name =~ /\|/){
	my ($isolate,$gene) = split /\|/,$name;
	$gene = $gene . "-RA" if $gene !~ /-RA/;
	my $pro = $name2protein{$isolate}{$gene};
	my $cds = $name2cds{$isolate}{$gene};
	print PRO ">$name\n$pro\n";
	print CDS ">$name\n$cds\n";
	}
	else {
		
		my $isolate = substr($name,0,5);
		my $gene = $name;
		print $isolate,",$gene\n";
		$gene = $gene . "-RA" if $gene !~ /-RA/;
		print $isolate, ",$gene\n";
		my $pro = $name2protein{$isolate}{$gene};
        	my $cds = $name2cds{$isolate}{$gene};
       	 	print PRO ">$name\n$pro\n";
        	print CDS ">$name\n$cds\n";
	}
}
close PRO;
close CDS;
