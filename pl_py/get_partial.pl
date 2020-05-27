#!/usr/bin/perl
use strict;
use warnings;

my $input = "Bfra_CDS.pep.fa";
open R, "$input";
my @lines = <R>;
chomp @lines;
close R;

my $name='';
my %name2seq;
my $sign = "";
my %name2sign;
foreach my $line (@lines){
	if ($line =~ />(.*)/){
		$name = $1;
	}
	else {
		$name2seq{$name} = $line;
	#	print "$sign\n";
		$sign = "complete" if ($line =~ /\*$/ && $line =~ /^M/);
		$sign = "noM" if $line !~ /^M/;
		$sign = "noS" if $line !~ /\*$/;	
		$sign = "noSM" if ($line !~ /\*$/ && $line !~ /^M/);
	#	print "$line\n$sign\n\n"; #sleep 1;
		$name2sign{$name}=$sign;
	}
}


open W, ">Bfra.partial.list";
foreach my $gene (sort {$a cmp $b} keys %name2seq){
	my $sign = $name2sign{$gene};
	print W "$gene\t$sign\n";
}
close W;



