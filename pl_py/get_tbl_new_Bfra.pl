#!/usr/bin/perl
use warnings;
use strict;
###############################################
###########GET low quality genes###############

my %gene2remove;
open R, "remove.list"; #genes with short intron <10nt
my @lines_rm = <R>;
close R;
chomp @lines_rm;

foreach my $line (@lines_rm){
	$line =~ s/Bfra/BFRA/;  
	$gene2remove{$line} ++;
	#print "$line\n"; sleep 2;
}


###############################################
###########GET LIST OF PARTIAL GENES##############
my %name2partial;
my @part_files = glob "*.partial.list";
foreach my $part_file (@part_files){
        open R, $part_file;
        my @lines = <R>;
        close R;
        chomp @lines;
        #shift @lines;
        foreach my $line (@lines){
                my ($name,$sign) = split /\t/,$line;
                $name2partial{$name}=$sign;
        }
}

###############################################
###########GET INFO FROM ANNO FILE###############
#my %name2complete;
my %name2desc;
my %name2se;
my %name2sp;
my @anno_tab = glob "*eggnog.annotations.short";
my @gff_file = glob "*all_no_seq.modified.sorted.gff";
print join ";", @anno_tab, "\n";
print join ";", @gff_file, "\n";
#foreach my $file (@anno_tab) {
#        print $file, "\n";
#        open R, $file;
#        my @lines = <R>;
#        close R;
#        chomp @lines;
#
#	print scalar @lines, "\n";
#        #shift @lines;
# 
#        foreach my $line (@lines){
#		if ($line !~ /\#/){
#			my @eles = split /\t/,$line;
#			my $name = $eles[0];
#			my $desc = $eles[21];
#			
#		}
#        }
#}

###############################################
###########GET DESC FROM ANNO FILE##############
foreach my $file2 (@anno_tab){
        open R, "$file2";
        my @lines = <R>;
        chomp @lines;
#        shift @lines;
        foreach my $line (@lines){
		if ($line !~ /\#/){
	                my @eles = split (/\t/,$line,-1);
        	        my $name = $eles[0];
#               		 print $name, "\n"; sleep 1;
                	my $desc = $eles[2];
#			print $line, "\n", $desc, "\n\n";
	                $desc = "hypothetical protein" if $desc eq '';
			#print $line, "\n", $desc, "\n\n";
	       #         $desc =~ s/\(.*\)// if $desc =~ /\(.*\)/;
	                $desc =~ s/ partial//g if $desc =~ / partial/;
        	        $desc =~ s/partial//g if $desc =~ /partial/;
	                $desc = "Function unknown protein" if $desc eq "DUF";
#	                $desc = "hypothetical protein" if $desc eq "family";
 	               $desc =~ s/family family/family/ if $desc =~ /family family/;
        	        $desc = "signal peptide protein" if $desc eq "signal peptide";
 	               $desc =~ s/uncharacterized/putative/ if $desc =~ /uncharacterized/;
        	        $desc =~ s/zinc finger/zinc finger protein/ if $desc =~ /zinc finger$/;
			$desc =~ s/repeats/protein/ if $desc =~ /repeats\.$/;
			$desc =~ s/domain\./domain containing protein/ if $desc =~ /domain\.$/;
			$desc =~ s/domain/domain containing protein/ if $desc =~ /domain$/;
			$desc =~ s/domains\./domain containing protein/ if $desc =~ /domains\.$/;
			$desc = "PhoX homologous protein" if $desc =~ /p47phox/;
			$desc =~ s/motif\./motif containing protein\./ if $desc =~ /motif\.$/;
			$desc =~ s/motifs\./motif containing protein\./ if $desc =~ /motifs\.$/;
			$desc ="Zinc-binding domain containing protein" if $desc =~ /Zinc-binding domain/;
			$desc =~ s/homologues\./homologues protein\./ if $desc =~ /homologues\.$/;
			$desc =~ s/family\./family protein\./ if $desc =~ /family\.$/;	
			if ($desc =~ /domain in/ || $desc =~ /Domain present in/ || $desc =~ /Domain involved in/ || $desc =~ /Domain in/){
				$desc =~ s/domain.*in//;
				$desc =~ s/Domain.*in//;
				$desc =~ s/Domain.*in//;
				$desc = $desc . "-like protein";
			}	
			$desc = "zinc-finger like protein" if $desc =~ /C4HC3/;
			$desc = "SCP family protein" if $desc =~ /SCP/ && $desc =~ /Ag5/;
			$desc = "hypothetical protein" if $desc =~/in LPS-induced tumor/;
			$desc =~ s/enzyme\./enzyme/ if $desc =~ /enzyme\.$/;
			$desc = "Putative catecholamine-binding domain containing protein" if $desc =~ /catecholamine\-binding/;	
			$desc =~ s/ with conserved.*// if $desc =~ /with conserved/;
			$desc = "cupin metalloenzyme superfamily protein" if $desc =~ /A domain family that is part/;
			$desc =~ s/containing/containing protein/ if $desc =~ /containing$/;
			$desc = "hypothetical protein" if $desc =~ /Four repeated domains in the Fasciclin/;
#	                $desc = "hypothetical protein" if $desc eq "amino acid";
#      	         	$desc = "hypothetical protein" if $desc eq "like superfamily";
#             	 	$desc = "hypothetical protein" if $desc eq "membrane";
#                	$desc =  "hypothetical protein" if $desc =~ /^related to/;
#                	$desc =  "hypothetical protein" if $desc =~ /unnamed protein/;
#                	$desc =  "hypothetical protein" if $desc eq "precursor";
#                	$desc = "hypothetical protein" if $desc eq "like helicase 2";
#                	$desc = "hypothetical protein" if $desc eq "like methyltransferase";
#                	$desc = "hypothetical protein" if $desc eq "like mitochondrial";
#                	$desc =~ s/probable/putative/ if $desc =~ /probable/;
#                	$desc = "hypothetical protein" if $desc eq "DNA" ;
#                	$desc = "hypothetical protein" if $desc eq "like subfamily B member 12";
#                	$desc = "hypothetical protein" if $desc eq "predicted protein";
#                	$desc =~ s/predicted/putative/ if $desc =~ /^predicted/;
                	$desc = $desc . "protein" if $desc =~ /\-$/;
                	$desc =~ s/\-// if $desc =~ /^\-/;
                	$desc = "hypothetical protein" if $desc =~ /hypothetical protein.*/;
#                	$desc =~ s/\(.*\)// if $desc =~ /\(.*\)/;
                	$desc =~ s/partial//g if $desc =~ /partial/;
			$desc =~ s/\(// if $desc =~ /^\(/;	
                	$desc =~ s/domain/domain containing protein/ if $desc =~ /domain$/;
			$desc = "zinc finger protein" if $desc =~ /zinc finger binding/;
			$desc =~ s/activity/activity protein/ if $desc =~ /activity$/;
			print $desc, "\n" if $desc =~ /\[.*\]/;
                	if ($desc =~ /\[.*\]/) {
                		$desc =~ s/\[.*\]//g;
                		print $desc, "\n";
	        	}
			$desc =~ s/\, $// if $desc =~ /\, $/;
                	my $sign = $name2partial{$name};
			$desc = $desc . " partial sequence" if $sign =~ /no/;
                	$desc =~ s/\.$// if $desc =~ /\.$/;
			$name2desc{$name}= $desc;
		
			}

        }
}

###############################################
###########GET FEATURE FROM GFF################
my %name2feature;
my %name2strand;
my %contig2gene;
foreach my $file3 (@gff_file){
        print "$file3\n"; sleep 2;
        my $isolate = 'Bfra'; #substr ($file3,0,5);
        open R, "$file3";
        my @lines = <R>;
        chomp @lines;
        foreach my $line (@lines){
                my ($contig,undef,$symbol,$s,$e,undef,$strand,undef,$left)=split /\t/,$line;
                if (!$s){
                        print $line, "\n"; sleep 1;
                }
                $s = 1 if $s ==0;
                $e = 1 if $e ==0;
                my @eles = split /\;/,$left;
                my $id = $eles[0];
                $id =~ s/\:.*//g;
                $id =~ s/ID\=//g;
                $id =~ s/\-RA//g;
                my $name = $id; #$isolate . "|" . $id;
		if (!$gene2remove{$name}){
		  #               print $id, "\n"; sleep 1;
		  $contig2gene{$isolate}{$contig}{$name} ++ ;#if $name2partial{$name};
		  #               print $name, "\n";
		  #               sleep 1;
		  $name2strand{$name} = $strand;
		  my $se = $s . "\t" . $e;
		  $se = $e . "\t" .$s if $strand eq "-";
		  if ($symbol eq "gene" || $symbol eq "mRNA" || $symbol eq "exon" || $symbol eq "CDS"){
		    $name2feature{$id}{$symbol}{$s} = $e if $s <= $e;
		    $name2feature{$id}{$symbol}{$e} = $s if $e < $s;
		    #                       print "$symbol, $s, $e\n" if $id =~ /UMSG2_06368/;
		  }
		}
        }
}


###############################################
###########CHANGE TAG AND PRINT TBL############
my %locus_tag;
$locus_tag{Bfra}="Bfra";
#$locus_tag{UMSG1}="GcM1";
#$locus_tag{UMSG2}="OnM2";
#$locus_tag{UMSG3}="GcM3";

my @syms = qw /gene mRNA exon CDS/;
my @strains = qw /Bfra/;
foreach my $isolate (@strains){
my $tag=$locus_tag{$isolate};
my $new_file = $isolate . ".tbl";
open W, ">$new_file";
foreach my $contig (sort { $a cmp $b} keys %{$contig2gene{$isolate}}){
        my $contig_new = $contig;
        #$contig_new =~ s/tig/Bfra/ if $contig =~ /tig/;
#        $contig_new =~ s/Gc_UMSG1/GcM1_contig/ if $contig =~ /UMSG1/;
#        $contig_new =~ s/Gc_UMSG3/GcM3/ if $contig =~ /UMSG3/;
#        $contig_new =~ s/Ol_UMSG2/OnM2/ if $contig =~ /UMSG2/;
        print W ">Feature $contig_new\n";

        foreach my $name (sort {$a cmp $b} keys %{$contig2gene{$isolate}{$contig}}){
               # my ($isolate,$id)=split /\|/,$name;
                my $id = $name;
		my $prefix = $locus_tag{$isolate};
                my $new_id = $id; $new_id =~ s/BFRA/$prefix/;
                my @eles = split /\_/,$new_id;
                if ((scalar @eles) > 2){
                        #print $new_id, "\n";
                        $new_id = $eles[0] . "_c". $eles[1] . "o". $eles[2];
                        #print $new_id , "\n";
                }
                my $comp = $name2partial{$name};
                my $desc = $name2desc{$name};
                $desc = "hypothetical protein" if !$name2desc{$name};
                $desc = $desc . "protein" if $desc =~ /\-$/;
                my @genese = sort {$a <=> $b} keys %{$name2feature{$id}{gene}};
                my @mrnase = sort {$a <=> $b} keys %{$name2feature{$id}{mRNA}};
                my @exonse = sort {$a <=> $b} keys %{$name2feature{$id}{exon}};
                my @cdsse = sort {$a <=>$b} keys %{$name2feature{$id}{CDS}};
                my $strand = $name2strand{$name};
                #print  "$isolate $id $strand\n";sleep 1;
		############################################################################################
		#######################PRINT GENE INFOR#####################################################
                my $genesfirst = $genese[0]; my $gs = $genesfirst;
                my $genesfirstend = $name2feature{$id}{gene}{$genesfirst}; my $ge = $genesfirstend;
                my $genesefirst = $genesfirst . "\t" . $genesfirstend;
                $genesefirst = $genesfirstend. "\t" . $genesfirst if $strand eq "-";
                if ($comp eq "complete"){
		  print W $genesefirst, "\tgene\n";
		#  print W "\t\t\tgene\t$desc\n";
		  print W "\t\t\tlocus_tag\t$new_id\n";
	        }
		elsif ($comp eq "noM"){
		   print W "<",$genesefirst, "\tgene\n";
		 #  print W "\t\t\tgene\t$desc\n";
		   print W "\t\t\tlocus_tag\t$new_id\n";
		 }
		elsif ($comp eq "noS"){
		  my ($one,$two) = split /\t/,$genesefirst;
		  $genesefirst = "$one\t".">" .$two;
		   print W $genesefirst, "\tgene\n";
		  # print W "\t\t\tgene\t$desc\n";
		   print W "\t\t\tlocus_tag\t$new_id\n";
		 }
		elsif ($comp eq "noSM"){
		  my ($one,$two) = split /\t/,$genesefirst;
		  $genesefirst = "$one\t".">" .$two;
		   print W "<",$genesefirst, "\tgene\n";
		   #print W "\t\t\tgene\t$desc\n";
		   print W "\t\t\tlocus_tag\t$new_id\n";
		 }
		############################################################################################
		#######################PRINT EXON INFOR#####################################################
		my $exonsfirst = $exonse[0];
                $exonsfirst = $exonse[-1] if $strand eq "-";
                my $exonsfirstend = $name2feature{$id}{exon}{$exonsfirst};
                #print W $exonfirst, "\n" if $id =~ /UMSG2_06368/;
                my $exonsefirst = $exonsfirst . "\t" . $exonsfirstend;
                $exonsefirst = $exonsfirstend. "\t" .  $exonsfirst if $strand eq "-";
		if  ((scalar @exonse) ==1 ){
		  my ($one,$two) = split /\t/,$exonsefirst;
		  my $new_exonsefirst = "$one\t".">" .$two;
		  print W $exonsefirst, "\tmRNA\n"  if ($comp eq "complete");
		  print W "<",$exonsefirst, "\tmRNA\n"  if ($comp eq "noM");
		  print W $new_exonsefirst, "\tmRNA\n"  if ($comp eq "noS");
		  print W "<",$new_exonsefirst, "\tmRNA\n"  if ($comp eq "noSM");
		}
                elsif ((scalar @exonse) >1 ){
		  print W $exonsefirst, "\tmRNA\n"  if ($comp eq "complete");
		   print W $exonsefirst, "\tmRNA\n"  if ($comp eq "noS");
		   print W "<",$exonsefirst, "\tmRNA\n"  if ($comp eq "noM");
		   print W "<",$exonsefirst, "\tmRNA\n"  if ($comp eq "noSM");
                        if ($strand eq "+"){
                                foreach my $i (1 .. ((scalar @exonse)-1)){
                                        my $s = $exonse[$i];
                                        my $e = $name2feature{$id}{exon}{$s};
                                        my $se = $s . "\t" . $e;
					if ($i < ((scalar @exonse)-1)){
					  print W "$se\n";
					}
					elsif ($i == ((scalar @exonse)-1)){
					   my ($one,$two) = split /\t/,$se;
					   my $new_se = "$one\t".">" .$two;
					  print W "$se\n" if  ($comp eq "complete");
					  print W "$se\n" if  ($comp eq "noM");
					  print W "$new_se\n" if  ($comp eq "noS");
					  print W "$new_se\n" if  ($comp eq "noSM");
					}
                                }
                        }
                        else {
                                @exonse = sort {$b <=> $a} @exonse;
                                foreach my $i (1 .. ((scalar @exonse)-1) ){
                                        my $s = $exonse[$i];
                                        my $e = $name2feature{$id}{exon}{$s};         
                                        my $se =$e . "\t" . $s;
                                       	if ($i < ((scalar @exonse)-1)){
					  print W "$se\n";
					}
					elsif ($i == ((scalar @exonse)-1)){
					   my ($one,$two) = split /\t/,$se;
					   my $new_se = "$one\t".">" .$two;
					   print W "$se\n" if  ($comp eq "complete");
					  print W "$se\n" if  ($comp eq "noM");
					  print W "$new_se\n" if  ($comp eq "noS");
					  print W "$new_se\n" if  ($comp eq "noSM");
					}
                                }
                        }
		      }
		print W "\t\t\tproduct\t$desc\n";
                print W "\t\t\tprotein_id\tgnl\|HuUMD\|$new_id\n";
                print W "\t\t\ttranscript_id\tgnl\|HuUMD\|$new_id","T\n";
                print W "\t\t\tnote\tpartial gene\n" if ($comp ne "complete"); ;
        #       if ($name2shortintron{$new_id}){
        #               print W "\t\t\tnote\tnonfunctional due to frameshift\n";
        #       }
		#       else {
		############################################################################################
		#######################PRINT EXON INFOR#####################################################
                my $cdssfirst = $cdsse[0];
                $cdssfirst = $cdsse[-1] if $strand eq "-";
               # print $cdssfirst, "\n" if $id =~ /UMSG2_06368/;
                my $cdssefirst = $cdssfirst . "\t" . $name2feature{$id}{CDS}{$cdssfirst};
                $cdssefirst = $name2feature{$id}{CDS}{$cdssfirst} . "\t" .  $cdssfirst if $strand eq "-";
		if ((scalar @cdsse)==1){
		   my ($one,$two) = split /\t/,$cdssefirst;
		  my $new_cdssefirst = "$one\t".">" .$two;
		  print W $cdssefirst, "\tCDS\n"  if ($comp eq "complete");
		  print W "<",$cdssefirst, "\tCDS\n"  if ($comp eq "noM");
		  print W $new_cdssefirst, "\tCDS\n"  if ($comp eq "noS");
		  print W "<",$new_cdssefirst, "\tCDS\n"  if ($comp eq "noSM");
		}
                if ((scalar @cdsse)>1){
		  print W $cdssefirst, "\tCDS\n"  if ($comp eq "complete");
		  print W $cdssefirst, "\tCDS\n"  if ($comp eq "noS");
		  print W "<",$cdssefirst, "\tCDS\n"  if ($comp eq "noM");
		  print W "<",$cdssefirst, "\tCDS\n"  if ($comp eq "noSM");
                        if ($strand eq "+"){
                                foreach my $i (1 .. ((scalar @cdsse)-1)){
                                        my $s = $cdsse[$i];
                                        my $se = $s . "\t" . $name2feature{$id}{CDS}{$s};
                                    	if ($i < ((scalar @cdsse)-1)){
					  print W "$se\n";
					}
					elsif ($i == ((scalar @cdsse)-1)){
					  my ($one,$two) = split /\t/,$se;
					   my $new_se = "$one\t".">" .$two;
					  print W "$se\n" if  ($comp eq "complete");
					   print W "$se\n" if  ($comp eq "noM");
					  print W "$new_se\n" if  ($comp eq "noS");
					  print W "$new_se\n" if  ($comp eq "noSM");
					}
                               }
                         }
                        else {
                                @cdsse =sort {$b <=> $a} @cdsse;
                                foreach my $i (1 .. ((scalar @cdsse)-1)){
                                        my $s = $cdsse[$i];
                                        my $se = $name2feature{$id}{CDS}{$s} . "\t" . $s;
                                        if ($i < ((scalar @cdsse)-1)){
					  print W "$se\n";
					}
					elsif ($i == ((scalar @cdsse)-1)){
					  my ($one,$two) = split /\t/,$se;
					   my $new_se = "$one\t".">" .$two;
					  print W "$se\n" if  ($comp eq "complete");
					  print W "$se\n" if  ($comp eq "noM");
					  print W "$new_se\n" if  ($comp eq "noS");
					  print W "$new_se\n" if  ($comp eq "noSM");
					}
                                }
                        }
                }
                print W "\t\t\tproduct\t$desc\n";
		print W "\t\t\tcodon_start\t1\n" if ($comp eq "noM" || $comp eq "noSM") ;
                print W "\t\t\tprotein_id\tgnl\|HuUMD\|$new_id\n";
                print W "\t\t\ttranscript_id\tgnl\|HuUMD\|$new_id","T\n";
                print W "\t\t\tnote\tpartial gene\n" if $comp ne "complete";
        }
}
        close W;
        print $new_file, "\n";

}
