#!/usr/bin/perl -ws
# Usage:
# perl get_intron_for_LaSSO.pl annotation/gencode.v28.annotation.gtf

my $gtfFile = shift;
my @lines = get_elements($gtfFile);
my %HASH_exon = ();
my %geneid2chr = ();
my %geneid2strand = ();
foreach my $tl (@lines) {
    unless ($tl =~ /^#/){
        my ($chr, $source, $type, $s, $e, $dot1, $strand, $dot2, $other) = split "\t", $tl;
        if ($type eq "exon") {
            my @temp = split "gene_id \"", $other;
            my @temp2 = split "\"", $temp[1];
            my $geneid = $temp2[0];

            @temp = split "transcript_id \"", $other;
            @temp2 = split "\"", $temp[1];
            my $transid = $temp2[0];

            $geneid2chr{$geneid} = $chr;
            $geneid2strand{$geneid} = $strand;
            $HASH_exon{$geneid}{$transid}{$s} = $s . "." . $e;
        }
    }
}


my %already_printed = ();
foreach my $geneid (keys %HASH_exon) {
    my ($max, $min) = (0,1000000000000000);
    my $chr = $geneid2chr{$geneid};
    my $strand = $geneid2strand{$geneid};
    my $strand_digit = 0;
    if ($strand eq "+"){
        $strand_digit = 1;
    } else {
        $strand_digit = -1;
    }
    foreach my $transid (keys %{$HASH_exon{$geneid}}) {
        my $last_exonid = "NA";
        if ($strand eq "+"){
            my $last_end = -1;
            foreach my $exonid (sort {$HASH_exon{$geneid}{$transid}{$a} <=> $HASH_exon{$geneid}{$transid}{$b}} keys %{$HASH_exon{$geneid}{$transid}}) {
                my $se = $HASH_exon{$geneid}{$transid}{$exonid};
                my ($s, $e) = split /\./, $se;
                if ($last_exonid ne "NA"){
                    my $id = $chr . ";" . $strand_digit . ";" . $last_end . ";" . ($s - 1);
                    print $chr, "\t", $last_end, "\t", $s - 1, "\t", $transid, "__", $chr, "__", $last_end, "__", $s - 1, "__", $strand, "\t1\t", $strand, "\n";
                    $already_printed{$id} = 1;
                }
                $last_exonid = $exonid;
                $last_end = $e
            }
        }
        else {
            my $last_start = -1;
            foreach my $exonid (sort {$HASH_exon{$geneid}{$transid}{$b} <=> $HASH_exon{$geneid}{$transid}{$a}} keys %{$HASH_exon{$geneid}{$transid}}) {
                my $se = $HASH_exon{$geneid}{$transid}{$exonid};
                my ($s, $e) = split /\./, $se;
                if ($last_exonid ne "NA"){
                    my $id = $chr . ";" . $strand_digit . ";" . $e . ";" . ($last_start - 1);
                    print $chr, "\t", $e, "\t", $last_start - 1, "\t", $transid, "__", $chr, "__", $e, "__", $last_start - 1, "__", $strand, "\t1\t", $strand, "\n";
                        $already_printed{$id} = 1;
               }
                $last_exonid = $exonid;
                $last_start = $s
            }
        }
        #print "\n\n\n";
    }
}


sub get_elements {
    my $filename = shift;
    open R, "$filename" or die "cannot open $filename $!";
    my @lines = <R>;
    close R;
    chomp $_ foreach @lines;
    return @lines;
}
