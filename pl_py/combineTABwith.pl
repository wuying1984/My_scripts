#!/usr/bin/perl
my $tab1 = shift;
my $col1 = shift;
my $tab2 = shift;
my $col2 = shift;

my %HASH = ();
### TAB2 is the supplementary table, and we get hash from it
my @lines = get_lines("$tab2");
foreach my $tl (@lines) {
    chomp $tl;
    my @temp = split "\t", $tl;
    $HASH{$temp[$col2-1]}=\@temp;
}


### TAB1 is the main table, that means we add new stuff at the end of each line 
@lines = ();
@lines = get_lines("$tab1");
#open W, ">$tab1.added";
foreach my $tl (@lines) {
    chomp $tl;
    my @temp = split "\t", $tl;
    my $key_element = $temp[$col1 - 1];
    if (my $sup = $HASH{$key_element}) {
        print $tl, "\t", join "\t", @$sup;
	print "\n";
    }
}
#close W;


sub get_lines {
    my $filename = shift;
    open R, "$filename";
    my @lines = <R>;
    close R;
    return @lines;
}

