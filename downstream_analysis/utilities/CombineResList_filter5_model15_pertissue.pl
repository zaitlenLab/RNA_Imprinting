##combine imprinting results by chr and tissue

#arguments: original_file_name tissue_list outfile_prefix

use strict;
use warnings;
use Scalar::Util qw(looks_like_number);

open(TIS, "$ARGV[1]") or die "can't open tissue list\n";
my @tis = <TIS>;
chomp @tis;

my $h;
my $data={};
foreach my $t (@tis) {
    for(my $chr=1;$chr<=22;$chr++) {
#    for(my $chr=21;$chr<=22;$chr++) {
        open(IMPR, join("", "impglr/", "$ARGV[0]", ".", $t, ".", $chr, ".filter5.model15.impglr.txt")) or die "can't open impr res\n";
        print STDERR "processing ", join(".", "$ARGV[0]", $t, $chr), "\n";
        $h = <IMPR>;
        while(<IMPR>) {
            chomp;
            my @ts = split(/[\s]+/);
			if (looks_like_number($ts[0])) {
				splice @ts, 0, 0, '';
			}
			if (looks_like_number($ts[3])) {
				splice @ts, 3, 0, 'NA';
			}
#			print join("\t", @ts[12..19]), "\n";
			if (looks_like_number($ts[18])) {
				splice @ts, 18, 0, 'NA';
				splice @ts, 19, 0, 'NA';
			}
	        for (my $i=3;$i<(scalar(@ts)+3);$i++) {						
		        $data->{$t}->{$ts[2]}->{$i} = $ts[$i];
            }
        }
        close IMPR;
    }
}

chomp $h;
my @header = split(/[\s]+/, $h);
$header[0] = "tissue";
$header[1] = "gene";

foreach my $t (keys(%{$data})) { #for tissues
	my $outfile = join("", ">", $ARGV[2], ".", $t);
	open(OUT, $outfile) or die "can't open outfile\n";
	print OUT join("\t", @header), "\n";

	foreach my $g (keys(%{$data->{$t}})) { #for genes
#    	my $eqtls = join(",", keys(%{$data->{$g}->{"EQTL"}}));
#    	if ($eqtls eq "0.000") {
#        	$eqtls = "-";
#	    }
#    	print OUT join("\t", $g, $eqtls, $data->{$g}->{"KNOWN"}); 
    	print OUT join("\t", $t, $g); 
	    for (my $i=3;$i<scalar(keys(%{$data->{$t}->{$g}}));$i++) {
#			print $i;
	        print OUT "\t" . $data->{$t}->{$g}->{$i};
       	 }
	   	print OUT "\n";   	    
	}
	close OUT;
}

        
        





