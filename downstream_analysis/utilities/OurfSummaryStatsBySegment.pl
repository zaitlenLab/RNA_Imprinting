##summarize methylation value data per annotation
#arguments: ourf_Y_file bed_intersect_file >output_file

#gencodeV10_geneBodies_2kbEndGene.gencodeV10_geneBodies_2kbEndGene.bed
#chr16	53468111	53468112	cg00000029	.	-1	-1	.	-1	.	.	0
#chr3	37459205	37459206	cg00000108	chr3	37429760	37476988	ENSG00000198590.7	2	+	protein_coding	1
#chr3	171916036	171916037	cg00000109	chr3	171759418	172119455	ENSG00000075420.8	2	+	protein_coding	1

use strict;
use lib '/data/NYGC/software/perl/perls/perl-5.18.0//lib/site_perl/5.18.0/';
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Statistics::Basic qw(:all);

    
open(METH, "$ARGV[0]") or die "can't open methylation file\n";
open(BEDINT, "$ARGV[1]") or die "can't open bed file\n";


my $bed = {};
while(<BEDINT>) {
    chomp;
    my @ts = split(/\t/);
    if ($ts[7] eq ".") {next}
    $bed->{$ts[3]} = $ts[7];
#    print STDERR $ts[7], "\n";
#    push(@{$bed->{ts[8]}}, $ts[4]);
}
close BEDINT;

my $h = <METH>;
chomp $h;
my @header = split(/\t/, $h);

my $data = {};
while(<METH>) {
    chomp;
    my @ts = split(/\t/);        
    if (exists($bed->{$ts[0]})) {
        for(my $i=4;$i<scalar(@ts);$i++) {
            push(@{$data->{$bed->{$ts[0]}}->{$header[$i]}}, $ts[$i]);
        }
    }
}
close METH;

print join("\t", "GENE", @header[4..(scalar(@header)-1)]), "\n";
foreach my $gene (keys(%{$data})) {
    print $gene;
    foreach my $sample (@header[4..(scalar(@header)-1)]) {
#	print STDERR join("\t", ( @{$data->{$gene}->{$sample}})), "\n";
        my $mean = mean( @{$data->{$gene}->{$sample}});
        my $median = median( @{$data->{$gene}->{$sample}});
        my $var = variance( @{$data->{$gene}->{$sample}});
        my $totcount = scalar(@{$data->{$gene}->{$sample}});
        my $semicount = 0;
        foreach my $m (@{$data->{$gene}->{$sample}}) {
            if ($m<=0.7 && $m>=0.3) {
                $semicount++;
            }
        }
        print "\t" . join(";", $totcount, $semicount, $mean, $median, $var); 
    }
    print "\n";
}


    
    
    



