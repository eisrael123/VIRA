=o
in this script, it will:
1. from the files including peaks information corresponding to the gene names extracted from annotation hg38 file, 
for each isoform: 
+ if muitple peaks are around 100 basepair, it will pick the highst peak(s) with highest summation depth (column 5)
for the same gene (if two or more peaks have the same highest depth from the same cluster (within 100bp), keep all of them)
+ if the highest peaks of different isoforms are at least 200 basepair apart for the same gene, keep them 
    *use the FIRST peak position of the first isoform as a standard to compare >= 200 basepair of other isoforms
for example:
cluster1 [100,100]
cluster2 [400]
cluster3 [500]
use position 100 of the first peak on cluster1 to compare with peak 400 and peak 500 from cluster 2 and 3
#usage: perl script.pl @files $cluster_distance $distinct_distance
#for example: perl choosing_main_peaks.pl file1.bed file2.bed 100 200

=cut

use strict;
use warnings;
use File::Basename;

my ($cluster_distance, $distinct_distance) = splice(@ARGV, -2);
my @files = @ARGV;

foreach my $file (@files) {

    my $filename = basename($file);
    my $output_file = "Output_main_PEAKS_$filename.txt";

    open(my $FH, "<", $file) or die "Couldn't open $file: $!";
    open(my $OUT, ">", $output_file) or die "Couldn't open $output_file: $!";

    my %genes;

    # ---- read input ----
    while (my $line = <$FH>) {
        chomp $line;
        my @f = split /\t/, $line;

        my @peak_info = split /;/, $f[3];
        my $peak_pos  = $peak_info[1];

        my ($chr, $start, $end, $depth, $strand, $gene) =
           ($f[0], $f[1], $f[2], $f[4], $f[5], $f[6]);

        push @{ $genes{$gene} }, {
            line     => $line,
            peak_pos => $peak_pos,
            depth    => $depth,
            chr      => $chr,
            strand   => $strand,
        };
    }
    close $FH;

    my @final_peaks;

   
    foreach my $gene (keys %genes) {

        my @peaks =
            sort { $a->{peak_pos} <=> $b->{peak_pos} }
            @{ $genes{$gene} };

        my @clusters;
        my @current_cluster = ($peaks[0]);

        for (my $i = 1; $i < @peaks; $i++) {
            if (abs($peaks[$i]{peak_pos}
                  - $current_cluster[-1]{peak_pos}) <= $cluster_distance) {

                push @current_cluster, $peaks[$i];
            }
            else {
                push @clusters, [@current_cluster];
                @current_cluster = ($peaks[$i]);
            }
        }
        push @clusters, [@current_cluster];

        my @kept;
        my $last_cluster_rep_pos;

        foreach my $cluster (@clusters) {

            # find max depth in this cluster
            my $max_depth;
            foreach my $p (@$cluster) {
                if (!defined $max_depth || $p->{depth} > $max_depth) {
                    $max_depth = $p->{depth};
                }
            }

            # collect all strongest peaks in this cluster
            my @cluster_strongest =
                grep { $_->{depth} == $max_depth } @$cluster;

            my $cluster_rep_pos = $cluster_strongest[0]{peak_pos};

            if (!defined $last_cluster_rep_pos
                || abs($cluster_rep_pos - $last_cluster_rep_pos)
                   >= $distinct_distance) {

                push @kept, @cluster_strongest;
                $last_cluster_rep_pos = $cluster_rep_pos;
            }
        }

        push @final_peaks, map { $_->{line} } @kept;
    }
    @final_peaks = sort {
    my @a = split /\t/, $a;
    my @b = split /\t/, $b;

    $a[0] cmp $b[0]        # column 1: chr
        ||
    $a[1] <=> $b[1]       # column 2: start
} @final_peaks;


    print $OUT "$_\n" for @final_peaks;
    close $OUT;

    print "✔ Written $output_file\n";
}


