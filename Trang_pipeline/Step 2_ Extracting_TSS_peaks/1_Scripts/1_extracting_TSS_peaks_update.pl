
=o
in this script, it will extract the peaks in the window of : exon 1 start (thick start) (col 7) + 50 bp from positions of annotaton file, and TSS (col 2)
minus 100 bases (for + strand) or plus 100 bases (for - strand), then look to see if there are peaks (col 2 and 3) from peak calling file that fall in those windows. if yes, print:
1. the gene name and peak information
2. if same peak matches two the same gene, even though they have different transcript ID,  count as one
3. if same peak matches two different gene, print different peaks 

# Usage:
#   perl match_peaks_to_genes.pl peaks.bed annotations.bed 100 50 
# Note: if you are in different directory, you have to drag the whole path of the script file, peak calling file and annotation file
# Output:
#   Matched_100_bases_from_TSS.bed
#   Unmatched_100_bases_from_TSS.bed
##this one removes peaks that have the same genes even though they have different transcript ID
##keep the peaks that share between two different genes
=cut
use strict;
use warnings;
use File::Basename;

my $base_dir = ".";
my ($peak_calling_file, $annotation_file, $TSS_cutoff, $TS_cutoff) = @ARGV;
my $filename = basename $peak_calling_file;

my $matched_peak   = "$base_dir/Matched_TSS_${TSS_cutoff}_TS_${TS_cutoff}_$filename";
my $unmatched_peak = "$base_dir/Unmatched_TSS_${TSS_cutoff}_TS_${TS_cutoff}_$filename";


open(my $FH1, "<", $peak_calling_file) or die "Couldn't open the CAGE peak calling file: $!";
open(my $FH2, "<", $annotation_file)   or die "Couldn't open the hg38 annotation bed file: $!";
open(my $OUT1, ">", $matched_peak)     or die "Couldn't open the matched output file: $!";
open(my $OUT2, ">", $unmatched_peak)   or die "Couldn't open the unmatched output file: $!";


my @annotations;
while (my $line2 = <$FH2>) {
    chomp $line2;
    my @split_line2 = split(/\t/, $line2);
    next unless @split_line2 > 6;
    push @annotations, \@split_line2;
}
close $FH2;


while (my $line1 = <$FH1>) {
    chomp $line1;
    my @split_line1 = split(/\t/, $line1);
    next unless @split_line1 > 4;

    my ($chr, $start_peak, $end_peak, $peak_strand) = @split_line1[0,1,2,5];
    my $matched = 0;

    # Track which genes have been printed for this peak
    my %seen_gene;

    foreach my $anno_ref (@annotations) {
        my @split_line2 = @$anno_ref;
        next unless $chr eq $split_line2[0];          
        next unless $peak_strand eq $split_line2[5];   

        # Extract gene base name (remove transcript or isoform suffix)
        my $gene_full = $split_line2[3];
        $gene_full =~ s/\s+//g;                       
        (my $gene_base = $gene_full) =~ s/[_-].*$//;  

        my ($TSS, $upstream_TSS, $thick_start);
        if ($split_line2[5] eq "+") {
            $TSS          = $split_line2[1];
            $upstream_TSS = $TSS - $TSS_cutoff;
            $thick_start  = $split_line2[6] + $TS_cutoff;

            if ($start_peak > $upstream_TSS && $end_peak < $thick_start) {
                next if $seen_gene{$gene_base}++;
                print $OUT1 join("\t", @split_line1, $gene_base), "\n";
                $matched = 1;
            }

        } elsif ($split_line2[5] eq "-") {
            $TSS          = $split_line2[2];
            $upstream_TSS = $TSS + $TSS_cutoff;
            $thick_start  = $split_line2[7] - $TS_cutoff;

            if ($start_peak > $thick_start && $end_peak < $upstream_TSS) {
                next if $seen_gene{$gene_base}++;
                print $OUT1 join("\t", @split_line1, $gene_base), "\n";
                $matched = 1;
            }
        }
    }

    # If no match was found for this peak, print to unmatched file
    if (!$matched) {
        print $OUT2 join("\t", @split_line1), "\n";
    }
}

close $FH1;
close $OUT1;
close $OUT2;

print "Done!\nMatched peaks:   $matched_peak\nUnmatched peaks: $unmatched_peak\n";
