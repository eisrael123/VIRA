=o
in this script, it will:

2. open the deseq2 file, select for the gene significantly induced/supressed/non-sig by BRRF1 in DG75, 
found those genes and extract seed peak from the main peak files and extract the tpm values (average tmp)
usage: perl script.pl $log2FC $padj deseq2.xlsx main_peaks.bed...
=cut


use warnings;
use strict;
use File::Basename;
use Spreadsheet::ParseXLSX;
use Text::CSV;
use Scalar::Util qw(looks_like_number);
use Spreadsheet::WriteExcel;

my ($log2FC, $padj, $deseq2, @main_peaks) = @ARGV;
my $base_dir = ".";
my $filename = basename $deseq2;
$filename =~ s/_test_vs_cntl_deseq2_results_genes.xlsx//;

my $tsv_file = "$base_dir/$filename.tsv";
my $induced_gene = "$base_dir/Induced_genes_$filename.txt";
my $suppressed_gene = "$base_dir/Suppressed_genes_$filename.txt";
my $non_sig = "$base_dir/Non-sig_genes_$filename.txt";
my $parser = Spreadsheet::ParseXLSX->new();
my $workbook = $parser->parse($deseq2) or die "Cannot open Excel file: $!";
my $tsv = Text::CSV->new({ sep_char => "\t" }) or die "Cannot create TSV object";

open(my $OUT1, ">", $tsv_file) or die "Cannot open tsv file: $!";

for my $worksheet ($workbook->worksheets()) {
    my ($row_min, $row_max) = $worksheet->row_range();
    my ($col_min, $col_max) = $worksheet->col_range();

    for my $row ($row_min .. $row_max) {
        my @row_data;
        for my $col ($col_min .. $col_max) {
            my $cell = $worksheet->get_cell($row, $col);
            my $cell_value = $cell ? $cell->value() : "";
            $cell_value = "" if $cell_value eq 'NA';
            push @row_data, $cell_value;
        }
        $tsv->print($OUT1, \@row_data);
        print $OUT1 "\n";
    }
}
close($OUT1);

print "Conversion completed......!!!\n";
print "Preparing the tsv file.....\n";

open(my $FH2, "<", $tsv_file) or die "Cannot open the tsv file: $!";
open(my $OUT2, ">", $induced_gene) or die "Cannot open the output file for induced genes: $!";
open(my $OUT3, ">", $suppressed_gene) or die "Cannot open the output file for suppressed genes: $!";
open(my $OUT6, ">", $non_sig) or die "Cannot open the output file for non-sig genes: $!";
print "Processing TSV file...\n\n";

my $log2FC_index = -1;
my $padj_index   = -1;
my @cntl_indices = ();
my @test_indices = ();

my $header_line = <$FH2>;
chomp($header_line);
my @split_header = split("\t", $header_line);

for my $i (0..$#split_header) {
    $log2FC_index = $i if $split_header[$i] eq 'log2FoldChange';
    $padj_index   = $i if $split_header[$i] eq 'padj';

    if ($split_header[$i] =~ /_cntl\d+_/i) {
        push @cntl_indices, $i;
    }
    elsif ($split_header[$i] =~ /_test\d+_/i) {
        push @test_indices, $i;
    }
}

die "Missing required columns!"
    if $log2FC_index == -1 || $padj_index == -1 || !@cntl_indices || !@test_indices;

#print $OUT2 "$header_line\n";
#print $OUT3 "$header_line\n";

while (my $line = <$FH2>) {
    chomp($line);
    my @split_line = split("\t", $line);

    my @test_values;
    for my $test_idx (@test_indices) {
        push @test_values, $split_line[$test_idx]
            if looks_like_number($split_line[$test_idx]);
    }
    next unless @test_values;

    my $test_sum = 0;
    $test_sum += $_ for @test_values;
    my $test_avg_tpm = $test_sum / scalar(@test_values);

    my @cntl_values;
    for my $idx (@cntl_indices) {
        push @cntl_values, $split_line[$idx]
            if looks_like_number($split_line[$idx]);
    }
    next unless @cntl_values;

    my $sum = 0;
    $sum += $_ for @cntl_values;
    my $cntl_avg_tpm = $sum / scalar(@cntl_values);

    next if $cntl_avg_tpm < 3;

    my $logFC_value = $split_line[$log2FC_index];
    my $padj_value  = $split_line[$padj_index];
    next if !defined $padj_value || $padj_value eq '' || $padj_value =~ /^NA$/i;

    if (looks_like_number($logFC_value) && looks_like_number($padj_value)) {
        if ($logFC_value < -$log2FC && $padj_value < $padj) {
            print $OUT3 "$split_line[0]\t$test_avg_tpm\n";
        }
        elsif ($logFC_value > $log2FC && $padj_value < $padj) {
            print $OUT2 "$split_line[0]\t$test_avg_tpm\n";
        }
        else{
            print $OUT6 "$split_line[0]\t$test_avg_tpm\n";
        }
    }
    else {
        die "logFC or padj is not numeric in line: $line\n";
    }
}

close($FH2);
close($OUT2);
close($OUT3);
close($OUT6);

print "Done filtering.\n";

foreach my $file (@main_peaks){

    open(my $FH3, "<", $induced_gene) or die "Cannot open the induced genes file for reading: $!";
    open(my $FH4, "<", $suppressed_gene) or die "Cannot open the suppressed gene file for reading: $!";
    open(my $FH6, "<", $non_sig) or die "Cannot open the non-sig gene file for reading: $!";
    open(my $FH5, "<", $file) or die "Cannot open the main peaks file for reading: $!";
    my $induced_filename = basename $induced_gene;
    my $suppressed_filename = basename $suppressed_gene;
    my $NonSig_filename = basename $non_sig;
    my $filename1 = basename $file;
    $filename1 =~ s/(.*from_TSS_).*/$1/;
    my $output1 = "$base_dir/${induced_filename}_$filename1.txt";
    my $output2 = "$base_dir/${suppressed_filename}_$filename1.txt";
    my $output3 = "$base_dir/${NonSig_filename}_$filename1.txt";
    open (my $OUT4, ">", $output1) or die "Cannot open the file: $output1";
    open (my $OUT5, ">", $output2) or die "Cannot open the file: $output2";
    open (my $OUT7, ">", $output3) or die "Cannot open the file: $output3";


my %induced;
my %suppressed;
my %non_sig;

while (my $l = <$FH3>) {
    chomp $l;
    my @f = split "\t", $l;
    $induced{$f[0]} = $f[1];   # gene => test_avg
}

while (my $l = <$FH4>) {
    chomp $l;
    my @f = split "\t", $l;
    $suppressed{$f[0]} = $f[1];  # gene => test_avg
}
while (my $l = <$FH6>) {
    chomp $l;
    my @f = split "\t", $l;
    $non_sig{$f[0]} = $f[1];  # gene => test_avg
}
    while (my $line5 = <$FH5>) {
    chomp $line5;
    my @split_line5 = split "\t", $line5;

    my $peak_info = $split_line5[3];
    my @peaks = split ";", $peak_info;
    my $seed_peak = join "\t", $peaks[1], $peaks[2];

    my $gene = $split_line5[6];

    if (exists $induced{$gene}) {
        print $OUT4 join("\t",
            $split_line5[0],     # chr
            $seed_peak,          # peak
            $split_line5[5],     # strand
            $gene,               # gene
            $induced{$gene}      # test_avg
        ), "\n";
    }
    elsif (exists $suppressed{$gene}) {
        print $OUT5 join("\t",
            $split_line5[0],
            $seed_peak,
            $split_line5[5],
            $gene,
            $suppressed{$gene}
        ), "\n";
    }
    elsif (exists $non_sig{$gene}) {
        print $OUT7 join("\t",
            $split_line5[0],
            $seed_peak,
            $split_line5[5],
            $gene,
            $non_sig{$gene}
        ), "\n";
    }
}

    

close $FH3;
close $FH4;
close $FH5;
close $OUT4;
close $OUT5;
close $OUT7;


}
