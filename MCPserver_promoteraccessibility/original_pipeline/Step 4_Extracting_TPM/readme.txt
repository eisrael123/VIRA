In this step:
- extract the genes that with log2fC < -0.5 (suppressed) or > 0.5 (induced) with padj < 0.005. The rest are non-significant, obtain the gene names and average tpm values
- then use the gene names, found those in main peak files and extract the seed peak (1 basepair) information
- once again,  I use the main peaks in BRRF1-transfected sample for the this step (it does not matter, you can use main peaks for control sample too)

run:
perl /Users/trangnguyen/Desktop/untitled\ folder/Step\ 4_Extracting_TPM/1_Scripts/3_extracting_tpm_for_atac_seq.pl 0.5 0.005 /Users/trangnguyen/Desktop/untitled\ folder/Step\ 3_Choosing_main_peaks/2_Output/Output_main_PEAKS_Matched_TSS_100_TS_50_BRRF1_PEAKS_max_width_16.fraction_max_peaks_0.2.min_single_pos_depth_10.strand_MERGED.bed.txt