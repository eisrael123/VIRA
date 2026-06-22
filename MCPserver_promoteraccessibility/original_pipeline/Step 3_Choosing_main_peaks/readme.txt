1. Use the output of step 2 to run this
2. In this step, I want to select for main peaks because peaks from step 2 are clusters and containing multiple peaks
- for the same isoform, if the peaks are within 100 basepairs, pick the peak with highest summation depth
- for the same gene but different isoforms, if the peak with the highest summation depth chosen above are >= 200 basepair from the next peak with the highest summation depth  of the next isoform, keep two of them. if they are < 200, keep the first one then use it to compare to the next new one









For example:
perl /Users/trangnguyen/Desktop/untitled\ folder/Step\ 3_Choosing_main_peaks/1_Scripts/2_choosing_main_peaks.pl Users/trangnguyen/Desktop/untitled\ folder/Step\ 2_\ Extracting_TSS_peaks/2_Output/Matched_TSS_100_TS_50_Cntl_PEAKS_max_width_16.fraction_max_peaks_0.2.min_single_pos_depth_10.strand_MERGED.bed /Users/trangnguyen/Desktop/untitled\ folder/Step\ 3_Choosing_main_peaks/1_Scripts/2_choosing_main_peaks.pl /Users/trangnguyen/Desktop/untitled\ folder/Step\ 2_\ Extracting_TSS_peaks/2_Output/Matched_TSS_100_TS_50_BRRF1_PEAKS_max_width_16.fraction_max_peaks_0.2.min_single_pos_depth_10.strand_MERGED.bed 100 200