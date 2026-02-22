# Overview of inferred loci from SPICE

This directory contains the reference loci set inferred by SPICE for the associated publication.
If you are interested in the top loci, either sort the loci by p-value or use the manually curated list of top 100 loci.


## Files
* `all_460_loci`: Complete loci set inferred by SPICE from 5,966 pan-cancer TCGA samples. This loci set was used for the analysis in the associated publication.
* `manual_curation_top_100_loci`: Manually curated list of the top 100 loci based on how pronounced the triangular shape is across all length scales.

## Columns
* `chrom`: Chromosome
* `pos`: Genomic position
* `type`: Locus class where `OG` indicates oncogene-like gain selection and `TSG` indicates tumor suppressor-like loss selection.
* `start_ci`: Start of the locus confidence interval
* `end_ci`: End of the locus confidence interval
* `nr_events_in_data`: Number of events in the input dataset supporting the locus signal
* `p_value`: Significance value for the locus in the SPICE selection analysis
* `genes`: Comma-separated genes overlapping the inferred locus interval
