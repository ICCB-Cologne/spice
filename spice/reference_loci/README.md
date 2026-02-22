# Overview of inferred loci from SPICE

This directory contains the reference loci set inferred by SPICE for the associated publication.

* `all_460_loci`: Complete loci set inferred by SPICE from 5,966 pan-cancer TCGA samples.

## Columns

* `chrom`: Chromosome
* `pos`: Genomic position
* `type`: Locus class where `OG` indicates oncogene-like gain selection and `TSG` indicates tumor suppressor-like loss selection.
* `start_ci`: Start of the locus confidence interval
* `end_ci`: End of the locus confidence interval
* `nr_events_in_data`: Number of events in the input dataset supporting the locus signal
* `p_value`: Significance value for the locus in the SPICE selection analysis
* `genes`: Comma-separated genes overlapping the inferred locus interval
