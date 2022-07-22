# Understanding SARS-CoV-2 transmission: Linking traditional contact tracing with genomic surveillance of a SARS-CoV-2 student cluster on a large urban campus
Bioinformatic methods for iScience paper "Understanding SARS-CoV-2 transmission: Linking traditional contact tracing with genomic surveillance of a SARS-CoV-2 student cluster on a large urban campus"

## Methods overview
Sequencing libraries were generated using the ARTIC v4 primer sets and the Illumina COVIDSeq assay according to the manufacturer's protocol. Library pools were sequenced on a NextSeq 500 platform, generating 75x75 bp paired reads.

### Processing pipeline
Software and pipeline are outlined below. For full details, please see `pipeline.sh`.

#### Software used

Processing was performed on the [BU Supercomputing Cluster (SCC)](https://www.bu.edu/tech/support/research/computing-resources/scc/). See the wrapper script, `pipeline.sh` for full details. Consensus sequences were assigned lineages using the [pangolin web portal](https://pangolin.cog-uk.io/).

| Package    | Version |
| :--------: | :-----: |
| bowtie2    | 2.3.4.1 |
| samtools   | 1.15.1  |
| lofreq     | 2.1.3.1 |
| R          | 4.0.2   |
| pangolin   | 4.1.1   |
| pango data | 1.11    |

#### Pipeline overview

1. FASTQ files were aligned to the Wuhan Hu-1 RefSeq ([MN908947.3](https://www.ncbi.nlm.nih.gov/nuccore/MN908947)) using `bowtie2`
2. Primers were softclipped from mapped reads using `samtools ampliconclip`
3. The alignment BAM file was sorted via `samtools sort`
4. Viral genome coverage was calculated using the sorted BAM and `samtools depth`
5. Indels were assessed with `lofreq indelqual` and the resulting BAM was indexed with `samtools index`
6. SNVs were quantified with `lofreq call` and annotated with a custom R script

### Clustering analysis

Clusters were assigned using a pairwise SNV difference matrix calculated from the coverage and VCF file. Any SNVs with <10X read depth were removed, as were SNVs in the 3' and 5' UTRs (bases <265 and >29675, respectively). Only consensus SNVs (>50%) were considered. Numeric cluster IDs were assigned using bins of 0, 1, and 2 SNV differences.

| Tool      | Version  |
| :-------: | :------: |
| R         | 4.0.2    |
| tidyverse | 1.3.0    |
| argparse  | 2.1.5    |


