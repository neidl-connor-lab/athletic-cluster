# Understanding SARS-CoV-2 transmission: Linking traditional contact tracing with genomic surveillance of a SARS-CoV-2 student cluster on a large urban campus
Bioinformatic methods for iScience paper "Understanding SARS-CoV-2 transmission: Linking traditional contact tracing with genomic surveillance of a SARS-CoV-2 student cluster on a large urban campus"

## Methods overview
Sequencing libraries were generated using the ARTIC v4 primer sets and the Illumina COVIDSeq assay according to the manufacturer's protocol. Library pools were sequenced on a NextSeq 500 platform, generating 75x75 bp paired reads.

### Processing pipeline
Software and pipeline are outlined below. For full details, please see `pipeline.sh`.

#### Software packages used

| Package     | Version |
| :---------: | :-----: |
| trimmomatic | 0.36    |
| bowtie2     | 2.3.4.1 |
| samtools    | 1.10    |
| lofreq      | 2.1.3.1 |

#### Pipeline overview

1. FASTQ files were trimmed using `trimmomatic`
2. Trimmed FASTQ files were aligned to the Wuhan Hu-1 RefSeq ([MN908947.3](https://www.ncbi.nlm.nih.gov/nuccore/MN908947)) using `bowtie2`
3. The alignment BAM file was sorted via `samtools sort`
4. Viral genome coverage was calculated using the sorted BAM and `samtools depth`
5. Indels were assessed with `lofreq indelqual` and the resulting BAM was indexed with `samtools index`
6. SNVs were quantified with `lofreq call` and annotated with a custom R script

### Clustering analysis
To be continued...

