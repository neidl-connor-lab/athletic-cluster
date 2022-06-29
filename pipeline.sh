#!/bin/bash -l

# qsub options
#$ -P lasvchal
#$ -l h_rt=24:00:00
#$ -pe omp 4
#$ -j y
#$ -o log-$JOB_NAME.qlog

# job info
echo "=========================================================="
echo "Start date: $(date)"
echo "Running on node: $(hostname)"
echo "Current directory: $(pwd)"
echo "Job name: $JOB_NAME"
echo "Job ID: $JOB_ID"
echo "=========================================================="
echo ""

## helper functions --------------------------------------------
mesg () { echo -e "[MSG] $@"; }
err () { echo -e "[ERR] $@"; exit 1; }
checkcmd () {
  if [ $? -eq 0 ]
  then
    mesg "$@ succeeded"
  else
    err "$@ failed"
  fi
}
checkmodule () {
  # using 2> to keep error quiet
  # check if it's already available
  if [ -z "$(which $@ 2> /dev/null)" ]
  then
    # try to load module on cluster
    module load $@ 2> /dev/null
    if [ -z "$(which $@ 2> /dev/null)" ]
    then
      err "$@ not found"
    fi
  fi
}
checkfile () {
  if [ ! -f "$@" ]
  then
    err "File not found: $@"
  fi
}

## pre-set variables --------------------------------------------
ANNOT="pipeline/annotate-snvs.r" # SNV annotation script
ASSEMBLE="pipeline/assemble-consensus.r" # consensus assembly
LOFREQ="pipeline/lofreq/lofreq" # LoFreq executable
IDX="pipeline/bowtie2/MN908947.3" # bowtie2 index
FASTA="pipeline/sars-cov2-MN908947.3.fa" # CoV2 genome
CSV="pipeline/sars-cov2-MN908947.3.csv" # CoV2 annotation

## default values and help -------------------------------------
ODIR="data"
BED="pipeline/artic-v4.bed"

HELP="usage: qsub -N JOBNAME $(basename "$0") [options] ID FQ1 FQ2

positional arguments:
  ID  sample ID 
  FQ1 forward FASTQ file
  FQ2 [optional] reverse FASTQ file

options (default):
  -o output directory ($ODIR)
  -b primer BED file ($BED)
  -h print this message and exit
"

# parsing arguments
while getopts ":ho:b:" opt 
do 
  case ${opt} in 
    o ) ODIR="${OPTARG}"
      ;;
    b ) BED="${OPTARG}"
      ;;
    h ) echo "$HELP" && exit 0
      ;;
    \? ) err "Invalid option ${opt}\n${HELP}"
      ;;
  esac
done
shift $((OPTIND -1))

# set positional arguments
ID="$1"
FQ1="$2"
FQ2="$3"

## check modules and resources ---------------------------------
mesg "Step 0/6: Check inputs & dependencies"

# modules
checkmodule "bowtie2"
checkmodule "samtools"
checkmodule "R"

# scripts & resources
checkfile $ANNOT
checkfile $ASSEMBLE
checkfile $LOFREQ
checkfile $FASTA
checkfile $CSV
checkfile $BED

# special check for bowtie2 index: 6 files
if [ "$(ls -1 ${IDX}.* | wc -l)" -ne 6 ]
then
  err "Invalid bowtie2 index: $IDX"
fi

# sample ID
if [ -z "$ID" ]
then
  err "No sample ID provided\n${HELP}"
else
  mesg "Sample ID: $ID"
fi

# input FASTQ(s)
checkfile $FQ1
mesg "R1 FASTQ file: $FQ1"
# if FQ2, using paired reads
if [ ! -z "$FQ2" ]
then
  PAIRED=true
  checkfile $FQ2
  mesg "R2 FASTQ file: $FQ2"
else
  PAIRED=false
  mesg "No R2 file detected; assuming unpaired reads"
fi

# output directory
if [ ! -d "$ODIR" ]
then
  mesg "Creating output directory: $ODIR"
  mkdir -p "$ODIR"
else
  mesg "Valid output directory: $ODIR"
fi

# done checking dependencies & inputs
VAR="$ODIR/$ID"
echo ""

## alignment ---------------------------------------------------
mesg "Step 1/6: Align to viral genome"

# paired vs. unpaired syntax
if $PAIRED
then
  CMD="bowtie2 --threads 4 -x '$IDX' -1 '$FQ1' -2 '$FQ2' > '${VAR}.sam'"
else
  CMD="bowtie2 --threads 4 -x '$IDX' -U '$FQ1' > '${VAR}.sam'"
fi

# run alignment
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Alignment"

# compress SAM to BAM
CMD="samtools view --threads 4 -b -h '${VAR}.sam' > '${VAR}-raw.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Compression"

# clean up SAM
rm "${VAR}.sam"
echo ""

## soft clip primers -------------------------------------------
mesg "Step 2/6: Soft-clip primers"
# TODO: make primer BED an input option
CMD="samtools ampliconclip --threads 4 -b '$BED' '${VAR}-raw.bam' -o '${VAR}-clipped.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Primer clipping"

# sort BAM; use max compression
CMD="samtools sort --threads 4 --output-fmt BAM -l 9 '${VAR}-clipped.bam' > '${VAR}.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Sorting"

# clean up raw and clipped bam file
rm "${VAR}-raw.bam" "${VAR}-clipped.bam"
echo ""

## calculate coverage ------------------------------------------
mesg "Step 3/6: Calculate coverage"
COVERAGE="${VAR}.tsv"
CMD="samtools depth --threads 4 -a -H '${VAR}.bam' > '$COVERAGE'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Coverage"
echo ""

## quantify SNVs -----------------------------------------------
mesg "Step 4/6: Quantify SNVs"
# score indels
CMD="$LOFREQ indelqual --dindel --ref '$FASTA' '${VAR}.bam' > '${VAR}-indel.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Indelqual"

# index indel BAM
CMD="samtools index '${VAR}-indel.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Indexing"

# run lofreq
VCF="${VAR}.vcf"
CMD="$LOFREQ call-parallel --pp-threads 4 --call-indels --min-cov 10 --ref '$FASTA' '${VAR}-indel.bam' > '$VCF'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "LoFreq"

# clean up
rm "${VAR}-indel.bam" "${VAR}-indel.bam.bai"
echo ""

## annotate SNVs -----------------------------------------------
mesg "Step 5/6: Annotate SNVs"
SNV="${VAR}.csv"
CMD="Rscript $ANNOT -c '$CSV' -v '$VCF' -o '$SNV'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Annotation"
echo ""

## assemble consensus ------------------------------------------
mesg "Step 6/6: Assemble consensus"
CONSENSUS="${VAR}.fa.gz"
CMD="Rscript $ASSEMBLE -r '$FASTA' -v '$VCF' -c '$COVERAGE' -o '$CONSENSUS'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Consensus"
echo ""

## package version ---------------------------------------------
mesg "Pipeline complete! Printing package versions..."
module list
bowtie2 --version
echo ""
samtools --version
echo ""
Rscript --version 
echo ""
$LOFREQ --version
echo ""
