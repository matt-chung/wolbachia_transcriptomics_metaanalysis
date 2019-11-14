# A meta-analysis of Wolbachia transcriptomics reveals a stage-specific Wolbachia transcriptional response shared across different hosts
 
# Table of Contents
1. [Set software, input, and output paths](#)


1. [Download, align, and quantify Wolbachia sequencing data](#)
	1. Create directories

# Set software and directory paths <a name="setpaths"></a>

For rerunning analyses, all paths in this section must be set by the user.

## Software <a name="setpaths_software"></a>

```{bash, eval = F}
HISAT2_BIN_DIR=/usr/local/packages/hisat2-2.1.0
SAMTOOLS_BIN_DIR=/usr/local/packages/samtools-1.9/bin
SRATOOLKIT_BIN_DIR=/usr/local/packages/sratoolkit-2.9.0/bin
```

## Directories <a name="setpaths_directories"></a>

```{bash, eval = F}
INPUTS_DIR=
READS_DIR=/local/aberdeen2rw/julie/Matt_dir/wolbachia_metatranscriptome_analysis/
REFERENCES_DIR=/local/aberdeen2rw/julie/Matt_dir/wolbachia_metatranscriptome_analysis/references/
WORKING_DIR=
```

# Download, align, and quantify Wolbachia sequencing data

## Create directories (2019-05-30)
```{bash, eval = F}
mkdir -p "$REFERENCES_DIR"

mkdir -p "$READS_DIR"/wOo_darby_2012/bam
mkdir -p "$READS_DIR"/wMel_darby_2014/bam
mkdir -p "$READS_DIR"/wDi_luck_2014/bam
mkdir -p "$READS_DIR"/wDi_luck_2015/bam
mkdir -p "$READS_DIR"/wMel_gutzwiller_2015/bam
mkdir -p "$READS_DIR"/wBm_grote_2017/bam
mkdir -p "$READS_DIR"/wBm_chung_2019/bam
```

## Download reference files and create combined references (2019-05-30)

The gff reference for wDi must be downloaded from this link: http://rast.nmpdr.org/rast.cgi?page=JobDetails&job=53659 (login:guest, password:guest). All subsequent code will make use of this wDi.gff at this path: "$REFERENCES_DIR"/wDi.gff

#### Commands:
```{bash, eval = F}
## wBm + B. malayi
wget -O "$REFERENCES_DIR"/wBm.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/385/GCF_000008385.1_ASM838v1/GCF_000008385.1_ASM838v1_genomic.fna.gz
wget -O "$REFERENCES_DIR"/bmalayi.fna.gz ftp://ftp.wormbase.org/pub/wormbase/releases/WS270/species/b_malayi/PRJNA10729/b_malayi.PRJNA10729.WS270.genomic.fa.gz
gunzip "$REFERENCES_DIR"/wBm.fna.gz
gunzip "$REFERENCES_DIR"/bmalayi.fna.gz
cat "$REFERENCES_DIR"/wBm.fna "$REFERENCES_DIR"/bmalayi.fna > "$REFERENCES_DIR"/wBm_bmalayi_combined.fna

wget -O "$REFERENCES_DIR"/wBm.gff.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/385/GCF_000008385.1_ASM838v1/GCF_000008385.1_ASM838v1_genomic.gff.gz
gunzip "$REFERENCES_DIR"/wBm.gff.gz

## wDi + D. immitis
wget -O "$REFERENCES_DIR"/wDi.fna.gz http://nematodes.org/downloads/959nematodegenomes/blast/db2/Dirofilaria_immitis_wolbachia_2.2.fna.gz
wget -O "$REFERENCES_DIR"/dimmitis.fna.gz http://nematodes.org/downloads/959nematodegenomes/blast/db2/Dirofilaria_immitis_2.2.fna.gz
gunzip "$REFERENCES_DIR"/wDi.fna.gz
gunzip "$REFERENCES_DIR"/dimmitis.fna.gz
cat "$REFERENCES_DIR"/wDi.fna "$REFERENCES_DIR"/dimmitis.fna > "$REFERENCES_DIR"/wDi_dimmitis_combined.fna

wget -O "$REFERENCES_DIR"/wBm.gff.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/385/GCF_000008385.1_ASM838v1/GCF_000008385.1_ASM838v1_genomic.gff.gz
gunzip "$REFERENCES_DIR"/wBm.gff.gz

## wMel + D. melanogaster
wget -O "$REFERENCES_DIR"/wMel.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/025/GCF_000008025.1_ASM802v1/GCF_000008025.1_ASM802v1_genomic.fna.gz
wget -O "$REFERENCES_DIR"/dmelanogaster.fna.gz ftp://ftp.flybase.net/releases/FB2019_02/dmel_r6.27/fasta/dmel-all-chromosome-r6.27.fasta.gz
gunzip "$REFERENCES_DIR"/wMel.fna.gz
gunzip "$REFERENCES_DIR"/dmelanogaster.fna.gz
cat "$REFERENCES_DIR"/wMel.fna "$REFERENCES_DIR"/dmelanogaster.fna > "$REFERENCES_DIR"/wMel_dmelanogaster_combined.fna

wget -O "$REFERENCES_DIR"/wMel.gff.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/025/GCF_000008025.1_ASM802v1/GCF_000008025.1_ASM802v1_genomic.gff.gz
gunzip "$REFERENCES_DIR"/wMel.gff.gz

## wOo + O. ochengi
wget -O "$REFERENCES_DIR"/wOo.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/306/885/GCF_000306885.1_ASM30688v1/GCF_000306885.1_ASM30688v1_genomic.fna.gz
wget -O "$REFERENCES_DIR"/oochengi.fna.gz http://www.nematodes.org/downloads/959nematodegenomes/blast/db2/Onchocerca_ochengi_nuclear_assembly_nOo.2.0.fna.gz
gunzip "$REFERENCES_DIR"/wOo.fna.gz
gunzip "$REFERENCES_DIR"/oochengi.fna.gz
cat "$REFERENCES_DIR"/wOo.fna "$REFERENCES_DIR"/oochengi.fna > "$REFERENCES_DIR"/wOo_oochengi_combined.fna

wget -O "$REFERENCES_DIR"/wOol.gff.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/306/885/GCF_000306885.1_ASM30688v1/GCF_000306885.1_ASM30688v1_genomic.gff.gz
gunzip "$REFERENCES_DIR"/wOo.gff.gz
```

## Download FASTQ files from the SRA (2019-05-30)

#### Input Sets:
```{bash, eval = F}
## wOo Darby 2012
SRR_ID_LIST="$INPUTS_DIR"/wOo_darby_2012_srr.list
OUTPUT_DIR="$READS_DIR"/wOo_darby_2012

## wMel__DIR"/wMel_darby_2014_srr.list
OUTPUT_DIR="$READS_DIR"/wMel_darby_2014

## wBm Grote 2017
SRR_ID_LIST="$INPUTS_DIR"/wDi_luck_2014_srr.list
OUTPUT_DIR="$READS_DIR"/wBm_grote_2017

## wBm Grote 2017
SRR_ID_LIST="$INPUTS_DIR"/wDi_luck_2015_srr.list
OUTPUT_DIR="$READS_DIR"/wBm_grote_2017

## wBm Grote 2017
SRR_ID_LIST="$INPUTS_DIR"/wBm_grote_2017_srr.list
OUTPUT_DIR="$READS_DIR"/wBm_grote_2017

## wBm Grote 2017
SRR_ID_LIST="$INPUTS_DIR"/wBm_grote_2017_srr.list
OUTPUT_DIR="$READS_DIR"/wBm_grote_2017

## wBm Chung 2019
SRR_ID_LIST="$INPUTS_DIR"/wBm_chung_2019_srr.list
OUTPUT_DIR="$READS_DIR"/wBm_chung_2019
```

#### Commands:
```{bash, eval = F}
while read Line
do

SRR_ID="$(echo $Line | awk -F "\t" '{print $1}')"

qsub -P jdhotopp-lab -l mem_free=2G -N fastq_dump -wd "$OUTPUT_DIR" -b y "$SRATOOLKIT_BIN_DIR"/fastq-dump --split-files "$SRR_ID" -O "$OUTPUT_DIR"

done < "$SRR_ID_LIST"
```

## Align FASTQ files to combined host and Wolbachia reference while allowing for splicing (2019-05-30)

#### Input Sets:
```{bash, eval = F}
THREADS=4

## wBm Chung 2019
REF_FNA="$REFERENCES_DIR"/wBm_bmalayi_combined.fna
SRR_ID_LIST="$INPUTS_DIR"/wBm_chung_2019_srr.list
OUTPUT_DIR="$READS_DIR"/wBm_chung_2019/bam
STRANDEDNESS="--rna-strandness RF"
```

```{bash, eval = F}
while read srr
do

echo ""$HISAT2_BIN_DIR"/hisat2 -p "$THREADS" -x "$REF_FNA" "$STRANDEDNESS" -1 "$READS_DIR"/"$SRR"_1.fastq  -2 "$READS_DIR"/"$SRR"_2.fastq | "$SAMTOOLS_BIN_DIR"/samtools view -bhSo "$OUTPUT_DIR"/"$SRR".bam -" | qsub -q threaded.q  -pe thread "$THREADS" -P jdhotopp-lab -l mem_free=20G -N hisat2 -wd "$OUTPUT_DIR"

done < "$SRR_ID_LIST"
```

## Identify Wolbachia-mapping reads in all BAM files (2019-05-30)


## Create subsets of all FASTQ files that contain only Wolbachia-mapping reads (2019-05-30)
