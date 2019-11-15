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
SEQTK_BIN_DIR=/usr/local/packages/seqtk-1.2/bin
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

## Download FASTQ files (2019-05-30)

All FASTQ files from 6 of the studies were downloaded from the SRA. Because the Wolbachia reads were not included in the FASTQ files uploaded for the Luck 2015 study, NEB uploaded them onto the IGS FTP site (ftp://nsflgt@ftp.igs.umaryland.edu; password: drosophila) for download to this directory: "$READS_DIR"/wDi_luck_2015

### SRA Datasets
#### Input Sets:
```{bash, eval = F}
## wOo Darby 2012
SRR_ID_LIST="$INPUTS_DIR"/wOo_darby_2012_srr.list
OUTPUT_DIR="$READS_DIR"/wOo_darby_2012

## wMel Darby 2014
SRR_ID_LIST=/wMel_darby_2014_srr.list
OUTPUT_DIR="$READS_DIR"/wMel_darby_2014

## wDi Luck 2014
SRR_ID_LIST="$INPUTS_DIR"/wDi_luck_2014_srr.list
OUTPUT_DIR="$READS_DIR"/wDi_luck_2014

## wMel Gutzwiller 2015
SRR_ID_LIST="$INPUTS_DIR"/wMel_gutzwiller_2015_srr.list
OUTPUT_DIR="$READS_DIR"/wMel_gutzwiller_2015

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

### Luck et al 2015 Dataset

#### Input Sets:
```{bash, eval = F}
## wDi Luck 2015
SRR_ID_LIST="$INPUTS_DIR"/wDi_luck_2015_srr.list
OUTPUT_DIR="$READS_DIR"/wDi_luck_2015
```

#### Commands:
```{bash, eval = F}
zcat "$READS_DIR"/wDi_luck_2015/Galaxy1-[Di_Female_Intestine_2].fastqsanger.zip > "$READS_DIR"/wDi_luck_2015/SRR1974179.fastq
zcat "$READS_DIR"/wDi_luck_2015/Galaxy2-[Di_Female_Intestine_1].fastqsanger.zip > "$READS_DIR"/wDi_luck_2015/SRR1974178.fastq
zcat "$READS_DIR"/wDi_luck_2015/Galaxy3-[Di_Female_Body_Wall_2].fastqsanger.zip > "$READS_DIR"/wDi_luck_2015/SRR1974175.fastq
zcat "$READS_DIR"/wDi_luck_2015/Galaxy4-[Di_Female_Uterus_1].fastqsanger.zip > "$READS_DIR"/wDi_luck_2015/SRR1974180.fastq
zcat "$READS_DIR"/wDi_luck_2015/Galaxy5-[Di_Female_Head_2].fastqsanger.zip > "$READS_DIR"/wDi_luck_2015/SRR1974177.fastq
zcat "$READS_DIR"/wDi_luck_2015/Galaxy6-[Di_Female_Uterus_2].fastqsanger.zip > "$READS_DIR"/wDi_luck_2015/SRR1974181.fastq
zcat "$READS_DIR"/wDi_luck_2015/Galaxy7-[Di_Female_Head_1].fastqsanger.zip > "$READS_DIR"/wDi_luck_2015/SRR1974176.fastq
zcat "$READS_DIR"/wDi_luck_2015/Galaxy8-[Di_Female_Body_Wall_1].fastqsanger.zip > "$READS_DIR"/wDi_luck_2015/SRR1974174.fastq
zcat "$READS_DIR"/wDi_luck_2015/Galaxy9-[Di_Male_Testes].fastq.zip > "$READS_DIR"/wDi_luck_2015/SRR1974185.fastq
zcat "$READS_DIR"/wDi_luck_2015/Galaxy10-[Di_Male_Body_Wall_1].fastq.zip > "$READS_DIR"/wDi_luck_2015/SRR1974182.fastq
zcat "$READS_DIR"/wDi_luck_2015/Galaxy11-[Di_Male_Body_Wall_2] (1).fastq.zip > "$READS_DIR"/wDi_luck_2015/SRR1974183.fastq
zcat "$READS_DIR"/wDi_luck_2015/Galaxy12-[Di_Male_Intestine] (1).fastq.zip > "$READS_DIR"/wDi_luck_2015/SRR1974184.fastq
```

## Align FASTQ files to combined host and Wolbachia reference while allowing for splicing (2019-05-30)

### Single-End FASTQs
#### Input Sets:
```{bash, eval = F}
THREADS=4

## wOo Darby 2012
REF_FNA="$REFERENCES_DIR"/wOo_oochengi_combined.fna
SRR_ID_LIST="$INPUTS_DIR"/wOo_oochengi_2012_srr.list
OUTPUT_DIR="$READS_DIR"/wOo_oochengi_2012/bam
STRANDEDNESS=""
PREFIX=wOo_darby_2012

## wDi Luck 2014
REF_FNA="$REFERENCES_DIR"/wDi_dimmitis_combined.fna
SRR_ID_LIST="$INPUTS_DIR"/wDi_luck_2014_srr.list
OUTPUT_DIR="$READS_DIR"/wDi_luck_2014/bam
STRANDEDNESS=""
PREFIX=wDi_luck_2014

## wDi Luck 2015
REF_FNA="$REFERENCES_DIR"/wDi_dimmitis_combined.fna
SRR_ID_LIST="$INPUTS_DIR"/wDi_luck_2015_srr.list
OUTPUT_DIR="$READS_DIR"/wDi_luck_2015/bam
STRANDEDNESS=""
PREFIX=wDi_luck_2015
```

#### Commands:
```{bash, eval = F}
while read SRR
do
	echo ""$HISAT2_BIN_DIR"/hisat2 -p "$THREADS" -x "$REF_FNA" "$STRANDEDNESS" -U "$READS_DIR"/"$PREFIX"/$SRR"_1.fastq | "$SAMTOOLS_BIN_DIR"/samtools view -bhSo "$OUTPUT_DIR"/"$SRR".bam -" | qsub -q threaded.q  -pe thread "$THREADS" -P jdhotopp-lab -l mem_free=20G -N hisat2 -wd "$OUTPUT_DIR"
done < "$SRR_ID_LIST"
```

### Paired-End FASTQs
#### Input Sets:
```{bash, eval = F}
THREADS=4

## wMel Darby 2014
REF_FNA="$REFERENCES_DIR"/wMel_dmelanogaster_combined.fna
SRR_ID_LIST="$INPUTS_DIR"/wMel_darby_2014_srr.list
OUTPUT_DIR="$READS_DIR"/wMel_darby_2014/bam
STRANDEDNESS="--rna-strandness FR"
PREFIX=wMel_darby_2014

## wMel Gutzwiller 2015
REF_FNA="$REFERENCES_DIR"/wMel_dmelanogaster_combined.fna
SRR_ID_LIST="$INPUTS_DIR"/wMel_gutzwiller_2015_srr.list
OUTPUT_DIR="$READS_DIR"/wMel_gutzwiller_2015/bam
STRANDEDNESS="--rna-strandness RF"
PREFIX=wMel_gutzwiller_2015

## wBm Grote 2017
REF_FNA="$REFERENCES_DIR"/wBm_bmalayi_combined.fna
SRR_ID_LIST="$INPUTS_DIR"/wBm_grote_2017_srr.list
OUTPUT_DIR="$READS_DIR"/wBm_grote_2017/bam
STRANDEDNESS=""
PREFIX=wBm_grote_2017

## wBm Chung 2019
REF_FNA="$REFERENCES_DIR"/wBm_bmalayi_combined.fna
SRR_ID_LIST="$INPUTS_DIR"/wBm_chung_2019_srr.list
OUTPUT_DIR="$READS_DIR"/wBm_chung_2019/bam
STRANDEDNESS="--rna-strandness RF"
PREFIX=wBm_chung_2019
```

#### Commands:
```{bash, eval = F}
while read srr
do
	echo ""$HISAT2_BIN_DIR"/hisat2 -p "$THREADS" -x "$REF_FNA" "$STRANDEDNESS" -1 "$READS_DIR"/"$PREFIX"/"$SRR"_1.fastq  -2 "$READS_DIR"/"$PREFIX"/"$SRR"_2.fastq | "$SAMTOOLS_BIN_DIR"/samtools view -bhSo "$OUTPUT_DIR"/"$SRR".bam -" | qsub -q threaded.q  -pe thread "$THREADS" -P jdhotopp-lab -l mem_free=20G -N hisat2 -wd "$OUTPUT_DIR"
done < "$SRR_ID_LIST"
```

## Create a list of Wolbachia-mapping reads in all BAM files (2019-05-30)

#### Input Sets:
```{bash, eval = F}
## wOo Darby 2012
BAM_DIR="$READS_DIR"/wOo_darby_2012/bam
SRR_ID_LIST="$INPUTS_DIR"/wOo_darby_2012_srr.list
CONTIG=NC_018267.1

## wMel Darby 2014
BAM_DIR="$READS_DIR"/wMel_darby_2014/bam
SRR_ID_LIST=/wMel_darby_2014_srr.list
CONTIG=NC_002978.6

## wDi Luck 2014
BAM_DIR="$READS_DIR"/wDi_luck_2014/bam
SRR_ID_LIST="$INPUTS_DIR"/wDi_luck_2014_srr.list
CONTIG=wDi22.scaf1

## wDi Luck 2015
BAM_DIR="$READS_DIR"/wDi_luck_2015/bam
SRR_ID_LIST="$INPUTS_DIR"/wDi_luck_2015_srr.list
CONTIG=wDi22.scaf1

## wMel Gutzwiller 2015
BAM_DIR="$READS_DIR"/wMel_gutzwiller_2015/bam
SRR_ID_LIST="$INPUTS_DIR"/wMel_gutzwiller_2015_srr.list
CONTIG=NC_002978.6

## wBm Grote 2017
BAM_DIR="$READS_DIR"/wBm_grote_2017/bam
SRR_ID_LIST="$INPUTS_DIR"/wBm_grote_2017_srr.list
CONTIG=NC_006833.1

## wBm Chung 2019BAM_DIR=="$READS_DIR"/wBm_chung_2019/bam
BAM_DIR="$READS_DIR"/wBm_chung_2019/bam
SRR_ID_LIST="$INPUTS_DIR"/wBm_chung_2019_srr.list
CONTIG=NC_006833.1
```

#### Commands:
```{bash, eval = F}
while read SRR
do

	"#SAMTOOLS_BIN_DIR"/samtools view "$BAM_DIR"/"SRR".bam | grep "$PREFIX" | awk '{print $1}' | sort -n | uniq > "$BAM_DIR"/"SRR".target_reads.list

done < "$SRR_ID_LIST"
```

## Create subsets of all FASTQ files that contain only Wolbachia-mapping reads (2019-05-30)

### Single-End FASTQs

#### Input Sets:
```{bash, eval = F}
## wOo Darby 2012
SRR_ID_LIST="$INPUTS_DIR"/wOo_oochengi_2012_srr.list
OUTPUT_DIR="$READS_DIR"/wOo_oochengi_2012

## wDi Luck 2014
SRR_ID_LIST="$INPUTS_DIR"/wDi_luck_2014_srr.list
OUTPUT_DIR="$READS_DIR"//wDi_luck_2014

## wDi Luck 2015
SRR_ID_LIST="$INPUTS_DIR"/wDi_luck_2015_srr.list
OUTPUT_DIR="$READS_DIR"/wDi_luck_2015
```

#### Commands:
```{bash, eval = F}
while read SRR
do

	echo -e ""$SEQTK_BIN_DIR"/seqtk subseq "$OUTPUT_DIR"/"$SRR"_1.fastq "$OUTPUT_DIR"/"$SRR".target_reads.list > "$OUTPUT_DIR"/"$SRR"_1.subset.fastq | qsub -P jdhotopp-lab -l mem_free=2G -N seqtk -wd "$OUTPUT_DIR"

done < "$SRR_ID_LIST"
```

### Paired-End FASTQs

#### Input Sets:
```{bash, eval = F}
## wMel Darby 2014
SRR_ID_LIST=/wMel_darby_2014_srr.list
OUTPUT_DIR="$READS_DIR"/wMel_darby_2014

## wMel Gutzwiller 2015
SRR_ID_LIST="$INPUTS_DIR"/wMel_gutzwiller_2015_srr.list
OUTPUT_DIR="$READS_DIR"/wMel_gutzwiller_2015

## wBm Grote 2017
SRR_ID_LIST="$INPUTS_DIR"/wBm_grote_2017_srr.list
OUTPUT_DIR="$READS_DIR"/wBm_grote_2017

## wBm Chung 2019BAM_DIR=="$READS_DIR"/wBm_chung_2019/bam
SRR_ID_LIST="$INPUTS_DIR"/wBm_chung_2019_srr.list
OUTPUT_DIR="$READS_DIR"/wBm_chung_2019
```

#### Commands:
```{bash, eval = F}
while read SRR
do

	echo -e ""$SEQTK_BIN_DIR"/seqtk subseq "$OUTPUT_DIR"/"$SRR"_1.fastq "$OUTPUT_DIR"/"$SRR".target_reads.list > "$OUTPUT_DIR"/"$SRR"_1.subset.fastq | qsub -P jdhotopp-lab -l mem_free=2G -N seqtk -wd "$OUTPUT_DIR"

	echo -e ""$SEQTK_BIN_DIR"/seqtk subseq "$OUTPUT_DIR"/"$SRR"_2.fastq "$OUTPUT_DIR"/"$SRR".target_reads.list > "$OUTPUT_DIR"/"$SRR"_2.subset.fastq | qsub -P jdhotopp-lab -l mem_free=2G -N seqtk -wd "$OUTPUT_DIR"

done < "$SRR_ID_LIST"
```

## Align FASTQ files containing only Wolbachia-mapping reads to combined host and Wolbachia reference while disallowing splicing (2019-06-01)

### Single-End FASTQs
#### Input Sets:
```{bash, eval = F}
THREADS=4

## wOo Darby 2012
REF_FNA="$REFERENCES_DIR"/wOo_oochengi_combined.fna
SRR_ID_LIST="$INPUTS_DIR"/wOo_oochengi_2012_srr.list
OUTPUT_DIR="$READS_DIR"/wOo_oochengi_2012/bam
STRANDEDNESS=""
PREFIX=wOo_darby_2012

## wDi Luck 2014
REF_FNA="$REFERENCES_DIR"/wDi_dimmitis_combined.fna
SRR_ID_LIST="$INPUTS_DIR"/wDi_luck_2014_srr.list
OUTPUT_DIR="$READS_DIR"/wDi_luck_2014/bam
STRANDEDNESS=""
PREFIX=wDi_luck_2014

## wDi Luck 2015
REF_FNA="$REFERENCES_DIR"/wDi_dimmitis_combined.fna
SRR_ID_LIST="$INPUTS_DIR"/wDi_luck_2015_srr.list
OUTPUT_DIR="$READS_DIR"/wDi_luck_2015/bam
STRANDEDNESS=""
PREFIX=wDi_luck_2015
```

#### Commands:
```{bash, eval = F}
while read SRR
do
	echo -e ""$HISAT2_BIN_DIR"/hisat2 -p "$THREADS" -x "$REF_FNA" "$STRANDEDNESS" -U "$READS_DIR"/"$PREFIX"/$SRR"_1.subset.fastq --no-spliced-alignment -X 1000 | "$SAMTOOLS_BIN_DIR"/samtools view -bhSo "$OUTPUT_DIR"/"$SRR".bam -" | qsub -q threaded.q  -pe thread "$THREADS" -P jdhotopp-lab -l mem_free=20G -N hisat2 -wd "$OUTPUT_DIR"
done < "$SRR_ID_LIST"
```

### Paired-End FASTQs
#### Input Sets:
```{bash, eval = F}
THREADS=4

## wMel Darby 2014
REF_FNA="$REFERENCES_DIR"/wMel_dmelanogaster_combined.fna
SRR_ID_LIST="$INPUTS_DIR"/wMel_darby_2014_srr.list
OUTPUT_DIR="$READS_DIR"/wMel_darby_2014/bam
STRANDEDNESS="--rna-strandness FR"
PREFIX=wMel_darby_2014

## wMel Gutzwiller 2015
REF_FNA="$REFERENCES_DIR"/wMel_dmelanogaster_combined.fna
SRR_ID_LIST="$INPUTS_DIR"/wMel_gutzwiller_2015_srr.list
OUTPUT_DIR="$READS_DIR"/wMel_gutzwiller_2015/bam
STRANDEDNESS="--rna-strandness RF"
PREFIX=wMel_gutzwiller_2015

## wBm Grote 2017
REF_FNA="$REFERENCES_DIR"/wBm_bmalayi_combined.fna
SRR_ID_LIST="$INPUTS_DIR"/wBm_grote_2017_srr.list
OUTPUT_DIR="$READS_DIR"/wBm_grote_2017/bam
STRANDEDNESS=""
PREFIX=wBm_grote_2017

## wBm Chung 2019
REF_FNA="$REFERENCES_DIR"/wBm_bmalayi_combined.fna
SRR_ID_LIST="$INPUTS_DIR"/wBm_chung_2019_srr.list
OUTPUT_DIR="$READS_DIR"/wBm_chung_2019/bam
STRANDEDNESS="--rna-strandness RF"
PREFIX=wBm_chung_2019
```

#### Commands:
```{bash, eval = F}
while read SRR
do
	echo -e ""$HISAT2_BIN_DIR"/hisat2 -p "$THREADS" -x "$REF_FNA" "$STRANDEDNESS" -1 "$READS_DIR"/"$PREFIX"/"$SRR"_1.subset.fastq  -2 "$READS_DIR"/"$PREFIX"/"$SRR"_2.subset.fastq --no-spliced-alignment -X 1000 | "$SAMTOOLS_BIN_DIR"/samtools view -bhSo "$OUTPUT_DIR"/"$SRR".bam -" | qsub -q threaded.q  -pe thread "$THREADS" -P jdhotopp-lab -l mem_free=20G -N hisat2 -wd "$OUTPUT_DIR"
done < "$SRR_ID_LIST"
```

## Sort BAM files (2019-06-01)

#### Input Sets:
```{bash, eval = F}
THREADS=4

## wOo Darby 2012
BAM_DIR="$READS_DIR"/wOo_darby_2012/bam

## wMel Darby 2014
BAM_DIR="$READS_DIR"/wMel_darby_2014/bam

## wDi Luck 2014
BAM_DIR="$READS_DIR"/wDi_luck_2014/bam

## wDi Luck 2015
BAM_DIR="$READS_DIR"/wDi_luck_2015/bam

## wMel Gutzwiller 2015
BAM_DIR="$READS_DIR"/wMel_gutzwiller_2015/bam

## wBm Grote 2017
BAM_DIR="$READS_DIR"/wBm_grote_2017/bam

## wBm Chung 2019
BAM_DIR="$READS_DIR"/wBm_chung_2019/bam
```

#### Commands:
```{bash, eval = F}
while read SRR
do
	echo -e ""$SAMTOOLS_BIN_DIR"/samtools sort -@ "$THREADS" -o "$BAM_DIR"/"$SRR".sortedbyposition.bam "$BAM_DIR"/"$SRR".bam" | qsub -q threaded.q -pe thread "$THREADS" -P jdhotopp-lab -l mem_free=2G -N samtools -wd "$BAM_DIR"
done < "$SRR_ID_LIST"
```

## Index BAM files (2019-06-01)

#### Input Sets:
```{bash, eval = F}
THREADS=4

## wOo Darby 2012
BAM_DIR="$READS_DIR"/wOo_darby_2012/bam

## wMel Darby 2014
BAM_DIR="$READS_DIR"/wMel_darby_2014/bam

## wDi Luck 2014
BAM_DIR="$READS_DIR"/wDi_luck_2014/bam

## wDi Luck 2015
BAM_DIR="$READS_DIR"/wDi_luck_2015/bam

## wMel Gutzwiller 2015
BAM_DIR="$READS_DIR"/wMel_gutzwiller_2015/bam

## wBm Grote 2017
BAM_DIR="$READS_DIR"/wBm_grote_2017/bam

## wBm Chung 2019
BAM_DIR="$READS_DIR"/wBm_chung_2019/bam
```

#### Commands:
```{bash, eval = F}
while read SRR
do
	echo -e ""$SAMTOOLS_BIN_DIR"/samtools index -@ "$THREADS" "$BAM_DIR"/"$SRR".sortedbyposition.bam" | qsub -q threaded.q -pe thread "$THREADS" -P jdhotopp-lab -l mem_free=2G -N samtools -wd "$BAM_DIR"
done < "$SRR_ID_LIST"
```