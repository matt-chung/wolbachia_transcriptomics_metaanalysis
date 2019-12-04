# A meta-analysis of Wolbachia transcriptomics reveals a stage-specific Wolbachia transcriptional response shared across different hosts
 
# Table of Contents


1. [Set software, input, and output directory paths](#setpaths)
	1. [Software](#setpaths_software)
	2. [Directories](#setpaths_directories)
2. [Download, align, and quantify Wolbachia sequencing data](#bashanalysis)
	1. [Create directories](#bashanalysis_createdirectories)
	2. [Download reference files and create combined references](#bashanalysis_getreferences)
	3. [Download FASTQ files](#bashanalysis_getfastq)
		1. [SRA Datasets](#bashanalysis_getfastq_sra)
		2. [Luck et al 2015 Dataset](#bashanalysis_getfastq_luck2015)
	4. [Align FASTQ files to combined host and Wolbachia reference while allowing for splicing](#bashanalysis_align)
		1. [Single-End FASTQs](#bashanalysis_align_se)
		2. [Paired-End FASTQs](#bashanalysis_align_pe)
	5. [Create a list of Wolbachia-mapping reads in all BAM files](#bashanalysis_wolbachiareads)
	6. [Create subsets of all FASTQ files that contain only Wolbachia-mapping reads](#bashanalysis_subsetfastq)
		1. [Single-End FASTQs](#bashanalysis_subsetfastq_se)
		2. [Paired-End FASTQs](#bashanalysis_subsetfastq_pe)
	7. [Align FASTQ files containing only Wolbachia-mapping reads to combined host and Wolbachia reference while disallowing splicing](#bashanalysis_realign)
		1. [Single-End FASTQs](#bashanalysis_realign_se)
		2. [Paired-End FASTQs](#bashanalysis_realign_pe)
	8. [Sort BAM files](#bashanalysis_realign_sortbam)
	9. [Index BAM files](#bashanalysis_indexbam)
	10. [Quantify Wolbachia genes](#bashanalysis_quant)
	11. [Assign InterPro descriptions and GO terms for Wolbachia transcripts]
3. [Transcriptomics meta-analysis](#ranalysis)

# Set software and directory paths <a name="setpaths"></a>

For rerunning analyses, all paths in this section must be set by the user.

## Software <a name="setpaths_software"></a>

```{bash, eval = F}
JULIA_DEPOT_PATH=/home/mattchung/.julia
PYTHON_LIB_PATH=/usr/local/packages/python-3.5/lib

JULIA_BIN_DIR=/usr/local/bin
PYTHON_BIN_DIR=/usr/local/bin

FADU_BIN_DIR=/home/mattchung/scripts/FADU
HISAT2_BIN_DIR=/usr/local/packages/hisat2-2.1.0
INTERPROSCAN_BIN_DIR=/local/aberdeen2rw/julie/Matt_dir/packages/interproscan-5.34-73.0
SAMTOOLS_BIN_DIR=/usr/local/packages/samtools-1.9/bin
SEQTK_BIN_DIR=/usr/local/packages/seqtk-1.2/bin
SRATOOLKIT_BIN_DIR=/usr/local/packages/sratoolkit-2.9.0/bin
```

## Directories <a name="setpaths_directories"></a>

```{bash, eval = F}
INPUTS_DIR=
SCRIPTS_DIR=/home/mattchung/scripts
READS_DIR=/local/aberdeen2rw/julie/Matt_dir/wolbachia_metatranscriptome_analysis/
REFERENCES_DIR=/local/aberdeen2rw/julie/Matt_dir/wolbachia_metatranscriptome_analysis/references/
WORKING_DIR=/local/projects-t3/EBMAL/mchung_dir/wolbachia_metatranscriptome_analysis/
```

# Download, align, and quantify Wolbachia sequencing data <a name="bashanalysis"></a>

## Create directories (2019-05-30) <a name="bashanalysis_createdirectories"></a>
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

## Download reference files and create combined references (2019-05-30) <a name="bashanalysis_getreferences"></a>

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

wget -O "$REFERENCES_DIR"/wOo.gff.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/306/885/GCF_000306885.1_ASM30688v1/GCF_000306885.1_ASM30688v1_genomic.gff.gz
gunzip "$REFERENCES_DIR"/wOo.gff.gz
```

## Download FASTQ files (2019-05-30) <a name="bashanalysis_getfastq"></a>

All FASTQ files from 6 of the studies were downloaded from the SRA. Because the Wolbachia reads were not included in the FASTQ files uploaded for the Luck 2015 study, NEB uploaded them onto the IGS FTP site (ftp://nsflgt@ftp.igs.umaryland.edu; password: drosophila) for download to this directory: "$READS_DIR"/wDi_luck_2015

### SRA Datasets <a name="bashanalysis_getfastq_sra"></a>
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

### Luck et al 2015 Dataset <a name="bashanalysis_getfastq_luck2015"></a>

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

## Align FASTQ files to combined host and Wolbachia reference while allowing for splicing (2019-05-30)  <a name="bashanalysis_align"></a>

### Single-End FASTQs <a name="bashanalysis_align_se"></a>
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
	echo ""$HISAT2_BIN_DIR"/hisat2 -p "$THREADS" -x "$REF_FNA" "$STRANDEDNESS" -U "$READS_DIR"/"$PREFIX"/"$SRR"_1.fastq | "$SAMTOOLS_BIN_DIR"/samtools view -bhSo "$OUTPUT_DIR"/"$SRR".bam -" | qsub -q threaded.q  -pe thread "$THREADS" -P jdhotopp-lab -l mem_free=20G -N hisat2 -wd "$OUTPUT_DIR"
done < "$SRR_ID_LIST"
```

### Paired-End FASTQs <a name="bashanalysis_align_pe"></a>
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

## Create a list of Wolbachia-mapping reads in all BAM files (2019-05-30) <a name="bashanalysis_wolbachiareads"></a>

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

## wBm Chung 2019
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

## Create subsets of all FASTQ files that contain only Wolbachia-mapping reads (2019-05-30) <a name="bashanalysis_subsetfastq"></a>

### Single-End FASTQs <a name="bashanalysis_subsetfastq_se"></a>

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
	echo -e ""$SEQTK_BIN_DIR"/seqtk subseq "$OUTPUT_DIR"/"$SRR"_1.fastq "$OUTPUT_DIR"/"$SRR".target_reads.list > "$OUTPUT_DIR"/"$SRR"_1.subset.fastq" | qsub -P jdhotopp-lab -l mem_free=2G -N seqtk -wd "$OUTPUT_DIR"
done < "$SRR_ID_LIST"
```

### Paired-End FASTQs <a name="bashanalysis_subsetfastq_pe"></a>

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

## wBm Chung 2019
SRR_ID_LIST="$INPUTS_DIR"/wBm_chung_2019_srr.list
OUTPUT_DIR="$READS_DIR"/wBm_chung_2019
```

#### Commands:
```{bash, eval = F}
while read SRR
do

	echo -e ""$SEQTK_BIN_DIR"/seqtk subseq "$OUTPUT_DIR"/"$SRR"_1.fastq "$OUTPUT_DIR"/"$SRR".target_reads.list > "$OUTPUT_DIR"/"$SRR"_1.subset.fastq" | qsub -P jdhotopp-lab -l mem_free=2G -N seqtk -wd "$OUTPUT_DIR"

	echo -e ""$SEQTK_BIN_DIR"/seqtk subseq "$OUTPUT_DIR"/"$SRR"_2.fastq "$OUTPUT_DIR"/"$SRR".target_reads.list > "$OUTPUT_DIR"/"$SRR"_2.subset.fastq" | qsub -P jdhotopp-lab -l mem_free=2G -N seqtk -wd "$OUTPUT_DIR"

done < "$SRR_ID_LIST"
```

## Align FASTQ files containing only Wolbachia-mapping reads to combined host and Wolbachia reference while disallowing splicing (2019-06-01) <a name="bashanalysis_realign"></a>

### Single-End FASTQs <a name="bashanalysis_realign_se"></a>
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
	echo -e ""$HISAT2_BIN_DIR"/hisat2 -p "$THREADS" -x "$REF_FNA" "$STRANDEDNESS" -U "$READS_DIR"/"$PREFIX"/"$SRR"_1.subset.fastq --no-spliced-alignment -X 1000 | "$SAMTOOLS_BIN_DIR"/samtools view -bhSo "$OUTPUT_DIR"/"$SRR".bam -" | qsub -q threaded.q  -pe thread "$THREADS" -P jdhotopp-lab -l mem_free=20G -N hisat2 -wd "$OUTPUT_DIR"
done < "$SRR_ID_LIST"
```

### Paired-End FASTQs <a name="bashanalysis_realign_pe"></a>
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

## Sort BAM files (2019-06-01) <a name="bashanalysis_sortbam"></a>

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

## Index BAM files (2019-06-01) <a name="bashanalysis_indexbam"></a>

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

## Quantify Wolbachia genes (2019-06-01) <a name="bashanalysis_quant"></a>

#### Input Sets:
```{bash, eval = F}
FEAT_TYPE="CDS"
ATTR_TYPE="ID"

## wOo Darby 2012
BAM_DIR="$READS_DIR"/wOo_darby_2012/bam
GFF3="$REFERENCES_DIR"/wOo.gff
STRANDED="no"

## wMel Darby 2014
BAM_DIR="$READS_DIR"/wMel_darby_2014/bam
GFF3="$REFERENCES_DIR"/wMel.gff
STRANDED="yes"

## wDi Luck 2014
BAM_DIR="$READS_DIR"/wDi_luck_2014/bam
GFF3="$REFERENCES_DIR"/wDi.gff
STRANDED="no"

## wDi Luck 2015
BAM_DIR="$READS_DIR"/wDi_luck_2015/bam
GFF3="$REFERENCES_DIR"/wDi.gff
STRANDED="no"

## wMel Gutzwiller 2015
BAM_DIR="$READS_DIR"/wMel_gutzwiller_2015/bam
GFF3="$REFERENCES_DIR"/wMel.gff
STRANDED="reverse"

## wBm Grote 2017
BAM_DIR="$READS_DIR"/wBm_grote_2017/bam
GFF3="$REFERENCES_DIR"/wBm.gff
STRANDED="no"

## wBm Chung 2019
BAM_DIR="$READS_DIR"/wBm_chung_2019/bam
GFF3="$REFERENCES_DIR"/wBm.gff
STRANDED="reverse"
```

#### Commands:
```{bash, eval = F}
mkdir -p "$WORKING_DIR"/fadu
for BAM in $(find "$BAM_DIR" -name "*sortedbyposition.bam")
do

	echo -e "export JULIA_DEPOT_PATH="$JULIA_DEPOT_PATH"\n"$JULIA_BIN_DIR"/julia "$FADU_BIN_DIR"/fadu.jl -g "$GFF3" -b "$BAM" -o "$WORKING_DIR"/fadu -s "$STRANDED" -f "$FEAT_TYPE" -a "$ATTR_TYPE"" | qsub -P jdhotopp-lab -l mem_free=5G -N fadu -wd "$WORKING_DIR"/fadu

done
```

## Create nucleotide coding sequence fasta files (2019-06-07) <a name="bashanalysis_createcdsfasta"></a>

#### Input Sets:
```{bash, eval = F}
FEAT_TYPE="CDS"
ATTR_TYPE="ID"

## wBm
FNA="$REFERENCES_DIR"/wBm.fna
GFF3="$REFERENCES_DIR"/wBm.gff

## wDi
FNA="$REFERENCES_DIR"/wDi.fna
GFF3="$REFERENCES_DIR"/wDi.gff

## wMel
FNA="$REFERENCES_DIR"/wMel.fna
GFF3="$REFERENCES_DIR"/wMel.gff

## wOo
FNA="$REFERENCES_DIR"/wOo.fna
GFF3="$REFERENCES_DIR"/wOo.gff
```

#### Commands:
```{bash, eval = F}
"$SCRIPTS_DIR"/createnuccdsfasta.sh -n "$FNA" -g "$GFF3" -f "$FEAT_TYPE" -i "$ATTR_TYPE" > "$(echo -e "$FNA" | sed "s/[.]fna$/.cds.fna/g")"
```

## Identify InterPro descriptions and GO terms for Wolbachia transcripts (2019-06-07) <a name="bashanalysis_interproscan"></a>

#### Input Sets:
```{bash, eval = F}
SEQ_TYPE="n"
THREADS=4

## wBm
CDS_FNA="$REFERENCES_DIR"/wBm.cds.fna
GFF3="$REFERENCES_DIR"/wBm.gff

## wDi
CDS_FNA="$REFERENCES_DIR"/wDi.cds.fna
GFF3="$REFERENCES_DIR"/wDi.gff

## wMel
CDS_FNA="$REFERENCES_DIR"/wMel.cds.fna
GFF3="$REFERENCES_DIR"/wMel.gff

## wOo
CDS_FNA="$REFERENCES_DIR"/wOo.cds.fna
GFF3="$REFERENCES_DIR"/wOo.gff
```

#### Commands:
```{bash, eval = F}
echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_PATH":"$LD_LIBRARY_PATH"\n"$INTERPROSCAN_BIN_DIR"/interproscan.sh -i "$CDS_FNA" -f tsv -o "$CDS_FNA".interproscan.tsv --seqtype "$SEQ_TYPE" --goterms --iprlookup" | qsub -P jdhotopp-lab -q threaded.q  -pe thread "$THREADS" -l mem_free=20G -N interproscan -wd "$(dirname "$CDS_FNA")"
```


## Convert InterProScan output to geneinfo file (2019-06-07) <a name="bashanalysis_iprscan2geneinfo"></a>


# Transcriptomics meta-analysis <a name="ranalysis"></a>

## Load R packages and view sessionInfo <a name="ranalysis_sessioninfo"></a>
```{R, eval = T}
library(cowplot)
library(edgeR)
library(emdbook)
library(FactoMineR)
library(dendextend)
library(ggdendro)
library(ggplot2)
library(gplots)
library(gridExtra)
library(gtools)
library(pvclust)
library(reshape2)
library(vegan)
library(WGCNA)

sessionInfo()
```

```{R, eval = F}
R version 3.5.1 (2018-07-02)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 7 x64 (build 7601) Service Pack 1

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] WGCNA_1.66            fastcluster_1.1.25    dynamicTreeCut_1.63-1 vegan_2.5-4           lattice_0.20-35       permute_0.9-4        
 [7] reshape2_1.4.3        pvclust_2.0-0         gtools_3.8.1          gridExtra_2.3         gplots_3.0.1.1        ggdendro_0.1-20      
[13] dendextend_1.9.0      FactoMineR_1.41       emdbook_1.3.11        edgeR_3.24.3          limma_3.38.3          cowplot_0.9.4        
[19] ggplot2_3.2.1        

loaded via a namespace (and not attached):
  [1] colorspace_1.4-1      class_7.3-14          modeltools_0.2-22     mclust_5.4.3          htmlTable_1.13.1      base64enc_0.1-3      
  [7] rstudioapi_0.10       flexmix_2.3-15        bit64_0.9-7           AnnotationDbi_1.44.0  mvtnorm_1.0-10        codetools_0.2-15     
 [13] splines_3.5.1         leaps_3.0             doParallel_1.0.14     impute_1.56.0         robustbase_0.93-4     knitr_1.22           
 [19] Formula_1.2-3         cluster_2.0.7-1       kernlab_0.9-27        GO.db_3.7.0           rrcov_1.4-7           compiler_3.5.1       
 [25] backports_1.1.3       assertthat_0.2.1      Matrix_1.2-14         lazyeval_0.2.2        acepack_1.4.1         htmltools_0.3.6      
 [31] tools_3.5.1           coda_0.19-2           gtable_0.3.0          glue_1.3.1            dplyr_0.8.3           Rcpp_1.0.1           
 [37] bbmle_1.0.20          Biobase_2.42.0        trimcluster_0.1-2.1   gdata_2.18.0          preprocessCore_1.44.0 nlme_3.1-137         
 [43] iterators_1.0.10      fpc_2.1-11.1          xfun_0.6              stringr_1.4.0         DEoptimR_1.0-8        MASS_7.3-50          
 [49] scales_1.0.0          parallel_3.5.1        RColorBrewer_1.1-2    yaml_2.2.0            memoise_1.1.0         rpart_4.1-13         
 [55] latticeExtra_0.6-28   stringi_1.4.3         RSQLite_2.1.1         S4Vectors_0.20.1      pcaPP_1.9-73          foreach_1.4.4        
 [61] checkmate_1.9.1       caTools_1.17.1.2      BiocGenerics_0.28.0   rlang_0.4.0           pkgconfig_2.0.2       prabclus_2.2-7       
 [67] bitops_1.0-6          matrixStats_0.54.0    purrr_0.3.2           htmlwidgets_1.3       bit_1.1-14            tidyselect_0.2.5     
 [73] robust_0.4-18         plyr_1.8.4            magrittr_1.5          R6_2.4.0              IRanges_2.16.0        Hmisc_4.2-0          
 [79] fit.models_0.5-14     DBI_1.0.0             pillar_1.3.1          whisker_0.3-2         foreign_0.8-70        withr_2.1.2          
 [85] mgcv_1.8-24           survival_2.42-3       scatterplot3d_0.3-41  nnet_7.3-12           tibble_2.1.1          crayon_1.3.4         
 [91] KernSmooth_2.23-15    viridis_0.5.1         locfit_1.5-9.1        grid_3.5.1            data.table_1.12.2     blob_1.1.1           
 [97] digest_0.6.18         diptest_0.75-7        flashClust_1.01-2     numDeriv_2016.8-1     stats4_3.5.1          munsell_0.5.0        
[103] viridisLite_0.3.0
```

## Load R functions <a name="ranalysis_loadfunctions"></a>
```{R, eval = T}
g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)
} 

get_dendro_structure <- function(result){
  structure <- hang.dendrogram(as.dendrogram(result$hclust))
  structure <- capture.output(str(structure))
  structure <- structure[grepl("leaf", structure)]
  structure <- as.numeric(as.character(substr(structure, regexpr("h=", structure ) + 3, regexpr("  )", structure))))
  return(structure)
}
get_dendro_data <- function(result){
  dendro.data <- dendro_data(result$hclust)
  dendro.data <- dendro.data$segments[which(dendro.data$segments$y == dendro.data$segments$yend),]
  for(i in 1:nrow(dendro.data)){
    dendro.data$minx[i] <- min(c(dendro.data$x[i], dendro.data$xend[i]))
  }
  dendro.data <- dendro.data[order(as.numeric(as.character(dendro.data$y)), as.numeric(as.character(dendro.data$minx))),]
  return(dendro.data)
}
get_dendro_bootstraps <- function(dendro_data){
  bootstrap.positions <- as.data.frame(matrix(nrow = length(dendro_data$y[duplicated(dendro_data$y)]),
                                              ncol = 2))
  for(i in 1:length(dendro_data$y[duplicated(dendro_data$y)])){
    dendro_data.subset <- dendro_data[which(dendro_data$y == dendro_data$y[duplicated(dendro_data$y)][i]),]
    bootstrap.positions[i,1] <- unique(dendro_data.subset$x)
    bootstrap.positions[i,2] <- unique(dendro_data.subset$y)
  }
  return(bootstrap.positions)
}
get_heatmap_separators <- function(vector){
  sep <- c()
  for(i in 2:length(unique(vector))){
    sep[length(sep) + 1] <- min(which(vector == unique(vector)[i])) - 1
  }
  return(sep)
}
formalize_study_names <- function(vector){
  vector <- gsub("wOo_darby_2012", "Darby et al., 2012",vector)
  vector <- gsub("wMel_darby_2014", "Darby et al., 2014",vector)
  vector <- gsub("wDi_luck_2014", "Luck et al., 2014",vector)
  vector <- gsub("wDi_luck_2015", "Luck et al., 2015",vector)
  vector <- gsub("wMel_gutzwiller_2015", "Gutzwiller et al., 2015",vector)
  vector <- gsub("wBm_grote_2017", "Grote et al., 2017",vector)
  vector <- gsub("wBm_chung_2019", "Chung et al., 2019",vector)
  vector <- factor(vector,levels=c("Darby et al., 2012",
                                   "Darby et al., 2014",
                                   "Luck et al., 2014",
                                   "Luck et al., 2015",
                                   "Gutzwiller et al., 2015",
                                   "Grote et al., 2017",
                                   "Chung et al., 2019"))
  return(vector)
}
find_soft_power <- function(sft){
  df <- as.data.frame(cbind(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2]))
  y <- -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
  dy <- diff(y) 
  softpower <- which(abs(diff(y)) < 0.05)[1]
  return(softpower)
}
eigengene_invert_id <- function(tpm.de.wgcna, mergedColors, mergedMEs){
  tpm.de.wgcna$invert <- T
  tpm.de.wgcna$module <- mergedColors
  for(i in 1:nrow(tpm.de.wgcna)){
    if(cor(t(tpm.de.wgcna[i,1:(ncol(tpm.de.wgcna) - 2)]), mergedMEs[,which(colnames(mergedMEs) == paste0("ME",tpm.de.wgcna$module[i]))], method = "pearson") > 0){
      tpm.de.wgcna$invert[i] <- F
    }
  }
  return(tpm.de.wgcna)
}
wgcna_heatmap_reorder <- function(tpm.de.wgcna){
  clusters <- as.data.frame(table(tpm.de.wgcna$module))
  clusters <- clusters[order(-clusters[,2]),1]
  
  tpm.de.wgcna.reordered <- as.data.frame(matrix(nrow = 0,
                                                 ncol = ncol(tpm.de.wgcna)))
  for(i in 1:length(clusters)){
    tpm.de.wgcna.reordered <- as.data.frame(rbind(tpm.de.wgcna.reordered,
                                                  tpm.de.wgcna[tpm.de.wgcna$module == clusters[i] & tpm.de.wgcna$invert == F,],
                                                  tpm.de.wgcna[tpm.de.wgcna$module == clusters[i] & tpm.de.wgcna$invert == T,]))
  }
  return(tpm.de.wgcna.reordered)
}
functionaltermenrichment <- function(genes, geneinfo){
  for(i in 1:ncol(geneinfo)){geneinfo[,i] <- as.character(geneinfo[,i])}
  geneinfo$interpro_description[which(is.na(geneinfo$interpro_description))] <- "No InterPro entry"
  geneinfo$go_biologicalprocess[which(is.na(geneinfo$go_biologicalprocess))] <- "No GO terms for biological process"
  geneinfo$go_cellularcomponent[which(is.na(geneinfo$go_cellularcomponent))] <- "No GO terms for cellular component"
  geneinfo$go_molecularfunction[which(is.na(geneinfo$go_molecularfunction))] <- "No GO terms for molecular function"
  
  functionalterms.list <- list(ipr=as.data.frame(table(unlist(strsplit(paste(geneinfo$interpro_description, collapse = "|"),  split = "[|]")))),
                               gobio=as.data.frame(table(unlist(strsplit(paste(geneinfo$go_biologicalprocess, collapse = "|"),  split = "[|]")))),
                               gocell=as.data.frame(table(unlist(strsplit(paste(geneinfo$go_cellularcomponent, collapse = "|"),  split = "[|]")))),
                               gomol=as.data.frame(table(unlist(strsplit(paste(geneinfo$go_molecularfunction, collapse = "|"),  split = "[|]")))))
  
  geneinfo.subset <- geneinfo[geneinfo$gene %in% genes,]
  term <- c()
  clusteroccurences <- c()
  genomeoccurences <- c()
  pvalue <- c()
  correctedpvalue <- c()
  oddsratio <- c()

  functionalterms.list.subset <- list(ipr=as.data.frame(table(unlist(strsplit(paste(geneinfo.subset$interpro_description, collapse = "|"),  split = "[|]")))),
                                      gobio=as.data.frame(table(unlist(strsplit(paste(geneinfo.subset$go_biologicalprocess, collapse = "|"),  split = "[|]")))),
                                      gocell=as.data.frame(table(unlist(strsplit(paste(geneinfo.subset$go_cellularcomponent, collapse = "|"),  split = "[|]")))),
                                      gomol=as.data.frame(table(unlist(strsplit(paste(geneinfo.subset$go_molecularfunction, collapse = "|"),  split = "[|]")))))
  
  for(i in 1:length(functionalterms.list)){
    for(j in 1:nrow(functionalterms.list[[i]])){
      freq.all <- functionalterms.list[[i]][j,2]
      freq.subset <- ifelse(functionalterms.list[[i]][j,1] %in% functionalterms.list.subset[[i]][,1],
                            functionalterms.list.subset[[i]][functionalterms.list.subset[[i]][,1] == as.character(functionalterms.list[[i]][j,1]),2],
                            0)
      genes.all <- nrow(geneinfo)
      genes.subset <- nrow(geneinfo.subset)

      fisherexact.matrix <- matrix(c(freq.subset, freq.all - freq.subset,
                                     genes.subset - freq.subset, genes.all - genes.subset - freq.all + freq.subset),
                                   nrow = 2,
                                   ncol = 2)
      fisher.test <- fisher.test(fisherexact.matrix)
      
      term[length(term) + 1] <- as.character(functionalterms.list[[i]][j,1])
      clusteroccurences[length(clusteroccurences) + 1] <- as.numeric(as.character(freq.subset))
      genomeoccurences[length(genomeoccurences) + 1] <- as.numeric(as.character(freq.all))
      pvalue[length(pvalue) + 1] <- as.numeric(as.character(fisher.test$p.value))
      correctedpvalue[length(correctedpvalue) + 1] <- p.adjust(as.numeric(as.character(fisher.test$p.value)), method = "fdr", n = nrow(functionalterms.list[[i]]))
      oddsratio[length(oddsratio) + 1] <- as.numeric(as.character(fisher.test$estimate))
    }
  }
  
  terms.df <- as.data.frame(cbind(term,
                                  clusteroccurences,
                                  genomeoccurences,
                                  pvalue,
                                  correctedpvalue,
                                  oddsratio))
  terms.df <- terms.df[order(as.numeric(as.character(terms.df$pvalue))),]
  return(terms.df)
}
```

## Set R inputs <a name="ranalysis_setinputs"></a>

R inputs from bash should be:

RAW_COUNTS.DIR = "$WORKING_DIR"/fadu
GFF3MAP.DIR = "$REFERENCES_DIR"

```{R}
RAW_COUNTS.DIR <- "Z:/EBMAL/mchung_dir/wolbachia_metatranscriptome_analysis/fadu/"
GFF3MAP.DIR <- 


maps.dir <- "Z:/EBMAL/mchung_dir/wolbachia_metatranscriptome_analysis/maps/"
sample_map.path <- "Z:/EBMAL/mchung_dir/wolbachia_metatranscriptome_analysis/study_sample_map.tsv.txt"
geneinfo.dir <- "Z:/EBMAL/mchung_dir/wolbachia_metatranscriptome_analysis/interproscan"
gff3map.dir <- "Z:/EBMAL/mchung_dir/wolbachia_metatranscriptome_analysis/references"
output.dir="Z:/EBMAL/mchung_dir/wolbachia_metatranscriptome_analysis/"


srr2wolbachiamappedreads.path <- "Z:/EBMAL/mchung_dir/wolbachia_metatranscriptome_analysis/srr_wolbachiamappedreads.tsv"
```

## Create counts and TPM tables <a name="ranalysis_counts"></a>

### Combine count files into a single table for each study <a name="ranalysis_counts_maketable"></a>

### Remove non-protein-coding genes from counts table <a name="ranalysis_counts_removenonproteincoding"></a>

### Remove non-protein-coding genes from counts table <a name="ranalysis_counts_removenonproteincoding"></a>
