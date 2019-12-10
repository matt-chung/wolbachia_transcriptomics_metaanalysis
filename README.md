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
	11. [Convert GFF to a mapping table](#bashanalysis_gff2map)
	12. [Create nucleotide coding sequence fasta files](#bashanalysis_createcdsfasta)
	13. [Identify InterPro descriptions and GO terms for Wolbachia transcripts](#bashanalysis_interproscan)
3. [Transcriptomics meta-analysis](#ranalysis)
	1. [Create counts table](#ranalysis_counts)
		1. [Set R inputs](#ranalysis_counts_setinputs)
		2. [View sessionInfo](#ranalysis_counts_sessioninfo)
		3. [Combine count files into a single table for each study](#ranalysis_counts_maketable)
	2. [Convert InterProScan output to geneinfo file](#ranalysis_iprscan2geneinfo)
	3. [Differential expression and WGCNA analyses](#ranalysis_de)
		1. [Set R inputs](#ranalysis_de_setinputs)
		2. [Load R packages and view sessionInfo](#ranalysis_de_sessioninfo)
		3. [Load R functions](#ranalysis_de_loadfunctions)
		4. [Process counts files and calculate TPM values](#ranalysis_de_counts_readcounts)
			1. [Read counts table into a list](#ranalysis_de_counts_readcounts)
			2. [Exclude non-protein-coding genes from counts table](#ranalysis_de_counts_removenonproteingenes)
			3. [Exclude genes without Agilent SureSelect probes from counts table in Chung et al 2019 study](#ranalysis_de_counts_removenoprobegenes)
			4. [Print the number of reads mapping to protein-coding genes for each sample in each study](#ranalysis_de_counts_countscolsum)
			5. [Create TPM table](#ranalysis_de_counts_tpm)
		5. [Create keys](#ranalysis_de_keys)
			1. [Sample keys](#ranalysis_de_keys_samples)
			2. [Size keys](#ranalysis_de_keys_size)
			3. [Heatmap log2TPM key](#ranalysis_de_keys_hm1)
			4. [Heatmap z-score log2TPM key](#ranalysis_de_keys_hm2)
		6. [Assess if samples have been sequenced to saturation using saturation curves](#ranalysis_de_keys_saturation)
		7. [Exclude samples and/or studies from differential expression analyses due to inadequate sequencing depth](#ranalysis_de_exclude)
		8. [Identify the number of genes in each study that meets the minimum CPM threshold](#ranalysis_de_exclude)
	9. [Identify differentially expressed genes (FDR < 0.05) in each study](#ranalysis_de_de)
	10. [Conduct principal component analyses](#ranalysis_de_pca)
	11. [Conduct hierarchical clustering analyses](#ranalysis_de_hc)
	12. [Create TPM and log2TPM heatmaps](#ranalysis_de_heatmap)
	13. [Choose studies for basic pairwise analyses or WGCNA](#ranalysis_de_selectanalysis)
	14. [Conduct pairwise analyses on studies with two biological groups](#ranalysis_de_pairwise)
		1. [Create z-score of log2TPM heatmaps with a column bar that shows up- and down-regulated genes](#ranalysis_de_pairwise_hm)
		2. [Conduct functional term enrichment analysis on up- and down-regulated gene subsets](#ranalysis_de_pairwise_fxnterm)
	15. [Conduct WGCNA on time-course studies](#ranalysis_de_wgcna)
		1. [Identify WGCNA soft power values](#ranalysis_de_wgcna_softpower)
		2. [Create WGCNA expression modules and merge similar modules](#ranalysis_de_wgcna_createmodules)
		3. [Plot WGCNA-derived expression modules](#ranalysis_de_wgcna_heatmap)
		4. [Conduct functional term enrichment analysis on WGCNA expression modules](#ranalysis_de_wgcna_fxnterm)



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

## Convert GFF to a mapping table (2019-06-07) <a name="bashanalysis_gff2map"></a>

#### Input Sets:
```{bash, eval = F}
## wBm
GFF3="$REFERENCES_DIR"/wBm.gff

## wDi
GFF3="$REFERENCES_DIR"/wDi.gff

## wMel
GFF3="$REFERENCES_DIR"/wMel.gff

## wOo
GFF3="$REFERENCES_DIR"/wOo.gff
```

#### Commands:
```{bash, eval = F}
"$R_BIN_DIR"/Rscript "$SCRIPTS_DIR"/gff3_to_map.R "$GFF_PATH"
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
SEQ_TYPE=n
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

# Transcriptomics meta-analysis <a name="ranalysis"></a>

## Create counts table <a name="ranalysis_counts"></a>

### Set R inputs <a name="ranalysis_setinputs"></a>

R inputs from bash should be:

SAMPLE_MAP.PATH = "$INPUTS_DIR"/study_sample_map.tsv.txt  
RAW_COUNTS.DIR = "$WORKING_DIR"/fadu  
OUTPUT.DIR = "$WORKING_DIR"  

```{R}
SAMPLE_MAP.PATH  <- "Z:/EBMAL/mchung_dir/wolbachia_metatranscriptome_analysis/study_sample_map.tsv.txt"
RAW_COUNTS.DIR <- "Z:/EBMAL/mchung_dir/wolbachia_metatranscriptome_analysis/fadu/"
OUTPUT.DIR <- "Z:/EBMAL/mchung_dir/wolbachia_metatranscriptome_analysis/"
```

### View sessionInfo <a name="ranalysis_sessioninfo"></a>

```{R, eval = T}
sessionInfo()
```

### Combine count files into a single table for each study <a name="ranalysis_counts_maketable"></a>

```{R}
sample_map <- read.delim(SAMPLE_MAP.PATH,header = T,sep = "\t")
counts.files <- list.files(RAW_COUNTS.DIR,
                           full.names = T,
                           pattern=".*counts.txt")

counts.list <- list(wOo_darby_2012="",
                    wMel_darby_2014="",
                    wDi_luck_2014="",
                    wDi_luck_2015="",
                    wMel_gutzwiller_2015="",
                    wBm_grote_2017="",
                    wBm_chung_2019="")


for(i in 1:length(counts.list)){
  sample_map.subset <- sample_map[sample_map$study == names(counts.list)[i],]
  samples <- unique(sample_map.subset$sample_identifier)
  
  genes <- read.delim(grep(sample_map.subset$sra_id[i],counts.files,value = T),header = T)[,1]
  counts.list[[i]] <- as.data.frame(matrix(0,
                                           nrow = length(genes),
                                           ncol = length(samples)))
  rownames(counts.list[[i]]) <- genes
  colnames(counts.list[[i]]) <- samples
                                    
  for(j in 1:ncol(counts.list[[i]])){
    srr <- sample_map.subset$sra_id[sample_map.subset$sample_identifier == colnames(counts.list[[i]])[j]]
    while(length(srr) > 0){
      raw_counts.file <- read.delim(paste0(RAW_COUNTS.DIR,"/",srr[1],".sortedbyposition.counts.txt"),
                                    header = T)
      raw_counts.file <- raw_counts.file[match(rownames(counts.list[[i]]),raw_counts.file[,1]),]
      counts.list[[i]][,j] <- counts.list[[i]][,j] + raw_counts.file$counts
      
      srr <- srr[-1]
    }
  }
  write.table(counts.list[[i]],
              paste0(OUTPUT.DIR,"/",names(counts.list)[i],"_counts.tsv"),
              row.names = T,
              col.names = T,
              quote = F,
              sep = "\t")
}
```

## Convert InterProScan output to geneinfo file (2019-06-07) <a name="ranalysis_iprscan2geneinfo"></a>

#### Input Sets:
```{bash, eval = F}
SEQ_TYPE=n
THREADS=4

## wBm
INTERPROSCAN_OUTPUT="$REFERENCES_DIR"/wBm.cds.fna.interproscan.tsv
COUNTS="$WORKING_DIR"/wBm_chung_2019_counts.tsv

## wDi
INTERPROSCAN_OUTPUT="$REFERENCES_DIR"/wDi.cds.fna.interproscan.tsv
COUNTS="$WORKING_DIR"/wDi_luck_2015_counts.tsv

## wMel
INTERPROSCAN_OUTPUT="$REFERENCES_DIR"/wMel.cds.fna.interproscan.tsv
COUNTS="$WORKING_DIR"/wMel_darby_2014_counts.tsv

## wOo
INTERPROSCAN_OUTPUT="$REFERENCES_DIR"/wOo.cds.fna.interproscan.tsv
COUNTS="$WORKING_DIR"/wOo_darby_2012_counts.tsv
```

#### Commands:
```{bash, eval = F}
"$R_BIN_DIR"/Rscript "$SCRIPTS_DIR"/interproscan2geneinfo.R "$INTERPROSCAN_OUTPUT" "$COUNTS"
```

## Differential expression and WGCNA analyses (2019-06-07) <a name="ranalysis_iprscan2geneinfo"></a>

### Set R inputs <a name="ranalysis_de_setinputs"></a>

R inputs from bash should be:

SAMPLE_MAP.PATH = "$INPUTS_DIR"/study_sample_map.tsv.txt  
COUNTS.DIR = "$WORKING_DIR"  
RAW_COUNTS.DIR = "$WORKING_DIR"/fadu  
GFF3MAP.DIR = "$REFERENCES_DIR"  
OUTPUT.DIR = "$WORKING_DIR"  

```{R}
SAMPLE_MAP.PATH  <- "Z:/EBMAL/mchung_dir/wolbachia_metatranscriptome_analysis/study_sample_map.tsv.txt"
COUNTS.DIR <- "Z:/EBMAL/mchung_dir/wolbachia_metatranscriptome_analysis/"
GENEINFO.DIR <- "Z:/EBMAL/mchung_dir/wolbachia_metatranscriptome_analysis/references"
GFF3MAP.DIR <- "Z:/EBMAL/mchung_dir/wolbachia_metatranscriptome_analysis/references"
RAW_COUNTS.DIR <- "Z:/EBMAL/mchung_dir/wolbachia_metatranscriptome_analysis/fadu"
OUTPUT.DIR <- "Z:/EBMAL/mchung_dir/wolbachia_metatranscriptome_analysis/"
```

### Load R packages and view sessionInfo <a name="ranalysis_de_sessioninfo"></a>
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

### Load R functions <a name="ranalysis_de_loadfunctions"></a>
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
### Process counts files and calculate TPM values  <a name="ranalysis_de_counts_readcounts"></a>
#### Read counts table into a list <a name="ranalysis_de_counts_readcounts"></a>

```{R, eval = T}
sample_map <- read.delim(SAMPLE_MAP.PATH)

counts.list <- list(wOo_darby_2012="",
                    wMel_darby_2014="",
                    wDi_luck_2014="",
                    wDi_luck_2015="",
                    wMel_gutzwiller_2015="",
                    wBm_grote_2017="",
                    wBm_chung_2019="")
for(i in 1:length(counts.list)){
  counts.list[[i]] <- read.delim(paste0(COUNTS.DIR,"/",names(counts.list)[i],"_counts.tsv"),
                                 row.names = 1)
}
```

#### Exclude non-protein-coding genes from counts table <a name="ranalysis_de_counts_removenonproteingenes"></a>

```{R, eval = T}
for(i in 1:length(counts.list)){
  species <- substr(names(counts.list)[i],1,regexpr("_",names(counts.list[i]))-1)
  if(species != "wDi"){
    map <- read.delim(paste0(GFF3MAP.DIR,"/",species,".gff.map"),
                      header = T,
                      sep = "\t")
    map <- map[map$gene_biotype == "protein_coding" & !is.na(map$gene_biotype),]
    map$ID <- gsub(".*[|]cds","cds",map$ID)
    counts.list[[i]] <- counts.list[[i]][rownames(counts.list[[i]]) %in% map$ID,]
  }
}
```
#### Exclude genes without Agilent SureSelect probes from counts table in Chung et al 2019 study <a name="ranalysis_de_counts_removenoprobegenes"></a>

```{R, eval = T}
exclude.genes <- c("cds28","cds31","cds63","cds64","cds74","cds99","cds118","cds127","cds135","cds168","cds169","cds170","cds171","cds185","cds240","cds245","cds246","cds307","cds309","cds321","cds322","cds324","cds328","cds354","cds359","cds360","cds363","cds371","cds372","cds380","cds426","cds427","cds431","cds444","cds452","cds453","cds466","cds490","cds491","cds492","cds493","cds500","cds543","cds573","cds581","cds607","cds608","cds614","cds622","cds635","cds636","cds638","cds643","cds657","cds658","cds733","cds734","cds737","cds746","cds750","cds774","cds775","cds784","cds785","cds787","cds805","cds816","cds817","cds822","cds828","cds851","cds853","cds870","cds871","cds872","cds874","cds886","cds909","cds910","cds932","cds957","cds963","cds965","cds978","cds979")

counts.list$wBm_chung_2019 <- counts.list$wBm_chung_2019[!(rownames(counts.list$wBm_chung_2019) %in% exclude.genes),]
```

#### Print the number of reads mapping to protein-coding genes for each sample in each study <a name="ranalysis_de_counts_countscolsum"></a>
```{R, eval = T}
for(i in 1:length(counts.list)){
  print(names(counts.list)[i])
  print(colSums(counts.list[[i]]))
}
```

```{R, eval = F}
[1] "wOo_darby_2012"
  darby_woo_adultmale_a   darby_woo_adultmale_b darby_woo_adultfemale_a darby_woo_adultfemale_b 
               26067.08                83692.20               146642.65                80723.69 
[1] "wMel_darby_2014"
  darby_wmel_doxy_a   darby_wmel_doxy_b   darby_wmel_doxy_c darby_wmel_nodoxy_c darby_wmel_nodoxy_a darby_wmel_nodoxy_b 
          120344.75            72768.46            73974.65            64291.14            42960.91            82329.35 
[1] "wDi_luck_2014"
         luck_wdi_L3_a          luck_wdi_L4_a   luck_wdi_adultmale_a luck_wdi_adultfemale_a          luck_wdi_MF_a 
                495.32                1238.56                3788.02               12840.55               27118.30 
[1] "wDi_luck_2015"
 luck_wdi_adultfemale_bodywall_a  luck_wdi_adultfemale_bodywall_b      luck_wdi_adultfemale_head_a      luck_wdi_adultfemale_head_b 
                         6691.13                          8436.58                          4045.55                          6019.27 
luck_wdi_adultfemale_intestine_a luck_wdi_adultfemale_intestine_b    luck_wdi_adultfemale_uterus_a    luck_wdi_adultfemale_uterus_b 
                           14.34                            11.18                          1785.09                           188.19 
   luck_wdi_adultmale_bodywall_a    luck_wdi_adultmale_bodywall_b   luck_wdi_adultmale_intestine_a      luck_wdi_adultmale_testis_a 
                        14444.09                         10506.72                            36.70                            13.34 
[1] "wMel_gutzwiller_2015"
                 gutzwiller_wmel_embryo0to2hr_a                  gutzwiller_wmel_embryo0to2hr_b                  gutzwiller_wmel_embryo2to4hr_a 
                                      114790.66                                        85945.65                                       110702.61 
                 gutzwiller_wmel_embryo2to4hr_b                gutzwiller_wmel_embryo10to12hr_a                gutzwiller_wmel_embryo10to12hr_b 
                                       53392.96                                       188273.14                                       134967.66 
               gutzwiller_wmel_embryo12to14hr_a                gutzwiller_wmel_embryo12to14hr_b                gutzwiller_wmel_embryo18to20hr_a 
                                      137773.04                                        40542.10                                       177066.53 
               gutzwiller_wmel_embryo18to20hr_b                gutzwiller_wmel_embryo22to24hr_a                gutzwiller_wmel_embryo22to24hr_b 
                                      148100.78                                       192603.87                                        18862.07 
                           gutzwiller_wmel_L1_a                            gutzwiller_wmel_L1_b                            gutzwiller_wmel_L2_a 
                                      153643.74                                       128799.29                                        65018.60 
                           gutzwiller_wmel_L2_b               gutzwiller_wmel_L3_12hrpostmolt_a               gutzwiller_wmel_L3_12hrpostmolt_b 
                                       55155.07                                        72087.91                                        45411.64 
        gutzwiller_wmel_L3_darkbluegutPS_1to2_a         gutzwiller_wmel_L3_darkbluegutPS_1to2_b        gutzwiller_wmel_L3_lightbluegutPS_3to6_a 
                                      122280.34                                        90161.56                                       152077.10 
       gutzwiller_wmel_L3_lightbluegutPS_3to6_b            gutzwiller_wmel_L3_cleargutPS_7to9_a            gutzwiller_wmel_L3_cleargutPS_7to9_b 
                                      112893.91                                       241950.76                                        76163.54 
                          gutzwiller_wmel_WPP_a                           gutzwiller_wmel_WPP_b             gutzwiller_wmel_pupae_12hrpostWPP_a 
                                      178456.20                                       115149.63                                       162560.29 
            gutzwiller_wmel_pupae_12hrpostWPP_b             gutzwiller_wmel_pupae_24hrpostWPP_a             gutzwiller_wmel_pupae_24hrpostWPP_b 
                                       92488.90                                       112776.08                                        48500.17 
           gutzwiller_wmel_pupae_2dayspostWPP_a            gutzwiller_wmel_pupae_2dayspostWPP_b            gutzwiller_wmel_pupae_3dayspostWPP_a 
                                      106883.82                                        60639.46                                       203627.16 
           gutzwiller_wmel_pupae_3dayspostWPP_b            gutzwiller_wmel_pupae_4dayspostWPP_a            gutzwiller_wmel_pupae_4dayspostWPP_b 
                                       76254.69                                       105137.01                                        95285.08 
 gutzwiller_wmel_adultfemale_1dayposteclosion_a  gutzwiller_wmel_adultfemale_1dayposteclosion_b  gutzwiller_wmel_adultfemale_5dayposteclosion_a 
                                      254985.87                                       140491.93                                       221664.24 
 gutzwiller_wmel_adultfemale_5dayposteclosion_b  gutzwiller_wmel_adultfemale_5dayposteclosion_c  gutzwiller_wmel_adultfemale_5dayposteclosion_d 
                                      111503.67                                       103197.92                                       180473.06 
 gutzwiller_wmel_adultfemale_5dayposteclosion_f gutzwiller_wmel_adultfemale_30dayposteclosion_a gutzwiller_wmel_adultfemale_30dayposteclosion_b 
                                      195493.19                                       341687.94                                       244695.15 
   gutzwiller_wmel_adultmale_1dayposteclosion_a    gutzwiller_wmel_adultmale_1dayposteclosion_b    gutzwiller_wmel_adultmale_5dayposteclosion_a 
                                      225874.52                                       137847.49                                       417013.04 
   gutzwiller_wmel_adultmale_5dayposteclosion_b    gutzwiller_wmel_adultmale_5dayposteclosion_c    gutzwiller_wmel_adultmale_5dayposteclosion_d 
                                      136930.51                                       328285.55                                       222262.97 
   gutzwiller_wmel_adultmale_5dayposteclosion_i    gutzwiller_wmel_adultmale_5dayposteclosion_j   gutzwiller_wmel_adultmale_30dayposteclosion_a 
                                       88320.71                                        93761.47                                       325975.98 
  gutzwiller_wmel_adultmale_30dayposteclosion_b 
                                      321449.11 
[1] "wBm_grote_2017"
                  grote_wbm_L4_a                   grote_wbm_L4_b   grote_wbm_30dpi_immaturemale_a   grote_wbm_30dpi_immaturemale_b 
                        22374.59                          6178.84                         12843.80                          9650.02 
grote_wbm_30dpi_immaturefemale_a grote_wbm_30dpi_immaturefemale_b   grote_wbm_42dpi_immaturemale_a   grote_wbm_42dpi_immaturemale_b 
                         9080.46                         14064.82                         13019.66                          8636.27 
grote_wbm_42dpi_immaturefemale_a grote_wbm_42dpi_immaturefemale_b     grote_wbm_120dpi_adultmale_a     grote_wbm_120dpi_adultmale_b 
                         6721.62                         13460.39                         16643.58                          8373.13 
  grote_wbm_120dpi_adultfemale_a   grote_wbm_120dpi_adultfemale_b 
                        31289.42                         19151.90 
[1] "wBm_chung_2019"
               chung_wbm_vector18hpi_a                chung_wbm_vector18hpi_b                 chung_wbm_vector4dpi_a 
                             481351.25                              625028.83                              240693.73 
                chung_wbm_vector4dpi_b                 chung_wbm_vector8dpi_a                 chung_wbm_vector8dpi_b 
                             777934.24                             3218498.46                             3987678.06 
                chung_wbm_mammal1dpi_a                 chung_wbm_mammal1dpi_b                 chung_wbm_mammal2dpi_a 
                             159809.14                              103367.39                               53272.07 
                chung_wbm_mammal2dpi_b                 chung_wbm_mammal3dpi_a                 chung_wbm_mammal3dpi_b 
                             123961.25                              109225.39                              122560.27 
                chung_wbm_mammal4dpi_a                 chung_wbm_mammal4dpi_b                 chung_wbm_mammal8dpi_a 
                             120910.97                              121522.47                              284943.40 
                chung_wbm_mammal8dpi_b   chung_wbm_mammal20dpi_immaturemale_a   chung_wbm_mammal20dpi_immaturemale_b 
                             221120.86                              108117.16                               78688.09 
chung_wbm_mammal24dpi_immaturefemale_a chung_wbm_mammal24dpi_immaturefemale_b                  chung_wbm_adultmale_a 
                           35276890.78                              192206.48                               23790.09 
                 chung_wbm_adultmale_b                chung_wbm_adultfemale_a                chung_wbm_adultfemale_b 
                              43446.49                            18840207.11                               10561.16 
               chung_wbm_adultfemale_c                chung_wbm_adultfemale_d                     chung_wbm_embryo_a 
                              13947.67                            24661817.10                                8525.69 
                    chung_wbm_embryo_b                 chung_wbm_immatureMF_a                 chung_wbm_immatureMF_b 
                              14704.98                               65069.68                               57976.58 
                  chung_wbm_matureMF_a                   chung_wbm_matureMF_b 
                              46554.09                               51856.71 

```

#### Create TPM table <a name="ranalysis_de_counts_tpm"></a>

```{R, eval = T}
tpm.list <- counts.list

for(i in 1:length(tpm.list)){
  sample_map.subset <- sample_map[sample_map$study == names(counts.list)[i] & sample_map$sample_identifier %in% colnames(counts.list[[i]]),]
  srr <- sample_map.subset$sra_id[sample_map.subset$sample_identifier == colnames(counts.list[[i]])[1]]
  raw_counts.file <- read.delim(paste0(RAW_COUNTS.DIR,"/",srr[1],".sortedbyposition.counts.txt"),
                                      header = T)
  genelength <- raw_counts.file$uniq_len[match(rownames(tpm.list[[i]]),raw_counts.file$featureID)]
  
  for(j in 1:ncol(tpm.list[[i]])){
    tpm.list[[i]][,j] <- tpm.list[[i]][,j]/genelength
    tpm.list[[i]][,j] <- tpm.list[[i]][,j]/(sum(tpm.list[[i]][,j])/1000000)
  }
  
  write.table(tpm.list[[i]],
            paste0(OUTPUT.DIR,"/",names(tpm.list)[i],"_tpm.tsv"),
            row.names = T,
            col.names = T,
            quote = F,
            sep = "\t")
}
```

### Create keys <a name="ranalysis_de_keys"></a>

#### Sample keys <a name="ranalysis_de_keys_samples"></a>
```{R,fig.height=2,fig.width=10}
for(i in 1:length(counts.list)){
  sample_map.subset <- sample_map[sample_map$study == names(counts.list)[i] & sample_map$sample_identifier %in% colnames(counts.list[[i]]),]
  groups <- as.data.frame(cbind(as.character(sample_map.subset$sample_identifier),
                                as.character(sample_map.subset$sample_group),
                                as.character(sample_map.subset$color)))
  groups <- unique(groups)
  
  groups[,1] <- factor(groups[,1],levels=groups[,1])
  groups[,2] <- factor(groups[,2],levels=unique(groups[,2]))
  groups[,3] <- factor(groups[,3],levels=unique(groups[,3]))
  
  legend.plot <- ggplot()+
    geom_point(aes(x=groups[,2], y=seq(1,length(groups[,2]),1), color = groups[,2]), show.legend = T, size = 2)+
    scale_color_manual(values = as.character(unique(groups[,3])))+
    guides(colour = guide_legend(title = "Samples", title.position = "top", nrow = 6))+
    theme_bw()+
    theme(legend.position="top",legend.title.align=0.5)
  
  sample.legend <- g_legend(legend.plot)
  
  pdf(paste0(OUTPUT.DIR,"/",names(counts.list)[i],"_samplekey.pdf"),
      height=2,
      width=10)
  grid.arrange(sample.legend)
  dev.off()
  
  grid.arrange(sample.legend)
}
```

##### wOo, Darby et al 2012
![Image description](/images/wOo_darby_2012_samplekey.png)

##### wMel, Darby et al 2014
![Image description](/images/wMel_darby_2014_samplekey.png)

##### wDi, Luck et al 2014
![Image description](/images/wDi_luck_2014_samplekey.png)

##### wDi, Luck et al 2015
![Image description](/images/wDi_luck_2015_samplekey.png)

##### wMel, Gutzwiller et al 2015
![Image description](/images/wMel_gutzwiller_2015_samplekey.png)

##### wBm, Grote et al 2017
![Image description](/images/wBm_grote_2017_samplekey.png)

##### wBm, Chung et al 2019
![Image description](/images/wBm_chung_2019_samplekey.png)

#### Size key  <a name="ranalysis_de_keys_size"></a>
```{R,fig.height=2,fig.width=10}
max_lib_sizes <- c()
min_lib_sizes <- c()
for(i in 1:length(counts.list)){
  max_lib_sizes <- c(max_lib_sizes,max(colSums(counts.list[[i]])))
  min_lib_sizes <- c(min_lib_sizes,min(colSums(counts.list[[i]])))
}
size_limits <- c(1,10^ceiling(log10(max(max_lib_sizes))))
size_breaks <- lseq(10^floor(log10(min(min_lib_sizes))),10^ceiling(log10(max(max_lib_sizes))),8)
size.plot <- ggplot()+
  geom_point(aes(x=seq(1,length(colSums(counts.list[[i]]))), y=colSums(counts.list[[i]]),size=colSums(counts.list[[i]])))+
  scale_size_continuous(limits=size_limits,breaks=size_breaks ,trans = "log")+
  guides(size = guide_legend(title = "Reads Mapping to Wolbachia Genes", title.position = "top", nrow = 2))+
  theme_bw()+
  theme(legend.position="top",legend.title.align=0.5)
  
size.legend <- g_legend(size.plot)

pdf(paste0(OUTPUT.DIR,"/sizekey.pdf"),
    height=2,
    width=10)
grid.arrange(size.legend)
dev.off()

png(paste0(OUTPUT.DIR,"/sizekey.png"),
    height=2,
    width=10,
    units = "in",res=300)
grid.arrange(size.legend)
dev.off()

grid.arrange(size.legend)
```
![Image description](/images/sizekey.png)

#### Heatmap log2TPM key <a name="ranalysis_de_keys_hm1"></a>
```{R,fig.height=2,fig.width=7}
hmcol1 <- colorRampPalette(c("#FFFFFF","#FDE0D2","#FABAA1","#F69274", "#F16B4E","#EF3D2D","#590A16","#A51E23","#650C16"))(20)
hmcol2 <- colorRampPalette(c("navyblue","white","firebrick3"))(12)

hmcolor1.plot <- ggplot() + 
  geom_raster(aes(x=seq(0,10,0.5), y=seq(0,10,0.5), fill = seq(0,10,0.5)))+
  scale_fill_gradientn(name = "log2TPM",
                       colours=hmcol1,
                       breaks=c(0,5,10))+
  theme(legend.position="bottom")+
  guides(fill = guide_colorbar(title.position = "top"))

sample.legend <- g_legend(hmcolor1.plot)

pdf(paste0(OUTPUT.DIR,"/hmcolor1key.pdf"),
    height=2,
    width=7)
grid.arrange(sample.legend)
dev.off()

png(paste0(OUTPUT.DIR,"/hmcolor1key.png"),
    height=2,
    width=7,
    units = "in",res=300)
grid.arrange(sample.legend)
dev.off()

grid.arrange(sample.legend)
```
![Image description](/images/hmcolor1key.png)

#### Heatmap z-score log2TPM key <a name="ranalysis_de_keys_hm2"></a>
```{R,fig.height=2,fig.width=7}
hmcolor2.plot <- ggplot() + 
  geom_raster(aes(x=seq(-3,3,0.5), y=seq(-3,3,0.5), fill = seq(-3,3,0.5)))+
  scale_fill_gradientn(name = "z-score of log2TPM",
                       colours=hmcol2,
                       breaks=c(-3,0,3))+
  theme(legend.position="bottom")+
  guides(fill = guide_colorbar(title.position = "top"))

heat.legend <- g_legend(hmcolor2.plot)

pdf(paste0(OUTPUT.DIR,"/hmcolor2key.pdf"),
    height=2,
    width=7)
grid.arrange(heat.legend)
dev.off()

png(paste0(OUTPUT.DIR,"/hmcolor2key.png"),
    height=2,
    width=7,
    units = "in",res=300)
grid.arrange(heat.legend)
dev.off()

grid.arrange(heat.legend)
```
![Image description](/images/hmcolor2key.png)

### Assess if samples have been sequenced to saturation using saturation curves <a name="ranalysis_de_keys_saturation"></a>
```{R, fig.height = 5, fig.width = 5}
for(i in 1:length(tpm.list)){
  rarefy.counts <- ceiling(counts.list[[i]])
  rarefy.counts <- rarefy.counts[,colSums(rarefy.counts) != 0,]
  
  sample_map.subset <- sample_map[sample_map$study == names(counts.list)[i] & sample_map$sample_identifier %in% colnames(counts.list[[i]]),]
  
  groups <- as.data.frame(cbind(as.character(sample_map.subset$sample_identifier),
                                as.character(sample_map.subset$sample_group),
                                as.character(sample_map.subset$color)))
  groups <- unique(groups)
  
  groups[,1] <- factor(groups[,1],levels=groups[,1])
  groups[,2] <- factor(groups[,2],levels=unique(groups[,2]))
  groups[,3] <- factor(groups[,3],levels=unique(groups[,3]))
  

  raremax <- round(min(rowSums(t(rarefy.counts))[rowSums(t(rarefy.counts)) != 0]),0)
  srare <- rarefy(t(rarefy.counts),raremax)
  raremax <- ifelse(raremax < 10, 10, raremax)
  
  rarefy.raw.df <- rarecurve(t(rarefy.counts), step = round(raremax/10,0), sample = raremax)
  
  rarefy.df <- as.data.frame(matrix(nrow = 0,
                                    ncol = 5))
  rarefy.points.df <- rarefy.df
  for(j in 1:length(rarefy.raw.df)){
    steps <- as.numeric(gsub("N","",names(rarefy.raw.df[[j]])))
    detected_genes <- as.numeric(rarefy.raw.df[[j]])
    rarefy.df <- as.data.frame(rbind(rarefy.df,
                                     cbind(as.numeric(steps),
                                           as.numeric(detected_genes),
                                           as.character(groups[j,1]),
                                           as.character(groups[j,2]),
                                           groups[j,3])))
    rarefy.points.df <- as.data.frame(rbind(rarefy.points.df,
                                     cbind(as.numeric(max(steps)),
                                           as.numeric(max(detected_genes)),
                                           as.character(groups[j,1]),
                                           as.character(groups[j,2]),
                                           groups[j,3])))
    
  }
  rarefy.plot <- ggplot()+
    geom_line(mapping=aes(x=as.numeric(as.character(rarefy.df[,1])), y=as.numeric(as.character(rarefy.df[,2])),group=rarefy.df[,3],color=rarefy.df[,4]))+
    #geom_point(mapping=aes(x=as.numeric(as.character(rarefy.df[,1])), y=as.numeric(as.character(rarefy.df[,2])),group=rarefy.df[,3],color=rarefy.df[,4]))+
    geom_point(mapping=aes(x=as.numeric(as.character(rarefy.points.df[,1])), y=as.numeric(as.character(rarefy.points.df[,2])),group=rarefy.points.df[,3],color=rarefy.points.df[,4]),size = 3)+
    guides(colour = F,shape = F)+
    scale_color_manual(values = levels(groups[,3]))+
    labs(x="reads mapping to protein-coding genes", y="genes detected", color = "Sample")+
    #coord_cartesian(xlim=c(0,100000))+
    theme_bw()
  
  pdf(paste0(OUTPUT.DIR,"/",names(counts.list)[i],"_rarefaction.pdf"),
      height=5,
      width=5)
  print(rarefy.plot)
  dev.off()
  
  png(paste0(OUTPUT.DIR,"/",names(counts.list)[i],"_rarefaction.png"),
      height=5,
      width=5,
      units = "in",res=300)
  print(rarefy.plot)
  dev.off()
  
  print(rarefy.plot)
}
```

##### wOo, Darby et al 2012
![Image description](/images/wOo_darby_2012_rarefaction.png)

##### wMel, Darby et al 2014
![Image description](/images/wMel_darby_2014_rarefaction.png)

##### wDi, Luck et al 2014
![Image description](/images/wDi_luck_2014_rarefaction.png)

##### wDi, Luck et al 2015
![Image description](/images/wDi_luck_2015_rarefaction.png)

##### wMel, Gutzwiller et al 2015
![Image description](/images/wMel_gutzwiller_2015_rarefaction.png)

##### wBm, Grote et al 2017
![Image description](/images/wBm_grote_2017_rarefaction.png)

##### wBm, Chung et al 2019
![Image description](/images/wBm_chung_2019_rarefaction.png)

### Exclude samples and/or studies from differential expression analyses due to inadequate sequencing depth <a name="ranalysis_de_exclude"></a>

All samples from Luck et al 2014 were excluded due to having inadequate replicates and sequencing depth for differential expression analyses. All samples with <2000 gene counts from Luck et al 2015 are excluded due to have inadequate sequencing depth (refer to saturation curves).

```{R}
counts.kept.list <- counts.list
tpm.kept.list <- tpm.list

counts.edgeR.list <- counts.list
tpm.edgeR.list <- tpm.list

counts.edgeR.list[[which(names(counts.edgeR.list) == "wDi_luck_2014")]] <- NULL
tpm.edgeR.list[[which(names(tpm.edgeR.list) == "wDi_luck_2014")]] <- NULL

counts.kept.list[[which(names(counts.kept.list) == "wDi_luck_2015")]] <- counts.kept.list[[which(names(counts.kept.list) == "wDi_luck_2015")]][colSums(counts.kept.list[[which(names(counts.kept.list) == "wDi_luck_2015")]]) > 2000]
tpm.kept.list[[which(names(tpm.kept.list) == "wDi_luck_2015")]] <- tpm.kept.list[[which(names(tpm.kept.list) == "wDi_luck_2015")]][,colnames(tpm.kept.list[[which(names(tpm.kept.list) == "wDi_luck_2015")]]) %in% colnames(counts.kept.list[[which(names(counts.kept.list) == "wDi_luck_2015")]])]

counts.edgeR.list[[which(names(counts.edgeR.list) == "wDi_luck_2015")]] <- counts.edgeR.list[[which(names(counts.edgeR.list) == "wDi_luck_2015")]][colSums(counts.edgeR.list[[which(names(counts.edgeR.list) == "wDi_luck_2015")]]) > 2000]
tpm.edgeR.list[[which(names(tpm.edgeR.list) == "wDi_luck_2015")]] <- tpm.edgeR.list[[which(names(tpm.edgeR.list) == "wDi_luck_2015")]][,colnames(tpm.edgeR.list[[which(names(tpm.edgeR.list) == "wDi_luck_2015")]]) %in% colnames(counts.edgeR.list[[which(names(counts.edgeR.list) == "wDi_luck_2015")]])]

qlf.list <- counts.edgeR.list
```

### Identify the number of genes in each study that meets the minimum CPM threshold <a name="ranalysis_de_exclude"></a>

```{R}
for(i in 1:length(counts.kept.list)){
  sample_map.subset <- sample_map[sample_map$study == names(counts.list)[i] & sample_map$sample_identifier %in% colnames(counts.list[[i]]),]
  groups <- as.data.frame(cbind(as.character(sample_map.subset$sample_identifier),
                                as.character(sample_map.subset$sample_group),
                                as.character(sample_map.subset$color)))
  groups <- groups[groups[,1] %in% colnames(counts.kept.list[[i]]),]
  groups <- unique(groups)

  groups[,1] <- factor(groups[,1],levels=groups[,1])
  groups[,2] <- factor(groups[,2],levels=unique(groups[,2]))
  groups[,3] <- factor(groups[,3],levels=unique(groups[,3]))

  cpm.cutoff <- 5/min(colSums(counts.kept.list[[i]])) * 1000000

  y <- DGEList(counts = counts.kept.list[[i]], group = groups[,2])
  y <- calcNormFactors(y)
  keep <- rowSums(cpm(y) >= cpm.cutoff) >= min(table(groups[,2]))

  keep.df <- as.data.frame(table(keep))
  print(paste0(names(counts.kept.list)[i],": ",keep.df[keep.df[,1] == F,2]," genes excluded using edgeR filter, ",keep.df[keep.df[,1] == T,2], " genes kept" ))

  counts.kept.list[[which(names(counts.kept.list) == names(counts.list)[i])]] <- counts.kept.list[[which(names(counts.kept.list) == names(counts.list)[i])]][keep,]
  tpm.kept.list[[which(names(tpm.kept.list) == names(tpm.list)[i])]] <- tpm.kept.list[[which(names(tpm.kept.list) == names(tpm.list)[i])]][keep,]
}
```

```{R, eval = F}
[1] "wOo_darby_2012: 91 genes excluded using edgeR filter, 560 genes kept"
[1] "wMel_darby_2014: 219 genes excluded using edgeR filter, 867 genes kept"
[1] "wDi_luck_2014: 853 genes excluded using edgeR filter, 18 genes kept"
[1] "wDi_luck_2015: 546 genes excluded using edgeR filter, 325 genes kept"
[1] "wMel_gutzwiller_2015: 305 genes excluded using edgeR filter, 781 genes kept"
[1] "wBm_grote_2017: 414 genes excluded using edgeR filter, 425 genes kept"
[1] "wBm_chung_2019: 206 genes excluded using edgeR filter, 548 genes kept"
```

### Identify differentially expressed genes (FDR < 0.05) in each study <a name="ranalysis_de_de"></a>

```{R}
FDRcutoff <- 0.05

for(i in 1:length(counts.edgeR.list)){
  sample_map.subset <- sample_map[sample_map$study == names(counts.edgeR.list)[i] & sample_map$sample_identifier %in% colnames(counts.edgeR.list[[i]]),]
  groups <- as.data.frame(cbind(as.character(sample_map.subset$sample_identifier),
                                as.character(sample_map.subset$sample_group),
                                as.character(sample_map.subset$color)))
  groups <- unique(groups)

  groups[,1] <- factor(groups[,1],levels=groups[,1])
  groups[,2] <- factor(groups[,2],levels=unique(groups[,2]))
  groups[,3] <- factor(groups[,3],levels=unique(groups[,3]))

  cpm.cutoff <- 5/min(colSums(counts.edgeR.list[[i]])) * 1000000

  y <- DGEList(counts = counts.edgeR.list[[i]], group = groups[,2])
  y <- calcNormFactors(y)
  keep <- rowSums(cpm(y) >= cpm.cutoff) >= min(table(groups[,2]))

  y <- y[keep, , keep.lib.sizes = F]
  design <- model.matrix(~groups[,2])
  y <- estimateDisp(y , design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = 2:ncol(fit))

  FDR <- as.data.frame(p.adjust(qlf$table$PValue, method="BH"))
  rownames(FDR) <- rownames(qlf$table)
  FDR.genes <- rownames(FDR[FDR[,1] < FDRcutoff, , drop = F])
  print(paste0(names(counts.edgeR.list)[i],": ",length(FDR.genes)," DE genes" ))
  
  qlf.list[[i]] <- as.data.frame(qlf$table)
  counts.edgeR.list[[i]] <- counts.edgeR.list[[i]][rownames(counts.edgeR.list[[i]]) %in% FDR.genes,]
  tpm.edgeR.list[[i]] <- tpm.edgeR.list[[i]][rownames(tpm.edgeR.list[[i]]) %in% FDR.genes,]
  
  write.table(tpm.edgeR.list[[i]],
            paste0(OUTPUT.DIR,"/",names(tpm.edgeR.list)[i],"_de_tpm.tsv"),
            row.names = T,
            col.names = T,
            quote = F,
            sep = "\t")
}
```

```{R, eval = F}
[1] "wOo_darby_2012: 0 DE genes"
[1] "wMel_darby_2014: 120 DE genes"
[1] "wDi_luck_2015: 36 DE genes"
[1] "wMel_gutzwiller_2015: 473 DE genes"
[1] "wBm_grote_2017: 94 DE genes"
[1] "wBm_chung_2019: 373 DE genes"
```

### Conduct principal component analyses <a name="ranalysis_de_pca"></a>
#### All genes <a name="ranalysis_de_pca_all"></a>

```{R, fig.height=5,fig.width=5}
for(i in 1:length(tpm.list)){ 
	sample_map.subset <- sample_map[sample_map$study == names(counts.list)[i] & sample_map$sample_identifier %in% colnames(counts.list[[i]]),]

	groups <- as.data.frame(cbind(as.character(sample_map.subset$sample_identifier),
	                            as.character(sample_map.subset$sample_group),
	                            as.character(sample_map.subset$color)))
	groups <- unique(groups)

	groups[,1] <- factor(groups[,1],levels=groups[,1])
	groups[,2] <- factor(groups[,2],levels=unique(groups[,2]))
	groups[,3] <- factor(groups[,3],levels=unique(groups[,3]))

	pca.df <- t(scale(t(log2(tpm.list[[i]] + 1))))
	pca.df[is.na(pca.df)] <- 0
	pca.df <- pca.df[rowSums(pca.df == 0) != ncol(pca.df),]
	pca <- PCA(as.data.frame(scale(t(pca.df))), graph = FALSE, ncp = ncol(counts.list[[i]]) - 1)

	pca.plot <- ggplot()+
	  geom_point(aes(x=pca$ind$coord[,1], y=pca$ind$coord[,2], color = groups[,2],size = colSums(counts.list[[i]])))+
	  labs(col = "Samples", size = "Reads Mapped\nto Features", 
	       x = paste0("PC1 (", round(pca$eig[1,2],1), "%)"), 
	       y = paste0("PC2 (", round(pca$eig[2,2],1), "%)"))+
	  guides(color = F,size = F)+
	  # guides(colour = guide_legend(ncol = 2))+
	  scale_color_manual(values = levels(groups[,3]))+
	  scale_size_continuous(limits=size_limits,breaks=size_breaks,trans = "log")+
	  theme_bw()

	pdf(paste0(OUTPUT.DIR,"/",names(counts.list)[i],"_pca_allgenes.pdf"),
	  height=5,
	  width=5)
	print(pca.plot)
	dev.off()
	
	png(paste0(OUTPUT.DIR,"/",names(counts.list)[i],"_pca_allgenes.png"),
	  height=5,
	  width=5,
	  units = "in",res=300)
	print(pca.plot)
	dev.off()

	print(pca.plot)
}
```

##### wOo, Darby et al 2012
![Image description](/images/wOo_darby_2012_pca_allgenes.png)

##### wMel, Darby et al 2014
![Image description](/images/wMel_darby_2014_pca_allgenes.png)

##### wDi, Luck et al 2014
![Image description](/images/wDi_luck_2014_pca_allgenes.png)

##### wDi, Luck et al 2015
![Image description](/images/wDi_luck_2015_pca_allgenes.png)

##### wMel, Gutzwiller et al 2015
![Image description](/images/wMel_gutzwiller_2015_pca_allgenes.png)

##### wBm, Grote et al 2017
![Image description](/images/wBm_grote_2017_pca_allgenes.png)

##### wBm, Chung et al 2019
![Image description](/images/wBm_chung_2019_pca_allgenes.png)

#### Kept genes <a name="ranalysis_de_pca_kept"></a>

```{R, fig.height=5,fig.width=5}
for(i in 1:length(tpm.kept.list)){
  sample_map.subset <- sample_map[sample_map$study == names(tpm.kept.list)[i] & sample_map$sample_identifier %in% colnames(tpm.kept.list[[i]]),]
  groups <- as.data.frame(cbind(as.character(sample_map.subset$sample_identifier),
                                as.character(sample_map.subset$sample_group),
                                as.character(sample_map.subset$color)))
  groups <- unique(groups)
  
  groups[,1] <- factor(groups[,1],levels=groups[,1])
  groups[,2] <- factor(groups[,2],levels=unique(groups[,2]))
  groups[,3] <- factor(groups[,3],levels=unique(groups[,3]))
  
  if(nrow(tpm.kept.list[[i]]) > 0){
    pca.df <- t(scale(t(log2(tpm.kept.list[[i]] + 1))))
    pca.df <- pca.df[rowSums(pca.df == 0) != ncol(pca.df),]
    pca <- PCA(as.data.frame(scale(t(pca.df))), graph = FALSE, ncp = ncol(tpm.kept.list[[i]]) - 1)
    
    pca.plot <- ggplot()+
      geom_point(aes(x=pca$ind$coord[,1], y=pca$ind$coord[,2], color = groups[,2],size = colSums(counts.kept.list[[i]])))+
      labs(col = "Samples", size = "Reads Mapped\nto Features", 
           x = paste0("PC1 (", round(pca$eig[1,2],1), "%)"), 
           y = paste0("PC2 (", round(pca$eig[2,2],1), "%)"))+
      guides(color = F,size = F)+
      # guides(colour = guide_legend(ncol = 2))+
      scale_color_manual(values = levels(groups[,3]))+
      scale_size_continuous(limits=size_limits,breaks=size_breaks,trans = "log")+
      theme_bw()
    
    pdf(paste0(OUTPUT.DIR,"/",names(counts.list)[i],"_pca_keptgenes.pdf"),
          height=5,
          width=5)
      print(pca.plot)
    dev.off()
      
    png(paste0(OUTPUT.DIR,"/",names(counts.list)[i],"_pca_keptgenes.png"),
  	    height=5,
  	    width=5,
  	    units = "in",res=300)
  	  print(pca.plot)
  	dev.off()
    
    print(pca.plot)
  }
}
```

##### wOo, Darby et al 2012
![Image description](/images/wOo_darby_2012_pca_keptgenes.png)

##### wMel, Darby et al 2014
![Image description](/images/wMel_darby_2014_pca_keptgenes.png)

##### wDi, Luck et al 2014
![Image description](/images/wDi_luck_2014_pca_keptgenes.png)

##### wDi, Luck et al 2015
![Image description](/images/wDi_luck_2015_pca_keptgenes.png)

##### wMel, Gutzwiller et al 2015
![Image description](/images/wMel_gutzwiller_2015_pca_keptgenes.png)

##### wBm, Grote et al 2017
![Image description](/images/wBm_grote_2017_pca_keptgenes.png)

##### wBm, Chung et al 2019
![Image description](/images/wBm_chung_2019_pca_keptgenes.png)

#### Differentially expressed genes <a name="ranalysis_de_pca_de"></a>

```{R, fig.height=5,fig.width=5}
for(i in 1:length(tpm.edgeR.list)){
  sample_map.subset <- sample_map[sample_map$study == names(tpm.edgeR.list)[i]& sample_map$sample_identifier %in% colnames(tpm.edgeR.list[[i]]),]
  groups <- as.data.frame(cbind(as.character(sample_map.subset$sample_identifier),
                                as.character(sample_map.subset$sample_group),
                                as.character(sample_map.subset$color)))
  groups <- unique(groups)
  
  groups[,1] <- factor(groups[,1],levels=groups[,1])
  groups[,2] <- factor(groups[,2],levels=unique(groups[,2]))
  groups[,3] <- factor(groups[,3],levels=unique(groups[,3]))
  
  if(nrow(tpm.edgeR.list[[i]]) > 0){
    pca.df <- t(scale(t(log2(tpm.edgeR.list[[i]] + 1))))
    pca.df <- pca.df[rowSums(pca.df == 0) != ncol(pca.df),]
    pca <- PCA(as.data.frame(scale(t(pca.df))), graph = FALSE, ncp = ncol(tpm.edgeR.list[[i]]) - 1)
    
    pca.plot <- ggplot()+
      geom_point(aes(x=pca$ind$coord[,1], y=pca$ind$coord[,2], color = groups[,2],size = colSums(counts.edgeR.list[[i]])))+
      labs(col = "Samples", size = "Reads Mapped\nto Features", 
           x = paste0("PC1 (", round(pca$eig[1,2],1), "%)"), 
           y = paste0("PC2 (", round(pca$eig[2,2],1), "%)"))+
      guides(color = F,size = F)+
      # guides(colour = guide_legend(ncol = 2))+
      scale_color_manual(values = levels(groups[,3]))+
      scale_size_continuous(limits=size_limits,breaks=size_breaks,trans = "log")+
      theme_bw()
    
    pdf(paste0(OUTPUT.DIR,"/",names(counts.list)[i],"_pca_degenes.pdf"),
        height=5,
        width=5)
      print(pca.plot)
    dev.off()
    
    png(paste0(OUTPUT.DIR,"/",names(counts.list)[i],"_pca_degenes.png"),
        height=5,
        width=5,
        units = "in",res=300)
      print(pca.plot)
    dev.off()
    
    print(pca.plot)
  }
}
```

##### wMel, Darby et al 2014
![Image description](/images/wMel_darby_2014_pca_degenes.png)

##### wDi, Luck et al 2015
![Image description](/images/wDi_luck_2015_pca_degenes.png)

##### wMel, Gutzwiller et al 2015
![Image description](/images/wMel_gutzwiller_2015_pca_degenes.png)

##### wBm, Grote et al 2017
![Image description](/images/wBm_grote_2017_pca_degenes.png)

##### wBm, Chung et al 2019
![Image description](/images/wBm_chung_2019_pca_degenes.png)

### Conduct hierarchical clustering analyses <a name="ranalysis_de_hc"></a>
#### All genes <a name="ranalysis_de_hc_all"></a>

```{R, fig.height=8, fig.width=11}
for(i in 1:length(tpm.list)){
  sample_map.subset <- sample_map[sample_map$study == names(tpm.list)[i]& sample_map$sample_identifier %in% colnames(tpm.list[[i]]),]
  groups <- as.data.frame(cbind(as.character(sample_map.subset$sample_identifier),
                                as.character(sample_map.subset$sample_group),
                                as.character(sample_map.subset$color)))
  groups <- unique(groups)
  
  groups[,1] <- factor(groups[,1],levels=groups[,1])
  groups[,2] <- factor(groups[,2],levels=unique(groups[,2]))
  groups[,3] <- factor(groups[,3],levels=unique(groups[,3]))
  
  dendrogram <- as.data.frame(t(scale(t(log2(tpm.list[[i]] + 1)))))
  dendrogram <- dendrogram[rowSums(is.na(dendrogram)) == 0,]
  result <- pvclust(dendrogram, method.dist="cor", method.hclust="average", nboot=100, quiet=T)
  
  structure <- get_dendro_structure(result)
  dendro.data <- get_dendro_data(result)
  bootstrap.positions <- get_dendro_bootstraps(dendro.data)
    
  points.df <- as.data.frame(cbind(seq(1,length(structure),1),
                                   structure))
  dendrogroups <- groups[,2][result$hclust$order]
  dendrocol <- groups[,3][result$hclust$order]
  dendrosize <- colSums(counts.list[[i]])[result$hclust$order]
  
  for(j in 1:length(unique(groups[,1]))){
    result$hclust$labels <- gsub(groups[,1][j],groups[,2][j],result$hclust$labels)
  }
  
  dendrogram.plot <- ggdendrogram(hang.dendrogram(as.dendrogram(result$hclust)), theme_dendro = T)+
    geom_point(aes(x=seq(1,length(structure)), y = structure, color = dendrogroups, size = dendrosize))+
    scale_color_manual(values = levels(groups[,3]))+
    labs(x = "", y = "", col = "Samples", size = "Reads Mapped\nto Features")+
    guides(colour = guide_legend(ncol = 2))+
    theme_minimal()+      
    scale_size_continuous(limits=size_limits,breaks=size_breaks,trans = "log")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1),
          axis.text.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position="none")
  
  for(j in 1:length(result$edges$bp)){
    text <- round(result$edges$bp[j] * 100,0)
    dendrogram.plot <- dendrogram.plot + annotate("text", label = text, x=bootstrap.positions[j,1] + 0.4, y=bootstrap.positions[j,2] + 0.03, size = 3)
  }
  
  pdf(paste0(OUTPUT.DIR,"/",names(counts.list)[i],"_dendrogram_allgenes.pdf"),
      width = 11,
      height = 8)
    print(dendrogram.plot)
  dev.off()
  
  png(paste0(OUTPUT.DIR,"/",names(counts.list)[i],"_dendrogram_allgenes.png"),
      width = 11,
      height = 8,
      units = "in",res=300)
    print(dendrogram.plot)
  dev.off()
  
  print(dendrogram.plot)
}
```

##### wOo, Darby et al 2012
![Image description](/images/wOo_darby_2012_dendrogram_allgenes.png)

##### wMel, Darby et al 2014
![Image description](/images/wMel_darby_2014_dendrogram_allgenes.png)

##### wDi, Luck et al 2014
![Image description](/images/wDi_luck_2014_dendrogram_allgenes.png)

##### wDi, Luck et al 2015
![Image description](/images/wDi_luck_2015_dendrogram_allgenes.png)

##### wMel, Gutzwiller et al 2015
![Image description](/images/wMel_gutzwiller_2015_dendrogram_allgenes.png)

##### wBm, Grote et al 2017
![Image description](/images/wBm_grote_2017_dendrogram_allgenes.png)

##### wBm, Chung et al 2019
![Image description](/images/wBm_chung_2019_dendrogram_allgenes.png)

#### Kept genes <a name="ranalysis_de_hc_kept"></a>

```{R, fig.height=8,fig.width=11}
for(i in 1:length(tpm.kept.list)){
  sample_map.subset <- sample_map[sample_map$study == names(tpm.kept.list)[i] & sample_map$sample_identifier %in% colnames(tpm.kept.list[[i]]),]
  groups <- as.data.frame(cbind(as.character(sample_map.subset$sample_identifier),
                                as.character(sample_map.subset$sample_group),
                                as.character(sample_map.subset$color)))
  groups <- unique(groups)
  
  groups[,1] <- factor(groups[,1],levels=groups[,1])
  groups[,2] <- factor(groups[,2],levels=unique(groups[,2]))
  groups[,3] <- factor(groups[,3],levels=unique(groups[,3]))
  
  dendrogram <- as.data.frame(t(scale(t(log2(tpm.kept.list[[i]]) + 1))))
  if(nrow(dendrogram) > 0){
    result <- pvclust(dendrogram, method.dist="cor", method.hclust="average", nboot=100, quiet=T)
    
    structure <- get_dendro_structure(result)
    dendro.data <- get_dendro_data(result)
    bootstrap.positions <- get_dendro_bootstraps(dendro.data)
      
    points.df <- as.data.frame(cbind(seq(1,length(structure),1),
                                     structure))
    dendrogroups <- groups[,2][result$hclust$order]
    dendrocol <- groups[,3][result$hclust$order]
    dendrosize <- colSums(counts.kept.list[[i]])[result$hclust$order]
    
    for(j in 1:length(unique(groups[,1]))){
      result$hclust$labels <- gsub(groups[,1][j],groups[,2][j],result$hclust$labels)
    }
    
    dendrogram.plot <- ggdendrogram(hang.dendrogram(as.dendrogram(result$hclust)), theme_dendro = T)+
      geom_point(aes(x=seq(1,length(structure)), y = structure, color = dendrogroups, size = dendrosize))+
      scale_color_manual(values = levels(groups[,3]))+
      labs(x = "", y = "", col = "Samples", size = "Reads Mapped\nto Features")+
      guides(colour = guide_legend(ncol = 2))+
      theme_minimal()+
      scale_size_continuous(limits=size_limits,breaks=size_breaks,trans = "log")+
      theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1),
            axis.text.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position="none")
    
    for(j in 1:length(result$edges$bp)){
      text <- round(result$edges$bp[j] * 100,0)
      dendrogram.plot <- dendrogram.plot + annotate("text", label = text, x=bootstrap.positions[j,1] + 0.4, y=bootstrap.positions[j,2] + 0.03, size = 3)
    }
    pdf(paste0(OUTPUT.DIR,"/",names(counts.kept.list)[i],"_dendrogram_keptgenes.pdf"),
        width = 11,
        height = 8)
      print(dendrogram.plot)
    dev.off()
    
    png(paste0(OUTPUT.DIR,"/",names(counts.kept.list)[i],"_dendrogram_keptgenes.png"),
        width = 11,
        height = 8,
        units = "in",res=300)
      print(dendrogram.plot)
    dev.off()
    
    print(dendrogram.plot)
  }
}
```

##### wOo, Darby et al 2012
![Image description](/images/wOo_darby_2012_dendrogram_keptgenes.png)

##### wMel, Darby et al 2014
![Image description](/images/wMel_darby_2014_dendrogram_keptgenes.png)

##### wDi, Luck et al 2014
![Image description](/images/wDi_luck_2014_dendrogram_keptgenes.png)

##### wDi, Luck et al 2015
![Image description](/images/wDi_luck_2015_dendrogram_keptgenes.png)

##### wMel, Gutzwiller et al 2015
![Image description](/images/wMel_gutzwiller_2015_dendrogram_keptgenes.png)

##### wBm, Grote et al 2017
![Image description](/images/wBm_grote_2017_dendrogram_keptgenes.png)

##### wBm, Chung et al 2019
![Image description](/images/wBm_chung_2019_dendrogram_keptgenes.png)

#### Differentially expressed genes <a name="ranalysis_de_hc_de"></a>

```{R,fig.height=8,fig.width=11}
for(i in 1:length(tpm.edgeR.list)){
  sample_map.subset <- sample_map[sample_map$study == names(tpm.edgeR.list)[i] & sample_map$sample_identifier %in% colnames(tpm.edgeR.list[[i]]),]
  groups <- as.data.frame(cbind(as.character(sample_map.subset$sample_identifier),
                                as.character(sample_map.subset$sample_group),
                                as.character(sample_map.subset$color)))
  groups <- unique(groups)
  
  groups[,1] <- factor(groups[,1],levels=groups[,1])
  groups[,2] <- factor(groups[,2],levels=unique(groups[,2]))
  groups[,3] <- factor(groups[,3],levels=unique(groups[,3]))
  
  dendrogram <- as.data.frame(t(scale(t(log2(tpm.edgeR.list[[i]]) + 1))))
  if(nrow(dendrogram) > 0){
    result <- pvclust(dendrogram, method.dist="cor", method.hclust="average", nboot=100, quiet=T)
    
    structure <- get_dendro_structure(result)
    dendro.data <- get_dendro_data(result)
    bootstrap.positions <- get_dendro_bootstraps(dendro.data)
      
    points.df <- as.data.frame(cbind(seq(1,length(structure),1),
                                     structure))
    dendrogroups <- groups[,2][result$hclust$order]
    dendrocol <- groups[,3][result$hclust$order]
    dendrosize <- colSums(counts.edgeR.list[[i]])[result$hclust$order]
    
    for(j in 1:length(unique(groups[,1]))){
      result$hclust$labels <- gsub(groups[,1][j],groups[,2][j],result$hclust$labels)
    }
    
    dendrogram.plot <- ggdendrogram(hang.dendrogram(as.dendrogram(result$hclust)), theme_dendro = T)+
      geom_point(aes(x=seq(1,length(structure)), y = structure, color = dendrogroups, size = dendrosize))+
      scale_color_manual(values = levels(groups[,3]))+
      labs(x = "", y = "", col = "Samples", size = "Reads Mapped\nto Features")+
      guides(colour = guide_legend(ncol = 2))+
      theme_minimal()+
      scale_size_continuous(limits=size_limits,breaks=size_breaks,trans = "log")+
      theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1),
            axis.text.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position="none")
    
    for(j in 1:length(result$edges$bp)){
      text <- round(result$edges$bp[j] * 100,0)
      dendrogram.plot <- dendrogram.plot + annotate("text", label = text, x=bootstrap.positions[j,1] + 0.4, y=bootstrap.positions[j,2] + 0.03, size = 3)
    }
    
    pdf(paste0(OUTPUT.DIR,"/",names(counts.edgeR.list)[i],"_dendrogram_degenes.pdf"),
        width = 11,
        height = 8)
      print(dendrogram.plot)
    dev.off()
    
    png(paste0(OUTPUT.DIR,"/",names(counts.edgeR.list)[i],"_dendrogram_degenes.png"),
        width = 11,
        height = 8,
        units = "in",res=300)
      print(dendrogram.plot)
    dev.off()
    
    print(dendrogram.plot)
  }
}
```

##### wMel, Darby et al 2014
![Image description](/images/wMel_darby_2014_dendrogram_keptgenes.png)

##### wDi, Luck et al 2015
![Image description](/images/wDi_luck_2015_dendrogram_keptgenes.png)

##### wMel, Gutzwiller et al 2015
![Image description](/images/wMel_gutzwiller_2015_dendrogram_keptgenes.png)

##### wBm, Grote et al 2017
![Image description](/images/wBm_grote_2017_dendrogram_keptgenes.png)

##### wBm, Chung et al 2019
![Image description](/images/wBm_chung_2019_dendrogram_keptgenes.png)

### Create TPM and log2TPM heatmaps <a name="ranalysis_de_heatmap"></a>
#### All genes <a name="ranalysis_de_heatmap_all"></a>

```{R,fig.height=10,fig.width=8}
for(i in 1:length(tpm.list)){
  sample_map.subset <- sample_map[sample_map$study == names(tpm.list)[i] & sample_map$sample_identifier %in% colnames(tpm.list[[i]]),]
  groups <- as.data.frame(cbind(as.character(sample_map.subset$sample_identifier),
                                as.character(sample_map.subset$sample_group),
                                as.character(sample_map.subset$color)))
  groups <- unique(groups)
  
  groups[,1] <- factor(groups[,1],levels=groups[,1])
  groups[,2] <- factor(groups[,2],levels=unique(groups[,2]))
  groups[,3] <- factor(groups[,3],levels=unique(groups[,3]))
  
  heatmap1.df <- as.data.frame(log2(tpm.list[[i]]+1))
  heatmap2.df <- as.data.frame(t(scale(t(log2(tpm.list[[i]] + 1)))))
  
  colcol <- as.character(groups[,3])
  
  #rowdendro <- pvclust(t(heatmap1.df), method.dist="cor", method.hclust="average", nboot=100)
  coldendro <- pvclust(heatmap2.df, method.dist="cor", method.hclust="average", nboot=100, quiet=T)
  
  heatmap2.df[is.na(heatmap2.df)] <- 0
  
  pdf(paste0(OUTPUT.DIR,"/",names(tpm.list)[i],"_heatmap_allgenes.pdf"),
      width = 8,
      height = 10)
  heatmap.2(as.matrix(heatmap1.df),
            col=hmcol1,
            trace="none",
            labRow=vector(mode = "character", length = nrow(heatmap1.df)),
            labCol=groups[,2],
            #Rowv = as.dendrogram(rowdendro),
            Colv = as.dendrogram(coldendro),
            ColSideColors=colcol,
            lwid=c(1,4),
            lhei = c(1,5),
            breaks = seq(0,10,by=.5),
            symkey=F,
            dendrogram = "both")
  
  heatmap.2(as.matrix(heatmap2.df),
            col=hmcol2,
            trace="none",
            labRow=vector(mode = "character", length = nrow(heatmap2.df)),
            labCol=groups[,2],
            #Rowv = as.dendrogram(rowdendro),
            Colv = as.dendrogram(coldendro),
            ColSideColors=colcol,
            lwid=c(1,4),
            lhei = c(1,5),
            breaks = seq(-3,3,by=.5),
            symkey=F,
            dendrogram = "both")
  dev.off()
  
  png(paste0(OUTPUT.DIR,"/",names(tpm.list)[i],"_heatmap_allgenes.png"),
      width = 8,
      height = 10,
      units = "in",res=300)
  heatmap.2(as.matrix(heatmap1.df),
            col=hmcol1,
            trace="none",
            labRow=vector(mode = "character", length = nrow(heatmap1.df)),
            labCol=groups[,2],
            #Rowv = as.dendrogram(rowdendro),
            Colv = as.dendrogram(coldendro),
            ColSideColors=colcol,
            lwid=c(1,4),
            lhei = c(1,5),
            breaks = seq(0,10,by=.5),
            symkey=F,
            dendrogram = "both")
  
  heatmap.2(as.matrix(heatmap2.df),
            col=hmcol2,
            trace="none",
            labRow=vector(mode = "character", length = nrow(heatmap2.df)),
            labCol=groups[,2],
            #Rowv = as.dendrogram(rowdendro),
            Colv = as.dendrogram(coldendro),
            ColSideColors=colcol,
            lwid=c(1,4),
            lhei = c(1,5),
            breaks = seq(-3,3,by=.5),
            symkey=F,
            dendrogram = "both")
  dev.off()
  
  heatmap.2(as.matrix(heatmap1.df),
            col=hmcol1,
            trace="none",
            labRow=vector(mode = "character", length = nrow(heatmap1.df)),
            labCol=groups[,2],
           #Rowv = as.dendrogram(rowdendro),
            Colv = as.dendrogram(coldendro),
            ColSideColors=colcol,
            lwid=c(1,4),
            lhei = c(1,5),
            breaks = seq(0,10,by=.5),
            symkey=F,
            dendrogram = "both")
  
  heatmap.2(as.matrix(heatmap2.df),
            col=hmcol2,
            trace="none",
            labRow=vector(mode = "character", length = nrow(heatmap2.df)),
            labCol=groups[,2],
            #Rowv = as.dendrogram(rowdendro),
            Colv = as.dendrogram(coldendro),
            ColSideColors=colcol,
            lwid=c(1,4),
            lhei = c(1,5),
            breaks = seq(-3,3,by=.5),
            symkey=F,
            dendrogram = "both")
}
```

##### wOo, Darby et al 2012
![Image description](/images/wOo_darby_2012_heatmap_allgenes.png)

##### wMel, Darby et al 2014
![Image description](/images/wMel_darby_2014_heatmap_allgenes.png)

##### wDi, Luck et al 2014
![Image description](/images/wDi_luck_2014_heatmap_allgenes.png)

##### wDi, Luck et al 2015
![Image description](/images/wDi_luck_2015_heatmap_allgenes.png)

##### wMel, Gutzwiller et al 2015
![Image description](/images/wMel_gutzwiller_2015_heatmap_allgenes.png)

##### wBm, Grote et al 2017
![Image description](/images/wBm_grote_2017_heatmap_allgenes.png)

##### wBm, Chung et al 2019
![Image description](/images/wBm_chung_2019_heatmap_allgenes.png)

#### Kept genes <a name="ranalysis_de_heatmap_kept"></a>

```{R, fig.height=10, fig.width=8}
for(i in 1:length(tpm.kept.list)){
  sample_map.subset <- sample_map[sample_map$study == names(tpm.kept.list)[i] & sample_map$sample_identifier %in% colnames(tpm.kept.list[[i]]),]
  groups <- as.data.frame(cbind(as.character(sample_map.subset$sample_identifier),
                                as.character(sample_map.subset$sample_group),
                                as.character(sample_map.subset$color)))
  groups <- unique(groups)
  
  groups[,1] <- factor(groups[,1],levels=groups[,1])
  groups[,2] <- factor(groups[,2],levels=unique(groups[,2]))
  groups[,3] <- factor(groups[,3],levels=unique(groups[,3]))
  
  heatmap1.df <- as.data.frame(log2(tpm.kept.list[[i]]+1))
  heatmap2.df <- as.data.frame(t(scale(t(log2(tpm.kept.list[[i]] + 1)))))
  
  colcol <- as.character(groups[,3])
  
  #rowdendro <- pvclust(t(heatmap1.df), method.dist="cor", method.hclust="average", nboot=100)
  coldendro <- pvclust(heatmap2.df, method.dist="cor", method.hclust="average", nboot=100, quiet=T)
  
  heatmap2.df[is.na(heatmap2.df)] <- 0
  
  pdf(paste0(OUTPUT.DIR,"/",names(tpm.kept.list)[i],"_heatmap_keptgenes.pdf"),
      width = 8,
      height = 10)
  heatmap.2(as.matrix(heatmap1.df),
            col=hmcol1,
            trace="none",
            labRow=vector(mode = "character", length = nrow(heatmap1.df)),
            labCol=groups[,2],
            #Rowv = as.dendrogram(rowdendro),
            Colv = as.dendrogram(coldendro),
            ColSideColors=colcol,
            lwid=c(1,4),
            lhei = c(1,5),
            breaks = seq(0,10,by=.5),
            symkey=F,
            dendrogram = "both")
  
  heatmap.2(as.matrix(heatmap2.df),
            col=hmcol2,
            trace="none",
            labRow=vector(mode = "character", length = nrow(heatmap2.df)),
            labCol=groups[,2],
            #Rowv = as.dendrogram(rowdendro),
            Colv = as.dendrogram(coldendro),
            ColSideColors=colcol,
            lwid=c(1,4),
            lhei = c(1,5),
            breaks = seq(-3,3,by=.5),
            symkey=F,
            dendrogram = "both")
  dev.off()
  
  png(paste0(OUTPUT.DIR,"/",names(tpm.kept.list)[i],"_heatmap_keptgenes.png"),
      width = 8,
      height = 10,
      units = "in",res=300)
  heatmap.2(as.matrix(heatmap1.df),
            col=hmcol1,
            trace="none",
            labRow=vector(mode = "character", length = nrow(heatmap1.df)),
            labCol=groups[,2],
            #Rowv = as.dendrogram(rowdendro),
            Colv = as.dendrogram(coldendro),
            ColSideColors=colcol,
            lwid=c(1,4),
            lhei = c(1,5),
            breaks = seq(0,10,by=.5),
            symkey=F,
            dendrogram = "both")
  
  heatmap.2(as.matrix(heatmap2.df),
            col=hmcol2,
            trace="none",
            labRow=vector(mode = "character", length = nrow(heatmap2.df)),
            labCol=groups[,2],
            #Rowv = as.dendrogram(rowdendro),
            Colv = as.dendrogram(coldendro),
            ColSideColors=colcol,
            lwid=c(1,4),
            lhei = c(1,5),
            breaks = seq(-3,3,by=.5),
            symkey=F,
            dendrogram = "both")
  dev.off()
  
  heatmap.2(as.matrix(heatmap1.df),
            col=hmcol1,
            trace="none",
            labRow=vector(mode = "character", length = nrow(heatmap1.df)),
            labCol=groups[,2],
           #Rowv = as.dendrogram(rowdendro),
            Colv = as.dendrogram(coldendro),
            ColSideColors=colcol,
            lwid=c(1,4),
            lhei = c(1,5),
            breaks = seq(0,10,by=.5),
            symkey=F,
            dendrogram = "both")
  
  heatmap.2(as.matrix(heatmap2.df),
            col=hmcol2,
            trace="none",
            labRow=vector(mode = "character", length = nrow(heatmap2.df)),
            labCol=groups[,2],
            #Rowv = as.dendrogram(rowdendro),
            Colv = as.dendrogram(coldendro),
            ColSideColors=colcol,
            lwid=c(1,4),
            lhei = c(1,5),
            breaks = seq(-3,3,by=.5),
            symkey=F,
            dendrogram = "both")
}
```

##### wOo, Darby et al 2012
![Image description](/images/wOo_darby_2012_heatmap_keptgenes.png)

##### wMel, Darby et al 2014
![Image description](/images/wMel_darby_2014_heatmap_keptgenes.png)

##### wDi, Luck et al 2014
![Image description](/images/wDi_luck_2014_heatmap_keptgenes.png)

##### wDi, Luck et al 2015
![Image description](/images/wDi_luck_2015_heatmap_keptgenes.png)

##### wMel, Gutzwiller et al 2015
![Image description](/images/wMel_gutzwiller_2015_heatmap_keptgenes.png)

##### wBm, Grote et al 2017
![Image description](/images/wBm_grote_2017_heatmap_keptgenes.png)

##### wBm, Chung et al 2019
![Image description](/images/wBm_chung_2019_heatmap_keptgenes.png)

#### Differentially expressed genes <a name="ranalysis_de_heatmap_de"></a>

```{R, fig.height=10, fig.width=8}
for(i in 1:length(tpm.edgeR.list)){
  sample_map.subset <- sample_map[sample_map$study == names(tpm.edgeR.list)[i] & sample_map$sample_identifier %in% colnames(tpm.edgeR.list[[i]]),]
  groups <- as.data.frame(cbind(as.character(sample_map.subset$sample_identifier),
                                as.character(sample_map.subset$sample_group),
                                as.character(sample_map.subset$color)))
  groups <- unique(groups)
  
  groups[,1] <- factor(groups[,1],levels=groups[,1])
  groups[,2] <- factor(groups[,2],levels=unique(groups[,2]))
  groups[,3] <- factor(groups[,3],levels=unique(groups[,3]))
  
  heatmap1.df <- as.data.frame(log2(tpm.edgeR.list[[i]]+1))
  heatmap2.df <- as.data.frame(t(scale(t(log2(tpm.edgeR.list[[i]] + 1)))))
  
  colcol <- as.character(groups[,3])
  
  if(nrow(tpm.edgeR.list[[i]]) > 0){
    #rowdendro <- pvclust(t(heatmap1.df), method.dist="cor", method.hclust="average", nboot=100)
    coldendro <- pvclust(heatmap2.df, method.dist="cor", method.hclust="average", nboot=100, quiet=T)
    
    heatmap2.df[is.na(heatmap2.df)] <- 0
    
    pdf(paste0(OUTPUT.DIR,"/",names(tpm.edgeR.list)[i],"_heatmap_degenes.pdf"),
        width = 8,
        height = 10)
    heatmap.2(as.matrix(heatmap1.df),
              col=hmcol1,
              trace="none",
              labRow=vector(mode = "character", length = nrow(heatmap1.df)),
              labCol=groups[,2],
              #Rowv = as.dendrogram(rowdendro),
              Colv = as.dendrogram(coldendro),
              ColSideColors=colcol,
              lwid=c(1,4),
              lhei = c(1,5),
              breaks = seq(0,10,by=.5),
              symkey=F,
              dendrogram = "both")
    
    heatmap.2(as.matrix(heatmap2.df),
              col=hmcol2,
              trace="none",
              labRow=vector(mode = "character", length = nrow(heatmap2.df)),
              labCol=groups[,2],
              #Rowv = as.dendrogram(rowdendro),
              Colv = as.dendrogram(coldendro),
              ColSideColors=colcol,
              lwid=c(1,4),
              lhei = c(1,5),
              breaks = seq(-3,3,by=.5),
              symkey=F,
              dendrogram = "both")
    dev.off()
    
    png(paste0(OUTPUT.DIR,"/",names(tpm.edgeR.list)[i],"_heatmap_keptgenes.png"),
        width = 8,
        height = 10,
        units = "in",res=300)
    heatmap.2(as.matrix(heatmap1.df),
              col=hmcol1,
              trace="none",
              labRow=vector(mode = "character", length = nrow(heatmap1.df)),
              labCol=groups[,2],
              #Rowv = as.dendrogram(rowdendro),
              Colv = as.dendrogram(coldendro),
              ColSideColors=colcol,
              lwid=c(1,4),
              lhei = c(1,5),
              breaks = seq(0,10,by=.5),
              symkey=F,
              dendrogram = "both")
    
    heatmap.2(as.matrix(heatmap2.df),
              col=hmcol2,
              trace="none",
              labRow=vector(mode = "character", length = nrow(heatmap2.df)),
              labCol=groups[,2],
              #Rowv = as.dendrogram(rowdendro),
              Colv = as.dendrogram(coldendro),
              ColSideColors=colcol,
              lwid=c(1,4),
              lhei = c(1,5),
              breaks = seq(-3,3,by=.5),
              symkey=F,
              dendrogram = "both")
    dev.off()
    
    heatmap.2(as.matrix(heatmap1.df),
              col=hmcol1,
              trace="none",
              labRow=vector(mode = "character", length = nrow(heatmap1.df)),
              labCol=groups[,2],
             #Rowv = as.dendrogram(rowdendro),
              Colv = as.dendrogram(coldendro),
              ColSideColors=colcol,
              lwid=c(1,4),
              lhei = c(1,5),
              breaks = seq(0,10,by=.5),
              symkey=F,
              dendrogram = "both")
    
    heatmap.2(as.matrix(heatmap2.df),
              col=hmcol2,
              trace="none",
              labRow=vector(mode = "character", length = nrow(heatmap2.df)),
              labCol=groups[,2],
              #Rowv = as.dendrogram(rowdendro),
              Colv = as.dendrogram(coldendro),
              ColSideColors=colcol,
              lwid=c(1,4),
              lhei = c(1,5),
              breaks = seq(-3,3,by=.5),
              symkey=F,
              dendrogram = "both")
  }
}
```

##### wMel, Darby et al 2014
![Image description](/images/wMel_darby_2014_heatmap_degenes.png)

##### wDi, Luck et al 2015
![Image description](/images/wDi_luck_2015_heatmap_degenes.png)

##### wMel, Gutzwiller et al 2015
![Image description](/images/wMel_gutzwiller_2015_heatmap_degenes.png)

##### wBm, Grote et al 2017
![Image description](/images/wBm_grote_2017_heatmap_degenes.png)

##### wBm, Chung et al 2019
![Image description](/images/wBm_chung_2019_heatmap_degenes.png)

### Choose studies for basic pairwise analyses or WGCNA <a name="ranalysis_de_selectanalysis"></a>

```{R}
pairwise.studies <- c("wMel_darby_2014")
wgcna.studies <- c("wDi_luck_2015","wMel_gutzwiller_2015","wBm_grote_2017","wBm_chung_2019")
```

### Conduct pairwise analyses on studies with two biological groups <a name="ranalysis_de_pairwise"></a>

#### Create z-score of log2TPM heatmaps with a column bar that shows up- and down-regulated genes <a name="ranalysis_de_pairwise_hm"></a>

```{R, fig.height = 10, fig.width = 5}
for(i in 1:length(pairwise.studies)){
  sample_map.subset <- sample_map[sample_map$study == pairwise.studies[i],]
  groups <- as.data.frame(cbind(as.character(sample_map.subset$sample_identifier),
                                as.character(sample_map.subset$sample_group),
                                as.character(sample_map.subset$color)))
  groups <- unique(groups)

  groups[,1] <- factor(groups[,1],levels=groups[,1])
  groups[,2] <- factor(groups[,2],levels=unique(groups[,2]))
  groups[,3] <- factor(groups[,3],levels=unique(groups[,3]))
  
  qlf <- as.data.frame(qlf.list[names(qlf.list) == pairwise.studies[i]])
  qlf <- qlf[order(qlf[,1],decreasing = T),]
  qlf <- qlf[rownames(qlf) %in% rownames(as.data.frame(tpm.edgeR.list[names(tpm.edgeR.list) == pairwise.studies[i]])),]
  
  log2tpm.de <- log2(as.data.frame(tpm.edgeR.list[names(tpm.edgeR.list) == pairwise.studies[i]])+1)
  zscore.log2tpm.de <- as.data.frame(t(scale(t(log2tpm.de))))
  zscore.log2tpm.de <- zscore.log2tpm.de[match(rownames(qlf),rownames(zscore.log2tpm.de)),]
  
  hmcol1 <- colorRampPalette(c("navyblue","white","firebrick3"))(12)
  hmcol2 <- colorRampPalette(c("#FFFFFF","#FDE0D2","#FABAA1","#F69274", "#F16B4E","#EF3D2D","#590A16","#A51E23","#650C16"))(30)
  rowcol1 <- ifelse(qlf[,1] > 0,"red","blue")
  colcol <- as.character(groups[,3])
  rowsep <- get_heatmap_separators(rowcol1)
  colsep <- get_heatmap_separators(colcol)
  
  coldendro <- pvclust(as.data.frame(tpm.edgeR.list[names(tpm.edgeR.list) == pairwise.studies[i]]), method.dist="cor", method.hclust="average", nboot=100, quiet=T)

  pdf(paste0(OUTPUT.DIR,"/",pairwise.studies[i],"_pairwise_zscorelog2tpm_heatmap.pdf"),
      width = 5,
      height = 10)
    
  heatmap.2(as.matrix(zscore.log2tpm.de),
              col=hmcol1,
              trace="none",
              labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
              Rowv = F,
              Colv = as.dendrogram(coldendro),
              RowSideColors=rowcol1,
              ColSideColors=colcol,
              lhei = c(2,8),
              breaks = seq(-3,3,by=.5),
              rowsep = rowsep,
              colsep = colsep,
              symkey=F,
              dendrogram = "column")
  dev.off()
  
  png(paste0(OUTPUT.DIR,"/",pairwise.studies[i],"_pairwise_zscorelog2tpm_heatmap.png"),
      width = 5,
      height = 10,
      units = "in",res=300)
    
  heatmap.2(as.matrix(zscore.log2tpm.de),
              col=hmcol1,
              trace="none",
              labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
              Rowv = F,
              Colv = as.dendrogram(coldendro),
              RowSideColors=rowcol1,
              ColSideColors=colcol,
              lhei = c(2,8),
              breaks = seq(-3,3,by=.5),
              rowsep = rowsep,
              colsep = colsep,
              symkey=F,
              dendrogram = "column")
  dev.off()
  
  heatmap.2(as.matrix(zscore.log2tpm.de),
              col=hmcol1,
              trace="none",
              labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
              Rowv = F,
              Colv = as.dendrogram(coldendro),
              RowSideColors=rowcol1,
              ColSideColors=colcol,
              lhei = c(2,8),
              breaks = seq(-3,3,by=.5),
              rowsep = rowsep,
              colsep = colsep,
              symkey=F,
              dendrogram = "column")
  
}
```

##### wMel, Darby et al 2014
![Image description](/images/wMel_darby_2014_pairwise_zscorelog2tpm_heatmap.pdf)

#### Conduct functional term enrichment analysis on up- and down-regulated gene subsets <a name="ranalysis_de_pairwise_fxnterm"></a>

```{R}
terms.pairwise.list <- list()
terms.sig.pairwise.list <- terms.pairwise.list

for(i in 1:length(pairwise.studies)){
  strain <- substr(pairwise.studies,1,regexpr("_",pairwise.studies)-1)
  geneinfo <- read.delim(paste0(GENEINFO.DIR,"/",strain,".cds.fna.interproscan.geneinfo.tsv"))
  
  qlf <- as.data.frame(qlf.list[names(qlf.list) == pairwise.studies[i]])
  qlf <- qlf[order(qlf[,1],decreasing = T),]
  qlf <- qlf[rownames(qlf) %in% rownames(as.data.frame(tpm.edgeR.list[names(tpm.edgeR.list) == pairwise.studies[i]])),]
  
  terms.colnames <- c("term","clusteroccurences","genomeoccurences","pvalue","correctedpvalue","oddsratio","group")
  terms.pairwise <- as.data.frame(matrix(nrow = 0,
                                         ncol = 7))
  colnames(terms.pairwise) <- terms.colnames
  
  terms.pairwise.up <- as.data.frame(cbind(functionaltermenrichment(rownames(qlf)[qlf[,1] > 0],geneinfo),
                                           "upregulated"))
  terms.pairwise.down <- as.data.frame(cbind(functionaltermenrichment(rownames(qlf)[qlf[,1] < 0],geneinfo),
                                       "downregulated"))
  
  colnames(terms.pairwise.up) <- terms.colnames
  colnames(terms.pairwise.down) <- terms.colnames
  terms.pairwise <- as.data.frame(rbind(terms.pairwise.up,
                                     terms.pairwise.down))
  
  for(j in 2:6){
    terms.pairwise[,j] <- as.numeric(as.character(terms.pairwise[,j]))
  }
  
  terms.pairwise.list[[i]] <- terms.pairwise
  terms.sig.pairwise.list[[i]] <- terms.pairwise[terms.pairwise$correctedpvalue < 0.05,]
  
  names(terms.pairwise.list)[i] <- pairwise.studies[i]
  
  write.table(terms.pairwise.list[[i]],
            paste0(OUTPUT.DIR,"/",names(terms.pairwise.list)[i],"_functionaltermenrichment.tsv"),
            row.names = T,
            col.names = T,
            quote = F,
            sep = "\t")
}
```

### Conduct WGCNA on time-course studies <a name="ranalysis_de_wgcna"></a>

#### Identify WGCNA soft power values <a name="ranalysis_de_wgcna_softpower"></a>
```{R, fig.height = 4, fig.width = 8}
tpm.wgcna.list <- tpm.edgeR.list
for(i in 1:length(tpm.wgcna.list)){
  tpm.wgcna.list[[i]] <- as.data.frame(t(tpm.edgeR.list[[i]]))
}
tpm.wgcna.list <- tpm.wgcna.list[names(tpm.wgcna.list) %in% wgcna.studies]
softpower.list <- tpm.wgcna.list

for(i in 1:length(tpm.wgcna.list)){
  powers <- c(1:20)
  sft <- pickSoftThreshold(tpm.wgcna.list[[i]], powerVector = powers, verbose = 5)
  softpower.list[[i]] <- find_soft_power(sft)
  
  text.color <- rep("black",length(sft$fitIndices[,1]))
  text.color[which(sft$fitIndices[,1] == softpower.list[[i]])] <- "red"
  
  scale_independence.plot <- ggplot()+
    geom_text(mapping = aes(x = sft$fitIndices[,1], y = -sign(sft$fitIndices[,3])*sft$fitIndices[,2], label=sft$fitIndices[,1]), color = text.color)+
    labs(title = "scale independence", x = "soft threshold (power)", y = "scale free topology model fit, signed R^2")+
    theme_bw()
  
  mean_connectivity.plot <- ggplot()+
    geom_text(mapping = aes(x = sft$fitIndices[,1], y = sft$fitIndices[,5], label=sft$fitIndices[,1]), color = text.color)+
    labs(title = "mean connectivity", x = "soft threshold (power)", y = "mean connectivity")+
    theme_bw()
  
  wgcna_soft_power_plots <- list(scale_independence.plot, mean_connectivity.plot)
  lay <- rbind(c(1,2))
  pdf(paste0(OUTPUT.DIR,"/",names(tpm.wgcna.list)[i],"_wgcna_soft_power_plot.pdf"),
      width = 8, 
      height = 4)
  grid.arrange(grobs = wgcna_soft_power_plots,
               widths = c(5,5),
               heights = c(5),
               layout_matrix = lay)
  dev.off()
  
  png(paste0(OUTPUT.DIR,"/",names(tpm.wgcna.list)[i],"_wgcna_soft_power_plot.png"),
      width = 8, 
      height = 4,
      units = "in",res=300)
  grid.arrange(grobs = wgcna_soft_power_plots,
               widths = c(5,5),
               heights = c(5),
               layout_matrix = lay)
  dev.off()
  
  grid.arrange(grobs = wgcna_soft_power_plots,
               widths = c(5,5),
               heights = c(5),
               layout_matrix = lay)
}
```

##### wDi, Luck et al 2015
![Image description](/images/wDi_luck_2015_wgcna_soft_power_plot.png)

##### wMel, Gutzwiller et al 2015
![Image description](/images/wMel_gutzwiller_2015_wgcna_soft_power_plot.png)

##### wBm, Grote et al 2017
![Image description](/images/wBm_grote_2017_wgcna_soft_power_plot.png)

##### wBm, Chung et al 2019
![Image description](/images/wBm_chung_2019_wgcna_soft_power_plot.png)

#### Create WGCNA expression modules and merge similar modules <a name="ranalysis_de_wgcna_createmodules"></a>

```{r wgcna_createmodules, fig.height = 4, fig.width = 8}
for(i in 1:length(tpm.wgcna.list)){
  adjacency <- adjacency(tpm.wgcna.list[[i]], power = softpower.list[[i]])
  TOM <- TOMsimilarity(adjacency)
  dissTOM <- 1-TOM
  geneTree <- hclust(as.dist(dissTOM), method = "average");
  
  minModuleSize <- 1
  dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                               deepSplit = 2, pamRespectsDendro = FALSE,
                               minClusterSize = minModuleSize)
  
  dynamicColors = labels2colors(dynamicMods)
  MEList = moduleEigengenes(tpm.wgcna.list[[i]], colors = dynamicColors)
  MEs = MEList$eigengenes
  MEDiss = 1-cor(MEs, use = "pairwise.complete.obs")
  METree = hclust(as.dist(MEDiss), method = "average")
  
  MEDissThres = 0.25
  
  pdf(paste0(OUTPUT.DIR,"/",names(tpm.wgcna.list)[i],"_wgcna_merge_eigengenes_plot.pdf"),
      width = 8,
      height = 5)
  plot(METree, main = "clustering of module eigengenes",
       xlab = "", sub = "")
  
  abline(h=MEDissThres, col = "red")
  dev.off()
  
  png(paste0(OUTPUT.DIR,"/",names(tpm.wgcna.list)[i],"_wgcna_merge_eigengenes_plot.png"),
      width = 8,
      height = 5,
      units = "in",res=300)
  plot(METree, main = "clustering of module eigengenes",
       xlab = "", sub = "")
  abline(h=MEDissThres, col = "red")
  dev.off()
  
  plot(METree, main = "clustering of module eigengenes",
       xlab = "", sub = "")
  abline(h=MEDissThres, col = "red")
  
  merge = mergeCloseModules(tpm.wgcna.list[[i]], dynamicColors, cutHeight = MEDissThres, verbose = 3)
  mergedColors = merge$colors
  mergedMEs = merge$newMEs
  
  pdf(paste0(OUTPUT.DIR,"/",names(tpm.wgcna.list)[i],"_wgcna_eigengene_dendrogram_plot.pdf"),
      width = 8,
      height = 5)
  plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("dynamic tree cut", "merged dynamic"),
                      dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
  dev.off()
  
  png(paste0(OUTPUT.DIR,"/",names(tpm.wgcna.list)[i],"_wgcna_eigengene_dendrogram_plot.png"),
      width = 8,
      height = 5,
      units = "in",res=300)
  plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("dynamic tree cut", "merged dynamic"),
                      dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
  dev.off()
  
  plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("dynamic tree cut", "merged dynamic"),
                      dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
  
  tpm.wgcna.list[[i]] <- eigengene_invert_id(as.data.frame(t(tpm.wgcna.list[[i]])), mergedColors, mergedMEs)
  tpm.wgcna.list[[i]] <- wgcna_heatmap_reorder(tpm.wgcna.list[[i]])
  
  write.table(tpm.wgcna.list[[i]],
              paste0(OUTPUT.DIR,"/",names(tpm.wgcna.list)[i],"_tpm_de_modules.tsv"),
              row.names = T,
              col.names = T,
              quote = F,
              sep = "\t")
}
```

##### wDi, Luck et al 2015
![Image description](/images/wDi_luck_2015_wgcna_eigengene_dendrogram_plot.png)
![Image description](/images/wDi_luck_2015_wgcna_merge_eigengenes_plot.png)

##### wMel, Gutzwiller et al 2015
![Image description](/images/wMel_gutzwiller_2015_wgcna_eigengene_dendrogram_plot.png)
![Image description](/images/wMel_gutzwiller_2015_wgcna_merge_eigengenes_plot.png)

##### wBm, Grote et al 2017
![Image description](/images/wBm_grote_2017_wgcna_eigengene_dendrogram_plot.png)
![Image description](/images/wBm_grote_2017_wgcna_merge_eigengenes_plot.png)

##### wBm, Chung et al 2019
![Image description](/images/wBm_chung_2019_wgcna_eigengene_dendrogram_plot.png)
![Image description](/images/wBm_chung_2019_wgcna_merge_eigengenes_plot.png)

#### Plot WGCNA-derived expression modules <a name="ranalysis_de_wgcna_heatmap"></a>

```{R, fig.height = 10, fig.width = 5}
for(i in 1:length(tpm.wgcna.list)){
  sample_map.subset <- sample_map[sample_map$study == names(tpm.wgcna.list)[i] & sample_map$sample_identifier %in% colnames(tpm.wgcna.list[[i]]),]
  groups <- as.data.frame(cbind(as.character(sample_map.subset$sample_identifier),
                                as.character(sample_map.subset$sample_group),
                                as.character(sample_map.subset$color)))
  groups <- unique(groups)

  groups[,1] <- factor(groups[,1],levels=groups[,1])
  groups[,2] <- factor(groups[,2],levels=unique(groups[,2]))
  groups[,3] <- factor(groups[,3],levels=unique(groups[,3]))
  
  log2tpm.de <- log2(tpm.wgcna.list[[i]][,1:(ncol(tpm.wgcna.list[[i]]) - 2)] + 1)
  zscore.log2tpm.de <- as.data.frame(t(scale(t(log2tpm.de))))
  
  hmcol1 <- colorRampPalette(c("navyblue","white","firebrick3"))(12)
  hmcol2 <- colorRampPalette(c("#FFFFFF","#FDE0D2","#FABAA1","#F69274", "#F16B4E","#EF3D2D","#590A16","#A51E23","#650C16"))(30)
  rowcol1 <- tpm.wgcna.list[[i]]$module
  colcol <- as.character(groups[,3])
  rowsep <- get_heatmap_separators(rowcol1)
  colsep <- get_heatmap_separators(colcol)
  
  rowcol2 <- tpm.wgcna.list[[i]]$invert
  rowcol2[rowcol2 == F] <- "grey"
  rowcol2[rowcol2 == T] <- "black"
  
  coldendro.df <- as.data.frame(tpm.kept.list[names(tpm.kept.list) == names(tpm.wgcna.list)[i]])
  coldendro.df <- t(scale(t(log2(coldendro.df + 1))))
  coldendro <- pvclust(coldendro.df, method.dist="cor", method.hclust="average", nboot=100, quiet=T)

  if(names(tpm.wgcna.list)[i] == "wMel_gutzwiller_2015"){
    coldendro <- rotate(as.dendrogram(coldendro), c(30:50,55:54,52:53,51,11:29,1:10))
  }

  pdf(paste0(OUTPUT.DIR,"/",names(tpm.wgcna.list)[i],"_wgcna_zscorelog2tpm_heatmap.pdf"),
      width = 5,
      height = 10)
  heatmap.2(as.matrix(zscore.log2tpm.de),
                col=hmcol1,
                trace="none",
                labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
                Rowv = F,
                Colv = as.dendrogram(coldendro),
                RowSideColors=rowcol1,
                ColSideColors=colcol,
                lhei = c(2,8),
                breaks = seq(-3,3,by=.5),
                rowsep = rowsep,
                colsep = colsep,
                dendrogram = "col")
  heatmap.2(as.matrix(zscore.log2tpm.de),
                col=hmcol1,
                trace="none",
                labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
                Rowv = F,
                Colv = as.dendrogram(coldendro),
                RowSideColors=rowcol2,
                ColSideColors=colcol,
                lhei = c(2,8),
                breaks = seq(-3,3,by=.5),
                rowsep = rowsep,
                colsep = colsep,
                dendrogram = "col")
  dev.off()
  
  png(paste0(OUTPUT.DIR,"/",names(tpm.wgcna.list)[i],"_wgcna_zscorelog2tpm_heatmap.png"),
      width = 5,
      height = 10,
      units = "in",res=300)
  heatmap.2(as.matrix(zscore.log2tpm.de),
                col=hmcol1,
                trace="none",
                labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
                Rowv = F,
                Colv = as.dendrogram(coldendro),
                RowSideColors=rowcol1,
                ColSideColors=colcol,
                lhei = c(2,8),
                breaks = seq(-3,3,by=.5),
                rowsep = rowsep,
                colsep = colsep,
                dendrogram = "col")
  heatmap.2(as.matrix(zscore.log2tpm.de),
                col=hmcol1,
                trace="none",
                labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
                Rowv = F,
                Colv = as.dendrogram(coldendro),
                RowSideColors=rowcol2,
                ColSideColors=colcol,
                lhei = c(2,8),
                breaks = seq(-3,3,by=.5),
                rowsep = rowsep,
                colsep = colsep,
                dendrogram = "col")
  dev.off()
  
  heatmap.2(as.matrix(zscore.log2tpm.de),
                col=hmcol1,
                trace="none",
                labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
                Rowv = F,
                Colv = as.dendrogram(coldendro),
                RowSideColors=rowcol1,
                ColSideColors=colcol,
                lhei = c(2,8),
                breaks = seq(-3,3,by=.5),
                rowsep = rowsep,
                colsep = colsep,
                dendrogram = "col")
  heatmap.2(as.matrix(zscore.log2tpm.de),
                col=hmcol1,
                trace="none",
                labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
                Rowv = F,
                Colv = as.dendrogram(coldendro),
                RowSideColors=rowcol2,
                ColSideColors=colcol,
                lhei = c(2,8),
                breaks = seq(-3,3,by=.5),
                rowsep = rowsep,
                colsep = colsep,
                dendrogram = "col")
}
```

##### wDi, Luck et al 2015
![Image description](/images/wDi_luck_2015_wgcna_zscorelog2tpm_heatmap.png)

##### wMel, Gutzwiller et al 2015
![Image description](/images/wMel_gutzwiller_2015_wgcna_zscorelog2tpm_heatmap.png)

##### wBm, Grote et al 2017
![Image description](/images/wBm_grote_2017_wgcna_zscorelog2tpm_heatmap.png)

##### wBm, Chung et al 2019
![Image description](/images/wBm_chung_2019_wgcna_zscorelog2tpm_heatmap.png)

#### Conduct functional term enrichment analysis on WGCNA expression modules <a name="ranalysis_de_wgcna_fxnterm"></a>

```{R}
terms.wgcna.list <- tpm.wgcna.list
terms.sig.wgcna.list <- tpm.wgcna.list

for(i in 1:length(tpm.wgcna.list)){
  strain <- substr(names(tpm.wgcna.list)[i],1,regexpr("_",names(tpm.wgcna.list)[i])-1)
  geneinfo <- read.delim(paste0(GENEINFO.DIR,"/",strain,".cds.fna.interproscan.geneinfo.tsv"))
  
  terms.colnames <- c("term","clusteroccurences","genomeoccurences","pvalue","correctedpvalue","oddsratio","module","invert")
  terms.wgcna <- as.data.frame(matrix(nrow = 0,
                                      ncol = 8))
  colnames(terms.wgcna) <- terms.colnames
  for(j in 1:length(unique(tpm.wgcna.list[[i]]$module))){
    terms.wgcna.f <- as.data.frame(cbind(functionaltermenrichment(rownames(tpm.wgcna.list[[i]])[tpm.wgcna.list[[i]]$module == unique(tpm.wgcna.list[[i]]$module)[j] & tpm.wgcna.list[[i]]$invert == F],geneinfo),unique(tpm.wgcna.list[[i]]$module)[j],F))
    terms.wgcna.t <- as.data.frame(cbind(functionaltermenrichment(rownames(tpm.wgcna.list[[i]])[tpm.wgcna.list[[i]]$module == unique(tpm.wgcna.list[[i]]$module)[j] & tpm.wgcna.list[[i]]$invert == T],geneinfo),unique(tpm.wgcna.list[[i]]$module)[j],T))
    
    colnames(terms.wgcna.f) <- terms.colnames
    colnames(terms.wgcna.t) <- terms.colnames
    terms.wgcna <- as.data.frame(rbind(terms.wgcna,
                                       terms.wgcna.f, 
                                       terms.wgcna.t))
  }
  
  for(j in 2:6){
    terms.wgcna[,j] <- as.numeric(as.character(terms.wgcna[,j]))
  }
  
  terms.wgcna.list[[i]] <- terms.wgcna
  terms.sig.wgcna.list[[i]] <- terms.wgcna[terms.wgcna$correctedpvalue < 0.05,]
  
  write.table(terms.wgcna.list[[i]],
            paste0(OUTPUT.DIR,"/",names(terms.wgcna.list)[i],"_functionaltermenrichment.tsv"),
            row.names = T,
            col.names = T,
            quote = F,
            sep = "\t")
}
```