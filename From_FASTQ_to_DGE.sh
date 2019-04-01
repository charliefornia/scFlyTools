#!/bin/bash

################################################################################
# TITLE: Convert barcoded reads to digital expression matrix.
# VERSION: 0.1.0 (beta)
# AUTHOR: Christoph Treiber, Waddell lab, University of Oxford
# DATE: 01/04/2019 (dd/mm/yyyy)
################################################################################

################################################################################
# This pipeline is largely based on the Macosko DropSeq Cookbook
# (see https://github.com/broadinstitute/Drop-seq)
# Please cite:
# "Highly Parallel Genome-wide Expression Profiling of Individual Cells Using
# Nanoliter Droplets"; Evan Z. Macosko, ..., Steven A. McCarroll; Cell 2015
################################################################################

################################################################################
# Make sure that DropSeqtools are installed (here v 2.1.0)
# Make sure that STAR aligner is installed and added to $PATH
# Make sure that java is installed and added to $PATH
################################################################################

################################################################################
# || VERY IMPORTANT :
# || Make sure that the cell- and UMI barcode lengths are correct !!!!
celltaglength=NUMBER_OF_READS
UMIlength=NUMBER_OF_READS
# || For 10x kit 3' - version2: celltaglength=16, UMIlength=10
# || For 10x kit 3' - version3: celltaglength=16, UMIlength=12
# || For Dropseq (Macosko, Chemgene beads): celltaglength=12, UMIlength=8
################################################################################

############################
# Setting the environment: #
############################

pathtoDROPSEQTOOLS=/PATH/TO/Drop-seq_tools-2.1.0
numberofcores=1
workingdirectory=/PATH/TO/WORKINGDIRECTORY/
pathtoREF=/PATH/TO/REF
genomename=NAME_OF_GENOME
# e.g. dmel625; match with file name in /PATH/TO/REF/
numberofcells=NUMBER_OF_CELLS
# e.g. 10000, this number determines the size of the digital expression matrix
# and therefore the size of the R object when using the Seurat package.

################################################################################

############################################
# UPDATE THESE PARAMETERS FOR EACH SAMPLE: #
############################################
fastafile1=/PATH/TO/FASTAFILE1.fasta.gz
fastafile2=/PATH/TO/FASTAFILE2.fasta.gz
samplename=NAME_OF_SAMPLE

################################################################################

cd ${workingdirectory}

mkdir ${samplename}_summaries

# STEP1: convert Fastq to SAM
java -jar ${pathtoDROPSEQTOOLS}/3rdParty/picard/picard.jar FastqToSam \
	F1=${fastafile1} \
	F2=${fastafile2} \
	O=./01_${samplename}_unmapped.bam \
	SM=${samplename} \
	&&

# STEP2: add flag containing cell barcode
${pathtoDROPSEQTOOLS}/TagBamWithReadSequenceExtended \
	SUMMARY=./${samplename}_summaries/01to02_summary.txt \
	BASE_RANGE=1-${celltaglength} \
	BASE_QUALITY=10 \
	BARCODED_READ=1 \
	DISCARD_READ=false \
	TAG_NAME=XC \
	NUM_BASES_BELOW_QUALITY=1 \
	INPUT=./01_${samplename}_unmapped.bam \
	OUTPUT=./02_${samplename}_unaligned_tagged_Cell.bam \
	&&

rm 01_${samplename}_unmapped.bam &&

# STEP3: add flag containing UMI barcode
${pathtoDROPSEQTOOLS}/TagBamWithReadSequenceExtended \
	SUMMARY=./${samplename}_summaries/02to03_summary.txt \
	BASE_RANGE=$((celltaglength + 1))-$((UMIlength + celltaglength)) \
	BASE_QUALITY=10 \
	BARCODED_READ=1 \
	DISCARD_READ=true \
	TAG_NAME=XM \
	NUM_BASES_BELOW_QUALITY=1 \
	INPUT=./02_${samplename}_unaligned_tagged_Cell.bam \
	OUTPUT=03_${samplename}_unaligned_tagged_CellMolecular.bam \
	&&

rm 02_${samplename}_unaligned_tagged_Cell.bam &&

# STEP4: remove reads where barcode has low quality score
${pathtoDROPSEQTOOLS}/FilterBAM \
	TAG_REJECT=XQ \
	INPUT=03_${samplename}_unaligned_tagged_CellMolecular.bam \
	OUTPUT=04_${samplename}_unaligned_tagged_trimmed_smart.bam \
	&&

rm 03_${samplename}_unaligned_tagged_CellMolecular.bam &&

# STEP5: trim poly-A sequences off reads, output file 05_ is kept and matched up
# 		 with aligned reads file 08_ in step 9.
${pathtoDROPSEQTOOLS}/PolyATrimmer \
	OUTPUT=05_${samplename}_unaligned_mc_tagged_polyA_filtered.bam \
	OUTPUT_SUMMARY=${samplename}_summaries/04to05_polyA_trimming_report.txt \
	MISMATCHES=0 \
	NUM_BASES=6 \
	INPUT=04_${samplename}_unaligned_tagged_trimmed_smart.bam \
	USE_NEW_TRIMMER=true \
	&&

rm 04_${samplename}_unaligned_tagged_trimmed_smart.bam &&

# STEP6: convert SAM file back to FASTQ
java -Xmx500m -jar ${pathtoDROPSEQTOOLS}/3rdParty/picard/picard.jar SamToFastq \
	INPUT=05_${samplename}_unaligned_mc_tagged_polyA_filtered.bam \
	FASTQ=06_${samplename}_unaligned_mc_tagged_polyA_filtered.fastq \
	&&

# align reads to reference genome
STAR --genomeDir ${pathtoREF} \
	--runThreadN ${numberofcores} \
	--outFileNamePrefix 07_${samplename}_${genomename}_star_ \
	--readFilesIn ./06_${samplename}_unaligned_mc_tagged_polyA_filtered.fastq \
	&&

mv 07_${samplename}_${genomename}_star_Log* ./${samplename}_summaries/ &&

rm 06_${samplename}_unaligned_mc_tagged_polyA_filtered.fastq &&

# sort aligned reads
java -Dsamjdk.buffer_size=131072 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 \
	-Xmx4000m -jar ${pathtoDROPSEQTOOLS}/3rdParty/picard/picard.jar SortSam \
	INPUT=07_${samplename}_${genomename}_star_Aligned.out.sam \
	OUTPUT=08_${samplename}_${genomename}_aligned_sorted.bam \
	SORT_ORDER=queryname \
	TMP_DIR=./ \
	&&

rm 07_${samplename}_${genomename}_star_Aligned.out.sam &&

rm 07_${samplename}_${genomename}_star_SJ.out.tab &&

# combine aligned reads with tags from previous file
java -Xmx4000m -jar ${pathtoDROPSEQTOOLS}/3rdParty/picard/picard.jar \
	MergeBamAlignment \
	REFERENCE_SEQUENCE=${pathtoREF}/${genomename}.fa \
	UNMAPPED_BAM=05_${samplename}_unaligned_mc_tagged_polyA_filtered.bam \
	ALIGNED_BAM=08_${samplename}_${genomename}_aligned_sorted.bam \
	INCLUDE_SECONDARY_ALIGNMENTS=false \
	PAIRED_RUN=false \
	OUTPUT=09_${samplename}_${genomename}_merged.bam \
	&&

rm 05_${samplename}_unaligned_mc_tagged_polyA_filtered.bam &&

rm 08_${samplename}_${genomename}_aligned_sorted.bam &&

# add flag to each read containing gene name (from refFlat) - strand-specific
${pathtoDROPSEQTOOLS}/TagReadWithGeneFunction \
	O=10_${samplename}_${genomename}_merged_gene_exon_tagged.bam \
	ANNOTATIONS_FILE=${pathtoREF}/${genomename}.refFlat \
	INPUT=09_${samplename}_${genomename}_merged.bam \
	&&
	
rm 09_${samplename}_${genomename}_merged.bam &&

# generate DGE. Here, reads that fall into intronic region are also counted
${pathtoDROPSEQTOOLS}/DigitalExpression \
	I=10_${samplename}_${genomename}_merged_gene_exon_tagged.bam \
	O=11_${samplename}_${genomename}_top${numberofcells}cells.dge.txt.gz \
	SUMMARY=./${samplename}_summaries/11_top${numberofcells}cells_summary.txt \
	NUM_CORE_BARCODES=${numberofcells} \
	LOCUS_FUNCTION_LIST=INTRONIC
