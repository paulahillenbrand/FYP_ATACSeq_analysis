** ATAC Analysis Protocol**

### Conda Environment
conda create -n atac
conda activate atac

conda install bioconda::fastqc=0.11.9
conda install bioconda::multiqc=1.0
conda install bioconda::trim-galore=0.6.10
conda install bioconda::bowtie2=2.3.4.3
conda install bioconda::samtools=1.10
conda install bioconda::bedtools=2.29.1
conda install bioconda::deeptools=2.0 
conda install bioconda::picard=2.18.29
conda install bioconda::sambamba=1.0.1
conda install bioconda::macs2=2.2.6

## FastQC
fastqc /atac/rawdata/{samp_name}

## MultiQC 
multiqc /atac/rawdata

## trim-galore
trim_galore \
--paired \
--o /atac/rawdata/trimmed \
/atac/rawdata/{samp_name_L1_1} \
/atac/rawdata/{samp_name_L1_2}


### Alignment
bowtie2 \
--very-sensitive \
-X 2000 \
-x /atac/genome/HSymV2.1.fasta \
-1 /atac/trimmed/{samp_name_L1_1} \
-2 /atac/trimmed/{samp_name_L1_2} \
-S /atac/alignment/{samp_name}_bowtie2.sam \
&> /atac/alignment/sam/{samp_name}_bowtie2.txt

## File Processing 
samtools view \
-Sb \
-F 0x04 \
alignment/sam/{samp_name}_bowtie2.bam \
-o alignment/bam/{samp_name}_bowtie2.bam

samtools sort \
/atac/alignment/bam/{samp_name}_bowtie2.bam \
-o /atac/alignment/bam/{samp_name}_bowtie2_sorted.bam

samtools index \
/atac/alignment/bam/{samp_name}_bowtie2_sorted.bam

### Remove Multimappers
sambamba view \
-h \
-t 2 \
-f bam \
-F "[XS] == null and not unmapped" \
/atac/alignment/bam/{samp_name}_bowtie2_sorted.bam \
-o /atac/alignment/bam/no_multimappers/{samp_name}_bowtie2_sorted_nomulti.bam

### Mark Duplicates
picard MarkDuplicates \
I=/atac/alignment/bam/no_multimappers/{samp_name}_bowtie2_sorted_nomulti.bam \
O=/atac/alignment/bam/no_multimappers/dupmarked/{samp_name}_bowtie2_sorted_nomulti_dupMarked.bam \
METRICS_FILE=/atac/alignment/bam/no_multimappers/dupmarked/{samp_name}_picard.dupMark.txt

### Shift alignements due to Tn5 strand bias
## index
samtools index \
/atac/alignment/bam/no_multimappers/dupmarked/{samp_name}_bowtie2_sorted_nomulti_dupMarked.bam

## shift alignments
alignmentSieve \
--ATACshift \
--bam /atac/alignment/bam/no_multimappers/dupmarked/{samp_name}_bowtie2_sorted_nomulti_dupMarked.bam \
-o /atac/alignment/bam/no_multimappers/dupmarked/{samp_name}_bowtie2_sorted_nomulti.dupMarked_shifted.bam

### Filter Fragment Size (0 - 247 bp)
samtools view \
-h /atac/alignment/bam/no_multimappers/dupmarked/{samp_name}_bowtie2_sorted_nomulti.dupMarked_shifted.bam | \
awk 'substr($0,1,1)=="@" || ($9>= 0 && $9<=247) || ($9<= 0 && $9>=-247)' | \
samtools view -b \
> /atac/alignment/bam/no_multimappers/dupmarked/shifted_filtered/{samp_name}_bowtie2_sorted_nomulti.dupMarked_shifted_filtered.bam

### index and sort
samtools sort \
/atac/alignment/bam/no_multimappers/dupmarked/shifted_filtered/{samp_name}_bowtie2_sorted_nomulti.dupMarked_shifted_filtered.bam \
-o /atac/alignment/bam/no_multimappers/dupmarked/shifted_filtered/sorted/{samp_name}_bowtie2_sorted_nomulti.dupMarked_shifted_filtered.bam

samtools index \
/atac/alignment/bam/no_multimappers/dupmarked/shifted_filtered/sorted/{samp_name}_bowtie2_sorted_nomulti.dupMarked_shifted_filtered.bam

### Convert to Bigwig
bamCoverage \
-b /atac/alignment/bam/no_multimappers/dupmarked/shifted_filtered/sorted/{samp_name}_bowtie2_sorted_nomulti.dupMarked_shifted_filtered.bam \
-o /atac/alignment/bigwig/{samp_name}_nomulti_dupMarked_shifted_filtered_norm.bw \
--normalizeUsing RPKM \
--exactScaling

### Query bam file for peak calling
samtools sort \
-n /atac/alignment/bam/no_multimappers/dupmarked/shifted_filtered/sorted/{samp_name}_bowtie2_sorted_nomulti.dupMarked_shifted_filtered.bam \
-o /Users/paula/Downloads/atac/alignment/query/{samp_name}_nomulti_dupMarked_shifted_filtered.bam

## MACS2 Peak Calling:
macs2 callpeak \
--broad \
-t /atac/alignment/query/{samp_name}_nomulti_dupMarked_shifted_filtered.bam \
-n {samp_name}_filtered \
-f BAMPE \
-g 185312851 \
-q 0.05 \
--nomodel \
--shift -37 
--extsize 73 \
--keep-dup all \
--outdir /Users/paula/Downloads/atac/MACS2

## Consensus peaks
bedtools intersect \
-a /atac/MACS2/{samp_name}_1_filtered_peaks.broadPeak \
-b /atac/MACS2/{samp_name}_2_filtered_peaks.broadPeak \
> nemato_1_intersect.bed

### Pygenometracks for marker genes
conda create -n pygenometracks -c bioconda -c conda-forge pygenometracks python=3.7
## create .ini files
## visualise tracks in pygenome tracks over gene regions of interest
    trichohyalin (0-200): NC_079883.1:6,733,000-6,740,695
    Ncol1 (0-200): NC_079886.1:21,987,446-21,990,374
    Piwi1 (0-450): NC_079887.1:16,491,611-16,514,937
    Piwi2 (0-450): NC_079889.1:11,087,216-11,108,520
    RFamide (0-150): NC_079883.1:20,780,393-20,787,079
    Pln (0-250): NC_079879.1:9,002,974-9,010,151

pyGenomeTracks --tracks tracks.ini --region {genomic_region} -o image.png 

#### TRANSCRIPTION START SITE ENRICHMENT 

### IN RSTUDIO
## Load packages
install.packages("BiocManager", repos = "https://cloud.r-project.org")
library(BiocManager)
BiocManager::install(c("ATACseqQC", "ChIPpeakAnno", "MotifDb", "GenomicAlignments"))
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("motifStack")
library(ATACseqQC)
library(ggplot2)
library(magrittr)
library(Rsamtools)
library(ChIPpeakAnno)

## Quality Control on bam files (pre shift)
bamQC("/atac/alignment/bam/no_multimappers/dupmarked/{samp_name}_bowtie2_sorted_nomulti_dupMarked.bam", outPath = NULL)

## Split Bam files into nucleosome free, mono- , di-, tri-nucleosome
outPath <- "/atac/splitBAM/{samp_name}"
if (dir.exists(outPath))
{
  unlink(outPath, recursive = TRUE, force = TRUE)
}
dir.create(outPath)

# bamfile tags to be read in
possibleTag <- combn(LETTERS, 2)
possibleTag <- c(paste0(possibleTag[1, ], possibleTag[2, ]),
                 paste0(possibleTag[2, ], possibleTag[1, ]))
bamTop100 <- scanBam(BamFile("/atac/alignment/bam/no_multimappers/dupmarked/{samp_name}_bowtie2_sorted_nomulti_dupMarked.bam", yieldSize = 100),
                     param = ScanBamParam(tag = possibleTag))[[1]]$tag
tags <- names(bamTop100)[lengths(bamTop100)>0]

# shift the coordinates of 5'ends of alignments in the bam file
gal <- readBamFile("/atac/alignment/bam/no_multimappers/dupmarked/{samp_name}_bowtie2_sorted_nomulti_dupMarked.bam", tag=tags,
                   asMates=TRUE, bigFile=TRUE)
shiftedBamFile <- file.path(outPath, "{samp_name}_shifted.bam")
gal1 <- shiftGAlignmentsList(gal, outbam=shiftedBamFile)

# split the reads into bins for reads derived from nucleosome-free regions, mononucleosome, dinucleosome and trinucleosome. And save the binned alignments into bam files.
objs <- splitGAlignmentsByCut(gal1, txs=txs, genome=genome, outPath = (outPath))
## list the files generated by splitGAlignmentsByCut.
dir(outPath)

#### IN COMMAND LINE
# Convert to bigwig
bamCoverage \
-b /atac/splitBAM/{samp_name}/NucleosomeFree.bam \
-o /atac/splitBAM/{samp_name}/bigwig/NucleosomeFree.bw 

bamCoverage \
-b /atac/splitBAM/{samp_name}/mononucleosome.bam \
-o /atac/splitBAM/{samp_name}/bigwig/mononucleosome.bw

## Heatmap over transcription units
# Nucleosome free
computeMatrix \
scale-regions \
-S /atac/splitBAM/{samp_name}/bigwig/NucleosomeFree.bw \
-R atac/Genome/HSymV2.1.tss.bed \
--beforeRegionStartLength 1000 \
--regionBodyLength 1 \
--afterRegionStartLength 1000 \
--binSize 1 \
--skipZeros \
-o /atac/TSS_coverage/{samp_name}/{samp_name}_NucleosomeFree_tss.mat.gz

plotHeatmap \
-m /atac/TSS_coverage/{samp_name}/{samp_name}_NucleosomeFree_tss.mat.gz \
-out /atac/TSS_coverage/{samp_name}/{samp_name}_NucleosomeFree \
--sortUsing sum \
--yMin 0 \
--yMax 5 \
--endLabel . 

# Mononucleosomes
computeMatrix \
scale-regions \
-S /atac/splitBAM/{samp_name}/bigwig/mononucleosome.bw  \
-R atac/Genome/HSymV2.1.tss.bed \
--beforeRegionStartLength 1000 \
--regionBodyLength 1 \
--afterRegionStartLength 1000 \
--binSize 1 \
--skipZeros \
-o /atac/TSS_coverage/{samp_name}/{samp_name}_mononucleosome_tss.mat.gz 

plotHeatmap \
-m /atac/TSS_coverage/{samp_name}/{samp_name}_mononucleosome_tss.mat.gz \
-out /atac/TSS_coverage/{samp_name}/{samp_name}_mononucleosome \
--sortUsing sum \
--yMin 0 \
--yMax 1.3 \
--endLabel . 


## Combined TSS enrichment
# nucleosome free
computeMatrix scale-regions \
-S /atac/splitBAM/{samp_name}/bigwig/NucleosomeFree.bw /atac/splitBAM/{samp_name}bigwig/NucleosomeFree.bw /atac/splitBAM/{samp_name}/bigwig/NucleosomeFree.bw /atac/splitBAM/{samp_name}bigwig/NucleosomeFree.bw \
-R /atac/Genome/HSymV2.1.tss.bed \
--beforeRegionStartLength 1000 \
--regionBodyLength 1 \
--afterRegionStartLength 1000 \
--binSize 1 \
--skipZeros \
-o /atac/TSS_coverage/NucleosomeFree_tss.mat.gz 

plotHeatmap \
-m /atac/TSS_coverage/NucleosomeFree_tss.mat.gz \
-out /atac/TSS_coverage/NucleosomeFree \
--sortUsing sum \
--endLabel . \
--yMin 0 

# mononucleosomes
computeMatrix scale-regions \
-S /atac/splitBAM/{samp_name}/bigwig/mononucleosome.bw /atac/splitBAM/{samp_name}bigwig/mononucleosome.bw /atac/splitBAM/{samp_name}/bigwig/mononucleosome.bw /atac/splitBAM/{samp_name}bigwig/mononucleosome.bw \
-R /atac/Genome/HSymV2.1.tss.bed \
--beforeRegionStartLength 1000 \
--regionBodyLength 1 \
--afterRegionStartLength 1000 \
--binSize 1 \
--skipZeros \
-o /atac/TSS_coverage/NucleosomeFree_tss.mat.gz 

plotHeatmap \
-m /atac/TSS_coverage/mononucleosome_tss.mat.gz \
-out /atac/TSS_coverage/mononucleosome \
--sortUsing sum \
--endLabel . \
--yMin 0 \
--yMax 1.5

## Coverage across genes
# nucleosome free
computeMatrix scale-regions \
-S /atac/splitBAM/{samp_name}/bigwig/NucleosomeFree.bw /atac/splitBAM/{samp_name}/bigwig/NucleosomeFree.bw /atac/splitBAM/{samp_name}/bigwig/NucleosomeFree.bw /atac/splitBAM/{samp_name}/bigwig/NucleosomeFree.bw \
-R /atac/Genome/HSymV2.1.gtf \
--beforeRegionStartLength 3000 \
--regionBodyLength 5000 \
--afterRegionStartLength 3000 \
--skipZeros \
-o /atac/gene_coverage/NucleosomeFree_genecoverage.mat.gz 

plotHeatmap \
-m /atac/gene_coverage/NucleosomeFree_genecoverage.mat.gz  \
-out /atac/gene_coverage/NucleosomeFree_genecoverage.png \
--sortUsing sum \
--yMin 0

# mononucleosome
computeMatrix scale-regions \
-S /atac/splitBAM/{samp_name}/bigwig/mononucleosome.bw /atac/splitBAM/{samp_name}/bigwig/mononucleosome.bw /atac/splitBAM/{samp_name}/bigwig/mononucleosome.bw /atac/splitBAM/{samp_name}/bigwig/mononucleosome.bw \
-R /atac/Genome/HSymV2.1.gtf \
--beforeRegionStartLength 3000 \
--regionBodyLength 5000 \
--afterRegionStartLength 3000 \
--skipZeros \
-o /atac/gene_coverage/mononucleosome_genecoverage.mat.gz 

plotHeatmap \
-m /atac/gene_coverage/mononucleosome_genecoverage.mat.gz \
 -out /atac/gene_coverage/mononucleosome_genecoverage.png \
 --sortUsing sum \
 --yMin 0

 ### HOMER MOTIF ENRICHMENT

## RStudio
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ChIPpeakAnno")
library(ChIPpeakAnno)

## load data of additional i-cell, neuron and sperm samples (from host lab)
# Peak overlap Venn diagram
icells <- toGRanges("/atac/motifs/icell_1_2_filtered_consensuspeaks.bed", format="broadPeak", header=FALSE)
Neurons <- toGRanges("/atac/motifs/RF_1_filtered_consensuspeaks.bed", format="broadPeak", header=FALSE)
Sperm <- toGRanges("/atac/motifs/sperm_consensuspeaks.bed", format="broadPeak", header=FALSE)
Nematocytes <- toGRanges("/atac/motifs/nemato_1_intersect.bed", format="broadPeak", header=FALSE)

ol <- findOverlapsOfPeaks(icells, Neurons, Sperm, Nematocytes)
ol <- addMetadata(ol, colNames="score", FUN=mean)

makeVennDiagram(ol, fill=c("#009E73", "#F0E442", "#ccccFF", "#ed67a6"), # circle fill color
                col=c("#D55E00", "#0072B2", "#2a52be", "#c071d0" ), #circle border color
                cat.col=c("#D55E00", "#0072B2", "#2a52be", "#c071d0")) # label color, keep same as circle border color




## Extract unique peaks that are present in nematocyte sample but not any other sample
# not in i-cells
bedtools intersect \
-a /atac/motifs/nemato_1_intersect.bed \
-b /atac/motifs/icell_1_2_filtered_consensuspeaks.bed \
-v >/atac/motifs/unique/icell_uniquepeaks.bed
# not in i-cells and neurons
bedtools intersect \
-a /atac/motifs/unique/icell_uniquepeaks.bed \
-b /atac/motifs/RF_1_filtered_consensuspeaks.bed \
-v >/atac/motifs/unique/icell_RF_uniquepeaks.bed
# not in i-cells, neurons and sperm
bedtools intersect \
-a /atac/motifs/unique/icell_RF_uniquepeaks.bed \
-b /atac/motifs/sperm_consensuspeaks.bed \
-v >/atac/motifs/unique/uniquepeaks.bed 


### Create new HOMER Conda environment 
conda create -n homer
conda activate homer
conda install bioconda::bioconductor-deseq2
conda install bioconda::bioconductor-edger
conda install anaconda::wget
conda install r::r-essentials
conda install bioconda::samtools

## install configurehomer.pl
chmod -x /atac/motifs/HOMER/configurehomer.pl
perl /atac/motifs/HOMER/configurehomer.pl -install
nano ~/.bash_profile
	PATH=$PATH:Users/paula/Downloads/atac/motifs/HOMER/bin/
source ~/.bash_profile

## find motifs
findMotifsGenome.pl \
/atac/motifs/unique/uniquepeaks.bed \
/atac/Genome/HSymV2.1.fasta \
/atac/motifs/HOMER/nemato_unique \
-size 200

### GENE ONTOLOGY ENRICHMENT

## RStudio: Peak annotation of unique nematocyte peaks (for GO enrichment)
BiocManager::install("GenomicFeatures")
BiocManager::install("ChIPseeker")
library("GenomicFeatures")
library(ChIPseeker)
library(ChIPpeakAnno)

TxDb <- makeTxDbFromGFF("/atac/Genome/HSymV2.1.gtf",
                format="gtf",
                dataSource=NA,
                organism=NA,
                taxonomyId=NA,
                circ_seqs=NULL,
                chrominfo=NULL,
                miRBaseBuild=NA,
                dbxrefTag=NULL)


# RStudio: Convert unique bed to GRanges
nemato_unique <- toGRanges("/atac/motifs/unique/uniquepeaks.bed", format="broadPeak", header=FALSE)
nemato_unique_peakanno <- annotatePeak(nemato_unique, TxDb = TxDb)

### Save geneIDs from annotation to tsv
nemato_unique_df <- as.data.frame(nemato_unique_peakanno)
write.table(nemato_unique_df$geneId, file='/atac/GO/nemato_genes.tsv', quote=FALSE, sep='\t', col.names = NA)

### Make GO_terms_ids file
install.packages("readr")
install.packages("GO.db")
library("readr")
library(GO.db)
columns(GO.db)
go <- keys(GO.db, keytype="GOID")
select(GO.db, columns=c("GOID","TERM","ONTOLOGY"), keys=go[1:5], keytype="GOID")
df <- select(GO.db, columns=c("GOID","TERM","ONTOLOGY"), keys=go, keytype="GOID")
write_tsv(df, "/Users/paula/Downloads/atac/Ome_NCBI/GO_database/GO_terms_ids")

### Make GAF file in genome folder
perl GAF_creation.pl GO_terms_ids Supplementary_File_2_GO.anno.tsv

### Run Ontologizer in folder with Ontologizer.jar folder on Linux
java -Xmx1G -cp swt.jar:OntologizerGui.jar ontologizer.gui.swt.Ontologizer

## RStudio: plot top 20 enriched GO terms
library(ggplot2)
nemato_GO <- read_tsv("atac/GO/nemato_go.txt")
nemato_GO

# Order by adjusted pvalue
nemato_GO=nemato_GO[order(nemato_GO$p.adjusted),]
nemato_GO$name=factor(nemato_GO$name,levels=nemato_GO$name)

# Extract top 20
nemato_GO_top20 <- head(nemato_GO, 20)
nemato_GO_top20

# Find class for each GO term and add column
Class_nemato <- c("cellular_component", "biological_process", "biological_process","biological_process", "biological_process", "molecular_function", "biological_process", "molecular_function", "biological_process", "biological_process", "biological_process", "biological_process", "molecular_function", "cellular_component", "biological_process", "biological_process","molecular_function", "biological_process", "biological_process", "biological_process")
nemato_GO_top20$Class <- Class_nemato

## Make plot
nemato_plot <- ggplot(nemato_GO_top20, aes(y = name, x = -log10(p.adjusted), fill=Class)) + geom_col(width = 0.7) + scale_y_discrete(limits=rev) + ggtitle("Nematocyte GO Enrichment") +labs(y= "GO term", x = "-log10(Padj)") + scale_fill_manual(values = c("#ff6666", "#8080ff", "#33cc00"))
nemato_plot
