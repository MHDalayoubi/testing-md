# Introduction:

The article “Common applications of next-generation sequencing technologies in genomics research” by Chien-Yueh Lee, Yu-Chiao Chiu, Liang-Bo Wang, Yu-Lun Kuo, Eric Y. Chuang, Liang-Chuan Lai, Mong-Hsun Tsai is a review of several next-generation sequencing technologies; which are the paradigm shifting and  new emerging technologies that revolutionized the sciences of gene-sequencing in terms of speed, effectiveness and accuracy, it presents a quick overview of the available technologies, a guide line for a pipeline data processing and it makes suggestions for the type of tools to be used in genomics transcriptomics and small RNA research.
The article gives a quick description of each one of the stages of gene-sequencing and list several recommended technologies suited for that particular stage along side a quick description of the underlying methods and algorithms that these technologies depends on. 

# DNA sequencing data analysis:

## Preprocessing Procedure:
The article starts by describing the DNA sequencing preprocessing procedure at this stage terabytes of raw image data generated from the sequencing  platform go through several preprocessing operations the general steps are: image analysis base-calling BCL conversion(optional) sample demultiplexing (optinal) and variant analysis.
during the image analysis procedure the software identify clusters and outputs positions and intensities and estimates the cluster noise, In the base calling step the software uses the cluster intensities and noise estimations to produce the base reading from each cluster and wither or not that base read can pass filters, the generated file are also converted to the fastaq format(if necessary) for further analysis and the samples are demultiplexed by using their index (called bar-codes) if they were multiplexed. programs like CASAVA and Bioscope do image preprocessing


## Read alignment:
The article mentions that this step is utilized by the read alignment tool to align several hundreds or thousands of reads back to an existing reference genome and it also lists few programs that can do this alignment such as MAQ,BFAST,Novoalign,Burrows-wheelr aligner and SOAP3.

## De novo assembly:
The de novo assembly step focuses on grouping short reads into into significant contigs and then grouping these contigs into scaffolds to reconstruct the original DNA, the article mentions few programs and algorithm to do de novo assembly like greedy graph approach used by VCAKE and the overlap/layout/consensus method and the reduction of the of path complexity of the brjuin graph that is used by the Velvet program	

## Signal nucleotide variant (SNV) detection:
Next in the analytical pipeline the detection of single nucleotide variation, programs that are mentioned in the article that do single nucleotide variation calling include GATK,SAMtolls,VarScan,SomaticSniper,jointSNVMix.
## Structural variation detection:
"structural variation ”generally implies a genetic change that is approximately 1kb to 3Mb in length,breakdancer ,variationhunter,SVDetect and PEMer are programs for detecting such structural variations
# RNA sequencing data analysis:
The next generation sequencing techniques has also been applied to study RNA transcripts typically called RNA-seq or transcriptome-seq.
RNA-seq goes through the same preprocessing steps as in DNA sequencing, RAN data can used for de novo assembly of transcriptome, expression profiling analysis and varint calling and transcriptomic epigenetics.

## de novo transcriptome assembly:
De novo assembly of RNA-seq data provides an overview and extracts clues to the “transcriptome”,some of the programs that uses brujin graph-based transcriptome assemblers include programs like TrinityiTrans-AbySS  and Oases, the long RNA contigs can be annotated with BLAST for further research and analysis.

## exprssion profiling analysis:
Currently the most famous use of RNA-seq is to profile gene expression levels and identify differently expressed transcripts among groups of samples, analysis of RNA-seq data typically includes mapping reads against a reference, per transcript counting and statistical testing for differential expression. Alignment programs include bowtie,bowtie2 and tophat.
 A python tool HTSeq extracts read counts for each transcript which can be then used for further statistical testing.
Many tools has been developed for normalization, bias correcting and statistical testing of RNA-seq read counts.
DESeq,bayseq and edgeR are three frequently used tools for detecting differential expression of transcripts among a set of samples, in addition to the mentioned tools Cuuflinks package can be used for assembling transcripts.

## Variant calling and transcriptomic epigenetics:
RNA sequencing provides a cost-effective way for detecting coding variants, in addition to variant calling with SAMtools, a mutation mapping pipeline has been proposed with allele frequency distance calculation signal processing and candidate SNP identification.
Research efforts has been directed to areas like transcription start site-associated RNAs,
promoter-associated RNAs, transcription-initiation RNAs (tiRNAs), and many others that may facilitate investigation into complex transcriptional regulation. 

# Small RNA sequencing:
Because of it’s highly accurate results Next generation small RNA sequencing technology
has become essential for sRNA discovery and profiling

## General workflow:
Although the sRNA sequencing workflow depends on the application and the platform one use some major steps are generally always followed, a library consisting of cDNA is obtained after sequencing,
reads containing sequences of adapters should be trimmed off by using tool kits such as flicker, FASTX clipper, scythe or cutadapt, next reads with low quality should be removed using tools like FASTQ quality filter or the NGS QC toolkit, after that tools like FASTQC or qrqc are used to check quality statistics finally the filtered reads should be validated by aligning to a reference genome database. 
Some of the databases commonly used include Rfam and miRBase. These discarded reads can be kept of disposed off depending on the purpose of the sequencing.

## Small RNA prediction:
The discovery of new sRNA is made easier by NGS technology there are many tools for miRNA discovery like miRTRAP,MIReNA,miRExplorer,miRAnalyzer,miRDeep/miRDeep2 and DSAP and for miRNA prediciton in plants miRDeep-p.

## miRNA characterization:
Currently micro arrays, quantitative real-time RT-PCR and sRNA-seq are all widely used for miRNA characterization and their attributes have been described in detail, in miRNA characterization
some reoccurring systematic biases towards different protocols can be noticed due to different usage of RNA ligase and this bias can be eliminated by using different adapters
	
	












