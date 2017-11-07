# ChIP-seq Hands-on
### Introduction
#### Goal
The aim is to :
  * have an understanding of the nature of ChIP-Seq data
  * perform a complete analysis workflow including quality check (QC), read mapping, visualization in a genome browser and peak-calling. Use command line and open source software for each each step of the workflow and feel the complexity of the task
  * have an overview of possible downstream analyses
  * perform a motif analysis with online web programs

#### Summary
This training gives an introduction to ChIP-seq data analysis, covering the processing steps starting from the reads to the peaks. Among all possible downstream analyses, the practical aspect will focus on motif analyses. A particular emphasis will be put on deciding which downstream analyses to perform depending on the biological question. This training does not cover all methods available today. It does not aim at bringing users to a professional NGS analyst level but provides enough information to allow biologists understand what DNA sequencing practically is and to communicate with NGS experts for more in-depth needs.

#### Dataset description
For this training, we will use a dataset produced by Myers et al (see on the right of the page for details on the publication) involved in the regulation of gene expression under anaerobic conditions in bacteria. We will focus on one factor: FNR.

## Downloading ChIP-seq reads from NCBI
**Goal**: Identify the dataset corresponding to the studied article and retrieve the data (reads as FASTQ file) corresponding to one experiment and the control.
**Related VIB training**: Downloading NGS data from NCBI

### 1 - Obtaining an identifier for a chosen dataset
Within an article of interest, search for a sentence mentioning the deposition of the data in a database. Here, the following sentence can be found at the end of the Materials and Methods section:
*"All genome-wide data from this publication have been deposited in NCBI’s Gene Expression Omnibus (**GSE41195**)."*
We will thus use the **GSE41195** identifier to retrieve the dataset from the **NCBI GEO** (Gene Expression Omnibus) database.

NGS datasets are (usually) made freely accessible for other scientists, by depositing these datasets into specialized databanks. Sequence Read Archive (SRA) located in USA hosted by NCBI, and its european equivalent European Nucleotide Archive (ENA) located in England hosted by EBI both contains **raw reads**.

Functional genomic datasets (transcriptomics, genome-wide binding such as ChIP-seq,...) are deposited in the databases Gene Expression Omnibus (GEO) or its European equivalent ArrayExpress.

### 2 - Accessing GSE41195 from GEO
1.  The GEO database hosts processed data files and many details related to the experiments. The SRA (Sequence Read Archive) stores the actual raw sequence data.
2. Search in google **GSE41195**. Click on the first link to directly access the correct page on the GEO database.
3. This GEO entry is a mixture of expression analysis and chip-seq. At the bottom of the page, click on the subseries related to the chip-seq datasets. (this subseries has its own identifier: **GSE41187**).
4. From this page, we will focus on the experiment **FNR IP ChIP-seq Anaerobic A**. At the bottom of the page, click on the link "GSM1010219 - FNR IP ChIP-seq Anaerobic A".
5. In the new page, go to the bottom to find the SRA identifier. This is the identifier of the raw dataset stored in the SRA database.
6. Copy the identifier **SRX189773** (do not click on the link that would take you to the SRA database, see below why)

### 3 - Downloading FASTQ file from the ENA database
Although direct access to the SRA database at the NCBI is doable, SRA does not store the sequence FASTQ format. In practice, it's simpler (and quicker!!) to download datasets from the ENA database (European Nucleotide Archive) hosted by EBI (European Bioinformatics Institute) in UK. ENA encompasses the data from SRA.

1. Go to the EBI website. Paste your SRA identifier (SRX189773) and click on the button "search".
2. Click on the first result. On the next page, there is a link to the FASTQ file. For efficiency, this file has already been downloaded and is available in the "data" folder (SRR576933.fastq)
