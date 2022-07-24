# EpiC Dog: Epigenome Catalog of the Dog

### Integrative mapping of the dog epigenome: reference annotation for comparative inter-tissue and cross-species studies

Keun Hong Son<sup>1,2,3†</sup>, Mark Borris Aldonza<sup>1,2,3†</sup>, A-reum Nam<sup>1,2,3†</sup>, Kang-Hoon Lee<sup>1,3</sup>, Jeong-Woon Lee<sup>1,2,3</sup>, Kyung-Ju Shin<sup>1,3</sup>, Keunsoo Kang<sup>4</sup>, and Je-Yoel Cho<sup>1,2,3*</sup>

<sup>1</sup> Department of Biochemistry, College of Veterinary Medicine, Seoul National University, Seoul, Korea<br>
<sup>2</sup> Comparative Medicine and Disease Research Center (CDRC), Science Research Center (SRC), Seoul National University, Seoul, Korea<br>
<sup>3</sup> BK21 PLUS Program for Creative Veterinary Science Research and Research Institute for Veterinary Science, Seoul National University, Seoul, Korea<br>
<sup>4</sup> Department of Microbiology, College of Natural Sciences, Dankook University, Cheonan, Korea<br>

<sup>†</sup> These authors contributed equally to this work as co-first authors: newhong@snu.ac.kr, borris@snu.ac.kr and arbjlvz@snu.ac.kr<br>
<sup>*</sup> Corresponding author: jeycho@snu.ac.kr<br><br>
<img src="https://user-images.githubusercontent.com/108702438/179219246-b79c13cb-7ca7-4de0-a1e1-56211ad574e8.jpg" width="700"><br>

## Introduction

&nbsp;The domestic dog has become a valuable model in exploring multifaceted diseases and biology important for human health. Large-scale dog genome projects produced high-quality draft references but still lack comprehensive annotation of encoded functional elements. Through the integrative next generation sequencing of transcriptomes and chromatin accessibility paired with DNA methylome profiling of 11 adult tissue types, implemented in a cross-species approach, we generated a reference epigenome of a domesticated dog. Using genome orthologues and synthenies, we deciphered the dog’s epigenetic code by defining distinct chromatin states, allowing for genome-wide, integratable data production. We then characterized somatic super-enhancer landscapes and showed that genes mapped on these regions are associated with a broad range of biological and disease traits and are traceable to their tissue-of-origin. Ultimately, we delineated conserved epigenomic changes at the tissue- and species-specific resolutions. Our study provides an epigenomic blueprint of the dog for comparative biology and medical research.

<img src="https://user-images.githubusercontent.com/108702438/180374875-61f24e55-a062-444e-9b27-5ec9e88e152e.jpg" width="700">

## Manuscript and Data availability
Paper (preprint): [bioRxiv](https://www.biorxiv.org/content/10.1101/2022.07.22.501075v1)<br>
Data repository: [Acession: GSE203107](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203107) (Openning soon)<br>
Data browser: [the UCSC_trackhub](http://genome.ucsc.edu/s/snu-cdrc/dog-reference-epigenome) (Optimizing due to SSL certificate problem in server)<br><br>

## Processed data repository
### 1. RNA contigs
&nbsp;Transcript contig means the region of the transcript covered by RNA-seq. Therefore, it is possible to study known and novel positions where transcripts are expressed without dependency on genome annotations. Strand-specific contig regions were defined using 11 tissues RNA-seq data through the approach described by Djebali et al.

Download: [RNA_Contigs](https://drive.google.com/drive/folders/1h362ZtVfhkabvMF7L2lGbicsTpQJOJNi?usp=sharing)<br>
#### Reference:
- A comparative encyclopedia of DNA elements in the mouse genome, Nature, 2014
- Landscape of transcription in human cells, Nature, 2012<br><br>



### 2. Comparative inter-tissue transcriptomics
&nbsp;<br>
<img src="https://user-images.githubusercontent.com/108702438/179224395-4ff8743f-2a80-44a1-bde2-0776b7cc8ace.jpg" width="600">

Download: [Categorized genes according to tissue specificity](https://drive.google.com/drive/folders/11QTuDknYqQmk_rjDw50PXp0VJN_TetEJ?usp=sharing)<br>
#### Reference:
- Tissue-based map of the human proteome, Science, 2015<br><br>



### 3. Conserved tissue-specific and species-specific genes
&nbsp;<br>
<img src="https://user-images.githubusercontent.com/108702438/177751739-30a13ea6-87ac-4e0f-b260-910bbe31fec3.jpg" width="600">

Download: [Conserved gene lists](https://drive.google.com/drive/folders/1w6Fh1ytLAoulALwxpQ3cFyFstP7vBWlR?usp=sharing)<br>
#### Reference:
- A comparative encyclopedia of DNA elements in the mouse genome, Nature, 2014
- Comparative transcriptomics in human and mouse, Nature Reviews Genetics, 2017<br><br>



### 4. Chromatin states
&nbsp;To advance the functional annotation of the dog genome, we produced integrated maps of histone modifications-informed, genome-wide 13-chromatin state model in 11 dog tissues. We defined the dog genome as having a core set of five histone H3 modification marks: H3K4me3, H3K4me1, H3K27ac, H3K27me3, and H3K9me3—marks well-known to have specific depositions on particular genomic regions and molecular signal associations (i.e., promoters, enhancers, heterochromatin, Polycomb repressive domains, etc).<br>
<img src="https://user-images.githubusercontent.com/108702438/179326386-5e1ea129-0e0c-4c05-a238-984e3d2ef548.jpg" width="300"><br>
<img src="https://user-images.githubusercontent.com/43947916/177108534-54b168a3-280b-49d2-8bc8-a33286eb2c8c.jpg" width="300">

Download: [Chromatin states](https://drive.google.com/drive/folders/1fPqttRt1x6f8RDYpmr8edwPoZyJF_fcv?usp=sharing)<br>
#### Reference:
- ChromHMM: automating chromatin-state discovery and characterization, Nature methods, 2012
- Integrative analysis of 111 reference human epigenomes, Nature, 2015<br>



### 5. Super-enhancer (SE)
&nbsp;To further probe tissue identity and function based on H3K27ac signals, a strong indicator of active promoter and enhancer states, we characterized super-enhancers (SEs) landscapes in the dog genome across multiple tissues.<br>
<img src="https://user-images.githubusercontent.com/108702438/179326388-c5e39a80-42c8-4034-9d8c-ce040ed304be.jpg" width="300"><br>

Download: [Super enhancer](https://drive.google.com/drive/folders/1DPk7Z39NHW1TTFaNz1Qief1taOxyDhbu?usp=sharing)<br>
#### Reference:
- Master Transcription Factors and Mediator Establish Super-Enhancers at Key Cell Identity Genes, Cell, 2013
- Selective Inhibition of Tumor Oncogenes by Disruption of Super-enhancers, Cell, 2013<br>



### 6. Common and tissue-specific differentially methylated regions
&nbsp;Methylation of cytosines in DNA is a prototypic, stable, nearly universal mechanism of the mammalian epigenome. In domestic dogs, DNA methylation studies have been performed yet still lack epigenome-scale resolution. So far, public resources of functionally annotated dog genomes (i.e., BarkBase and DoGA) do not include methylome data. To profile global DNA methylome landscape of the dog, we performed genome-wide MBD-seq experiments on 11 somatic tissues. In these experiments, captured and enriched genomic DNA fragments covering a CpG are used to assay the total amount of methylation for a locus about the size of the fragments, which dictate the resolution of association signals.<br>
<img src="https://user-images.githubusercontent.com/108702438/179228475-55caeda4-d3af-40c1-9ef2-fef25a8e33dc.jpg" width="400">

Download: [Location of CMRs and tsDMRs](https://drive.google.com/drive/folders/106a95qTG8qN84Jp8FyQGq9dkxK1vMoYH?usp=sharing)
#### Reference:
- MethylAction: detecting differentially methylated regions that distinguish biological subtypes, Nucleic Acids Res., 2016<br><br>
