# NMDtxDB: Data-driven identification and annotation of human NMD target transcripts

We employ a data-driven approach aimed at novel transcript discovery. Our primary objective is to identify isoforms that can activate the Nonsense-Mediated RNA Decay (NMD) machinery that are missing from reference databases and facilitate the access to this resource by a broad audience. Currently, the database offers evidence of premature stop codons (PTC) based on transcript structure. Accordingly, the feature (gene or transcript) abundances is also present. In that sense, NMDtxDB, in its current form, does not aim to replace the reference annotations, but offer an initial infrastructure to collect evidence of NMD activation and further characterize the rules of NDM. It's also important to present the resource on an easy-to-use platform, that can replicated by others. Thus, we hope this resource is used by a broad audience of researchers. This should help researchers to pinpoint relevant features from their genes of interest, its regulatory role in post-transcriptional regulation. 

## Introduction to NMDtxDB

### Content description

NMDtxDB is a database dedicated to exploring transcripts that may be targeted by the NMD pathway. This web-application currently collects the results obtained for four cell lines depleted for NMD factors, as further detailed in our below and in our research manuscript. 

The NMDtxDB is constructed using data from cells lines depleted for key NMD-factors, offering a comprehensive insight into these NMD depletion conditions. We've also incorporated Nanopore sequencing, which improves on the detection of splicing isoforms, along with three CDS sources, that allow us to infer coding sequence (CDS) that are currently not annotated on the reference database. The figure below has a simplified overview of our workflow:

![NMDtxDB Workflow](intro.jpeg)

Our work is open-source and available at GitHub. The Shiny App source is available at https://github.com/dieterich-lab/NMDtxDB. The workflow is available at https://github.com/dieterich-lab/nmd-wf. Both source codes are licensed under MIT.  We welcome contributions and feedback from the community to further improve  source code, and are considering new data sets and data modalities to enhance the database.   

A full description of the database and its construction can be found in our [manuscript](http://dx.doi.org/10.1261/rna.080066.124).

### Computational Workflow

We've developed a unique computational workflow that addresses the need for transcriptome construction and premature stop codon (PTC) annotation. This workflow enables users to compare PTC status with transcript expression, helping identify unannotated NMD targets and providing a more comprehensive understanding of NMD activation.

We highlight the challenges of building transcriptomes from RNA-seq dataset. The context of NMD makes this task even more challenging, as the pathway leads to the degradation of transcripts, which suppresses their detection under control conditions. The use of long read sequence alleviates some issues, but comes with it own limitations. Overall, we suggest users to be cautious when interpreting the database results, and taking the database evidences a form for prioritizing targets for experimental validation. 

The multiple source for CDS allow users to explore for further characterization. Of note, besides the canonical CDS source, imported from Ensembl, the other Ribo-seq and OpenProt CDS sources are from external high-throughput datasets, and so require further checks. The Ensembl source is an adaptation of Ensembl-annotated CDS to novel transcripts. 

Our approach is particularly useful in studying the impact of alternative splicing on NMD, as not all potential NMD substrates are degraded equally, and some may even escape degradation.

### Front-end Application

Our web application is designed with a strong focus on simplicity and usability, making the database accessible to a wider audience. The initial feedback from our prototypes, presented at prominent RNA conferences, has been extremely encouraging.

The application is built using Shiny and Golem and is hosted with ShinyProxy via docker containers. In principle, this infra-structure could be used elsewhere, although we currently offer little support to external servers. We recently improved the user experience with a tutorial to guide users through the web application. Additionally, we optimized the front-end initialization time and database query performance to ensure a seamless user journey.

Please give us some **Feedback** and help us improve further!
