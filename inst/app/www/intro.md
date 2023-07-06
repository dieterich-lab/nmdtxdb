# NMDtxDB

## Description

NMDtxDB is under development.

The NMDtxDB is a database of transcripts that may activate the Nonsense-Mediated RNA Decay (NMD) pathway. The primary goal of the database is to allow non-NMD specialist to check for any concrete evidence of NMD activation for their gene of interest, and so allowing for decision-making whether NMD could be a regulatory process in play in their studies.The database is constructed using data derived from NMD depleted cell lines, providing a new and comprehensive dataset for NMD depletion conditions. We have also incorporated Nanopore sequencing for a subset of the dataset, which improved the detection of the splicing isoforms and CDS from sources, which expand the coding sequence space. The CDS sources are the following:

  - CCDS from Ensembl
  - [Gao et al.](https://doi.org/10.1038/nmeth.3208)
  - [Zhang et al.](https://doi.org/10.1038/s41467-017-01981-8)
  - [OpenProt](https://doi.org/10.1093/nar/gkaa1036)
  - Transdecoder

## Computational workflow

The project encompasses a novel computational workflow that tackles the transcriptome construction and premature stop codon (PTC) annotation. This approach allow users to compare the PTC status in context to the transcript expression and identify unannotated NMD targets. Importantly, the computational workflow will be made open-source, enabling users to generate their own databases from in-house datasets, thereby promoting collaboration and further research in the field.

## Front-end application

This web application aims to enhance accessibility to the database. Its design prioritizes simplicity and ease of use. Initial feedback from prototypes at the "Complex life of RNA" conference 2022 in Heidelberg and the RNA Society conference 2023 in Singapore was encouraging. The application is based on Shiny and Golem, and hosted with ShinyProxy via docker containers. We are currently working on enhancing the user experience by developing a tutorial that explains how to effectively utilize the web application and interpret the provided data and plots. Additionally, efforts are underway to optimize the front-end initialization time and improve database query performance. 

Please submit **Feedback** by clickign on that button and filling the form.

## ABSTRACT
Nonsense-mediated RNA decay (NMD) is a cellular mechanism that detects and degrades mRNA molecules with premature stop codons (PTC). NMD is essential for normal cellular processes because it removes transcripts that would generate truncated proteins, causing deleterious functional impact. Genetic mutations resulting in stop codons are the most well-studied cause of PTC, but they are not the only cause, and the role of alternative splicing on NMD is poorly characterized. Most importantly, not all potential NMD substrates are degraded with the same efficiency, and some might completely escape the degradation pathway. Currently, in-depth, evidence-based information on a transcript's NMD-activation status is largely missing. Databases such as EnsEMBL provide annotations of PTC-containing transcripts that consider simple rules and ignore experimental data. This calls for the NMD-targeted transcriptome to be described in a systematic and multi-layered way.
In this work, we propose a workflow for the comprehensive identification of NMD substrates transcriptome-wide. We generated 4 cell lines depleted of key NMD factors, and obtained RNA libraries that were sequenced with Nanopore direct-RNA and Illumina sequencers. We combined both datasets to generate a de novo transcriptome containing 302,889. Furthermore, we use evidence of CDS to further annotate the transcripts, and 114,609 have a compatible CDS from the canonical CDS, OpenProt or experimental resources. Of these, 56,732 (50%) have Nanopore read support or 28,408 (25%) are unnanoted isoforms. We obtained a total of 38,487 (34%) isoforms that elicit the 50nt rule, of which 9,205 (8%) are novel. Taken together, our approach enables the robust identification of NMD substrates and the re-annotation of transcripts using multiple CDS sources, allowing users to query to analyze transcript expression in context with transcript PTC status. These results are available at our database https://shiny.dieterichlab.org/app/nmd_transcriptome, which will make these findings available to a wider audience, in addition to providing the source code for the workflow and web-application. Overall, the database can contribute to a better understanding of the basic biology of gene expression, as well as prioritizing new strategies to assess the role of the NMD pathway in disease.
