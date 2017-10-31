# MG-Cover #

## Overview: ## 

MG-Cover is a taxa specific Binning to Coverage Calculating Pipeline that is complimentary to the metagenomics annotation pipeline MG-RAST. It addresses the gap between open source and proprietary metagenomics pipelines by providing a binning and coverage pipeline. MG-Cover extracts both the annotated sequences and their pair-end sequences by utilizing BioPythons fastq parsing algorithms. Additionally, if no reference genome was provided we would use the taxa annotation to query NCBI for a reference genome. Burrow-Wheeler Aligner mem algorithm was utilized to align the pair ends of taxa specific sequences with the reference genome. Finally, the depth of coverage was calculated using Samtools Depth command. MG-Cover is a taxa specific Binning to Coverage Calculating Pipeline that is complimentary to MG-RAST. It addresses the gap between MG-RAST and proprietary metagenomics pipelines by providing a binning and coverage pipeline. 

## How to Install: ##
1. Check to see if Python 2.7 is installed
2. Copy or Download this github repo.
3. Run the command: **pip install -r requirements.txt**
    

## Command Line: ## 

    python mgcover.py --id [MG-RAST job id] --taxa [taxa being extracted] --key [user web key] --pair-end [fastq file containing pair end sequences]

MG-Cover is a command line program that was implemented in python 2.7. Here is an example on how to run MG-Cover from the command line currently. 

These are the flags available for users:
- (required) **--key**: users MG-RAST web key 
- (required) **--taxa**: names of the taxa(s) to extract. For multiple taxas separate with ":". –id
- (required) **--pair-end**: The file containing pair end sequences of the MG-RASTs annotated sequences.
- (required) **--id**: the MG-RAST job run id we are extracted our annotated sequences from.
- (optional) **--ref-genome**: Fasta file containing reference genome
- (optional) **--email**: Email account we report to NCBIs API 

## MG-Cover Workflow: ##

MG-Cover workflow is broken down into a four step process: Prepping, Processing, Aligning, and Coverage. 

- Prepping: The prepping phase of MG-Cover is when the program is fetching all of the externally needed files via interactions with the MG-RAST API. During these interactions with the API MG-Cover requests an annotated table of taxa specific sequences. The taxa being used is specified by User input.
- Processing: The processing phase of MG-Cover is extracting and processing DNA sequences from the annotated tables and pair-end fastq files. Additionally, if there are no reference genome MG-Cover searches for one from NCBI nucleotide database.
- Aligning: The aligning phase of MG-Cover is when taxa specific annotated and pair-end sequences are aligned with BWA mem algorithms. The alignment results are stored in files in the current working directory of the user.
- Coverage: The coverage phase of MG-Cover is responsible for calculating genome coverage depth of taxa specific sequences. The calculated coverage depth is stored in a tab-separated file ending with .coverage notation.

## MG-Cover Dependencies: ##

MG-Cover is dependent upon the MG-RAST API, BioPython, Burrow-Wheeler Aligner (BWA) and the Samtools kit. The MG-RAST API is used to extract an annotated sequence table based upon a specific Taxa a user specifies. BioPython is used to extract the DNA sequences from the annotated sequence table MG-RAST provided. BWA is utilized to align taxa specific sequences with their respective genome. Samtools is utilized to manipulate the results of BWA and then calculate the coverage depth of taxa specific sequences.

## Summary: ##

MG-Cover is a taxa specific Binning to Coverage Calculating Pipeline that is complimentary to MG-RAST. It addresses the gap between MG-RAST and proprietary metagenomics pipelines by providing a binning and coverage pipeline. 
