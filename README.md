# T2T Sequences

## Overview

To identify Telomere-2-Telomere (T2T) sequences based on:

- Presence of a user-defined motif on both ends of the sequence.
- Sequences that cover the full length of a reference chromosome and have the motif on both ends.

## Dependencies
- [Minigraph](https://github.com/lh3/minigraph)
- [SeqKit](https://github.com/shenwei356/seqkit)
- [tidk](https://github.com/tolkit/telomeric-identifier)

__#Note:__ In case provided pre-compiled tools do not work, kindly install these tools and import the path of executables to your working directory

## Install
```
git clone https://github.com/prasad693/T2T_Sequences.git
cd T2T_Sequences
chmod 777 *
```
Upon changing the permission, include the scripts path in your working directory
```
export PATH=$PATH:/folder_path/T2T_Sequences
```
## Usage
```
T2T_chromosomes.sh -a assembly.fasta -r reference.fasta -o output_prefix -m TTAGGG -t 10
```

Options:
```
          -a file         Assembly fasta file [required]
          -r file         Reference fasta file [optional]
          -o file         Output file prefix [default: T2T_sequences]
          -m string       Telomere motif to search [default: TTAGGG]
          -t int          Number of threads used for minigraph, seqkit [default: 10]
          -h,             Show this help message
```

## Output
- `<assembly.fasta>_motif_<output_prefix>.txt`     : Sequences with telomere motifs on both ends
  
In case reference fasta is provided:

- `<assembly.fasta>_alignment_<output_prefix>.txt` : Sequences covering reference chromosomes and contain telomere motifs on both ends

Output file contains:
- Sequence Name
- Sequence Length
- Number of N's
  
In case reference fasta is provided 

- Reference chromosome
- Chromosome length
  
