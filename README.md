# Data and scripts used for analyses of diatom mitochondrial genomes

Recurrent loss, horizontal transfer, and the obscure origins of mitochondrial introns in diatoms (Bacillariophyta)

Wilson X. Guillory, Anastasiia Onyshchenko, Elizabeth C. Ruck, Matthew Parks, Teofil Nakov, Norman J. Wickett, and Andrew J. Alverson

Data
======
- Genome sequences are available in GenBank under accessions MG148339, MG182051, MG271845, MG271846, and MG271847
- Alignments and workflows for phylogenetic trees are available in `intron_phylogenies/cox1` and `intron_phylogenies/rnl`

Scripts
======

## Drawing intron maps (Figs. 1 and 3)
### `map_introns.pl`

#### Call the program
`map_introns.pl [options] dna_alignment.fasta`

#### Program options
```
   --scale     - proportion to scale width of gene line (default: none)
   --bar_width - width of scale bar (default: 100 nt)
   --header    - print header at top of page (filename minus file extension) (default: yes)
   --ps2pdf    - call ps2pdf to convert postscript to pdf (default: yes)
```

#### Format of FASTA headers
```
>Citrullus 500 600 700m
ACGTACGT . . . 
>Cucumis 500 600 700
ACGTACGT . . . 
>Cucumis 500 600 700
ACGTACGT . . . 
```

#### Example usage
`map_introns.pl --scale=0.15 --bar=500 --noheader cox1.fa`

#### Output
A postscript file and, if Ghostscript is installed (see Notes below), a pdf file of the intron maps

#### Notes
- Introns are mapped as filled triangles and missing introns (e.g. 700m) are mapped as open triangles
- It is assumed that intron locations do not consider any gap characters that might exist in the alignment; that is, the script will adjust intron locations to account for gap characters ("-") that precede them in that sequence
- By default, each nucleotide is one pixel width, but a letter-sized palette is only 612 pixels wide, so the --scale option is almost always necessary

#### Dependencies
- Perl
- Optional: `ps2pdf` is part of the [Ghostscript](https://www.ghostscript.com/) package, which can be installed with [Homebrew](https://brew.sh/); if not installed, the postscript is still written out

## Compiling BLAST hits for intron phylogenies
This script was part of the workflow detailed in the README files in the `intron_phylogenies/cox1` and `intron_phylogenies/rnl` directories

### `blast_table_2_subject_fasta.py`

#### Call the program
`blast_table_2_subject_fasta.py -help`

#### Program options
```
usage: blast_table_2_subject_fasta.py [-h] [-e EMAX] [-q QCOVS] blast

This script reads a BLAST table and calls eutils to extract subject matches
from NCBI

positional arguments:
  blast                 BLAST output: -outfmt '6 std qcovs'

optional arguments:
  -h, --help            show this help message and exit
  -e EMAX, --emax EMAX  maximum e-value
  -q QCOVS, --qcovs QCOVS
                        maximum e-value
```

#### Example usage
`blast_table_2_subject_fasta.py Ulnaria_acus_intron_cox1_g2.2_A.blastn -q 50 -e 1e-12 > subject_hits.fa`

#### Output
- A FASTA file of the subject matches that meet the specified qcovs and e-value cutoffs
- If the same subject has multiple hsp's, those hsp's are combined into a single large subject sequence. This can have unintended consequences if, for example, there are many distantly spaced (i.e., different) hsp's rather than a single match that was split into multiple hsp's

#### Notes
- The script reads the following BLAST format: `-outfmt '6 std qcovs'`

#### Dependencies
- [E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK179288/) installed and in your PATH
- Python 3


