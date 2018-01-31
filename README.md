# Scripts used for analyses of diatom mitochondrial genomes

Recurrent loss, horizontal transfer, and the obscure origins of mitochondrial introns in diatoms (Bacillariophyta)

Wilson X. Guillory, Anastasiia Onyshchenko, Elizabeth C. Ruck, Matthew Parks, Teofil Nakov, Norman J. Wickett, and Andrew J. Alverson

### `map_introns.pl` usage

#### Call the program
`map_introns.pl [options] dna_alignment.fasta`

#### Program options
```
   --scale     - proportion to scale width of gene line (default: none)
   --bar_width - width of scale bar (default: 100 nt)
   --header    - print header at top of page (filename minus file extension) (default: yes)
   --ps2pdf    - call ps2pdf to convert postscript to pdf (default: yes)
```

#### Format of FASTA headers:
```
>Citrullus 500 600 700m
ACGTACGT . . . 
>Cucumis 500 600 700
ACGTACGT . . . 
>Cucumis 500 600 700
ACGTACGT . . . 
```

#### Example usage:
`map_introns.pl --scale=0.15 --bar=500 --noheader cox1.fa`

#### Notes:
>
- Introns are mapped as filled triangles and missing introns (e.g. 700m) are mapped as open triangles
- It is assumed that intron locations do not consider any gap characters that might exist in the alignment; that is, the script will adjust intron locations to account for gap characters ("-") that precede them in that sequence
- By default, each nucleotide is one pixel width, but a letter-sized palette is only 612 pixels wide, so the --scale option is almost always necessary
- ps2pdf is part of the [Ghostscript](https://www.ghostscript.com/) package, which can be installed with [Homebrew](https://brew.sh/)
