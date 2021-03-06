## Compiling BLASTN hits and alignments for each _cox1_ intron

1. BLAST each intron individually to GenBank
`for i in *.fa;do test -e ${i%.fa}.blastn && continue; blastn -word_size 11 -reward 2 -penalty -3 -gapopen 5 -gapextend 2 -query $i -remote -db nr -outfmt '6 std qcovs' -evalue 1e-10 > ${i%.fa}.blastn;done &`

1. Retrieve BLAST hits with e-value <= 1e-12 and qcovs > 50
`for i in *.blastn;do blast_table_2_subject_fasta.py $i -q 50 -e 1e-12 > ${i%.blastn}.subjects.fa;done`

1. Manually move files into group-specific intron folders (e.g., group_g2A, group_g2B, etc.)

1. Concatenate subject sequences for each intron cluster
`for f in group_g*;do cat $f/*.fa > $f/${f}_sbjcts.fa;done`

1. Concatenate query (i.e, diatom) sequences for each intron cluster (run within the group directory)
`for i in *.fa; do if [[ $i =~ subjects || $i =~ group ]];then continue;else cat $i;fi;done > group_g2E_queries.fa`

1. Run CD-HIT to remove redundant subject sequences (e.g., some queries hit to both the genome and a separate cox1 gene accession AND each sequence from an intron cluster was BLASTed individually)
`cd-hit-est -i group_g2E_sbjcts.fa -o group_g2E_sbjcts_filtered.fa -c 0.99 -d 20000`

1. Concatenate query and subject sequences
`cat group_g2E_sbjcts_filtered.fa group_g2E_queries.fa > group_g2E_all.fa`

1. Remove redundant sequences from the full set query and subject sequences
`cd-hit-est -i group_g2E_all.fa -o group_g2E_all_filtered.fa -c 1 -d 20000`

1. For **group_g2A** only, reduced cutoff to `-c 0.99` to remove a few highly similar sequences:
	- \>Chattonella_marina_intron_cox1_g2.1_A
	- \>Navicula_ramosissima_intron_cox1_g2.2_A

>Note: both of these species are still represented in the final alignment: `group_g2A_all_filtered_0.99.fa`

### Final sequence counts for each intron group
`for i in group_g*/*_all_filtered*.fa;do echo '| ' $i '|' | tr -d '\n' | perl -pwe 's/\s{2,}/ /';grep '>' $i | wc -l | perl -pwe 's/\s{2,}/ /' | tr -d '\n';echo " |";done`

|          Group           | Num sequences |
|-------------------------------------|----|
| group_g1/group_g1_all_filtered.fa   | 2 |
| group_g2A/group_g2A_all_filtered.fa | 16 |
| group_g2A/group_g2A_all_filtered_0.99.fa | 14 |
| group_g2B/group_g2B_all_filtered.fa | 11 |
| group_g2C/group_g2C_all_filtered.fa | 3 |
| group_g2D/group_g2D_all_filtered.fa | 7 |
| group_g2E/group_g2E_all_filtered.fa | 2 |
| group_g2F/group_g2F_all_filtered.fa | 6 |
| group_g2G/group_g2G_all_filtered.fa | 1 |
| group_g2H/group_g2H_all_filtered.fa | 3 |
| group_g2I/group_g2I_all_filtered.fa | 1 |
| group_g2J/group_g2J_all_filtered.fa | 2 |

### Phylogenetic analyses of groups g2A, g2B, g2D, and g2F

1. Align the set of diatom introns (queries) and non-redundant subject sequences
`muscle -in group_g2B_all_filtered.fa -out group_g2B_all_filtered.aligned.fa`

1. Use GBLOCKS to trim alignments
`gblocks.exe <ALIGNMENT> -t d -b5 h`

1. Use IQ-TREE to calculate the tree
`iqtree-omp.exe -s <TRIMMED_ALIGNMENT> -m MFP -bb 10000 -nt AUTO`
