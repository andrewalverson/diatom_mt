## Compiling BLASTN hits and alignments for each _rnl_ intron

1. BLAST each intron individually to GenBank
`for i in *.fa;do test -e ${i%.fa}.blastn && continue; blastn -word_size 11 -reward 2 -penalty -3 -gapopen 5 -gapextend 2 -query $i -remote -db nr -outfmt '6 std qcovs' > ${i%.fa}.blastn;done &`

1. Retrieve BLAST hits with e-value <= 1e-12 and qcovs > 50
`for i in *.blastn;do blast_table_2_subject_fasta.py $i -q 50 -e 1e-12 > ${i%.blastn}.subjects.fa;done`

1. Manually move files into group-specific intron folders (e.g., group_g2A, group_g2B, etc.)

1. Concatenate subject sequences for each intron cluster
`for f in group_g*;do cat $f/*.subjects.fa > $f/${f}_subjects.fa;done`

1. Concatenate query (i.e, diatom) sequences for each intron cluster (run within the group directory)
`for i in *.fa; do if [[ $i =~ subjects || $i =~ group ]];then continue;else cat $i;fi;done > group_g2E_queries.fa`

1. Run CD-HIT to remove redundant subject sequences (e.g., some queries hit to both the genome and a separate cox1 gene accession AND each sequence from an intron cluster was BLASTed individually)
`cd-hit-est -i group_g2E_sbjcts.fa -o group_g2E_subjects_filtered.fa -c 0.99 -d 20000`

1. Concatenate query and subject sequences
`cat group_g2E_queries.fa group_g2E_subjects_filtered.fa > group_g2E_all.fa`

1. Remove redundant sequences from the full set query and subject sequences
`cd-hit-est -i group_g2E_all.fa -o group_g2E_all_filtered.fa -c 1 -d 20000`

### Final sequence counts for each _rnl_ intron group
`for i in group_g*/*_all_filtered.fa;do echo '| ' $i '|' | tr -d '\n' | perl -pwe 's/\s{2,}/ /';grep '>' $i | wc -l | perl -pwe 's/\s{2,}/ /' | tr -d '\n';echo " |";done`

|          Group          | Num sequences |
|-------------------------------------|---|
| group_g1/group_g1_all_filtered.fa   | 1 |
| group_g2A/group_g2A_all_filtered.fa | 6 |
| group_g2B/group_g2B_all_filtered.fa | 4 |
| group_g2C/group_g2C_all_filtered.fa | 1 |
| group_g2D/group_g2D_all_filtered.fa | 1 |
| group_g2E/group_g2E_all_filtered.fa | 1 |

>Note: for group\_g2A, `Navicula_ramosissima_intron_rnl_g2_2.fa` had four hsp's to MF997424.1 (_Halamphora_ mt genome). The script `blast_table_2_subject_fasta.py` combined these into a single subject sequence, but they were clearly just 3 different matches. We extracted the region corresponding to the first two hsp's (66739-68234). The other two hits individually had qcovs < 50, so they were discarded. So, we removed the large MF997424.1 sequence from `group_g2A_all_filtered.fa` and replaced it with `MF997424.1_66739-68234.fa`.

### Phylogenetic analysis of group g2A

1. Align the set of diatom introns (queries) and non-redundant subject sequences
`muscle -in group_g2A_all_filtered.fa -out group_g2A_all_filtered.aligned.fa`

1. Use GBLOCKS to trim the alignment

1. Use IQ-TREE to calculate the tree

