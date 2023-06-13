####Before performing the following operations, you need to organize the motif's pwm matrix into a format that can be recognized by cluster-buster. 
####Here only the processing process of pwm matrix from cis-bp database is shown(cis_bp_pwm_matrix_merge.R)
####Only the key steps to getting a network are shown here
#### "selected_region_of_genome" means genome spanning 2 kb upstream, intron and 1 kb downstream regions, with the exons masked



./cbust -c 0 -f 5 -G 0 JASPAR.txt selected_region_of_genome.fa > motif_c_0_JASPAR.bed
./cbust -m 5 -f 5 -G 0 JASPAR.txt selected_region_of_genome.fa > motif_m_5_JASPAR.bed
./cbust -c 0 -f 5 -G 0 cisbp.txt selected_region_of_genome.fa > motif_c_5_cisbp.bed
./cbust -m 5 -f 5 -G 0 cisbp.txt selected_region_of_genome.fa > motif_m_5_cisbp.bed
###Extract motif location files from cluster-buster results(the example file merged.txt is as following)
###chrom	start	stop	TF_name	score	strand
#4	3886484	3886497	M01315_2.00	.	+
#4	3887162	3887175	M01315_2.00	.	+
#3	10432527	10432540	M01315_2.00	.	+
#10	23871352	23871365	M01315_2.00	.	+
#7	45722339	45722352	M01315_2.00	.	+
#10	54385569	54385582	M01315_2.00	.	-
#9	72742862	72742875	M01315_2.00	.	+
#1	298881545	298881559	M01162_2.00	.	-
#8	97474596	97474610	M01162_2.00	.	+


##Extract the gene location file from genome file and convert it to bed format 
awk '{if($3~/^gene$/)print}' file.gff > genes.gff
gff2bed <genes.gff> genes.bed
###The intersection of cluster-buster results and gene location files was used to obtain the correspondence between motif and target
bedtools intersect -nonamecheck -a genes.bed -b merged.txt -wa -wb -s | bedtools groupby -i - -g 1-15 -c 11 -o collapse > gene.tsv
###we can find the interactions between motif id and target in "gene.tsv" by choosing the column of it.
####we can get the interactions between motif_id and TF name in JASPAR and the file named "TF_Information_all_motifs.txt" from cis-bp database.
###Then the correspondence between TF and target is found through motif and TF, that is, the DNA motif network is captured
