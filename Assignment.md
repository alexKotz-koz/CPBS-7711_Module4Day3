# Module 4 Day 3 Assignment:

Implement the genetic algorithm from Tasan et al. and integrate with previous submissions to reimplement the Prix Fixe algorithm in the paper, in particular referencing the sections ‘Prix fixe subnetwork enrichment’ and ‘Prix fixe gene-scoring’ in the Online Methods section. However, the original algorithm assumes unweighted edges (i.e. binary indication of an edge) whereas your implementation will allow weighted edges. You must create and justify your method for incorporating weighted network edges in the subnetwork scoring metric.

As output, your program will provide a final relative ranking of each gene in each locus in a file titled ‘Day3_Output.gmt.’ The ranking file will be a modified format such that each gene name has appended to it the gene’s score (space between gene name and gene score, tab-delimited as in GMT between gene-score pairs). Use ‘NA’ as the score if the gene is not found in the network or otherwise not scored. Additional output files will specify the final top 10 scoring networks, each in a tab-delimited file in the format of ‘Day3_STRING.txt’, which itself is suitable for input into Cytoscape for visualization (http://manual.cytoscape.org/en/stable/Supported_Network_File_Formats.html).

The filename for each top-scoring network must follow the naming convention ‘Day3_Output_Network<DD>_pval<p-value>.txt’ for DD from 01 to 10 and the corresponding calculated p-value for the population. Example can be found in the image below. Also, you should provide a file or files containing statistics about the genetic algorithm (GA) search per generation.

1. Create 5000 Random FA Subnetworks

2. For each of the 5000 random fa subnetworks (initial population):

3. Mutation Step:

4. Mating Step:

5. Repeat 3 and 4 until average gene scores do not improve by more the 0.5%, this would result in the "final population"