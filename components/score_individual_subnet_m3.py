import time, os, json, concurrent.futures
import numpy as np
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import ProcessPoolExecutor
import concurrent.futures


class ScoreIndividualSubnet:
    def __init__(self, parentNetwork, loci):
        self.individualSubnet = []
        self.parentNetwork = parentNetwork  ###REFACTOR: to replace module1FaNetwork
        self.loci = loci
        self.canditdateGeneScores = {}

    # 8
    # Input: swapped subnets (a version of the individual subnet per gene in the locus {swapped gene from locus: subnet}), emptyLocusScore
    # Output: candidateGeneScores list (containing a list of individualGeneScores dictionaries, returned from process_subnet_count_edges)
    def count_edges(self, subnets, emptyLocusScore, batch_size=60):
        candidateGeneScores = []

        # create thread pool for each of the subnets (passed from candidate_gene_score) and execute in batches (batch_size of 60 seemed to be the optimal bin size for runtime)
        with concurrent.futures.ThreadPoolExecutor() as executor:
            batched_subnets = []

            for i in range(0, len(subnets), batch_size):
                batched_subnets.append(subnets[i : i + batch_size])

            futures = {
                executor.submit(self.process_subnet_count_edges, batch, emptyLocusScore)
                for batch in batched_subnets
            }

        # as the threads complete store each set of individualCandidtateGeneScores dictionaries in candidateGeneScores list
        for future in concurrent.futures.as_completed(futures):
            try:
                individualCandidateGeneScore = future.result()
                candidateGeneScores.extend(individualCandidateGeneScore)
            except Exception as exc:
                print(f"Count Edges - Generated an exception: {exc} {emptyLocusScore}")
        return candidateGeneScores

    # 9
    # Input: batch of swapped subnetworks, emptyLocusScore
    # Output: batch of individualCandidateGeneScore dictionaries (each containing the gene that was swapped and the gene score for that swapped subnet )
    ## geneScore in individualCandidateGeneScore represents the contribution "gene" makes to the network
    def process_subnet_count_edges(self, subnet, emptyLocusScore):
        batchScores = []
        for item in subnet:
            gene, subnetGenes = list(item.items())[0]

            # create a mask for the conditional (if either gene in the parent network data frame row is in the subnetwork) -> returns true or false
            mask = self.parentNetwork["gene1"].isin(subnetGenes) & self.parentNetwork[
                "gene2"
            ].isin(subnetGenes)

            selectedRows = self.parentNetwork[mask].copy()

            # sort the mask by row
            selectedRows["sorted_genes"] = np.sort(
                selectedRows[["gene1", "gene2"]], axis=1
            ).tolist()

            selectedRows.drop_duplicates(subset=["sorted_genes"], inplace=True)

            # count the number of times the mask conditional evaluated to true
            weightSum = (selectedRows["weight"].astype(float)).sum()

            individualCandidateGeneScore = {
                "gene": gene,
                "geneScore": weightSum - emptyLocusScore,
            }
            batchScores.append(individualCandidateGeneScore)
        return batchScores

    # 3
    # Input: geneLocus (locus returned from find_gene_locus), subnet
    # Output: edge count of empty locus case subnetwork
    def empty_locus_case(self, geneLocus, subnet):
        subnet = subnet.copy()
        subnet.remove(geneLocus)

        edgeCount = self.process_empty_locus_case(geneLocus, subnet)

        return edgeCount

    # 4
    # Input: gene in individual subnetwork, individual subnetwork
    # Output: empty locus case edge count -> process_gene_gene_score()
    def process_empty_locus_case(self, gene, subnet):
        # create a mask with the same logic as count_edges() -> process_subnet_count_edges()
        mask = self.parentNetwork["gene1"].isin(subnet) & self.parentNetwork[
            "gene2"
        ].isin(subnet)
        selectedRows = self.parentNetwork[mask].copy()

        selectedRows["sorted_genes"] = np.sort(
            selectedRows[["gene1", "gene2"]], axis=1
        ).tolist()
        selectedRows.drop_duplicates(subset=["sorted_genes"], inplace=True)

        weightSum = (selectedRows["weight"].astype(float)).sum()
        return weightSum

    # 5
    # Input: gene to find locus for, FA loci object
    # Output: locus of gene, locus -> process_gene_gene_score()
    def find_gene_locus(self, gene):
        for locus in self.loci:
            if gene in self.loci[locus]:
                return self.loci[locus], locus

    # 6
    # Input: locus for gene to swap (conatins all genes from locus), original gene to swap, individual subnetwork, empty locus score
    # Output: Dictionary of candidate gene scores: {{gene:gene, geneScore:geneScore}} -> process_gene_gene_score()
    def candidate_gene_score(self, locus, gene, subnet, emptyLocusScore):
        swappedSubnets = []
        # get the index of the gene in the subnet list
        geneIndexInSubnet = subnet.index(gene)

        # create a thread pool for each of the swapped subnets and execute process_locus_gene_candidate_gene_score for each of the swapped subnets
        with concurrent.futures.ThreadPoolExecutor() as executor:
            futures = {
                executor.submit(
                    self.process_locus_gene_candidate_gene_score,
                    locusGene,
                    geneIndexInSubnet,
                    subnet,
                )
                for locusGene in locus
            }

        # as the threads complete, add to swappedSubnets list for use in count_edges
        for future in concurrent.futures.as_completed(futures):
            try:
                result = future.result()
                swappedSubnets.append(result)
            except Exception as exc:
                print(f"Candidate Gene Score - Generated an exception: {exc}")

        candidateGeneScores = self.count_edges(swappedSubnets, emptyLocusScore)
        return candidateGeneScores

    # 7
    # Input: gene to swap, gene index in individual subnetwork, individual subnetwork
    # Output: dictionary that contains the locus specific gene to swap: swapped subnetwork -> candidate_gene_score()
    def process_locus_gene_candidate_gene_score(
        self, locusGene, geneIndexInSubnet, subnet
    ):
        subnet[geneIndexInSubnet] = locusGene
        tempSubnet = subnet.copy()
        return {locusGene: tempSubnet}

    # 2
    # Input: Gene from subnet and subnet
    # Output: Candidate Gene Scores object that contains:
    def process_gene_gene_score(self, gene, subnet):
        #####NOTE: each of the function calls are timed. IF DESIRED, uncomment the print statements after each of the end times to observe the total runtime per function call

        # get empty locus score for each gene
        estart = time.time()
        emptyLocusScore = self.empty_locus_case(gene, subnet)
        eend = time.time()
        # print(f"Empty Locus Time: {eend - estart}")

        # find the locus and return list of locus genes
        lstart = time.time()
        locus, locusNumber = self.find_gene_locus(gene)
        lend = time.time()
        # print(f"Find Locus Time: {lend-lstart}")

        # get candidate gene score for each gene in the subnet
        cstart = time.time()
        candidateGeneScores = self.candidate_gene_score(
            locus, gene, subnet, emptyLocusScore
        )
        cend = time.time()
        # print(f"Candidate Gene Score Time: {cend-cstart}")

        return locusNumber, candidateGeneScores

    # 1
    # Input: Individual FA subnetwork
    # Output: Average Gene Scores object -> main.py
    def gene_score(self, individualSubnet):
        self.individualSubnet = individualSubnet
        print(f"Subnet Scoring Initialized for subnet: {self.individualSubnet}")
        start = time.time()
        subnet = individualSubnet
        geneScores = {}

        # for each gene in the subnet, call process_gene which handles the intialization logic of scoring the genes within the genes' locus
        with ProcessPoolExecutor() as executor:
            futures = {
                executor.submit(self.process_gene_gene_score, gene, subnet): gene
                for gene in subnet
            }
            # gather results of candidate_gene_scores, as they complete and store in geneScores
            for future in concurrent.futures.as_completed(futures):
                locusNumber, candidateGeneScores = future.result()
                geneScores[locusNumber] = candidateGeneScores

        end = time.time()
        print(f"gene score time: {end-start}")
        currentDir = os.path.dirname(os.path.abspath(__file__))
        relativePath = os.path.join(
            currentDir, "..", "created_at_runtime", "scoring_geneScores.txt"
        )
        with open(relativePath, "a") as file:
            for key, value in geneScores.items():
                file.write(str(key) + ": " + str(value) + "\n")
        return geneScores
