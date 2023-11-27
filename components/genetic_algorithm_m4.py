import random
import os
import json
from components.fa_utilities import Fa_Utilities


class Genetic_Algorithm:
    def __init__(self, initialPopulation, faLoci, parentNetwork):
        self.initialPopulation = initialPopulation
        self.faLoci = faLoci
        self.parentNetwork = parentNetwork
        self.faUtilitiesInstance = Fa_Utilities(
            loci=faLoci, parentNetworkFile=parentNetwork
        )

    # enrichment():
    # 1) given the initial population, subject subnets to mutation step at 5% probability, add all subnets to stagedSubnetworks[].
    # 2) perform mating step on stagedSubnetworks[]
    # 3) repeat steps 1 and 2 until the generations average density fails to improve by 0.5%
    ## a) ADD GENERATION STATISTICS TO FILE
    ### i) Generation statistics: - average density, - number of subnetworks mutated, -

    def start_genetic_algorithm(self):
        initialPopulation = self.initialPopulation
        sumOfSelectionScoresList = []
        sumOfSelectionScores = 0
        generationXSelectionScores = {}
        generationXNormalizedDensityScores = {}
        # MUTATION
        generationXSubnets = self.mutate(initialPopulation=initialPopulation)

        # MATING
        # From Nourah: use a cubically transformed density score to determine the likelihood of a subnetwork being selected for mating
        # Ref: s*i =
        # Ref: Pairs of subnetworks were sampled (with replacement), where the probability of selecting a parent subnetwork i was equal to s*i
        for index, subnet in enumerate(generationXSubnets):
            selectionScore = self.calculate_selection_score(subnet)
            generationXSelectionScores[index] = {
                "selectionScore": selectionScore,
                "subnet": subnet,
            }
        for subnet in generationXSelectionScores.items():
            sumOfSelectionScoresList.append(subnet[1]["selectionScore"])
        sumOfSelectionScores = sum(sumOfSelectionScoresList)

        for index, subnet in enumerate(generationXSelectionScores.items()):
            subnetGenes = subnet[1]["subnet"]
            subnetSelectionScore = subnet[1]["selectionScore"]

            subnetProbabilityScore = self.calculate_normalized_subnet_probability_score(
                sumOfSelectionScores, subnetSelectionScore
            )
            generationXNormalizedDensityScores[index] = {
                "subnetProbabiltyScore": subnetProbabilityScore,
                "subnetSelectionScore": subnetSelectionScore,
                "subnet": subnetGenes,
            }
            # print(generationXNormalizedDensityScores[index])

        # print(f"probabilities: {generationXNormalizedDensityScores.items()}")

    def mutate(self, initialPopulation=None, generationX=None):
        swappedSubnets = []
        mutationStepSubnetworks = []

        if initialPopulation != None and generationX == None:
            for subnet in initialPopulation.items():
                swappedSubnets = self.mutate_create_swapped_subnets(subnet[1])
                mutationStepSubnetworks.append(swappedSubnets)

        currentDir = os.path.dirname(os.path.abspath(__file__))
        relativePath = os.path.join(currentDir, "..", "created_at_runtime", "temp.txt")

        with open(relativePath, "w") as outputFile:
            json.dump(mutationStepSubnetworks, outputFile)

        return mutationStepSubnetworks

    def mutate_create_swapped_subnets(self, subnetObj):
        subnet = subnetObj["subnet"]
        swappedSubnet = []
        mutationProbability = 0
        # print(f"og subnet: {subnet}")

        for gene in subnet:
            # Ref: "Each locus was mutated with a 5% probability"
            mutationProbability = random.randint(1, 100)

            if mutationProbability <= 5:
                # find the locus for the gene
                locus = self.faUtilitiesInstance.find_gene_locus(gene)[0]

                # find the index of a random gene from that locus
                randomGeneIndexFromLocus = random.randint(0, len(locus) - 1)

                # if the new gene does not equal the old gene, swap the original gene with a random gene from the locus ... achieved via random_gene_from_locus()
                # Ref: "replacement gene was chosen form the remaining available genes in that locus uniformly at random"
                newGene = self.random_gene_from_locus(gene, locus)
                if newGene == gene:
                    print(f"in: {gene} | {newGene}")
                swappedSubnet.append(newGene)

            else:
                swappedSubnet.append(gene)

        # print(f"new subnet: {swappedSubnet}")
        return swappedSubnet

    def random_gene_from_locus(self, gene, locus):
        randomGeneIndexFromLocus = random.randint(0, len(locus) - 1)

        if locus[randomGeneIndexFromLocus] != gene:
            return locus[randomGeneIndexFromLocus]
        else:
            self.random_gene_from_locus(gene, locus)

    # Input:
    def calculate_selection_score(self, subnet):
        # si = raw score for subnet (edge count ^ 3)

        edgeCount = self.faUtilitiesInstance.count_edges(
            subnetGenes=subnet, parentNetwork=self.parentNetwork
        )
        """if edgeCount == 0:
            print(f"edgeCount: {edgeCount}")"""

        return edgeCount**3

    def calculate_normalized_subnet_probability_score(
        self, sumOfSelectionScores, subnetSelectionScore
    ):
        # print(f"generationX: {generationXSelectionScores}")

        # QUESTION: disregard subnets that have a selection score of 0
        if subnetSelectionScore == 0:
            return "NA"
        else:
            return subnetSelectionScore / sumOfSelectionScores
