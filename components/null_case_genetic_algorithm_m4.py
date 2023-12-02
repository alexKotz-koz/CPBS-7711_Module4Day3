import random, os, time, json
from concurrent.futures import ThreadPoolExecutor
import pandas as pd
import numpy as np
from components.fa_utilities import Fa_Utilities


class Null_Case_Genetic_Algorithm:
    def __init__(self, initialPopulation, bins, faLoci, masterParentNetworkDF):
        self.initialPopulation = initialPopulation
        self.bins = bins

        self.parentNetwork = masterParentNetworkDF

        self.faUtilitiesInstance = Fa_Utilities(parentNetworkFile=self.parentNetwork)
        self.invertedBins = self.faUtilitiesInstance.genes_to_bins(bins=bins)
        self.generations = {}

    def start_genetic_algorithm(self):
        print("Starting Genetic Algorithm Optimization")
        initialPopulation = self.initialPopulation
        onestart = time.time()
        generationXSubnets = self.mutate(initialPopulation=initialPopulation)
        initialPopulationAverageDensity = self.calculate_average_density(
            generationXSubnets
        )
        print(f"init avg den: {initialPopulationAverageDensity}")
        currentDir = os.path.dirname(os.path.abspath(__file__))
        relativePath = os.path.join(currentDir, "..", "created_at_runtime", "TEST.json")

        with open(relativePath, "w") as outputFile:
            json.dump(generationXSubnets, outputFile)
        oneend = time.time()
        print(f"1 gen completed in: {oneend-onestart}")
        # print(f"initialPopulationAverageDensity:{initialPopulationAverageDensity}")

        # TEST MATING
        averageDensity = 0.1235346322
        secondGeneration = self.mating(generationXSubnets)

        """
        averageDensity, secondGeneration = self.mating(generationXSubnets)
        """

        self.generations[1] = {
            "averageDensity": averageDensity,
            "subnets": secondGeneration,
        }
        finalPopulation = self.run_optimization(secondGeneration, averageDensity)

        currentDir = os.path.dirname(os.path.abspath(__file__))
        relativePath = os.path.join(
            currentDir, "..", "created_at_runtime", "generations.json"
        )

        with open(relativePath, "w") as outputFile:
            json.dump(self.generations, outputFile)

        currentDir = os.path.dirname(os.path.abspath(__file__))
        relativePath = os.path.join(
            currentDir, "..", "created_at_runtime", "finalNfaPopulation.json"
        )

        with open(relativePath, "w") as outputFile:
            json.dump(finalPopulation, outputFile)
        return finalPopulation
        # MATING
        # Ref: s*i =
        # Ref: Pairs of subnetworks were sampled (with replacement), where the probability of selecting a parent subnetwork i was equal to s*i

    # generationX = any generation of subnetworks after the initial generation (population)
    def mutate(self, initialPopulation=None, generationX=None):
        swappedSubnets = []
        mutationStepSubnetworks = []

        if initialPopulation != None and generationX == None:
            print(f"Mutating Initial Population")
            for subnet in initialPopulation["subnetworks"]:
                swappedSubnets = self.mutate_create_swapped_subnets(subnet)
                mutationStepSubnetworks.append(swappedSubnets)

        elif initialPopulation == None and generationX != None:
            print(f"Mutating next generation population")
            for subnet in generationX:
                swappedSubnets = self.mutate_create_swapped_subnets(subnet)
                mutationStepSubnetworks.append(swappedSubnets)

        """currentDir = os.path.dirname(os.path.abspath(__file__))
        relativePath = os.path.join(currentDir, "..", "created_at_runtime", "temp.txt")

        with open(relativePath, "w") as outputFile:
            json.dump(mutationStepSubnetworks, outputFile)"""

        return mutationStepSubnetworks

    def mutate_create_swapped_subnets(self, subnet):
        # print(f"subnet: {subnet}")
        swappedSubnet = []
        mutationProbability = 0

        for gene in subnet:
            # print(gene)
            mutationProbability = random.randint(1, 100)
            if gene == None:
                print(f"gene is none before mutation, subnet is: {subnet}")

            if mutationProbability <= 5:
                if gene == "NA" or gene == "FAGENEROW":
                    swappedSubnet.append(gene)
                    continue

                genesBinName = self.invertedBins[gene]
                genesBin = self.bins[genesBinName]

                # if the new gene does not equal the old gene, swap the original gene with a random gene from the locus ... achieved via random_gene_from_locus()
                # Ref: "replacement gene was chosen form the remaining available genes in that locus uniformly at random"
                newGene = self.mutate_random_gene_from_bin(gene, genesBin)

                if newGene == None:
                    print(f"gene is none after mutation: {newGene}")

                swappedSubnet.append(newGene)

            else:
                if gene == None:
                    print(f"GENE is none: {gene}")
                swappedSubnet.append(gene)

        """print(f"og subnet:{subnet}")
        print(f"swapped subnet: {swappedSubnet}")"""
        return swappedSubnet

    def mutate_random_gene_from_bin(self, gene, bin):
        randomGeneIndexFromBin = random.randint(0, len(bin) - 1)
        # print(f"Start: rand -> {locus[randomGeneIndexFromLocus]} | gene: {gene}")
        if bin[randomGeneIndexFromBin] != gene:
            newGene = bin[randomGeneIndexFromBin]
            return newGene
        else:
            print(f"here")
            return self.mutate_random_gene_from_bin(gene, bin)

    def mating(self, generationXSubnets):
        sumOfSelectionScoresList = []
        sumOfSelectionScores = 0
        generationXSelectionScores = {}
        generationXNormalizedDensityScores = {}
        newGeneration = []

        # Create Selection Scores
        for index, subnet in enumerate(generationXSubnets):
            selectionScore = self.mating_calculate_selection_score(subnet)
            generationXSelectionScores[index] = {
                "selectionScore": selectionScore,
                "subnet": subnet,
            }
        sumOfSelectionScores = sum(
            subnet[1]["selectionScore"] for subnet in generationXSelectionScores.items()
        )

        # Calculate Probability Scores
        for index, subnet in enumerate(generationXSelectionScores.items()):
            subnetGenes = subnet[1]["subnet"]
            subnetSelectionScore = subnet[1]["selectionScore"]

            subnetProbabilityScore = (
                self.mating_calculate_normalized_subnet_probability_score(
                    sumOfSelectionScores, subnetSelectionScore
                )
            )
            generationXNormalizedDensityScores[index] = {
                "subnetProbabiltyScore": subnetProbabilityScore,
                "subnetSelectionScore": subnetSelectionScore,
                "subnet": subnetGenes,
            }

        # Create weights and generate 5000 mutated subnets
        generationXNormalizedDensityScoresList = list(
            generationXNormalizedDensityScores.values()
        )
        weights = [
            float(item["subnetProbabiltyScore"])
            for item in generationXNormalizedDensityScoresList
        ]
        i = 0
        while i < 5000:
            child = []
            parentSubnets = random.choices(
                generationXNormalizedDensityScoresList, weights=weights, k=2
            )
            child = self.mating_create_child_subnet(parentSubnets)
            faUtilitiesInstance = Fa_Utilities()
            edgeCount = faUtilitiesInstance.count_edges(child, self.parentNetwork)
            newGeneration.append(child)
            i += 1
        averageDensity = self.calculate_average_density(newGeneration)
        return averageDensity, newGeneration

    def mating_create_child_subnet(self, parentSubnets):
        child = []
        parent1 = parentSubnets[0]["subnet"]
        parent2 = parentSubnets[1]["subnet"]
        index = 0
        while index < 12:
            whichParent = random.randint(0, 1)

            if whichParent == 0:
                child.append(parent1[index])
            elif whichParent == 1:
                child.append(parent2[index])
            index += 1

        return child

        # Input:

    def mating_calculate_selection_score(self, subnet):
        edgeCount = self.faUtilitiesInstance.count_edges(
            subnetGenes=subnet, parentNetwork=self.parentNetwork
        )
        # return edgeCount**3
        return edgeCount

    def mating_calculate_normalized_subnet_probability_score(
        self, sumOfSelectionScores, subnetSelectionScore
    ):
        # print(f"generationX: {generationXSelectionScores}")

        # QUESTION: disregard subnets that have a selection score of 0
        return subnetSelectionScore / sumOfSelectionScores

    def count_edges_wrapper(self, subnet):
        return float(self.faUtilitiesInstance.count_edges(subnet, self.parentNetwork))

    def calculate_average_density(self, subnets):
        print("Calculating average density of mutated random non fa subnetworks")
        weights = []

        with ThreadPoolExecutor() as executor:
            weights = list(executor.map(self.count_edges_wrapper, subnets))

        return sum(weights) / len(subnets)

    def optimize(self, generation2Subnets, generation2AverageDensity):
        print("Starting Optimization")
        print(f"second gen average: {generation2AverageDensity}")

        # Mutate the second generation subnets
        nextGenerationMutated = self.mutate(generationX=generation2Subnets)

        # Perform mating on the mutated generation and get the average density
        nextGenAverageDensity, nextGenerationMated = self.mating(nextGenerationMutated)
        print(f"Next Generation: {nextGenAverageDensity}")

        # Get the last generation's average density
        lastIndex, lastValues = list(self.generations.items())[-1]
        previousAverageDensity = lastValues["averageDensity"]

        # Prepare the next generation's index
        nextIndex = int(lastIndex) + 1

        # Store the next generation's average density and subnets
        self.generations[nextIndex] = {
            "averageDensity": nextGenAverageDensity,
            "subnets": nextGenerationMated,
        }
        print("\n")
        print(f"lastAvg: {previousAverageDensity}\n")
        print(f"nextAvg: {nextGenAverageDensity}\n")

        # Calculate the improvement in density
        densityImprovement = (
            (nextGenAverageDensity - previousAverageDensity) / previousAverageDensity
        ) * 100

        print(f"Percentage: {densityImprovement}")
        # Return the next generation, its average density, and the improvement
        return nextGenerationMated, nextGenAverageDensity, densityImprovement

    def run_optimization(self, initialSubnets, initialAverageDensity):
        densityImprovement = 100  # Start with a large value

        # Continue optimizing until the improvement is less than 0.5
        while densityImprovement > 0.5:
            (
                nextGenerationMated,
                nextGenAverageDensity,
                densityImprovement,
            ) = self.optimize(initialSubnets, initialAverageDensity)

            # Update the initial subnets and average density for the next iteration
            initialSubnets, initialAverageDensity = (
                nextGenerationMated,
                nextGenAverageDensity,
            )
        finalPopulation = {
            "finalAverageDensity": nextGenAverageDensity,
            "finalSubnets": nextGenerationMated,
        }

        return finalPopulation
