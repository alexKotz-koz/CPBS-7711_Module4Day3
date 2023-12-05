import random, os, time, json, concurrent
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
import pandas as pd
import numpy as np
from components.fa_utilities import Fa_Utilities


class Null_Case_Genetic_Algorithm:
    def __init__(
        self, initialPopulation, bins, faLoci, masterParentNetworkDF, geneDict
    ):
        self.initialPopulation = initialPopulation
        self.bins = bins

        self.parentNetwork = masterParentNetworkDF

        self.faUtilitiesInstance = Fa_Utilities(parentNetworkFile=self.parentNetwork)
        self.invertedBins = self.faUtilitiesInstance.genes_to_bins(bins=bins)
        self.generations = {}
        self.geneDict = geneDict

    def start_genetic_algorithm(self):
        print("Starting Genetic Algorithm Optimization")
        initialPopulation = self.initialPopulation

        onestart = time.time()
        # mutate and mate initial population
        generationXSubnets = self.mutate(initialPopulation=initialPopulation)
        print(f"First mutated generation: {generationXSubnets}")
        densityStart = time.time()
        initialPopulationAverageDensity = self.calculate_average_density(
            generationXSubnets
        )
        densityEnd = time.time()
        print(f"Average Density calculated in {densityEnd-densityStart}")

        currentDir = os.path.dirname(os.path.abspath(__file__))
        relativePath = os.path.join(currentDir, "..", "created_at_runtime", "TEST.json")
        with open(relativePath, "w") as outputFile:
            json.dump(generationXSubnets, outputFile)
        oneend = time.time()
        print(f"1 gen completed in: {oneend-onestart}")

        generationXSubnets = {}
        with open("created_at_runtime/TEST.json", "r") as file:
            generationXSubnets = json.load(file)

        # TEST MATING
        """averageDensity = 0.1235346322
        secondGeneration = self.mating(generationXSubnets)"""

        averageDensity, secondGeneration = self.mating(generationXSubnets)

        self.generations[1] = {
            "averageDensity": averageDensity,
            "subnets": secondGeneration,
        }

        # start optimization routine
        finalPopulation = self.run_optimization(secondGeneration, averageDensity)

        currentDir = os.path.dirname(os.path.abspath(__file__))
        relativePath = os.path.join(
            currentDir, "..", "created_at_runtime", "null_generations.json"
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
        
    # Input: either the initial population or a secondary population (this function is called in the optimization function, this generationX population will be the next population to optimize)
    # Output: a population of 5000 mutated subnetworks
    def mutate(self, initialPopulation=None, generationX=None):
        swappedSubnets = []
        mutationStepSubnetworks = []

        # initial population case, the only difference is the data structure of the initial population vs generationX population
        if initialPopulation != None and generationX == None:
            print(f"Mutating Initial Population\n")
            for subnet in initialPopulation["subnetworks"]:
                swappedSubnets = self.mutate_create_swapped_subnets(subnet)
                mutationStepSubnetworks.append(swappedSubnets)

        # secondary population case (generationX population)
        elif initialPopulation == None and generationX != None:
            print(f"Mutating next generation population\n")
            for subnet in generationX:
                swappedSubnets = self.mutate_create_swapped_subnets(subnet)
                mutationStepSubnetworks.append(swappedSubnets)

        return mutationStepSubnetworks

    # Input: an individual subnetwork from the population of subnetworks
    # Output: a mutated subnetwork (swapped: using this term as an indicator of randomly replacing genes from the same bin)
    def mutate_create_swapped_subnets(self, subnet):
        swappedSubnet = []
        mutationProbability = 0

        for gene in subnet:
            # create a random number from 1 to 100 to use as a test (whether to swap this gene with another gene from the same bin or not)
            mutationProbability = random.randint(1, 100)

            # reference from paper: "Each locus was mutated with a 5% probablility" pg 161
            if mutationProbability <= 5:
                # check to see if the gene is one of the FA genes that only exist in the Input.gmt.txt file and not in the STRING 1.txt file
                if gene == "NA" or gene == "FAGENEROW":
                    swappedSubnet.append(gene)
                    continue

                # use invertedBins to find the bin number of the gene
                genesBinName = self.invertedBins[gene]
                # use the bin number from invertedBins search to find the bin (containing the genes)
                genesBin = self.bins[genesBinName]
                # call mutate_random_gene_from_bin to actually mutate the gene
                newGene = self.mutate_random_gene_from_bin(gene, genesBin)
                # append the mutated gene to the new subnetwork
                swappedSubnet.append(newGene)
            # if the random number was not equal to or less than 5, just add the original gene to the new subnetwork
            else:
                swappedSubnet.append(gene)
        return swappedSubnet

    def mutate_random_gene_from_bin(self, gene, bin):
        randomGeneIndexFromBin = random.randint(0, len(bin) - 1)
        # print(f"Start: rand -> {locus[randomGeneIndexFromLocus]} | gene: {gene}")
        if bin[randomGeneIndexFromBin] != gene:
            newGene = bin[randomGeneIndexFromBin]
            return newGene
        else:
            return self.mutate_random_gene_from_bin(gene, bin)

    def calculate_subnet_density_wrapper(self, subnet):
        start = time.time()

        result = float(
            self.faUtilitiesInstance.count_edges_null_case(
                subnet, self.parentNetwork, self.geneDict
            )
        )
        end = time.time()

        return result

    def mating_calculate_selection_score(self, generationXSubnets):
        generationXSelectionScores = {}
        with ThreadPoolExecutor() as executor:
            future_to_subnet = {
                executor.submit(self.calculate_subnet_density_wrapper, subnet): tuple(
                    subnet
                )
                for subnet in generationXSubnets
            }
            for future in concurrent.futures.as_completed(future_to_subnet):
                subnet = future_to_subnet[future]
                selectionScore = future.result()
                generationXSelectionScores[subnet] = {
                    "selectionScore": selectionScore,
                    "subnet": subnet,
                }

        sumOfSelectionScores = sum(
            subnet[1]["selectionScore"] for subnet in generationXSelectionScores.items()
        )

        return generationXSelectionScores, sumOfSelectionScores

    def mating(self, generationXSubnets):
        sumOfSelectionScores = 0
        generationXSelectionScores = {}
        generationXNormalizedDensityScores = {}
        newGeneration = []

        (
            generationXSelectionScores,
            sumOfSelectionScores,
        ) = self.mating_calculate_selection_score(generationXSubnets=generationXSubnets)

        for index, subnet in enumerate(generationXSelectionScores.items()):
            subnetGenes = list(subnet[1]["subnet"])
            subnetSelectionScore = subnet[1]["selectionScore"]

            subnetProbabilityScore = subnetSelectionScore / sumOfSelectionScores
            generationXNormalizedDensityScores[index] = {
                "subnetProbabiltyScore": subnetProbabilityScore,
                "subnetSelectionScore": subnetSelectionScore,
                "subnet": subnetGenes,
            }

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

    def calculate_average_density(self, subnets):
        print("Calculating average density of mutated random non fa subnetworks")
        weights = []

        with ThreadPoolExecutor() as executor:
            weights = list(executor.map(self.calculate_subnet_density_wrapper, subnets))

        return sum(weights) / len(subnets)

    def optimize(self, generation2Subnets, generation2AverageDensity):
        print("Starting Optimization")
        print(f"second gen average: {generation2AverageDensity}")

        nextGenerationMutated = self.mutate(generationX=generation2Subnets)

        nextGenAverageDensity, nextGenerationMated = self.mating(nextGenerationMutated)
        print(f"Next Generation: {nextGenAverageDensity}")

        lastIndex, lastValues = list(self.generations.items())[-1]
        previousAverageDensity = lastValues["averageDensity"]

        nextIndex = int(lastIndex) + 1

        self.generations[nextIndex] = {
            "averageDensity": nextGenAverageDensity,
            "subnets": nextGenerationMated,
        }
        print("\n")
        print(f"lastAvg: {previousAverageDensity}\n")
        print(f"nextAvg: {nextGenAverageDensity}\n")

        densityImprovement = (
            (nextGenAverageDensity - previousAverageDensity) / previousAverageDensity
        ) * 100

        print(f"Percentage: {densityImprovement}")
        return nextGenerationMated, nextGenAverageDensity, densityImprovement

    def run_optimization(self, initialSubnets, initialAverageDensity):
        densityImprovement = 100  

        while densityImprovement > 0.5:
            (
                nextGenerationMated,
                nextGenAverageDensity,
                densityImprovement,
            ) = self.optimize(initialSubnets, initialAverageDensity)

            initialSubnets, initialAverageDensity = (
                nextGenerationMated,
                nextGenAverageDensity,
            )
        finalPopulation = {
            "finalAverageDensity": nextGenAverageDensity,
            "finalSubnets": nextGenerationMated,
        }

        return finalPopulation
