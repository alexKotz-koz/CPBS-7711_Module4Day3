import random
import os
import time
import json
import concurrent
from concurrent.futures import ThreadPoolExecutor
from components.fa_utilities import Fa_Utilities


class Genetic_Algorithm:
    def __init__(self, initialPopulation, faLoci, parentNetwork):
        self.initialPopulation = initialPopulation
        self.faLoci = faLoci
        self.parentNetwork = parentNetwork
        self.faUtilitiesInstance = Fa_Utilities(
            loci=faLoci, parentNetworkFile=parentNetwork
        )
        self.generations = {}

    # Start for genetic algorithm
    # Input: initial fa population
    # Output: final optimized fa population
    def start_genetic_algorithm(self):
        print("Starting Genetic Algorithm Optimization")
        initialPopulation = self.initialPopulation

        #INITIAL POPULATION OPTIMIZATION ROUTINE
        generationXSubnets = self.mutate(initialPopulation=initialPopulation)
        initialPopulationAverageDensity = self.calculate_average_density(
            generationXSubnets
        )
        averageDensity, secondGeneration = self.mating(generationXSubnets)
        self.generations[1] = {
            "averageDensity": averageDensity,
            "subnets": secondGeneration,
        }

        # call optimization function to run genetic algorithm on first optimized generation
        finalPopulation = self.run_optimization(secondGeneration, averageDensity)

        currentDir = os.path.dirname(os.path.abspath(__file__))
        relativePath = os.path.join(
            currentDir, "..", "created_at_runtime", "generations.json"
        )

        with open(relativePath, "w") as outputFile:
            json.dump(self.generations, outputFile)

        currentDir = os.path.dirname(os.path.abspath(__file__))
        relativePath = os.path.join(
            currentDir, "..", "created_at_runtime", "finalFaPopulation.json"
        )

        with open(relativePath, "w") as outputFile:
            json.dump(finalPopulation, outputFile)
        return finalPopulation
        
    # Input: either initialPopulation or generationX (subsequent population to optimize)
    # Ouput: Mutated subnetwork population
    def mutate(self, initialPopulation=None, generationX=None):
        swappedSubnets = []
        mutationStepSubnetworks = []

        if initialPopulation != None and generationX == None:
            print(f"Mutating Initial Population")
            for subnet in initialPopulation.items():
                swappedSubnets = self.mutate_create_swapped_subnets(subnet[1])
                mutationStepSubnetworks.append(swappedSubnets)
        elif initialPopulation == None and generationX != None:
            print(f"Mutating next generation population")
            for subnet in generationX:
                swappedSubnets = self.mutate_create_swapped_subnets(subnet)
                mutationStepSubnetworks.append(swappedSubnets)

        return mutationStepSubnetworks

    # Input: subnetwork to mutate
    # Output: mutated subnetwork
    def mutate_create_swapped_subnets(self, subnetObj):
        if isinstance(subnetObj, dict):
            subnet = subnetObj["subnet"]
        elif isinstance(subnetObj, list):
            subnet = subnetObj

        swappedSubnet = []
        mutationProbability = 0

        for gene in subnet:

            # Paper reference: "Each locus was mutated with a 5% probability"
            mutationProbability = random.randint(1, 100)
            if mutationProbability <= 5:
                locus = self.faUtilitiesInstance.find_gene_locus(gene)[0]

                # Paper reference for function: "replacement gene was chosen form the remaining available genes in that locus uniformly at random"
                newGene = self.mutate_random_gene_from_locus(gene, locus)

                swappedSubnet.append(newGene)

            else:
                swappedSubnet.append(gene)

        return swappedSubnet

    # Input: gene and gene's locus
    # Ouput: new gene to add to mutated subnetwork
    def mutate_random_gene_from_locus(self, gene, locus):
        randomGeneIndexFromLocus = random.randint(0, len(locus) - 1)
        if locus[randomGeneIndexFromLocus] != gene:
            newGene = locus[randomGeneIndexFromLocus]
            return newGene
        else:
            return self.mutate_random_gene_from_locus(gene, locus)

    # Input: mutated subnetwork population
    # Output: average density of mated subnetwork population and mated subnetwork population
    def mating(self, generationXSubnets):
        sumOfSelectionScoresList = []
        sumOfSelectionScores = 0
        generationXSelectionScores = {}
        generationXNormalizedDensityScores = {}
        newGeneration = []

        # create selection scores
        for index, subnet in enumerate(generationXSubnets):
            selectionScore = self.mating_calculate_selection_score(subnet)
            generationXSelectionScores[index] = {
                "selectionScore": selectionScore,
                "subnet": subnet,
            }
        sumOfSelectionScores = sum(
            subnet[1]["selectionScore"] for subnet in generationXSelectionScores.items()
        )

        # create probability scores
        for index, subnet in enumerate(generationXSelectionScores.items()):
            subnetGenes = subnet[1]["subnet"]
            subnetSelectionScore = subnet[1]["selectionScore"]

            subnetProbabilityScore = subnetSelectionScore / sumOfSelectionScores
            generationXNormalizedDensityScores[index] = {
                "subnetProbabiltyScore": subnetProbabilityScore,
                "subnetSelectionScore": subnetSelectionScore,
                "subnet": subnetGenes,
            }

        # create weights and generate 5000 mutated subnets
        generationXNormalizedDensityScoresList = list(
            generationXNormalizedDensityScores.values()
        )
        weights = [
            float(item["subnetProbabiltyScore"])
            for item in generationXNormalizedDensityScoresList
        ]
        # create 5000 mated subnetworks
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

    # Input: two parent subnetworks to mate
    # Output: mated child subnetwork
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

    # Input: subnetwork
    # Output: selection score of subnetwork
    def mating_calculate_selection_score(self, subnet):
        weight = self.faUtilitiesInstance.count_edges(
            subnetGenes=subnet, parentNetwork=self.parentNetwork
        )

        return weight

    # wrapper function for thread pool
    def count_edges_wrapper(self, subnet):
        return float(self.faUtilitiesInstance.count_edges(subnet, self.parentNetwork))

    # Input: subnetwork population
    # Ouput: average density of population
    def calculate_average_density(self, subnets):
        print("Calculating average density of mutated random fa subnetworks")
        weights = []

        with ThreadPoolExecutor() as executor:
            weights = list(executor.map(self.count_edges_wrapper, subnets))

        return sum(weights) / len(subnets)

    # Input: generation x subnetwork and average density
    # Output: optimized subnetwork, average density, and percent difference between input population and output population
    def optimize(self, generation2Subnets, generation2AverageDensity):
        print("Starting Optimization")
        print(f"second gen average: {generation2AverageDensity}")
        start = time.time()
        # mutate the second generation subnets
        nextGenerationMutated = self.mutate(generationX=generation2Subnets)

        # perform mating on the mutated generation and get the average density
        nextGenAverageDensity, nextGenerationMated = self.mating(nextGenerationMutated)
        print(f"Next Generation: {nextGenAverageDensity}")

        # get the last generation's average density
        lastIndex, lastValues = list(self.generations.items())[-1]
        previousAverageDensity = lastValues["averageDensity"]

        # prepare the next generation's index
        nextIndex = int(lastIndex) + 1

        # store the next generation's average density and subnets
        print("\n")
        print(f"lastAvg: {previousAverageDensity}\n")
        print(f"nextAvg: {nextGenAverageDensity}\n")

        # calculate the improvement in density
        densityImprovement = (
            (nextGenAverageDensity - previousAverageDensity) / previousAverageDensity
        ) * 100
        end = time.time()
        self.generations[nextIndex] = {
            "averageDensity": nextGenAverageDensity,
            "timeToCompletion": end - start,
            "percentImproved": densityImprovement,
            "subnets": nextGenerationMated,
        }
        print(f"Percentage: {densityImprovement}")

        # return the next generation, its average density, and the improvement
        return nextGenerationMated, nextGenAverageDensity, densityImprovement

    # Input: initial subnetwork population
    # Output: final optimized subnetwork population
    def run_optimization(self, initialSubnets, initialAverageDensity):
        densityImprovement = 100  # Start with a large value

        # continue optimizing until the improvement is less than 0.5
        while densityImprovement > 0.5:
            (
                nextGenerationMated,
                nextGenAverageDensity,
                densityImprovement,
            ) = self.optimize(initialSubnets, initialAverageDensity)

            # update the initial subnets and average density for the next iteration
            initialSubnets, initialAverageDensity = (
                nextGenerationMated,
                nextGenAverageDensity,
            )
        subnetFinal = []
        for subnet in nextGenerationMated:
            weight = self.faUtilitiesInstance.count_edges(subnet, self.parentNetwork)
            subnetDict = {"averageDensity": weight, "subnet": subnet}
            subnetFinal.append(subnetDict)

        finalPopulation = {
            "finalAverageDensity": nextGenAverageDensity,
            "finalSubnets": subnetFinal,
        }

        return finalPopulation
