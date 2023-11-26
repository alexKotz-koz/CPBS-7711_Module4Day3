import random
import os
import json
from components.fa_utilities import Fa_Utilities


class Optimization_Algorithm:
    def __init__(self, initialPopulation, faLoci):
        self.initialPopulation = initialPopulation
        self.faLoci = faLoci
        self.faUtilitiesInstance = Fa_Utilities(loci=faLoci)

    # enrichment():
    # 1) given the initial population, subject subnets to mutation step at 5% probability, add all subnets to stagedSubnetworks[].
    # 2) perform mating step on stagedSubnetworks[]
    # 3) repeat steps 1 and 2 until the generations average density fails to improve by 0.5%
    ## a) ADD GENERATION STATISTICS TO FILE
    ### i) Generation statistics: - average density, - number of subnetworks mutated, -

    def start_optimization_algorithm(self):
        initialPopulation = self.initialPopulation

        statNumberOfMutatedSubnets, generationX = self.mutate(
            initialPopulation=initialPopulation
        )

    def mutate(self, initialPopulation=None, generationX=None):
        statNumberOfMutatedSubnets = 0
        swappedSubnets = []
        mutationStepSubnetworks = []

        if initialPopulation != None and generationX == None:
            for subnet in initialPopulation.items():
                # "Each locus was mutated with a 5% probability"
                mutationProbability = random.randint(1, 100)
                if mutationProbability <= 5:
                    swappedSubnets = self.create_swapped_subnets(subnet[1])
                    print(f"subnet to swap: {subnet[1]}")
                    statNumberOfMutatedSubnets += 1
                    mutationStepSubnetworks.extend(swappedSubnets)
                    # break
                else:
                    print(f"adding original subnetwork: {subnet[1]['subnet']}")
                    flatSubnet = []
                    for gene in subnet[1]['subnet']:
                        if isinstance(gene, list):
                            for item in gene:
                                if item not in subnet
                                flatSubnet.append(item)

                    # mutationStepSubnetworks.append(subnet)
        print(f"number of subnetworks after mutation: {len(mutationStepSubnetworks)}")
        currentDir = os.path.dirname(os.path.abspath(__file__))
        relativePath = os.path.join(currentDir, "..", "created_at_runtime", "temp.txt")

        with open(relativePath, "w") as outputFile:
            json.dump(mutationStepSubnetworks, outputFile)

        return statNumberOfMutatedSubnets, mutationStepSubnetworks

    def create_swapped_subnets(self, subnetObj):
        subnet = []
        swappedSubnets = []
        # subnet = [item for sublist in subnetObj.values() for item in sublist]

        # flatten subnet if sublist exists
        for sublist in subnetObj.values():
            for gene in sublist:
                if isinstance(gene, list):
                    for item in gene:
                        if item not in subnet:
                            subnet.append(item)
                if gene not in subnet and isinstance(gene, list) == False:
                    subnet.append(gene)

        for gene in subnet:
            geneIndex = subnet.index(gene)
            locus = self.faUtilitiesInstance.find_gene_locus(gene)[0]
            for locusGene in locus:
                tempSubnet = subnet.copy()
                # print(f"tempSubnet:  {tempSubnet}")
                tempSubnet[geneIndex] = locusGene
                swappedSubnets.append(tempSubnet)
                # print(f"after: {tempSubnet}")
        return swappedSubnets
