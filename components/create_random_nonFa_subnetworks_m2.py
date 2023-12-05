import os
import json
import pandas as pd
import numpy as np
import time
import random
import itertools
from concurrent.futures import ThreadPoolExecutor
from components.fa_utilities import Fa_Utilities


class Create_Random_Non_Fa_Subnetworks:
    def __init__(
        self,
        parentNetworkFile,
        faLoci,
        finalPopulationSubnets,
        parentFaNetworkDF,
        masterParentNetworkDF,
    ):
        self.parentNetworkFile = parentNetworkFile
        self.faLoci = faLoci
        self.finalPopulationSubnets = finalPopulationSubnets
        self.parentFaNetworkDF = parentFaNetworkDF
        self.nfaSubnetworks = {}
        self.bins = ""
        self.nfaGenes = ""
        self.parentNetworkDF = masterParentNetworkDF
        self.invertedBins = ""
        self.faUtilitiesInstance = Fa_Utilities()

    # Input: self (containing nessecary objects for non fa subnetwork creation)
    # Output: non fa subnetwork population with edge counts
    def create_non_fa_subnetworks(self):
        nfaSubnetworksToScore = []
        
        #extract non fa genes
        self.nfaGenes = self.extract_non_fa_genes()

        # prepare weight range for bins
        geneDensity = self.calculate_weights_for_bins(self.parentNetworkDF)

        # create bins based on average gene weight 
        self.bins = self.create_bins(geneDensity)

        # create gene index look up object to improve runtime of finding a genes bin
        self.invertedBins = self.genes_to_bins(self.bins)

        # create a thread pool to asynchronously create non fa subnetworks
        with ThreadPoolExecutor() as executor:
            for individualSubnet in executor.map(
                self.create_individual_non_fa_subnet, self.finalPopulationSubnets
            ):
                nfaSubnetworksToScore.append(individualSubnet)

        # add population average density and subnetworks to object to write to a file
        self.nfaSubnetworks = {
            "averageDensity": self.calculate_average_density(nfaSubnetworksToScore),
            "subnetworks": nfaSubnetworksToScore,
        }

        currentDir = os.path.dirname(os.path.abspath(__file__))
        relativePath = os.path.join(
            currentDir, "..", "created_at_runtime", "random_non_fa_subnetworks.json"
        )
        with open(relativePath, "a") as outputFile:
            json.dump(self.nfaSubnetworks, outputFile)
            outputFile.write(",")

    # Input: a subnetwork from the final fa population
    # Ouput: a non fa subnetwork
    def create_individual_non_fa_subnet(self, subnet):
        # create a dictionary: {density score: density score, subnet: subnet}
        individualSubnet = []
        subnet = subnet["subnet"]

        for gene in subnet:
            # find genes bin
            genesBinName = self.invertedBins[gene]
            genesBin = self.bins[genesBinName]

            # create random index based on the size of the genes bin
            randomIndex = random.randint(0, len(genesBin) - 1)
            # grab a random gene from the bin
            randomGene = genesBin[randomIndex]

            # make sure the bin contains a non fa gene, if it doesn't just add "NA" to the subnetwork, indicating the bin only contained FA genes
            if all(gene not in self.nfaGenes for gene in genesBin):
                individualSubnet.append("NA")

            # check to see if the new gene is a non fa gene and make sure its not already added to the subentwork
            elif randomGene in self.nfaGenes and randomGene not in individualSubnet:
                individualSubnet.append(randomGene)
            
            # catch all for conditions above
            else:
                # create a new random index, find a new gene from the original genes bin, check to see if new gene is a non fa gene
                # repeat ^ until all conditions are met
                while True:
                    randomIndex = random.randint(0, len(genesBin) - 1)
                    randomGene = genesBin[randomIndex]
                    if randomGene in self.nfaGenes:
                        individualSubnet.append(randomGene)
                        break

        return individualSubnet

    # Input: bins object
    # Ouput: a look up object that contains the parent network row indexs where each gene exists
    def genes_to_bins(self, bins):
        invertedBins = {
            gene: binName for binName, genes in bins.items() for gene in genes
        }
        return invertedBins

    # Input: master parent network DF
    # Output: bin incrament size (average weight for all genes)
    def calculate_weights_for_bins(self, parentNetworkDF):
        print("Creating bins for null case")
        start = time.time()

        # convert weight column to float
        parentNetworkDF["weight"] = parentNetworkDF["weight"].astype(float)

        # concatenate gene1 and gene2 columns and calculate mean for each gene
        geneDensity = (
            pd.concat(
                [
                    parentNetworkDF[["gene1", "weight"]],
                    parentNetworkDF[["gene2", "weight"]].rename(
                        columns={"gene2": "gene1"}
                    ),
                ]
            )
            .groupby("gene1")
            .mean()
            .to_dict()["weight"]
        )

        end = time.time()
        print(f"Calculated average edge weights for each gene in: {end - start}")

        return geneDensity
    

    # Input: gene scores object
    # Ouput: bins object containing all genes
    def create_bins(self, geneCounts):
        start = time.time()
        # get the minimum and maximum gene scores
        minGeneEdgeCount = min(geneCounts.values())
        maxGeneEdgeCount = max(geneCounts.values())

        # calculate range of gene scores
        edgeCountRange = maxGeneEdgeCount - minGeneEdgeCount
        # calculate bin size
        binSize = edgeCountRange / 128  

        # initialize bins
        bins = {i: [] for i in range(128)}  

        for gene, weight in geneCounts.items():
            # calculate which bin the weight falls into
            binIndex = int((weight - minGeneEdgeCount) / binSize)
            # ensure the maximum weight falls into the last bin
            binIndex = min(binIndex, 127)
            # add gene to bin
            bins[binIndex].append(gene)

        end = time.time()
        print(f"Bins created in: {end-start}")
        return bins

    # Input: self contains nessecary items for non fa gene extraction process
    # Output: a non fa genes object
    def extract_non_fa_genes(self):
        # extract fa genes from faLoci object
        faGenes = [string for sublist in self.faLoci.values() for string in sublist]
        nfaGenes = set()
        # read in parentNetworkFile, if both genes in each row are FA genes, add to faNetwork
        with open(self.parentNetworkFile, "r") as file:
            for line in file:
                line = line.strip().split("\t")
                if line[0] not in faGenes:
                    nfaGenes.add(line[0])
                if line[1] not in faGenes:
                    nfaGenes.add(line[1])

        return nfaGenes

    # Input: subnetwork
    # Output: weight of subnetwork
    # wrapper function for thread pool
    def count_edges_wrapper(self, subnet):
        return float(self.faUtilitiesInstance.count_edges(subnet, self.parentNetworkDF))

    # Input: non fa subnetwork population
    # Output: average weight of population 
    def calculate_average_density(self, subnets):
        print("Calculating average density of random non fa subnetworks")
        weights = []
        # REFACTOR: Change to process pool becuase count edges is more CPU bound
        with ThreadPoolExecutor() as executor:
            weights = list(executor.map(self.count_edges_wrapper, subnets))

        return sum(weights) / len(subnets)
