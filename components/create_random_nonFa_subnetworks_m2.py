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
        self, parentNetworkFile, faLoci, finalPopulationSubnets, parentFaNetworkDF
    ):
        self.parentNetworkFile = parentNetworkFile
        self.faLoci = faLoci
        self.finalPopulationSubnets = finalPopulationSubnets
        self.parentFaNetworkDF = parentFaNetworkDF
        self.nfaSubnetworks = {}
        self.bins = ""
        self.nfaGenes = ""
        self.parentNetworkDF = ""
        self.invertedBins = ""
        self.faUtilitiesInstance = Fa_Utilities()

    def create_non_fa_subnetworks(self):
        nfaSubnetworksToScore = []

        self.nfaGenes = self.extract_non_fa_genes()

        print(f"Number of unique nfa genes: {len(self.nfaGenes)}")

        # Create parent network data frame from String file
        parentNetworkDF = self.create_parent_network()
        self.parentNetworkDF = parentNetworkDF

        geneDensity = self.calculate_weights_for_bins(parentNetworkDF)

        print(f"Number of unique genes: {len(geneDensity)}")

        self.bins = self.create_bins(geneDensity)

        self.invertedBins = self.genes_to_bins(self.bins)

        """currentDir = os.path.dirname(os.path.abspath(__file__))
        relativePath = os.path.join(
            currentDir, "..", "created_at_runtime", "invertedbins.json"
        )
        with open(relativePath, "w") as outputFile:
            json.dump(self.invertedBins, outputFile)

        currentDir = os.path.dirname(os.path.abspath(__file__))
        relativePath = os.path.join(currentDir, "..", "created_at_runtime", "bins.json")
        with open(relativePath, "w") as outputFile:
            json.dump(self.bins, outputFile)"""

        # test
        """testData = list(itertools.islice(self.finalPopulationSubnets, 1))
        self.create_individual_non_fa_subnet(testData)"""

        # Thread pool for final population
        substart = time.time()
        with ThreadPoolExecutor() as executor:
            for individualSubnet in executor.map(
                self.create_individual_non_fa_subnet, self.finalPopulationSubnets
            ):
                nfaSubnetworksToScore.append(individualSubnet)
        subend = time.time()
        print(f"thread pool stopped: {subend-substart}")
        self.nfaSubnetworks = {
            "averageDensity": self.calculate_average_density(nfaSubnetworksToScore),
            "subnetworks": nfaSubnetworksToScore,
        }

        currentDir = os.path.dirname(os.path.abspath(__file__))
        relativePath = os.path.join(
            currentDir, "..", "created_at_runtime", "random_non_fa_subnetworks.json"
        )
        with open(relativePath, "w") as outputFile:
            json.dump(self.nfaSubnetworks, outputFile)

    def create_individual_non_fa_subnet(self, subnet):
        # create a dictionary: {density score: density score, subnet: subnet}
        individualSubnet = []

        for gene in subnet:
            genesBinName = self.invertedBins[gene]
            genesBin = self.bins[genesBinName]
            # added = False

            if all(gene not in self.nfaGenes for gene in genesBin):
                individualSubnet.append("NA")

            randomIndex = random.randint(0, len(genesBin) - 1)
            randomGene = genesBin[randomIndex]

            if randomGene in self.nfaGenes and randomGene not in individualSubnet:
                individualSubnet.append(randomGene)
                # added = True

        print(individualSubnet)
        return individualSubnet

    def genes_to_bins(self, bins):
        invertedBins = {
            gene: binName for binName, genes in bins.items() for gene in genes
        }
        return invertedBins

    def calculate_weights_for_bins(self, parentNetworkDF):
        print("Creating bins for null case")
        start = time.time()

        # Convert weight column to float
        parentNetworkDF["weight"] = parentNetworkDF["weight"].astype(float)

        # Concatenate gene1 and gene2 columns and calculate mean for each gene
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

    def create_bins(self, geneCounts):
        start = time.time()
        minGeneEdgeCount = min(geneCounts.values())
        maxGeneEdgeCount = max(geneCounts.values())
        print(f"min: {minGeneEdgeCount}")
        print(f"max:{maxGeneEdgeCount}")

        edgeCountRange = maxGeneEdgeCount - minGeneEdgeCount
        binSize = edgeCountRange / 128  # Calculate bin size

        bins = {i: [] for i in range(128)}  # Initialize bins

        for gene, weight in geneCounts.items():
            # Calculate which bin the weight falls into
            binIndex = int((weight - minGeneEdgeCount) / binSize)
            # Ensure the maximum weight falls into the last bin
            binIndex = min(binIndex, 127)

            bins[binIndex].append(gene)

        end = time.time()
        print(f"Bins created in: {end-start}")
        return bins

    def create_parent_network(self):
        start = time.time()
        # read in faNetwork and store the first two columns in a dataframe
        # note: replace "faNetwork.txt" with "STRING 1.txt" to create full parent network
        parentNetwork = pd.read_csv(
            "static/STRING 1.txt",
            sep="\t",
            header=None,
            names=["gene1", "gene2", "weight"],
            usecols=[0, 1, 2],
        )

        # convert gene1 and gene2 to strings
        parentNetwork["gene1"] = parentNetwork["gene1"].astype(str)
        parentNetwork["gene2"] = parentNetwork["gene2"].astype(str)
        parentNetwork["weight"] = parentNetwork["weight"].astype(str)

        # create a set of sorted gene pairs
        sortedGenePairs = set(
            map(tuple, np.sort(parentNetwork[["gene1", "gene2"]].values, axis=1))
        )

        # filter the data frame based on the set of sorted gene pairs
        parentNetwork = parentNetwork[
            parentNetwork.apply(
                lambda row: tuple(sorted([row["gene1"], row["gene2"]]))
                in sortedGenePairs,
                axis=1,
            )
        ]
        fagenes = [gene for sublist in self.faLoci.values() for gene in sublist]
        # Create a new DataFrame with genes from faLoci
        faLociDF = pd.DataFrame(fagenes, columns=["gene1"])

        # Add gene2 and weight columns with null and 0 values respectively
        faLociDF["gene2"] = "FAGENEROW"
        faLociDF["weight"] = 0

        # Concatenate faLoci_df and parentNetwork
        masterParentNetwork = pd.concat([parentNetwork, faLociDF])
        end = time.time()
        ex = end - start
        print(f"Parent Null Case Network finished in : {ex}")

        return masterParentNetwork

    def extract_non_fa_genes(self):
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

    def calculate_average_density(self, subnets):
        weights = []
        for subnet in subnets:
            print("score subnet")
            weights.append(
                float(
                    self.faUtilitiesInstance.count_edges(subnet, self.parentNetworkDF)
                )
            )
        return sum(weights) / len(subnets)
