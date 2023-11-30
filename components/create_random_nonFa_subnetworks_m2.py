import os
import json
import pandas as pd
import numpy as np
import time
import random
import itertools
import concurrent


class Create_Random_Non_Fa_Subnetworks:
    def __init__(self, parentNetworkFile, faLoci):
        self.parentNetworkFile = parentNetworkFile
        self.faLoci = faLoci
        self.parentNetwork = pd.DataFrame()
        self.nfaSubnetworks = {}

    def create_non_fa_subnetworks(self):
        nfaNetwork = self.extract_non_fa_genes()

        self.create_parent_network()

        geneCounts = self.count_edges()

        bins = self.create_bins(geneCounts)
        currentDir = os.path.dirname(os.path.abspath(__file__))
        relativePath = os.path.join(currentDir, "..", "created_at_runtime", "bins.json")

        with open(relativePath, "w") as outputFile:
            json.dump(bins, outputFile)
        initialSubnet = self.create_initial_individual_non_fa_subnetwork(nfaNetwork)

        self.nfaSubnetworks[1] = initialSubnet

        subnets = self.create_random_non_fa_subnetworks(bins)
        currentDir = os.path.dirname(os.path.abspath(__file__))
        relativePath = os.path.join(
            currentDir, "..", "created_at_runtime", "nfaSubnetworks.json"
        )

        with open(relativePath, "w") as outputFile:
            json.dump(self.nfaSubnetworks, outputFile)

    def create_random_non_fa_subnetworks(self, bins):
        geneToBin = {
            gene: self.find_bin(gene, bins)
            for subnet in self.nfaSubnetworks.values()
            for gene in subnet
        }

        new_subnetworks = []  # List to hold the new subnetworks

        for _ in range(5000):
            tempsubnet = set()
            for subnet in self.nfaSubnetworks.values():
                for gene in subnet:
                    bin = geneToBin[gene]

                    if len(bin) <= 1:
                        randomIndex = 0
                    else:
                        randomIndex = random.randint(0, len(bin) - 1)

                    newGene = bin[randomIndex]
                    if newGene not in tempsubnet:
                        tempsubnet.add(newGene)

            new_subnetworks.append(
                list(tempsubnet)
            )  # Add the new subnetwork to the list

        # Add all the new subnetworks to self.nfaSubnetworks at once
        start_key = list(self.nfaSubnetworks.keys())[-1] + 1
        self.nfaSubnetworks.update(
            {i + start_key: subnet for i, subnet in enumerate(new_subnetworks)}
        )

        return self.nfaSubnetworks

    def create_initial_individual_non_fa_subnetwork(self, nfaNetwork):
        nfaGenes = set(gene for row in nfaNetwork for gene in row[:2])
        initialNfaSubnet = []
        randomIndicies = []
        for _ in range(0, 12):
            randomIndex = random.randint(0, len(nfaGenes) - 1)
            if randomIndex not in randomIndicies:
                randomIndicies.append(randomIndex)
            else:
                randomIndex = random.randint(0, len(nfaGenes) / 2)
                randomIndicies.append(randomIndex)

        nfaGenes = list(nfaGenes)
        for item in randomIndicies:
            initialNfaSubnet.append(nfaGenes[item])
        return initialNfaSubnet

    def create_bins(self, geneCounts):
        minGeneEdgeCount = int(geneCounts.values.min())
        maxGeneEdgeCount = int(geneCounts.values.max())

        edgeCountRange = maxGeneEdgeCount - minGeneEdgeCount
        binSize = edgeCountRange / 128  # Calculate bin size

        bins = {}

        for gene, count in geneCounts.items():
            if count not in bins:
                bins[count] = []  # Create the key with an empty list as its value

            bins[count].append(gene)

        return bins

    def find_bin(self, gene, bins):
        for bin in bins.items():
            if gene in bin[1]:
                return bin[1]
        print(f"error finding bin: {gene}")

    def create_parent_network(self):
        start = time.time()
        # read in faNetwork and store the first two columns in a dataframe
        # note: replace "faNetwork.txt" with "STRING 1.txt" to create full parent network
        self.parentNetwork = pd.read_csv(
            "created_at_runtime/nonFaNetwork.txt",
            sep="\t",
            header=None,
            names=["gene1", "gene2", "weight"],
            usecols=[0, 1, 2],
        )

        # convert gene1 and gene2 to strings
        self.parentNetwork["gene1"] = self.parentNetwork["gene1"].astype(str)
        self.parentNetwork["gene2"] = self.parentNetwork["gene2"].astype(str)
        self.parentNetwork["weight"] = self.parentNetwork["weight"].astype(str)

        # create a set of sorted gene pairs
        sortedGenePairs = set(
            map(tuple, np.sort(self.parentNetwork[["gene1", "gene2"]].values, axis=1))
        )

        # filter the data frame based on the set of sorted gene pairs
        self.parentNetwork = self.parentNetwork[
            self.parentNetwork.apply(
                lambda row: tuple(sorted([row["gene1"], row["gene2"]]))
                in sortedGenePairs,
                axis=1,
            )
        ]

        end = time.time()
        ex = end - start
        print(f"Parent NFA Network finished in : {ex}")
        return self.parentNetwork

    def count_edges(self):
        parentNetwork = self.parentNetwork
        geneCounts = pd.concat(
            [self.parentNetwork["gene1"], self.parentNetwork["gene2"]]
        ).value_counts()
        return geneCounts

    def extract_non_fa_genes(self):
        faGenes = [string for sublist in self.faLoci.values() for string in sublist]
        nfaNetwork = []
        # read in parentNetworkFile, if both genes in each row are FA genes, add to faNetwork
        with open(self.parentNetworkFile, "r") as file:
            for line in file:
                line = line.strip().split("\t")
                if line[0] not in faGenes and line[1] not in faGenes:
                    nfaNetwork.append(line)

        currentDir = os.path.dirname(os.path.abspath(__file__))
        relativePath = os.path.join(
            currentDir, "..", "created_at_runtime", "nonFaNetwork.txt"
        )

        with open(relativePath, "w") as file:
            file.write("\n".join("\t".join(sublist) for sublist in nfaNetwork))
        return nfaNetwork
