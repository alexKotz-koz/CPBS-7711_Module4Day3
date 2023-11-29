import json
import random
import time
import os

from components.fa_utilities import Fa_Utilities


class Create_Random_Fa_Subnetworks:
    def __init__(
        self,
        prevFaSubnetworkFile,
        faInputFile,
        stringInputFile,
        faLoci,
        module1FaNetwork,
        parentNetwork,
    ):
        self.prevFaSubnetworkFile = prevFaSubnetworkFile
        self.faInputFile = faInputFile
        self.stringInputFile = stringInputFile
        self.faLoci = faLoci
        self.module1FaNetwork = module1FaNetwork
        self.module1FASubnetwork = []
        self.faGenes = []
        self.sortedDictionary = {}
        self.parentNetwork = parentNetwork
        self.test = {}

    # Input: loci dictionary (faGenes)
    def generate_12_genes(self):
        faGenes = self.faLoci

        genesForSubnetwork = set()
        # for each locus in the faGenes dictionary, extract one gene at random, from each locus and store in a new list
        for index, item in enumerate(faGenes):
            random_int = random.randint(0, len(faGenes[item]))
            try:
                genesForSubnetwork.add(faGenes[item][random_int])
            except IndexError:
                random_int = random.randrange(0, len(faGenes[item]))
                genesForSubnetwork.add(faGenes[item][random_int])
        return list(genesForSubnetwork)

    def extract_module1_fa_network(self):
        module1FASubnetwork = []
        with open(self.prevFaSubnetworkFile, "r") as file:
            for row in file:
                row = row.split("\t")
                module1FASubnetwork.append(row)
        return module1FASubnetwork

    """def create_individual_subnetwork(self, module1FASubnetwork):
        subnetworkToWrite = []
        flattenedSubnetwork = []
        geneSet12 = self.generate_12_genes()

        print(f"a: {geneSet12}")
        for gene in geneSet12:
            for row in module1FASubnetwork:
                if gene == row[0] and row[1] in geneSet12:
                    subnetworkToWrite.append(row[:2])
                elif gene == row[1] and row[0] in geneSet12:
                    subnetworkToWrite.append(row[:2])

        print(f"1: {subnetworkToWrite}")

        for item in subnetworkToWrite:
            for gene in item:
                flattenedSubnetwork.append(gene)
        print(f"2: {flattenedSubnetwork}")

        for gene in geneSet12:
            if gene not in flattenedSubnetwork:
                subnetworkToWrite.append(gene)

        print(f"3: {subnetworkToWrite}")

        subnetworkToWrite = [
            sublist
            for index, sublist in enumerate(subnetworkToWrite)
            if sublist not in subnetworkToWrite[:index]
        ]

        print(f"4: {subnetworkToWrite}")

        return subnetworkToWrite"""

    def check_individual_subnet_edge_count(self):
        individualSubnetwork = self.generate_12_genes()

        faUtilitiesInstance = Fa_Utilities()

        individualSubnetworkEdgeCount = faUtilitiesInstance.count_edges(
            subnetGenes=individualSubnetwork, parentNetwork=self.parentNetwork
        )

        return individualSubnetworkEdgeCount, individualSubnetwork

        # filter for subnetworks that only return > 0 edge count
        """if individualSubnetworkEdgeCount == 0:
            return self.check_individual_subnet_edge_count()
        elif individualSubnetworkEdgeCount > 0:
            return individualSubnetworkEdgeCount, individualSubnetwork
"""

    def create_random_subnetworks(self):
        print("Creating stage 1 random subnetworks...")
        module1FASubnetwork = self.module1FaNetwork

        finalList = []
        finalDictionary = {}

        count = 0
        while count < 5000:
            edgeCount, individualSubnetwork = self.check_individual_subnet_edge_count()
            finalList.append(individualSubnetwork)
            count += 1

        for index, item in enumerate(finalList):
            index = str(index)

            finalDictionary[index] = {"subnet": item}

        currentDir = os.path.dirname(os.path.abspath(__file__))
        relativePath = os.path.join(
            currentDir, "..", "created_at_runtime", "random_fa_subnetworks.json"
        )

        with open(relativePath, "w") as outputFile:
            json.dump(finalDictionary, outputFile)
        print("First 5,000 subnetworks created")

        return finalDictionary
