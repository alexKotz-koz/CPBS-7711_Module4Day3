import json
import random
import os


class Create_Random_Fa_Subnetworks:
    def __init__(
        self,
        prevFaSubnetworkFile,
        faInputFile,
        stringInputFile,
        faLoci,
        module1FaNetwork,
    ):
        self.prevFaSubnetworkFile = prevFaSubnetworkFile
        self.faInputFile = faInputFile
        self.stringInputFile = stringInputFile
        self.faLoci = faLoci
        self.module1FaNetwork = module1FaNetwork
        self.module1FASubnetwork = []
        self.faGenes = []
        self.sortedDictionary = {}

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

    def create_individual_subnetwork(self, module1FASubnetwork):
        subnetworkToWrite = []
        flattenedSubnetwork = []
        geneSet12 = self.generate_12_genes()

        for gene in geneSet12:
            for row in module1FASubnetwork:
                if gene == row[0] and row[1] in geneSet12:
                    subnetworkToWrite.append(row[:2])
                elif gene == row[1] and row[0] in geneSet12:
                    subnetworkToWrite.append(row[:2])

        for item in subnetworkToWrite:
            for gene in item:
                flattenedSubnetwork.append(gene)

        for gene in geneSet12:
            if gene not in flattenedSubnetwork:
                subnetworkToWrite.append(gene)

        subnetworkToWrite = [
            sublist
            for index, sublist in enumerate(subnetworkToWrite)
            if sublist not in subnetworkToWrite[:index]
        ]

        return subnetworkToWrite

    def create_random_subnetworks(self):
        print("Creating stage 1 random subnetworks...")
        module1FASubnetwork = self.module1FaNetwork

        finalList = []
        finalDictionary = {}

        count = 0
        while count < 5000:
            individualSubnetwork = []
            individualSubnetwork = self.create_individual_subnetwork(
                module1FASubnetwork
            )
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
