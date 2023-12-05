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
    # Ouput: a list of 12 random genes (one from each locus)
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

    # Input: module 1 Fa subnetwork
    # Output: a list of the subnetwork
    def extract_module1_fa_network(self):
        module1FASubnetwork = []
        with open(self.prevFaSubnetworkFile, "r") as file:
            for row in file:
                row = row.split("\t")
                module1FASubnetwork.append(row)
        return module1FASubnetwork

    # Input: individual subnetwork from generation_12_genes()
    # Ouput: weight of individual subnetwork and the subnetwork itself
    def check_individual_subnet_edge_count(self):
        individualSubnetwork = self.generate_12_genes()

        faUtilitiesInstance = Fa_Utilities()
        # create an instance of the fa utilities and call the weight calculation function (passing in the parentFaNetworkDF)
        individualSubnetworkEdgeCount = faUtilitiesInstance.count_edges(
            subnetGenes=individualSubnetwork, parentNetwork=self.parentNetwork
        )

        return individualSubnetworkEdgeCount, individualSubnetwork

    #Input: module1FaNetwork list, individual subnetwork and weight from check_individual_subnet_edge_count
    #Ouput: subnetwork population
    def create_random_subnetworks(self):
        print("Creating initial FA population...")
        module1FASubnetwork = self.module1FaNetwork

        finalList = []
        finalDictionary = {}

        # create 5000 individual subnetworks and weigh them
        count = 0
        while count < 5000:
            edgeCount, individualSubnetwork = self.check_individual_subnet_edge_count()
            finalList.append(individualSubnetwork)
            count += 1

        # write subnetwork population to file
        for index, item in enumerate(finalList):
            index = str(index)

            finalDictionary[index] = {"subnet": item}

        currentDir = os.path.dirname(os.path.abspath(__file__))
        relativePath = os.path.join(
            currentDir, "..", "created_at_runtime", "random_fa_subnetworks.json"
        )

        with open(relativePath, "w") as outputFile:
            json.dump(finalDictionary, outputFile)
        print("Initial FA subnetwork population created")

        return finalDictionary
