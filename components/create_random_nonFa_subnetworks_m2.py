import os
import json


class Create_Random_Non_Fa_Subnetworks:
    def __init__(self, parentNetworkFile, faLoci):
        self.parentNetworkFile = parentNetworkFile
        self.faLoci = faLoci

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
        return faGenes, nfaNetwork

    def create_bins(self):
        pass

    def create_individual_non_fa_subnetwork(self):
        pass

    def create_non_fa_subnetworks(self):
        faGenes, nfaNetwork = self.extract_non_fa_genes()
        for i in nfaNetwork:
            if i[0] in faGenes or i[1] in faGenes:
                print(f"here: {i}")
