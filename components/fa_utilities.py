import pandas as pd
import numpy as np
import os, json, time


class Fa_Utilities:
    def __init__(
        self,
        parentNetworkFile=None,
        individualSubnetwork=None,
        inputFile=None,
        module1FaNetworkFile=None,
        loci=None,
    ):
        self.individualSubnetwork = individualSubnetwork
        self.inputFile = inputFile
        self.module1FaNetworkFile = module1FaNetworkFile
        self.loci = loci

        # check to identify parentNetworkFile datatype
        # (useful for modularity, not complete abstration will need refactor to make completely uncoupled)
        if isinstance(parentNetworkFile, pd.DataFrame):
            self.parentNetwork = parentNetworkFile
        elif isinstance(parentNetworkFile, list):
            self.parentNetwork = parentNetworkFile
        elif isinstance(parentNetworkFile, str):
            self.parentNetworkFile = parentNetworkFile
            self.parentNetwork = pd.DataFrame()
        else:
            self.parentNetworkFile = parentNetworkFile

    # Input: parentNetworkFile (assuming .txt format)
    # Output: faNetwork list (contains a sublist for each FA-FA connection), creation of filtered faNetwork text file
    def filter_parent_network(self):
        # creating a filtered Parent Network: Only contains FA to FA connections. Limitation.
        # print("Filtering Parent Network for FA Genes")
        faLoci = self.extract_loci()
        faGenes = [string for sublist in faLoci.values() for string in sublist]
        faNetwork = []
        # read in parentNetworkFile, if both genes in each row are FA genes, add to faNetwork
        with open(self.parentNetworkFile, "r") as file:
            for line in file:
                line = line.strip().split("\t")
                if line[0] in faGenes:
                    if line[1] in faGenes:
                        faNetwork.append(line)
                if line[1] in faGenes:
                    if line[0] in faGenes:
                        faNetwork.append(line)

        currentDir = os.path.dirname(os.path.abspath(__file__))
        relativePath = os.path.join(
            currentDir, "..", "created_at_runtime", "faNetwork.txt"
        )

        with open(relativePath, "w") as file:
            file.write("\n".join("\t".join(sublist) for sublist in faNetwork))
        return faNetwork

    # Input: faNetwork.txt (generated from self.filter_parent_network())
    # Output: either filtered fa network dataframe or unfiltered parent network dataframe
    def create_parent_network(self):
        print("Creating Parent Network")
        start = time.time()

        # note: remove this line if desired result is a full parent network object
        self.filter_parent_network()

        # read in faNetwork and store the first two columns in a dataframe
        # note: replace "faNetwork.txt" with "STRING 1.txt" to create full parent network
        self.parentNetwork = pd.read_csv(
            "created_at_runtime/faNetwork.txt",
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
        print(f"Parent Network finished in : {ex}")
        return self.parentNetwork

    def create_master_parent_network(self, faLoci):
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
        fagenes = [gene for sublist in faLoci.values() for gene in sublist]
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

    # Input: Input.gmt.txt
    # Output: loci dictionary, containing one list per locus
    def extract_loci(self):
        loci = {}
        # print("Extracting FA Loci")
        with open(self.inputFile, "r") as file:
            for line in file:
                name = line.split()
                loci[name[3]] = line.strip().split("\t")[2:]
        # print("Loci Extracted")
        return loci

    def extract_module1_fa_network(self):
        module1FASubnetwork = []
        with open(self.module1FaNetworkFile, "r") as file:
            for row in file:
                row = row.split("\t")
                module1FASubnetwork.append(row)
        return module1FASubnetwork

    def find_gene_locus(self, gene):
        for locus in self.loci:
            if gene in self.loci[locus]:
                return self.loci[locus], locus

    def count_edges(self, subnetGenes, parentNetwork):
        mask = parentNetwork["gene1"].isin(subnetGenes) & parentNetwork["gene2"].isin(
            subnetGenes
        )

        selectedRows = parentNetwork[mask].copy()

        # sort the mask by row
        selectedRows["sorted_genes"] = np.sort(
            selectedRows[["gene1", "gene2"]], axis=1
        ).tolist()

        selectedRows.drop_duplicates(subset=["sorted_genes"], inplace=True)

        weightSum = (selectedRows["weight"].astype(float)).sum()
        # print(f"weight sum: {weightSum}")
        return weightSum

    def calculate_average_gene_score(self, faNetworkFile, geneScoresFile):
        faNetwork = []
        geneScoresFromFile = []
        scoresByGene = {}
        averageGeneScores = {}

        with open(faNetworkFile, "r") as file:
            for line in file:
                line = line.split()
                faNetwork.append(line)
            # create gene scores dictionary
        with open(geneScoresFile, "r") as file:
            for line in file:
                locusId = line.split()[0][:-1]
                locus_str = " ".join(line.split()[1:])
                locus_str = locus_str.replace("'", '"')
                locus = json.loads(locus_str)
                for gene in locus:
                    geneScoresFromFile.append({locusId: gene})

        # restructure geneScoresFromFile for use in calculating average gene scores
        for score in geneScoresFromFile:
            locusId = ",".join(score.keys())
            score = list(score.values())[0]
            gene = score["gene"]
            if gene not in scoresByGene:
                scoresByGene[gene] = {"locusId": locusId, "scores": []}
                scoresByGene[gene]["locusId"] = locusId
            scoresByGene[gene]["scores"].append(score["geneScore"])

        # calculate average gene scores and store in new dictionary
        for item in scoresByGene.items():
            gene = item[0]
            subitem = item[1]

            locusId = subitem["locusId"]
            scores = subitem["scores"]

            averageGeneScores[gene] = {
                "averageScore": sum(scores) / len(scores),
                "locusId": locusId,
            }

        # if a gene from averageGeneScores is not in the FA-FA network, mark genescore as "NA"
        for gene in averageGeneScores:
            if not any(gene in sublist for sublist in faNetwork):
                averageGeneScores[gene]["averageScore"] = "NA"

        return averageGeneScores

    def genes_to_bins(self, bins):
        invertedBins = {
            gene: binName for binName, genes in bins.items() for gene in genes
        }
        return invertedBins
