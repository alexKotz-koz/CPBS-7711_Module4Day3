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
        # creating a filtered Parent Network: Only contains FA to FA connections.
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
    # Output: filtered FA parent dataframe
    def create_parent_network(self):
        print("Creating Parent Network")
        start = time.time()

        self.filter_parent_network()

        # read in faNetwork and store each row into the dataframe
        self.parentNetwork = pd.read_csv(
            "created_at_runtime/faNetwork.txt",
            sep="\t",
            header=None,
            names=["gene1", "gene2", "weight"],
            usecols=[0, 1, 2],
        )

        # convert gene1, gene2, and weight to strings for sorting
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

    # Input: STRING 1.txt
    # Output: masterParentDataFrame | This df contains all unduplicated rows from STRING 1.txt and all genes from Input.gmt.txt
    def create_master_parent_network(self, faLoci):
        start = time.time()

        # read in STRING 1.txt and store each row into the dataframe
        parentNetwork = pd.read_csv(
            "static/STRING 1.txt",
            sep="\t",
            header=None,
            names=["gene1", "gene2", "weight"],
            usecols=[0, 1, 2],
        )

        # convert gene1, gene2, and weight to strings for sorting
        parentNetwork["gene1"] = parentNetwork["gene1"].astype(str)
        parentNetwork["gene2"] = parentNetwork["gene2"].astype(str)
        parentNetwork["weight"] = parentNetwork["weight"].astype(str)

        # create a set of sorted gene pairs
        sortedGenePairs = set(
            map(tuple, np.sort(parentNetwork[["gene1", "gene2"]].values, axis=1))
        )

        # filter the data frame based on the set of sorted gene pairs
        # note: axis=1 sorts on rows
        parentNetwork = parentNetwork[
            parentNetwork.apply(
                lambda row: tuple(sorted([row["gene1"], row["gene2"]]))
                in sortedGenePairs,
                axis=1,
            )
        ]

        # create a new DataFrame with genes from faLoci
        fagenes = [gene for sublist in faLoci.values() for gene in sublist]
        faLociDF = pd.DataFrame(fagenes, columns=["gene1"])

        # add gene2 and weight columns with "FAGENEROW" and 0 values respectively
        faLociDF["gene2"] = "FAGENEROW"
        faLociDF["weight"] = 0

        # concatenate faLoci_df and parentNetwork
        masterParentNetwork = pd.concat([parentNetwork, faLociDF])
        masterParentNetwork.drop_duplicates(
            subset=["gene1", "gene2", "weight"], inplace=True
        )

        # create a gene index dictionary that contains {gene: index of row where gene exists in masterParentNetwork}.
        # this is used in the subnet weight calculation to decrease runtime of searching for each of the 12 genes in the parentNetwork dataframe.
        self.geneIndexDict = {
            gene: list(indices)
            for gene, indices in masterParentNetwork.groupby("gene1").groups.items()
        }

        end = time.time()
        ex = end - start
        print(f"Parent Null Case Network finished in : {ex}")

        return masterParentNetwork, self.geneIndexDict

    # Input: Input.gmt.txt
    # Output: loci dictionary, containing one list per locus
    def extract_loci(self):
        loci = {}
        with open(self.inputFile, "r") as file:
            for line in file:
                name = line.split()
                loci[name[3]] = line.strip().split("\t")[2:]
        return loci

    # Input: module1FaNetworkFile | A network file containing FA-FA gene connections
    # Output: a list of FA-FA gene connections
    def extract_module1_fa_network(self):
        module1FASubnetwork = []
        with open(self.module1FaNetworkFile, "r") as file:
            for row in file:
                row = row.split("\t")
                module1FASubnetwork.append(row)
        return module1FASubnetwork

    # Input: fa loci dictionary, gene from subnetwork
    # Output: the input gene's locus
    def find_gene_locus(self, gene):
        for locus in self.loci:
            if gene in self.loci[locus]:
                return self.loci[locus], locus

    # Input: subnetGenes(subnetwork), filtered parentFANetwork dataframe
    # Ouput: sum of all weights, signifying the density of the subnetwork
    # USE: FA GENETIC ALGORITHM
    def count_edges(self, subnetGenes, parentNetwork):
        # create a mask to find all rows where both genes in the row are in the input subnetwork
        mask = parentNetwork["gene1"].isin(subnetGenes) & parentNetwork["gene2"].isin(
            subnetGenes
        )
        # create a copy of the mask for further data manipulation
        selectedRows = parentNetwork[mask].copy()

        # sort the mask by row
        selectedRows[["gene1", "gene2"]] = np.sort(
            selectedRows[["gene1", "gene2"]], axis=1
        )

        selectedRows.drop_duplicates(subset=["gene1", "gene2"], inplace=True)

        # calculate the sum of all weights found in the mask
        weightSum = selectedRows["weight"].astype(float).sum()
        # print(f"weightSum: {weightSum}")
        return weightSum

    # Input: subnetGenes (subnetwork), parentNetwork (masterParentNetwork dataframe), geneIndexDict (dictionary containing all genes and the index of the rows where they appear in the dataframe)
    # Output: a sum of the weights found in the function, signifying the density of the subnetwork
    # USE: NONFA (NULL CASE) GENETIC ALGORITHM
    # This is redundant and can be collapsed into the count_edges function
    def count_edges_null_case(self, subnetGenes, parentNetwork, geneIndexDict):
        # get the indices of the rows that contain the genes in subnetGenes
        indices = [
            index
            for gene in subnetGenes
            if gene in geneIndexDict
            for index in geneIndexDict[gene]
        ]
        # get the rows from the parentNetwork, based on the indices retrieved from geneIndexDict
        selectedRows = parentNetwork.loc[indices].copy()

        # sort the mask by row
        selectedRows[["gene1", "gene2"]] = np.sort(
            selectedRows[["gene1", "gene2"]], axis=1
        )

        selectedRows.drop_duplicates(subset=["gene1", "gene2"], inplace=True)

        # calculate the sum of all weights found in the mask
        weightSum = selectedRows["weight"].astype(float).sum()
        # print(f"weightSum: {weightSum}")
        return weightSum

    # Input: faNetworkFile (a file containing the FA-FA gene connections (from M1)), geneScoresFile (a file containing 5000 versions of [lociNumber: [{gene1:genescore}, ...]])
    # Ouput: an average of all gene scores
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
        print(f"genescoresfromfile:{geneScoresFromFile}")
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

    # Input: bins object
    # Ouput: a look up object that contains the parent network row indexs where each gene exists
    def genes_to_bins(self, bins):
        invertedBins = {
            gene: binName for binName, genes in bins.items() for gene in genes
        }
        return invertedBins

