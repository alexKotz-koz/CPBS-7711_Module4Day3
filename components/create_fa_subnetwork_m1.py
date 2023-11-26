import os


class Create_Fa_Subnetwork:
    def __init__(self, loci, inputFile):
        self.loci = loci
        self.inputFile = inputFile

    def create_subnetwork(self):
        inputFile = self.inputFile
        faLoci = []

        for locus in self.loci:
            faLoci.extend(gene for gene in self.loci[locus])
        print(faLoci)
        # Create sub-network with FA related genes, using STRING 1.txt

        # 1)Open STRING 1 file, iterate thru each row in file
        # 2)Cross check the genes in each row with the genes in faLoci
        # 3)If both genes are in the faLoci list add to results list
        with open(inputFile, "r") as file:
            results = []
            outliers = []  # (A, B, C, meausurement)
            for row in file:
                row = row.split("\t")

                # Error check format of STRING 1.txt
                if len(row) != 3:
                    print("Data format not consistent with STRING format")
                    break

                if len(results) == 0:
                    if row[0] in faLoci:
                        if row[1] in faLoci:
                            results.append(row)
                else:
                    if row[0] in faLoci:
                        if row[1] in faLoci:
                            for index in range(len(outliers)):
                                if (
                                    row[0] == outliers[index][0]
                                    or row[0] == outliers[index][1]
                                    or row[1] == outliers[index][0]
                                    or row[1] == outliers[index][1]
                                ):
                                    rowToRemove = row
                                    outliers.remove(rowToRemove)
                                    print("Row removed from outliers")

                            for index in range(len(results)):
                                if row[0] != results[index][0]:
                                    if row[0] != results[index][1]:
                                        if row[1] != results[index][0]:
                                            if row[1] != results[index][1]:
                                                results.append(row)
                                                print("Row added to results")
                                                for x in range(len(outliers)):
                                                    if (
                                                        row[0] != outliers[x][0]
                                                        or row[0] != outliers[x][1]
                                                        or row[1] != outliers[x][0]
                                                        or row[1] != outliers[x][1]
                                                    ):
                                                        outliers.append(row)
                                                        print("Row added to outliers")

        currentDir = os.path.dirname(os.path.abspath(__file__))
        relativePath = os.path.join(
            currentDir, "..", "created_at_runtime", "fa_subnetwork_m1.txt"
        )
        print(results)

        with open(relativePath, "w") as outputFile:
            for row in results:
                outputFile.write("\t".join(row))
        outputFile.close()
