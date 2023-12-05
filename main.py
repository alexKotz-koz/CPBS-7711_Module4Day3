import numpy as np
import pandas as pd
import itertools, time, json, concurrent, os

# Fa Utilities: create_parent_network, extract_loci
from components.fa_utilities import Fa_Utilities
from components.create_random_fa_subnetworks_m2 import Create_Random_Fa_Subnetworks
from components.genetic_algorithm_m4 import Genetic_Algorithm
from components.null_case_genetic_algorithm_m4 import Null_Case_Genetic_Algorithm
from components.create_random_nonFa_subnetworks_m2 import (
    Create_Random_Non_Fa_Subnetworks,
)
from components.score_individual_subnet_m3 import ScoreIndividualSubnet


# Input: final fa population, non fa populations
# Output: the 10 files containing the top ten most dense subnetworks (did this becuase of my non fa population situation) with the density score and pvalue
def statistical_test(finalFaPopulationFile, nonFaPopulationsFile=None):
    finalFaPopulation = {}
    nonFaPopulations = []

    finalFaSubnetScores = []
    nullPopulationScores = []

    with open(finalFaPopulationFile, "r") as file:
        finalFaPopulation = json.load(file)

    with open(nonFaPopulationsFile, "r") as file:
        for line in file:
            nonFaPopulations.append(json.loads(line))

    sortedFinalFaPopulation = sorted(
        finalFaPopulation["finalSubnets"],
        key=lambda x: x["averageDensity"],
        reverse=True,
    )

    topTenFaSubnets = sortedFinalFaPopulation[:10]

    finalFaPopulationAverageDensityScore = finalFaPopulation["finalAverageDensity"]
    for item in finalFaPopulation["finalSubnets"]:
        subnetScore = item["averageDensity"]
        finalFaSubnetScores.append(subnetScore)

    for item in nonFaPopulations[0]:
        nullPopulationScores.append(item["averageDensity"])

    observed_statistic = sum(finalFaSubnetScores)

    nullPopulationScores = [
        sum(item["averageDensity"] for item in population)
        for population in nonFaPopulations
    ]

    p_value = sum(
        1
        for null_statistic in nullPopulationScores
        if null_statistic >= observed_statistic
    ) / len(nullPopulationScores)

    print(f"P-value of permutation test: {p_value}")

    for index, subnet in enumerate(topTenFaSubnets):
        with open(
            f"supplemental_documents_for_submission/Day3_Output_Network<{index+1}>_pval<{p_value}>.txt",
            "w",
        ) as file:
            for gene in subnet["subnet"]:
                file.write(gene + "\t")
            file.write(str(subnet["averageDensity"]) + "\t")
            file.write(str(p_value))


def create_stats_file():
    generations = {}
    with open("created_at_runtime/generations.json") as file:
        generations = json.load(file)

    with open(
        "supplemental_documents_for_submission/generation_statistics.txt", "w"
    ) as file:
        for item in generations:
            if item == "1":
                print(generations[item]["averageDensity"])
                file.write("First Generation: \n")
                file.write(
                    f"   Average Density: {generations[item]['averageDensity']}\n"
                )
            else:
                file.write(f"Generation: {item}: \n")
                file.write(
                    f"   Average Density: {generations[item]['averageDensity']}\n"
                )
                file.write(
                    f"   Generation completed in: {generations[item]['timeToCompletion']} seconds\n"
                )
                file.write(
                    f"   Density improved by: {generations[item]['percentImproved']}%\n\n"
                )


def create_day3_output_file():
    # create supplemental documents for submission, adding this section becuase there is not enough time to refactor the existing code to write the files in the required format
    locus0 = []
    locus1 = []
    locus2 = []
    locus3 = []
    locus4 = []
    locus5 = []
    locus6 = []
    locus7 = []
    locus8 = []
    locus9 = []
    locus10 = []
    locus11 = []

    with open("created_at_runtime/averageGeneScores.json", "r") as file:
        scoresFromFile = json.load(file)
        for item in scoresFromFile.items():
            score = item[1]["locusId"]
            if score == "0":
                locus0.append(item)
            elif score == "1":
                locus1.append(item)
            elif score == "2":
                locus2.append(item)
            elif score == "3":
                locus3.append(item)
            elif score == "4":
                locus4.append(item)
            elif score == "5":
                locus5.append(item)
            elif score == "6":
                locus6.append(item)
            elif score == "7":
                locus7.append(item)
            elif score == "8":
                locus8.append(item)
            elif score == "9":
                locus9.append(item)
            elif score == "10":
                locus10.append(item)
            elif score == "11":
                locus11.append(item)

    with open("supplemental_documents_for_submission/Day3_Output.gmt", "a") as file:
        file.write("locus0" + "\t")
        for item in locus0:
            file.write(item[0] + " " + str(item[1]["averageScore"]) + "\t")
        file.write("\n")
        file.write("locus1" + "\t")
        for item in locus1:
            file.write(item[0] + " " + str(item[1]["averageScore"]) + "\t")
        file.write("\n")
        file.write("locus2" + "\t")
        for item in locus2:
            file.write(item[0] + " " + str(item[1]["averageScore"]) + "\t")
        file.write("\n")
        file.write("locus3" + "\t")
        for item in locus3:
            file.write(item[0] + " " + str(item[1]["averageScore"]) + "\t")
        file.write("\n")
        file.write("locus4" + "\t")
        for item in locus4:
            file.write(item[0] + " " + str(item[1]["averageScore"]) + "\t")
        file.write("\n")
        file.write("locus5" + "\t")
        for item in locus5:
            file.write(item[0] + " " + str(item[1]["averageScore"]) + "\t")
        file.write("\n")
        file.write("locus6" + "\t")
        for item in locus6:
            file.write(item[0] + " " + str(item[1]["averageScore"]) + "\t")
        file.write("\n")
        file.write("locus7" + "\t")
        for item in locus7:
            file.write(item[0] + " " + str(item[1]["averageScore"]) + "\t")
        file.write("\n")
        file.write("locus8" + "\t")
        for item in locus8:
            file.write(item[0] + " " + str(item[1]["averageScore"]) + "\t")
        file.write("\n")
        file.write("locus9" + "\t")
        for item in locus9:
            file.write(item[0] + " " + str(item[1]["averageScore"]) + "\t")
        file.write("\n")
        file.write("locus10" + "\t")
        for item in locus10:
            file.write(item[0] + " " + str(item[1]["averageScore"]) + "\t")
        file.write("\n")
        file.write("locus11" + "\t")
        for item in locus11:
            file.write(item[0] + " " + str(item[1]["averageScore"]) + "\t")


def main():
    start = time.time()
    # FA UTILITIES
    faStart = time.time()
    faUtilitiesInstance = Fa_Utilities(
        inputFile="static/Input.gmt.txt",
        parentNetworkFile="static/STRING 1.txt",
        module1FaNetworkFile="static/module1_fa_network.txt",
    )
    parentFaNetworkDF = faUtilitiesInstance.create_parent_network()

    faLoci = faUtilitiesInstance.extract_loci()
    masterParentNetworkDF, geneDict = faUtilitiesInstance.create_master_parent_network(
        faLoci=faLoci
    )
    module1FaNetwork = faUtilitiesInstance.extract_module1_fa_network()
    faEnd = time.time()
    print(
        f"faUtilities: parentFaNetworkDF, masterParentNetworkDF, geneDict, faLoci, and module1FaNetwork(from faUtilities-extract_...) created in: {faEnd - faStart}\n"
    )

    # M2 - 5000 Random FA Subnetworks
    crstart = time.time()
    createRandomFaSubnetworksInstance = Create_Random_Fa_Subnetworks(
        "static/module1_fa_network.txt",
        "static/Input.gmt.txt",
        "static/STRING 1.txt",
        faLoci,
        module1FaNetwork,
        parentFaNetworkDF,
    )
    randomFaSubnetworks = createRandomFaSubnetworksInstance.create_random_subnetworks()
    crend = time.time()
    print(f"M2: 5000 random fa subnetworks created in {crend-crstart}\n")

    # M4 - Genetic Algorithm on FA population
    gastart = time.time()
    geneticAlgorithmInstance = Genetic_Algorithm(
        initialPopulation=randomFaSubnetworks,
        faLoci=faLoci,
        parentNetwork=parentFaNetworkDF,
    )
    finalPopulation = geneticAlgorithmInstance.start_genetic_algorithm()
    gaend = time.time()
    print(f"M4: Optimization completed in {gaend-gastart}\n")

    finalPopulationFromJSON = {}
    with open("created_at_runtime/finalFaPopulation.json", "r") as file:
        finalPopulationFromJSON = json.load(file)

    # NULL CASE: Not currently being used, if desired please review the null_case_genetic_algorithm class to view the current status of this attempt
    """nfastart = time.time()
    createRandomNonFaSubnetworksInstance = Create_Random_Non_Fa_Subnetworks(
        "static/STRING 1.txt",
        faLoci,
        finalPopulationSubnets=finalPopulationFromJSON["finalSubnets"],
        parentFaNetworkDF=parentFaNetworkDF,
        masterParentNetworkDF=masterParentNetworkDF,
    )"""
    # change for number of null case non fa subnetwork populations desired to run !!!Warning changing this number will drastically increase runtime!!!
    for _ in range(0, 8):
        createRandomNonFaSubnetworksInstance.create_non_fa_subnetworks()
    nfaend = time.time()
    print(f"Non Fa Subnetworks created in: {nfaend-nfastart}\n")

    nullCasePopulationFromJSON = {}
    with open("created_at_runtime/random_non_fa_subnetworks.json", "r") as file:
        nullCasePopulationFromJSON = json.load(file)
    nullCaseBinsFromJSON = {}
    with open("created_at_runtime/bins.json", "r") as file:
        nullCaseBinsFromJSON = json.load(file)

    # NULL CASE GA
    nullCaseGeneticAlgorithmInstance = Null_Case_Genetic_Algorithm(
        initialPopulation=nullCasePopulationFromJSON,
        bins=nullCaseBinsFromJSON,
        faLoci=faLoci,
        masterParentNetworkDF=masterParentNetworkDF,
        geneDict=geneDict,
    )
    nullCaseGeneticAlgorithmInstance.start_genetic_algorithm()

    # M3 - Score Genes from Final Population

    # using a fraction of the final Subnets for gene scoring to save runtime
    finalPopulationFromJSON = {}
    with open("created_at_runtime/finalFaPopulation.json", "r") as file:
        finalPopulationFromJSON = json.load(file)

    testData = list(itertools.islice(finalPopulationFromJSON["finalSubnets"], 1000))

    scoreStart = time.time()
    print("Starting Gene Scoring")
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = []
        for item in testData:
            subnet = item["subnet"]
            print(f"subnet:{subnet}")

            scoreIndividualSubnetInstance = ScoreIndividualSubnet(
                individualSubnet=subnet,
                inputFile="Input.gmt.txt",
                parentNetwork=parentFaNetworkDF,
                loci=faLoci,
            )
            futures.append(executor.submit(scoreIndividualSubnetInstance.gene_score))

        for future in concurrent.futures.as_completed(futures):
            geneScores = future.result()

    averageGeneScores = faUtilitiesInstance.calculate_average_gene_score(
        faNetworkFile="created_at_runtime/faNetwork.txt",
        geneScoresFile="created_at_runtime/scoring_geneScores.txt",
    )
    scoreEnd = time.time()
    print(f"M3: Score genes from final population completed in {scoreEnd-scoreStart}")
    currentDir = os.path.dirname(os.path.abspath(__file__))
    relativePath = os.path.join(
        currentDir, "created_at_runtime", "averageGeneScores.json"
    )

    with open(relativePath, "w") as outputFile:
        json.dump(averageGeneScores, outputFile)

    end = time.time()
    print(f"Total Runtime: {end-start}")

    # FUNCTIONS BELOW ARE TO GENERATE REQUIRED OUTPUT FILES, DO NOT UNCOMMENT:

    # create_day3_output_file()
    # create_stats_file()
    statistical_test(
        finalFaPopulationFile="created_at_runtime/finalFaPopulation.json",
        nonFaPopulationsFile="created_at_runtime/random_non_fa_subnetworks.json",
    )


if __name__ == "__main__":
    main()
