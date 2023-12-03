import numpy as np
import pandas as pd
import itertools, time, json, concurrent, os

# Fa Utilities: create_parent_network, extract_loci
from components.fa_utilities import Fa_Utilities
from components.create_random_fa_subnetworks_m2 import Create_Random_Fa_Subnetworks
from components.create_fa_subnetwork_m1 import Create_Fa_Subnetwork
from components.genetic_algorithm_m4 import Genetic_Algorithm
from components.null_case_genetic_algorithm_m4 import Null_Case_Genetic_Algorithm
from components.create_random_nonFa_subnetworks_m2 import (
    Create_Random_Non_Fa_Subnetworks,
)
from components.score_individual_subnet_m3 import ScoreIndividualSubnet


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
    masterParentNetworkDF = faUtilitiesInstance.create_master_parent_network(
        faLoci=faLoci
    )
    module1FaNetwork = faUtilitiesInstance.extract_module1_fa_network()
    faEnd = time.time()
    print(
        f"faUtilities: parentFaNetworkDF, faLoci, and module1FaNetwork(from faUtilities-extract_...) created in: {faEnd - faStart}\n"
    )

    # M1 - FA Subnetwork
    """cstart = time.time()
    createFaSubnetworkM1Instance = Create_Fa_Subnetwork(faLoci, "static/STRING 1.txt")
    createFaSubnetworkM1Instance.create_subnetwork()
    cend = time.time()
    print(f"createFaSubnetworkM1: module 1 subnetwork created in: {cend-cstart}")"""

    # M2 - 5000 Random FA Subnetworks
    """crstart = time.time()
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
    print(f"M4: Optimization completed in {gaend-gastart}\n")"""

    """finalPopulationFromJSON = {}
    with open("created_at_runtime/finalFaPopulation.json", "r") as file:
        finalPopulationFromJSON = json.load(file)

    # NULL CASE
    nfastart = time.time()
    createRandomNonFaSubnetworksInstance = Create_Random_Non_Fa_Subnetworks(
        "static/STRING 1.txt",
        faLoci,
        finalPopulationSubnets=finalPopulationFromJSON["finalSubnets"],
        parentFaNetworkDF=parentFaNetworkDF,
        masterParentNetworkDF=masterParentNetworkDF,
    )
    createRandomNonFaSubnetworksInstance.create_non_fa_subnetworks()
    nfaend = time.time()
    print(f"Non Fa Subnetworks created in: {nfaend-nfastart}\n")"""

    """nullCasePopulationFromJSON = {}
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
    )
    nullCaseGeneticAlgorithmInstance.start_genetic_algorithm()"""

    # M3 - Score Genes from Final Population
    scoreIndividualSubnetInstance = ScoreIndividualSubnet(
        loci=faLoci, parentNetwork=parentFaNetworkDF
    )

    finalPopulationFromJSON = {}
    with open("created_at_runtime/finalFaPopulation.json", "r") as file:
        finalPopulationFromJSON = json.load(file)

    testData = list(itertools.islice(finalPopulationFromJSON["finalSubnets"], 3))

    scoreStart = time.time()
    # Process Pool to create gene scores for testData or stage1Subnetwork data
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = []
        for subnet in testData:
            futures.append(
                executor.submit(scoreIndividualSubnetInstance.gene_score, subnet)
            )

        # as each process completes store the results in geneScores
        for future in concurrent.futures.as_completed(futures):
            geneScores = future.result()

    # Calculate average gene scores from final population
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

    # Test Dataset: cut out the first 1000 randomly generated fa subnetworks to reduce runtime...
    # If desired, replace the 1000 on line 165 with the number of randomly generated fa subnetworks you wish to test ***MUST BE LESS THAN OR EQUAL TO 5000
    # THIS LINE TRUNCATES THE 5000 SUBNETWORKS TO REDUCE RUNTIME
    # TO TEST FULL FUNCTIONALITY: COMMENT THIS LINE AND REPLACE testData.items() WITH stage1Subnetworks.items() on line 170
    # testData = dict(itertools.islice(randomFaSubnetworks.items(), 1000))

    end = time.time()
    print(f"Total Runtime: {end-start}")


if __name__ == "__main__":
    main()
