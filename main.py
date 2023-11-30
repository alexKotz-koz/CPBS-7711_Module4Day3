import numpy as np
import pandas as pd
import itertools
import time

# Fa Utilities: create_parent_network, extract_loci
from components.fa_utilities import Fa_Utilities
from components.create_random_fa_subnetworks_m2 import Create_Random_Fa_Subnetworks
from components.create_fa_subnetwork_m1 import Create_Fa_Subnetwork
from components.genetic_algorithm_m4 import Genetic_Algorithm
from components.create_random_nonFa_subnetworks_m2 import (
    Create_Random_Non_Fa_Subnetworks,
)


def main():
    start = time.time()

    # FA UTILITIES
    faStart = time.time()
    faUtilitiesInstance = Fa_Utilities(
        inputFile="static/Input.gmt.txt",
        parentNetworkFile="static/STRING 1.txt",
        module1FaNetworkFile="static/module1_fa_network.txt",
    )
    parentNetworkDF = faUtilitiesInstance.create_parent_network()
    faLoci = faUtilitiesInstance.extract_loci()
    module1FaNetwork = faUtilitiesInstance.extract_module1_fa_network()
    faEnd = time.time()
    print(
        f"faUtilities: parentNetworkDF, faLoci, and module1FaNetwork(from faUtilities-extract_...) created in: {faEnd - faStart}\n"
    )

    # M1 - FA Subnetwork
    """cstart = time.time()
    createFaSubnetworkM1Instance = Create_Fa_Subnetwork(faLoci, "static/STRING 1.txt")
    createFaSubnetworkM1Instance.create_subnetwork()
    cend = time.time()
    print(f"createFaSubnetworkM1: module 1 subnetwork created in: {cend-cstart}")"""

    # M2 - 5000 Random FA Subnetworks
    crstart = time.time()
    createRandomFaSubnetworksInstance = Create_Random_Fa_Subnetworks(
        "static/module1_fa_network.txt",
        "static/Input.gmt.txt",
        "static/STRING 1.txt",
        faLoci,
        module1FaNetwork,
        parentNetworkDF,
    )
    randomFaSubnetworks = createRandomFaSubnetworksInstance.create_random_subnetworks()
    crend = time.time()
    print(f"5000 random fa subnetworks created in: {crend-crstart}\n")

    gastart = time.time()
    geneticAlgorithmInstance = Genetic_Algorithm(
        initialPopulation=randomFaSubnetworks,
        faLoci=faLoci,
        parentNetwork=parentNetworkDF,
    )
    geneticAlgorithmInstance.start_genetic_algorithm()
    gaend = time.time()
    print(f"Optimization completed in:{gaend-gastart}\n")

    # Create 1000 final (optimized) non fa subnetwork populations to simulate the null case
    nfastart = time.time()
    createRandomNonFaSubnetworksInstance = Create_Random_Non_Fa_Subnetworks(
        "static/STRING 1.txt", faLoci
    )
    createRandomNonFaSubnetworksInstance.create_non_fa_subnetworks()
    nfaend = time.time()
    print(f"Non Fa Subnetworks created in: {nfaend-nfastart}\n")

    # Test Dataset: cut out the first 1000 randomly generated fa subnetworks to reduce runtime...
    # If desired, replace the 1000 on line 165 with the number of randomly generated fa subnetworks you wish to test ***MUST BE LESS THAN OR EQUAL TO 5000
    # THIS LINE TRUNCATES THE 5000 SUBNETWORKS TO REDUCE RUNTIME
    # TO TEST FULL FUNCTIONALITY: COMMENT THIS LINE AND REPLACE testData.items() WITH stage1Subnetworks.items() on line 170
    # testData = dict(itertools.islice(randomFaSubnetworks.items(), 1000))

    end = time.time()
    print(f"Total Runtime: {end-start}")


if __name__ == "__main__":
    main()
