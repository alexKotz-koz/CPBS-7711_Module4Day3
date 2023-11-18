import numpy as np
import pandas as pd
import itertools

# Fa Utilities: create_parent_network, extract_loci
from components.fa_utilities import FaUtilities
from components.create_random_fa_subnetworks import Create_Random_Fa_Subnetworks


def main():
    faUtilitiesInstance = FaUtilities(
        inputFile="static/Input.gmt.txt",
        parentNetworkFile="static/STRING 1.txt",
        module1FaNetworkFile="static/module1_fa_network.txt",
    )
    parentNetworkDF = faUtilitiesInstance.create_parent_network()
    faLoci = faUtilitiesInstance.extract_loci()
    module1FaNetwork = faUtilitiesInstance.extract_module1_fa_network()

    createRandomFaSubnetworksInstance = Create_Random_Fa_Subnetworks(
        "static/module1_fa_network.txt",
        "static/Input.gmt.txt",
        "static/STRING 1.txt",
        faLoci,
        module1FaNetwork,
    )
    randomFaSubnetworks = createRandomFaSubnetworksInstance.create_random_subnetworks()

    # Test Dataset: cut out the first 1000 randomly generated fa subnetworks to reduce runtime...
    # If desired, replace the 1000 on line 165 with the number of randomly generated fa subnetworks you wish to test ***MUST BE LESS THAN OR EQUAL TO 5000
    # THIS LINE TRUNCATES THE 5000 SUBNETWORKS TO REDUCE RUNTIME
    # TO TEST FULL FUNCTIONALITY: COMMENT THIS LINE AND REPLACE testData.items() WITH stage1Subnetworks.items() on line 170
    testData = dict(itertools.islice(randomFaSubnetworks.items(), 1000))


if __name__ == "__main__":
    main()
