class Statistical_Test:
    def __init__(self, population1, population2):
        self.population1 = population1
        self.population2 = population2

    # set up function for p_test
    def statistic(x, y, axis):
        return np.mean(x, axis=axis) - np.mean(y, axis=axis)

    def p_test(stage1SubnetworksFile, stage2SubnetworksFile):
        stage1SubnetworksEdgeCount = []
        stage2SubnetworksEdgeCount = []
        alpha = 0.05

        # find the mean edge count of both stage 1 subnetworks and stage 2 subnetworks
        with open(stage1SubnetworksFile, "r") as file:
            stage1Subnetworks = json.load(file)
        for bin in stage1Subnetworks:
            stage1SubnetworksEdgeCount.append(stage1Subnetworks[bin]["edgeCount"])
        stage1SubnetworksMean = sum(stage1SubnetworksEdgeCount) / 5000
        print(f"stage1Mean: {stage1SubnetworksMean}")

        with open(stage2SubnetworksFile, "r") as file:
            stage2Subnetworks = json.load(file)
        for bin in stage2Subnetworks:
            stage2SubnetworksEdgeCount.append(stage2Subnetworks[bin]["edgeCount"])
        stage2SubnetworksMean = sum(stage2SubnetworksEdgeCount) / 5000
        print(f"stage2Mean: {stage2SubnetworksMean}")

        # calculate difference in means
        observedStatistic = stage2SubnetworksMean - stage1SubnetworksMean

        numPermutations = 10000

        permStatistic = np.zeros(numPermutations)

        # perform the permutation test
        for i in range(numPermutations):
            # combine the data and shuffle it
            combinedSubnetworks = np.concatenate(
                (stage1SubnetworksEdgeCount, stage2SubnetworksEdgeCount)
            )
            np.random.shuffle(combinedSubnetworks)

            # calculate the test statistic for this permutation
            permGroup1 = combinedSubnetworks[: len(stage1SubnetworksEdgeCount)]
            permGroup2 = combinedSubnetworks[len(stage1SubnetworksEdgeCount) :]
            permStat = np.mean(permGroup2) - np.mean(permGroup1)

            # ctore the permuted statistic
            permStatistic[i] = permStat

        # calculate the p-value by comparing the observed statistic to the permuted statistics
        pVal = (np.abs(permStatistic) >= np.abs(observedStatistic)).mean()

        return pVal
