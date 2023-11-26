class Visualize_Gene_Scores:
    def __init__(averageGeneScores, faNetwork):
        self.averageGeneScores = averageGeneScores

    # Input: averageGeneScores dictionary, fa network text file


# Output: creation of nextworkx network graph
def visualize_gene_scores(self):
    faNetwork = []
    averageGeneScores = self.averageGeneScores
    faNetwork = self.faNetwork

    ###REFACTOR: replace read code block logic to include the faNetwork object
    # get edges from filtered parent network
    with open(faNetworkFile, "r") as file:
        for row in file:
            row = row.strip().split("\t")[:2]
            faNetwork.append(row)

    G = nx.Graph()
    # add nodes to graph from averageGeneScores and add score and locusId as attributes
    for gene, data in averageGeneScores.items():
        for row in faNetwork:
            if gene in row and data["averageScore"] != "NA":
                G.add_node(
                    gene, averageScore=data["averageScore"], locusId=data["locusId"]
                )

    # add edges to graph object from the filtered parent network object
    nodesList = list(G.nodes)
    for edge in faNetwork:
        edgeOne = edge[0]
        edgeTwo = edge[1]
        if edgeOne in nodesList and edgeTwo in nodesList:
            G.add_edge(edgeOne, edgeTwo)

    # create community dictionary to group nodes by locus id
    commDict = {}
    for index, node in enumerate(list(G.nodes)):
        commDict[node] = G.nodes[node]["locusId"]

    color_map_dict = {
        "0": "blue",
        "1": "green",
        "2": "red",
        "3": "cyan",
        "4": "magenta",
        "5": "yellow",
        "6": "black",
        "7": "pink",
        "8": "brown",
        "9": "orange",
        "10": "purple",
        "11": "grey",
        "12": "olive",
    }

    pos = nx.circular_layout(G)
    # create node size and color maps for visualization
    nodeSize = {node: (G.nodes[node]["averageScore"] * 10) for node in list(G.nodes)}
    nodeColor = {node: color_map_dict[locusId] for node, locusId in commDict.items()}
    nodeLabels = {node: f"{node}" for node in G.nodes}
    fig, ax = plt.subplots(figsize=(10, 8))

    # trigger visualization
    Graph(
        G,
        node_color=nodeColor,
        node_size=nodeSize,
        node_edge_width=0.2,
        edge_width=0.1,
        edge_alpha=0.5,
        node_layout=pos,
        node_layout_kwargs=dict(node_to_community=commDict),
        node_alpha=0.75,
        node_labels=nodeLabels,
        node_label_fontdict={"size": 5, "color": "black", "weight": "bold"},
        node_label_offset=(0.05, 0.05),
        edge_layout="bundled",
        ax=ax,
    )

    plt.show()
