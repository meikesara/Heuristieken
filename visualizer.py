"""
# TODO: comment
"""

import matplotlib.pyplot as plt

# Global variables
colorDict = {"P": 'b', "H": 'r', "C": 'g'}
nameDict = {"P": "polair", "H": "hydrofoob", "C": "cysteine"}
xPoints = {"P": [], "H": [], "C": []}
yPoints = {"P": [], "H": [], "C": []}
xLines = []
yLines = []

def createPlot(protein):
    """
    # TODO: comment
    """

    preprocessData(protein)

    fig, ax = plt.subplots()

    ax.plot(xLines, yLines, "-k", zorder=-1)

    for i in colorDict:
        ax.scatter(xPoints[i], yPoints[i], color=colorDict[i], label=nameDict[i])

    plt.title("stability = " + str(protein.stability))

    ax.legend()
    ax.axis('equal')
    ax.grid(True)
    ax.axis('off')

    plt.show()


def preprocessData(protein):
    """
    # TODO: comment
    """

    for amino in protein.aminoList:
        xLines.append(amino.coordinate[0])
        yLines.append(amino.coordinate[1])

        xPoints[amino.type].append(amino.coordinate[0])
        yPoints[amino.type].append(amino.coordinate[1])
