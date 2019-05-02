"""
# TODO: comment
"""

import matplotlib.pyplot as plt

# Global variables


def plotProtein(protein):
    """
    Function for creating a visual representation of a folded protein
    """

    colorDict = {"P": "b", "H": "r", "C": "C6"}
    nameDict = {"P": "polair", "H": "hydrofoob", "C": "cysteine"}

    xPoints, yPoints, xLines, yLines = processData(protein)
    # plt.cla()
    # plt.clf()

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


def plotStability(stabilityList):
    """
    Function for plotting how stability changes over time
    """

    fig = plt.figure()
    # ax = plt.axes()
    plt.plot(stabilityList)
    plt.xlim(0, len(stabilityList))
    plt.ylim(min(stabilityList)-5, 0)
    plt.xlabel("Iteraties")
    plt.ylabel("Stabiliteit");
    plt.show()


def processData(protein):
    """
    # TODO: comment
    """

    xPoints = {"P": [], "H": [], "C": []}
    yPoints = {"P": [], "H": [], "C": []}
    xLines = []
    yLines = []

    for amino in protein.aminoList:
        xLines.append(amino.coordinate[0])
        yLines.append(amino.coordinate[1])

        xPoints[amino.type].append(amino.coordinate[0])
        yPoints[amino.type].append(amino.coordinate[1])

    return xPoints, yPoints, xLines, yLines
