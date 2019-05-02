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

    if not xPoints["C"]:
        del colorDict["C"]

    fig, ax = plt.subplots()

    ax.plot(xLines, yLines, "-k", zorder=-1)

    for amino in protein.aminoList:
        coordinates = getStabilityCo(amino, protein)
        for coordinate in coordinates:
            ax.plot([amino.coordinate[0], coordinate[0]],
                    [amino.coordinate[1], coordinate[1]], ":k", zorder=-1)

    for i in colorDict:
        ax.scatter(xPoints[i], yPoints[i], color=colorDict[i], label=nameDict[i])

    plt.title("stabiliteit = " + str(protein.stability))

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
    plt.ylim(min(stabilityList)-2, 0)
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


def getStabilityCo(amino, protein):
    """
    Function for determining if bonds are made with the amino acid
    """

    coordinates = []

    # Get type of current amino acid # TODO: bepaal of deze comment nodig is
    typeCo = amino.type

    if typeCo == "H" or typeCo == "C":
        # Get surrounding coordinates of given amino acid that are occupied
        aroundCos = protein.getSurroundCo(amino.coordinate, True)

        # Loop over occupied coordinates around current amino-acid # TODO: comment nodig?
        for aroundCo in aroundCos:
            # TODO: comment(?)
            idAround = protein.occupied.index(aroundCo)

            # Get the type of the amino-acid # TODO: comment nodig(?)
            nextCo = protein.aminoList[idAround].type

            if nextCo == "H" or nextCo == "C":
                # Check if amino is not connected in protein to amino
                if idAround != (amino.id + 1) and idAround != (amino.id - 1):
                    # Only add coordinates to list if aroundAmino is earlier in protein
                    if amino.id > idAround:
                        coordinates.append(aroundCo)

    return coordinates
