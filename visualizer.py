"""
Script for visualizing folded proteins and change in stability over a number of
iterations.

Meike Kortleve, Nicole Jansen
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plotProtein(protein):
    """
    Function for creating a visual representation of a folded protein.
    """

    # Dictionaries containing colors and labels for the types of amino acids
    colorDict = {"P": "b", "H": "r", "C": "C6"}
    nameDict = {"P": "polair", "H": "hydrofoob", "C": "cysteine"}

    if protein.plane == "2D":
        xPoints, yPoints, xLines, yLines = processData(protein)
    else:
        xPoints, yPoints, zPoints, xLines, yLines, zLines = processData(protein)

    # When protein only contains H's and P's, type C should not be in legend
    if not xPoints["C"]:
        del colorDict["C"]

    fig = plt.figure()
    if protein.plane == "2D":
        ax = plt.axes()
    elif protein.plane == "3D":
        ax = plt.axes(projection='3d')

    # Plot lines between amino acids
    if protein.plane == "2D":
        ax.plot(xLines, yLines, "-k", zorder=-1)
    else:
        ax.plot(xLines, yLines, zLines, "-k", zorder=-1)

    # Plot dotted lines to indicate bonds (HH, HC, CC)
    for amino in protein.aminoList:
        # Check if there are bonds for current amino acid
        coordinates = getStabilityCo(amino, protein)

        # Plot each dotted line seperately
        for coordinate in coordinates:
            if protein.plane == "2D":
                ax.plot([amino.coordinate[0], coordinate[0]],
                        [amino.coordinate[1], coordinate[1]], ":k", zorder=-1)
            else:
                ax.plot([amino.coordinate[0], coordinate[0]],
                        [amino.coordinate[1], coordinate[1]],
                        [amino.coordinate[2], coordinate[2]], ":k", zorder=-1)

    # Plot amino acids
    for i in colorDict:
        if protein.plane == "2D":
            ax.scatter(xPoints[i], yPoints[i], color=colorDict[i], label=nameDict[i])
        else:
            ax.scatter(xPoints[i], yPoints[i], zPoints[i], color=colorDict[i], label=nameDict[i])

    plt.title("stabiliteit = " + str(protein.stability))
    ax.legend()
    ax.axis('equal')
    # ax.set_aspect('equal')
    #ax.grid(True) # TODO: voor als wel grid willen, moet het op andere manier
    # ax.axis('off')
    plt.show()


def plotStability(stabilityList):
    """
    Function for plotting how stability changes over a number of iterations.
    """

    # Create plot
    # fig = plt.figure()
    # ax = plt.axes()
    plt.plot(stabilityList, "k")
    plt.xlim(0, len(stabilityList))
    plt.ylim((min(stabilityList) - 5), 0)
    plt.xlabel("Iteraties")
    plt.ylabel("Stabiliteit")
    plt.show()


def processData(protein):
    """
    Function for processing data needed to create plot of folded protein.
    """

    # Initialize dictionaries and lists for coordinates
    xPoints = {"P": [], "H": [], "C": []}
    yPoints = {"P": [], "H": [], "C": []}
    zPoints = {"P": [], "H": [], "C": []}
    xLines = []
    yLines = []
    zLines = []

    # Add coordinates of each amino acid to dictionaries and lists
    for amino in protein.aminoList:
        xLines.append(amino.coordinate[0])
        yLines.append(amino.coordinate[1])
        if protein.plane == "3D":
            zLines.append(amino.coordinate[2])

        xPoints[amino.type].append(amino.coordinate[0])
        yPoints[amino.type].append(amino.coordinate[1])
        if protein.plane == "3D":
            zPoints[amino.type].append(amino.coordinate[2])

    if protein.plane == "2D":
        return xPoints, yPoints, xLines, yLines
    else:
        return xPoints, yPoints, zPoints, xLines, yLines, zLines


def getStabilityCo(amino, protein):
    """
    Function for determining if bonds are made with the given amino acid.
    """

    coordinates = []

    # Get type of given amino acid
    typeCo = amino.type

    # Bonds only created between amino acids of type H and C
    if typeCo in {"H", "C"}:
        # Get surrounding coordinates of given amino acid that are occupied
        aroundCos = protein.getSurroundCo(amino.coordinate, True)

        # For each amino next to given amino, check if they create a bond
        for aroundCo in aroundCos:
            # Get id and type of amino acid next to given amino
            idAround = protein.occupied.index(aroundCo)
            aroundType = protein.aminoList[idAround].type

            if aroundType in {"H", "C"}:
                # Check if amino is not connected in protein to amino
                if idAround != (amino.id + 1) and idAround != (amino.id - 1):
                    # Bond only needs to be plotted once
                    if amino.id > idAround:
                        coordinates.append(aroundCo)

    return coordinates
