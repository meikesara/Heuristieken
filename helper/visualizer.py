"""
Script for visualizing folded proteins and change in stability over a
number of iterations.

Meike Kortleve, Nicole Jansen
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def plotProtein(protein):
    """
    Create a visual representation of a protein.

    Argument:
    protein -- object of class Protein
    """

    colorDict = {"P": "b", "H": "r", "C": "C6"}
    nameDict = {"P": "polar", "H": "hydrofobic", "C": "cysteine"}

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
        ax = plt.axes(projection="3d")

    # Plot lines between amino acids
    if protein.plane == "2D":
        ax.plot(xLines, yLines, "-k", zorder=-1)
    else:
        ax.plot(xLines, yLines, zLines, "-k", zorder=-1)

    # Plot dotted lines to indicate bonds (HH, HC, CC)
    for amino in protein.aminoList:
        coordinates = getStabilityCo(amino, protein)

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
            ax.scatter(xPoints[i], yPoints[i], color=colorDict[i],
                       label=nameDict[i])
        else:
            ax.scatter(xPoints[i], yPoints[i], zPoints[i], color=colorDict[i],
                       label=nameDict[i])

    plt.title("stability = " + str(protein.stability))
    ax.legend()
    if protein.plane == "2D":
        ax.axis('equal')
    else:
        # Create cubic grid
        minValue = min((xLines + yLines + zLines))
        maxValue = max((xLines + yLines + zLines))
        ax.set_xlim(minValue, maxValue)
        ax.set_ylim(minValue, maxValue)
        ax.set_zlim(minValue, maxValue)

    ax.axis("off")
    plt.show()


def plotStability(stabilityList):
    """
    Plot how stability changes over a number of iterations.

    Argument:
    stabilityList -- list with stability at each iteration
    """

    plt.plot(stabilityList, "k")
    plt.xlim(0, len(stabilityList))
    plt.ylim((min(stabilityList) - 5), 0)
    plt.xlabel("Iterations")
    plt.ylabel("Stability")
    plt.show()


def processData(protein):
    """
    Process data needed to create plot of protein.
    Return dictionaries of x, y, and z coordinates where types are classified,
    and lists of x, y, and z coordinates for creating lines between amino acids.

    Argument:
    protein -- object of class Protein
    """

    xPoints = {"P": [], "H": [], "C": []}
    yPoints = {"P": [], "H": [], "C": []}
    zPoints = {"P": [], "H": [], "C": []}
    xLines = []
    yLines = []
    zLines = []

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
    Determine if stability bonds are made with the given amino acid.
    Return list coordinates of aminos with which bonds are created.

    Arguments:
    amino -- object of class Amino
    protein -- object of class Protein
    """

    coordinates = []
    typeCo = amino.type

    # Bonds only created between amino acids of type H and C
    if typeCo in {"H", "C"}:
        aroundCos = protein.getSurroundCo(amino.coordinate, occupied=True)

        # For each amino next to given amino, check if they create a bond
        for aroundCo in aroundCos:
            idAround = protein.occupied.index(aroundCo)
            aroundType = protein.aminoList[idAround].type

            if aroundType in {"H", "C"}:
                # Check if given amino is not connected in protein to amino
                if idAround != (amino.id + 1) and idAround != (amino.id - 1):
                    # Bond only needs to be plotted once
                    if amino.id > idAround:
                        coordinates.append(aroundCo)

    return coordinates
