"""
This code contains a depth first search method to find the optimal stability
for a protein. It checks all possible ways to fold the protein, but the
mirror images and the rotationally symmetric folds.

Meike Kortleve, Nicole Jansen
"""

from protein import Protein
from amino import Amino
import copy
import visualizer
import matplotlib.pyplot as plt
import visualizer
import time
import re

optStability=[0]

def constructive(proteinString):
    """
    This function creates the protein, places the first two amino acids and starts
    the placing of the rest of the amino acids.

    The input (proteinString) of this function is a protein string consisting of the letters
    H, P and C.
    """

    # Initialize protein
    protein = Protein(proteinString, "2D")

    # Place first and second amino acid to prevent rotational symmetry
    for i in range(2):
        protein.aminoList.append(Amino(i, protein.proteinString[i]))
        protein.aminoList[i].addCoordinate([0, i])
        protein.occupied.append([0, i])

    # This starts the recursive function createFolded
    createFolded(protein, 2)

    # Print the results
    print("Protein", proteinString, "has an optimal stability of", optStability[0])


def createFolded(protein, idToMove):
    """
    This function recusrively places amino aminoacids. To create all
    possible ways to fold a protein.

    The stability of the full protein only needs to be checked if the protein
    folds right/left twice.

    This first argument (protein) is an object of the protein class.
    The second argument (idToMove) is the id of the aminoAcid that needs to be moved.
    """

    # Stop this function if idToMove exceeds the length of the protein
    if idToMove > (protein.proteinLength - 1):
        return

    # Delete the amino acids from the current amino acid
    del protein.aminoList[idToMove :]
    del protein.occupied[idToMove :]

    # Get the coordinate of the previous amino acid
    prevCo = protein.aminoList[(idToMove - 1)].coordinate

    # Get the non-occupied surrounding amino acids of the previous coordinates.
    possibleCos = protein.getSurroundCo(prevCo, False)

    # Place current amino acid in the aminoList
    protein.aminoList.append(Amino(idToMove, protein.proteinString[idToMove]))

    # Calculate the sum of all the x-coordinates
    xTotal = sum([xCo[0] for xCo in protein.occupied])

    if xTotal == 0 and len(possibleCos) == 3:

        # Remove the last coordinate to prevent mirror images.
        possibleCos.pop(1)

    for possibleCo in possibleCos:

        # Place the current amino acid on possibleCos
        protein.aminoList[idToMove].addCoordinate(possibleCo)

        # Add the current coordinate to the list of occupied coordinates.
        try:
            protein.occupied[idToMove] = possibleCo
        except:
            protein.occupied.append(possibleCo)


        # Calculate the sum of the x and y coordinates.
        xTotal = abs(sum([xCo[0] for xCo in protein.occupied]))
        yTotal = abs(sum([yCo[1] for yCo in protein.occupied]))

        # This total is two times the protein length if it is a straight line,
        # only folds once either to the left or right or only folds to the
        # left and right.
        total =  xTotal + yTotal

        # Check if all the amino acids have been placed
        if (idToMove == (protein.proteinLength - 1) and
            total != protein.proteinLength * 2):

            getStability(protein)

            # Update the optimalStability
            if protein.stability < optStability[0]:
                optStability[0] = protein.stability

        # Place the next amino acid
        createFolded(protein, (idToMove + 1))


def getStability(protein):
    """
    This function calculates the stability of a protein.

    This function needs an object from the protein class.
    """

    # Reset stability voor nieuwe vouwing
    protein.stability = 0

    for amino in protein.aminoList:

        typeCo = amino.type
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
                    id = amino.id
                    if idAround != (id + 1) and idAround != (id - 1):

                        # Bond only needs to be plotted once
                        if id > idAround:

                            # Stronger bond created when both aminos are type C
                            if typeCo == "C" and aroundType == "C":
                                protein.stability -= 5

                            # Weaker bond created when at least one amino is type H
                            else:
                                protein.stability -= 1
