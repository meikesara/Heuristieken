"""
This code contains a depth first search method to find the optimal stability
for a protein. It checks all possible ways to fold the protein, except for the
mirror images and the rotationally symmetric folds.

Meike Kortleve, Nicole Jansen
"""

import copy
import re
import time
import matplotlib.pyplot as plt
from amino import Amino
from protein import Protein
import visualizer

optStability = [0]

def constructive(proteinString):
    """
    This function creates the protein, places the first two amino acids and
    starts the placing of the rest of the amino acids.

    Argument:
    proteinString -- a string that contains the amino acids of the protein
    """

    protein = Protein(proteinString, "2D")

    # Place first and second amino acid to prevent rotational symmetry
    for i in range(2):
        protein.aminoList.append(Amino(i, protein.proteinString[i]))
        protein.aminoList[i].coordinate = [0, i]
        protein.occupied.append([0, i])

    # This starts the recursive function createFolded
    createFolded(protein, 2)

    print("Protein", proteinString, "has an optimal stability of",
          optStability[0])


def createFolded(protein, idToMove):
    """
    This function recursively places amino acids. To create all
    possible ways to fold a protein.

    Mirror images are prevented by using the fact that when the sum of the
    x-coordinates is zero and the next amino acid has three possible options,
    two of those options will be mirror images, so one of those will be removed.

    If the protein is not folded at least twice consecutively in the same
    direction, there will be no amino acids laying next to each other. This will
    result in the sum of the x- and y-coordinates to be two times the protein
    length - and therefore the stability will not to be checked.

    Arguments:
    protein -- object of class Protein
    idToMove -- positive integer, id of the amino acid that will be moved
    """

    # Stop this function if idToMove exceeds the length of the protein
    if idToMove > (protein.proteinLength - 1):
        return

    del protein.aminoList[idToMove:]
    del protein.occupied[idToMove:]

    # Get the coordinate of the previous amino acid
    prevCo = protein.aminoList[(idToMove - 1)].coordinate

    # Get the unoccupied surrounding amino acids of the previous coordinates
    possibleCos = protein.getSurroundCo(prevCo, occupied=False)

    protein.aminoList.append(Amino(idToMove, protein.proteinString[idToMove]))

    xTotal = sum([xCo[0] for xCo in protein.occupied])

    if xTotal == 0 and len(possibleCos) == 3:
        # Remove the last coordinate to prevent mirror images
        possibleCos.pop(1)

    for possibleCo in possibleCos:
        protein.aminoList[idToMove].coordinate = possibleCo

        try:
            protein.occupied[idToMove] = possibleCo
        except:
            protein.occupied.append(possibleCo)

        xTotal = abs(sum([xCo[0] for xCo in protein.occupied]))
        yTotal = abs(sum([yCo[1] for yCo in protein.occupied]))
        total =  xTotal + yTotal

        # Check if all the amino acids have been placed
        if (idToMove == (protein.proteinLength - 1) and
            total != protein.proteinLength * 2):
            getStability(protein)

            if protein.stability < optStability[0]:
                optStability[0] = protein.stability

        # Place the next amino acid
        createFolded(protein, (idToMove + 1))


def getStability(protein):
    """
    This function calculates the stability of a protein.

    Argument:
    protein -- object of class Protein
    """

    # Reset stability for new folding
    protein.stability = 0

    for amino in protein.aminoList:

        typeCo = amino.type
        if typeCo in {"H", "C"}:
            aroundCos = protein.getSurroundCo(amino.coordinate, occupied=True)

            for aroundCo in aroundCos:
                idAround = protein.occupied.index(aroundCo)
                aroundType = protein.aminoList[idAround].type

                if aroundType in {"H", "C"}:
                    id = amino.id

                    # Check if amino is not connected in protein to given amino
                    if idAround not in {(id + 1), (id - 1)}:
                        if id > idAround:
                            if typeCo == "C" and aroundType == "C":
                                protein.stability -= 5
                            else:
                                protein.stability -= 1
