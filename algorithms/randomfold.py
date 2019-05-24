"""
This script contains the algorithm to randomly fold a protein.

Meike Kortleve, Nicole Jansen
"""

import copy
import random
from classes.amino import Amino
from classes.protein import Protein


def randomFold(proteinString, plane):
    """
    This function uses the fucntion createRandom to fold protein until the
    protein is validly folded.
    Returns the folded protein.

    Arguments:
    proteinString -- this is a string containing the amino acids
    plane -- dimension in which protein is placed/ folded, either "2D" or "3D"
    """

    protein = Protein(proteinString, plane)

    # Fold the protein
    validFolding = createRandom(protein)

    while not validFolding:
        # Fold the protein again
        validFolding = createRandom(protein)

    return protein


def createRandom(protein):
    """
    Random folding of the protein. The first and second amino acids are
    always placed at coordinates (0, 0) and (0, 1) (or for 3D (0, 0, 0) and
    (0, 1, 0)) respectively. This excludes rotationally symmetric folds.
    Returns True if the protein is validly folded and False if the protein is
    not validly folded.

    Argument:
    protein -- an object of the Protein class
    """

    protein.occupied = []
    protein.aminoList = []

    for id in range(protein.proteinLength):
        protein.aminoList.append(Amino(id, protein.proteinString[id]))

        # Place the first and second amino acid
        if id in {0, 1}:
            thisCoordinate = [0, id]
            if protein.plane == "3D":
                thisCoordinate.append(0)
            protein.aminoList[id].coordinate = thisCoordinate
            protein.occupied.append(thisCoordinate)
        else:
            prevCo = protein.aminoList[(id - 1)].coordinate
            posCo = protein.getSurroundCo(prevCo, occupied=False)

            # If there are no surrounding coordinates available stop the folding
            if not posCo:
                protein.stability = 0
                return False

            coordinate = random.choice(posCo)
            protein.aminoList[id].coordinate = coordinate
            protein.occupied.append(coordinate)

            protein.stabilityUpdate(protein.aminoList[id])

    return True
