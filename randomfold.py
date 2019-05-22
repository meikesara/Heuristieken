"""
This script contains the algorithm to randomly fold a protein.

Meike Kortleve, Nicole Jansen
"""

from amino import Amino
from protein import Protein
import random
import copy

def randomFold(proteinString, plane):
    """
    Random folding of the protein. The first and second amino acids are
    always placed at coordinates (0, 0) and (0, 1) (or for 3D (0, 0, 0) and
    (0, 1, 0)) respectively. This excludes rotationally symmetric folds.
    Returns the randomly folded protein.

    Arguments:
    proteinString -- this is a string containing the amino aminoacids
    plane -- is either "2D" or "3D"
    """

    # Initialize new protein
    protein = Protein(proteinString, plane)

    validFolding = False

    while not validFolding:

        protein.aminoList = []
        protein.occupied = []

        for id in range(protein.proteinLength):
            protein.aminoList.append(Amino(id, protein.proteinString[id]))

            if id in {0, 1}:
                thisCoordinate = [0, id]
                if protein.plane == "3D":
                    thisCoordinate.append(0)
                protein.aminoList[id].addCoordinate(thisCoordinate)
                protein.occupied.append(thisCoordinate)

            # The remaining amino-acids are randomly placed
            else:
                # Get the coordinates of the previous amino-acid
                prevCo = protein.aminoList[(id - 1)].coordinate

                # Get the surrounding coordinates that are not occupied
                posCo = protein.getSurroundCo(prevCo, False)

                # If there are no surrounding coordinates available break from the loop
                if not posCo:
                    protein.stability = 0
                    break

                # Randomly choose one of the possible coordinates
                coordinate = random.choice(posCo)

                # Place the amino-acid on that coordinate
                protein.aminoList[id].addCoordinate(coordinate)

                # Add the coordinate to the list of occupied coordinates
                protein.occupied.append(coordinate)

                # Update the stability
                protein.stabilityUpdate(protein.aminoList[id])

        validFolding = True

    return protein
