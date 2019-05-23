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
    This method uses the method createRandom to fold protein until the protein
    is validly folded.
    Returns the folded protein.

    Arguments:
    proteinString -- this is a string containing the amino aminoacids
    plane -- is either "2D" or "3D"
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

        # Add amino acid to the aminoList
        protein.aminoList.append(Amino(id, protein.proteinString[id]))

        # Place the first and second amino-acid.
        if id in {0, 1}:
            thisCoordinate = [0, id]
            if protein.plane == "3D":
                thisCoordinate.append(0)
            protein.aminoList[id].coordinate = thisCoordinate
            protein.occupied.append(thisCoordinate)

        # The remaining amino acids are randomly placed
        else:

            # Get the coordinates of the previous amino-acid
            prevCo = protein.aminoList[(id - 1)].coordinate

            # Get the surrounding coordinates that are not occupied
            posCo = protein.getSurroundCo(prevCo, False)

            # If there are no surrounding coordinates available stop the folding
            if not posCo:
                protein.stability = 0
                return False

            # Randomly choose one of the possible coordinates
            coordinate = random.choice(posCo)

            # Place the amino-acid on that coordinate
            protein.aminoList[id].coordinate = coordinate

            # Add the coordinate to the list of occupied coordinates
            protein.occupied.append(coordinate)

            protein.stabilityUpdate(protein.aminoList[id])

    return True
