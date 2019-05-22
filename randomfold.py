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

    Arguments:
    proteinString -- this is a string containing the amino aminoacids
    plane -- is either "2D" or "3D"
    """

    # Initialize new protein
    protein = Protein(proteinString, plane)

    # Fold the protein
    validFolding = protein.createAminoList()

    while not validFolding:

        # Fold the protein again
        validFolding = protein.createAminoList()

    return protein.stability
