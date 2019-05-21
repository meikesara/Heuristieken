"""
# TODO: comment hier

Meike Kortleve, Nicole Jansen
"""

from amino import Amino
from protein import Protein
import random
import copy

def randomFold(protein, minStability):
    """
    Function for randomly creating new folds of protein.
    Because random, will not give same output each time.
    minStability is the stability it tries to reach
    """

    # Initialize stability
    stability = protein.stability

    # Counter for exiting while loop to avoid initinite while loop
    counter = 0

    # Loop while the stability is bigger than the minStability
    while stability > minStability:
        # Exit while loop if iterations exceeds one billion
        if counter > pow(10, 6):
            return protein

        # Initialize new protein
        newProtein = Protein(protein.proteinString, protein.plane)

        # Fold the protein again
        validFolding = newProtein.createAminoList()

        while not validFolding:
            newProtein = Protein(protein.proteinString, protein.plane)

            # Fold the protein again
            validFolding = newProtein.createAminoList()

        # If stability of newly folded protein is higher update stability and fold
        if newProtein.stability < stability:
            protein = newProtein
            stability = newProtein.stability
            print(stability)

        # Update counter
        counter += 1

    return protein
