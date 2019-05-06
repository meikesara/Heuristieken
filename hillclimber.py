"""
# TODO: comment hier

Meike Kortleve, Nicole Jansen
"""

from amino import Amino
from protein import Protein
import random
import copy

def hillClimber(protein, iterations):
    """
    Performs a hill climber
    protein is the protein that will be changed
    iterations is the amount of iterations the hill climber should be performed
    """

    # NOTE: ben niet helemaal zeker of het < of == moet zijn
    if iterations == 0:
        print("hier", protein)
        a = 1
        return a
        #protein
    else:
        newProtein = protein.pullMove()
        print(protein)
        print(newProtein)
        if newProtein.stability <= protein.stability:
            print("kwam hier!")
            protein = newProtein
        print(protein)
        print()
        hillClimber(protein, (iterations - 1))
