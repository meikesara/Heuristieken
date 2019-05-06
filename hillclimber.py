"""
# TODO: comment hier

Meike Kortleve, Nicole Jansen
"""

from amino import Amino
from protein import Protein
import random
import copy


def hillClimber(protein, iterations, stabilityChange=False):
    """
    Performs a hill climber
    protein is the protein that will be changed
    iterations is the amount of iterations the hill climber should be performed
    """

    if stabilityChange:
        stabilityList = []

    for i in range(iterations):
        if stabilityChange:
            stabilityList.append(protein.stability)
        newProtein = protein.pullMove()
        if newProtein.stability <= protein.stability:
            protein = newProtein

        if stabilityChange and (i == (iterations - 1)):
            stabilityList.append(protein.stability)

    if stabilityChange:
        return protein, stabilityList
    else:
        return protein
