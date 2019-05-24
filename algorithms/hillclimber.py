"""
Script to run a hill climber to fold proteins.
A new protein is accepted if the stability is the same or better than the
previous protein.

Meike Kortleve, Nicole Jansen
"""

import copy
import random
from classes.amino import Amino
from classes.protein import Protein
import helper.visualizer as visualizer 


def hillClimber(protein, iterations, runnings=False):
    """
    Runs a hill climber.
    When more runnings will be performed, the protein with the best stability
    will be returned, otherwise the protein will be visualized.

    Arguments:
    protein -- object of class Protein
    iterations -- positive integer, the amount of iterations of the algorithm

    Keyword argument:
    runnings -- boolean, True if multiple runnings of the algorithm are run,
                False if algorithm is run once (default)
    """

    stabilityList = []

    for i in range(iterations):
        stabilityList.append(protein.stability)

        # Create a new protein from the current protein
        newProtein = protein.pullMove()

        if newProtein.stability <= protein.stability:
            protein = newProtein

    stabilityList.append(protein.stability)

    if runnings:
        return protein
    else:
        visualizer.plotProtein(protein)
        visualizer.plotStability(stabilityList)
