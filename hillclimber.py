"""
Script to run a hill climber to fold proteins.
A new protein is accepted if the stability is the same or better than the
previous protein.

Meike Kortleve, Nicole Jansen
"""

from amino import Amino
from protein import Protein
import random
import copy
import visualizer


def hillClimber(protein, iterations, runnings = False):
    """
    Run a hill climber.

    Arguments:
    protein -- object of class Protein.
    iterations -- positive integer, the amount of iterations of the algorithm.
    runnings -- boolean, True if multiple runnings of the algorithm are run.
                False if algorithm is run once (default).
    """

    stabilityList = []

    for i in range(iterations):

        # Add current stability to stabilityList
        stabilityList.append(protein.stability)

        # Create a new protein from the current protein
        newProtein = protein.pullMove()

        # If protein if stability is the same or better
        if newProtein.stability <= protein.stability:
            protein = newProtein

    stabilityList.append(protein.stability)

    if runnings:
        return protein
    else:
        visualizer.plotProtein(protein)
        visualizer.plotStability(stabilityList)
