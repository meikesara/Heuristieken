"""
Script to run a simulated annealing algorithm with a logarithmic cooling
rate (temperature = D/ln(iteration + 1) - D/ln(1000000)).

Meike Kortleve, Nicole Jansen
"""

import copy
import math
import random
import matplotlib.pyplot as plt
from classes.amino import Amino
from classes.protein import Protein
import helper.visualizer as visualizer


def simulatedAnnealing(protein, D, iterations, runnings=False):
    """
    Run simulated annealing.
    Returns the protein the with best stability and also the stability if only
    one running is performed.

    Arguments:
    proteinString -- a string that contains the amino acids of the protein
    D -- positive integer, determines the cooling rate
    iterations -- positive integer, the amount of iterations of the algorithm
    runnings -- boolean, True if multiple runnings of the algorithm are run,
                False if algorithm is run once (default)
    """

    stabilityList = []
    tempList = []
    temperature = 100

    bestProtein = protein

    for k in range(iterations):
        stabilityList.append(protein.stability)

        # Create a new protein by moving the amino acids
        newProtein = protein.pullMove()

        # Update the temperature
        temperature = D/math.log(k + 2) - D/math.log(1000000)
        # temperature = D/math.log(k + 2)

        tempList.append(temperature)

        # Set acceptance to 1 if the stability of the new protein is better.
        if newProtein.stability < protein.stability:
            acceptance = 1
        else:
            acceptance = math.exp((protein.stability - newProtein.stability)
                                   / temperature)

        # Accept newProtein if the acceptance exceeds random number from 0 to 1.
        if (acceptance > random.random()):
            protein = newProtein

        # If the stability of newProtein is the best set as bestProtein
        if newProtein.stability < bestProtein.stability:
            bestProtein = newProtein

    if runnings:
        return bestProtein
        
    else:
        # Add final stability to stabilityList
        stabilityList.append(protein.stability)

        print("Final solution stability: ", protein.stability)
        print("Best stability: ", bestProtein.stability)

        visualizer.plotProtein(bestProtein)
        visualizer.plotStability(stabilityList)

        # Plot temperatuur change
        plt.plot(tempList)
        plt.title("Temperature")


        # Plot acceptance change
        plt.figure()
        plt.plot([math.exp(-1/elem) for elem in tempList])
        plt.title("acceptance")
        plt.show()
