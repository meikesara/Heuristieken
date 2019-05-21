"""
# TODO: comment hier

Meike Kortleve, Nicole Jansen
"""

from amino import Amino
from protein import Protein
import random
import copy
import math
import visualizer
import matplotlib.pyplot as plt
from randomfold import randomFold


def simulatedAnnealing(proteinString, beginTemp, iterations, runnings):
    """
    Performs simulated annealing with a logarithmic cooling rate.
    The first argument
    # TODO: ddd
    protein is the protein that will be changed
    iterations is the amount of iterations the hill climber should be performed
    """

    # Initialize variables
    stabilityList = []
    tempList = []
    temperature = beginTemp
    D = 6

    # TODO: dit in de main doen voor alle functies

    # Initialize protein
    protein = Protein(proteinString, "2D")

    # protein = randomFold
    validFolding = protein.createAminoList()

    # Check if protein is folded correctly
    while not validFolding:
        newProtein = Protein(proteinString, "2D")

        # Fold the protein again
        validFolding = newProtein.createAminoList()
        protein = newProtein

    # Inititalize bestProtein
    bestProtein = protein

    for k in range(iterations):

        # Append current stability to the stabilityList
        stabilityList.append(protein.stability)

        # Create a new protein by moving the amino acids
        newProtein = protein.pullMove()

        # Update the temperature
        temperature = D/math.log(k + 2 + 10)

        # Add temperature to the list of temperatures
        tempList.append(temperature)

        # Set acceptance to 1 if the stability of the new protein is better.
        if newProtein.stability < protein.stability:
            acceptance = 1
        else:
            # Calculate the acceptance
            acceptance = math.exp((protein.stability - newProtein.stability) / temperature)

        # Accept newProtein if the acceptance exceeds random number from 0 to 1.
        if (acceptance > random.random()):
            protein = newProtein

        # If the stability of newProtein is the best set as bestProtein
        if newProtein.stability < bestProtein.stability:
            bestProtein = newProtein

        # temperature = beginTemp/(1 + math.log(k))
        # temperature = beginTemp * pow(alfa,k)

    if runnings == 1:
        # Add final stability to stabilityList
        stabilityList.append(protein.stability)

        print("Final solution stability: ", protein.stability)
        print("Best stability: ", bestProtein.stability)

        visualizer.plotProtein(bestProtein)
        visualizer.plotStability(stabilityList)

        # Plot temperatuur verloop
        plt.plot(tempList)
        plt.title("Temperature")


        # Plot acceptance verloop
        plt.figure()
        plt.plot([math.exp(-1/elem) for elem in tempList])
        plt.title("acceptance")
        plt.show()
    else:
        return bestProtein.stability
