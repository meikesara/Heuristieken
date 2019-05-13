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

# bestProtein = 0


def simulatedAnnealing(protein, beginTemp, bestProtein):
    """
    Performs a hill climber
    protein is the protein that will be changed
    iterations is the amount of iterations the hill climber should be performed
    """

    stabilityList = []
    temperature = beginTemp
    alfa = 0.995
    k = 1

    while temperature > 0.0005:
        stabilityList.append(protein.stability)

        newProtein = protein.pullMove()

        if newProtein.stability <= protein.stability:
            acceptance = 1.0
        else:
            acceptance = math.exp((protein.stability - newProtein.stability) / temperature)

        # print("acceptance: ", acceptance)
        r = random.random()
        # print("random: ", r)

        if (acceptance > r):
            # print("newProtein is accepted")
            protein = newProtein
        # else:
            # print("newProtein is not accepted")

        if newProtein.stability < bestProtein.stability:
            bestProtein = newProtein

        # temperature = math.log(temperature)
        # temperature = beginTemp/(1 + math.log(k))
        k += 1
        # temperature *= alfa
        # print(temperature)
        temperature = beginTemp * pow(alfa,k)

    stabilityList.append(protein.stability)
    print("Final solution stability: ", protein.stability)
    print("Best stability: ", bestProtein.stability)
    visualizer.plotProtein(bestProtein)
    visualizer.plotStability(stabilityList)

            # print("ik ben slechter! i =", i, counter)

        #     counter = 0
        # else:
        #     counter += 1

        # if stabilityChange and (i == (iterations - 1)):
        #     stabilityList.append(protein.stability)

    # if stabilityChange:
    #     print("ik ben slechter! i = ", i, counter)
    #     return protein, stabilityList
    # else:
    #     return protein

if __name__ == "__main__":
    proteinString = "HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH"

    protein = Protein(proteinString)

    validFolding = protein.createAminoList()

    temperature = 10000

    # coolingRate = 0.0003

    while not validFolding:
        newProtein = Protein(proteinString)

        # Fold the protein again
        validFolding = newProtein.createAminoList()
        protein = newProtein

    bestProtein = protein

    print("Initial solution stability: ", protein.stability)

    simulatedAnnealing(protein, temperature, bestProtein)
