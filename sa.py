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

# bestProtein = 0


def simulatedAnnealing(protein, beginTemp):
    """
    Performs a hill climber
    protein is the protein that will be changed
    iterations is the amount of iterations the hill climber should be performed
    """

    stabilityList = []
    tempList = []
    temperature = beginTemp
    alfa = 0.99
    D = 6
    bestProtein = protein


    # NOTE: gebruik voor de log/ln een hogere temp om er eerder uit te komen;
        #Is alleen misschien niet helemaal precies wat je officieel zou willen doen,
        #zou ook D kleiner kunnen maken, dan komt de acceptance eerder bij een
        #lagere temp
    # while temperature > 0.18: #0.5:

    for k in range(5000):
        stabilityList.append(protein.stability)

        newProtein = protein.pullMove()

        if newProtein.stability < protein.stability:
            acceptance = 1.0
        else:
            acceptance = math.exp((protein.stability - newProtein.stability) / temperature)

        r = random.random()

        if (acceptance > r):
            protein = newProtein

        if newProtein.stability < bestProtein.stability:
            bestProtein = newProtein

        # temperature = math.log(temperature)
        # temperature = beginTemp/(1 + math.log(k))
        # k += 1

        # if k == 100:
        #     beginTemp -= 200
        #     k=0
        # temperature *= alfa
        # print(temperature)
        # temperature = beginTemp * pow(alfa,k)
        temperature = D/math.log(k + 2)
        tempList.append(temperature)

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
    return bestProtein


if __name__ == "__main__":
    proteinString = "PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP"

<<<<<<< HEAD
    temperature = 1
    finalStability = []
=======
    protein = Protein(proteinString, "2D")
>>>>>>> 0c062a5e7584bc53e09d77cefacc1782d938dd73

    for i in range(1):
        # print(i)
        protein = Protein(proteinString)
        validFolding = protein.createAminoList()

<<<<<<< HEAD
        while not validFolding:
            newProtein = Protein(proteinString)
                # Fold the protein again
            validFolding = newProtein.createAminoList()
            protein = newProtein
=======
    temperature = 10000

    # coolingRate = 0.0003

    while not validFolding:
        newProtein = Protein(proteinString, "2D")

        # Fold the protein again
        validFolding = newProtein.createAminoList()
        protein = newProtein

    bestProtein = protein
>>>>>>> 0c062a5e7584bc53e09d77cefacc1782d938dd73

        # print("Initial solution stability: ", protein.stability)

        protein = simulatedAnnealing(protein, temperature)
        finalStability.append(protein.stability)
    plt.hist(finalStability)
    print(finalStability)
    plt.show()
