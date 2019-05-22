"""
Main file
# TODO: dd

Meike Kortleve, Nicole Jansen
"""

from amino import Amino
from protein import Protein
import sys
import matplotlib.pyplot as plt
import visualizer
from hillclimber import hillClimber
from randomfold import randomFold
from constructive import constructive
from simulatedannealing import simulatedAnnealing

def checkInput():
    """
    This function checks the input into the main
    """

    # Check if the amount of arguments is 2
    if len(sys.argv) != 3:
        print("A proteinstring is needed and method is needed")
        exit(1)

    # Check if the second argument only contains H, P or C's
    for i in sys.argv[2]:
        if i not in {"H", "P", "C", "h", "p", "c"}:
            print("Protein should only contain H, P or C")
            exit(2)
    return sys.argv[2]


if __name__ == "__main__":

    # Check the input and save the protein string
    proteinString = checkInput()
    finalStability = []

    if sys.argv[1] == "constructive":
        constructive(proteinString)

    elif sys.argv[1] == "simulated":

        protein = randomFold(proteinString, "2D")

        D = int(input("Enter the D: "))
        runnings = int(input("Enter the amount of runnings: "))
        iterations = int(input("Enter the amount of iterations: "))

        if runnings != 1:
            for i in range(runnings):
                stability = simulatedAnnealing(protein, D, iterations, True)
                finalStability.append(stability)
            print(finalStability)
        else:
            simulatedAnnealing(protein, D, iterations, False)

    elif sys.argv[1] == "hillclimber":

        protein = randomFold(proteinString, "2D")

        runnings = int(input("Enter the amount of runnings: "))
        iterations = int(input("Enter the amount of iterations: "))

        if runnings != 1:
            finalStabilityList = []
            for i in range(runnings):
                stability = hillClimber(protein, iterations, True)
                finalStabilityList.append(stability)
            print(finalStabilityList)
        else:
            hillClimber(protein, iterations, False)

    elif sys.argv[1] == "random":

        randomList = []

        iterations = int(input("Enter the amount of iterations: "))

        for i in range(iterations):
            random = randomFold(proteinString, "2D")
            randomList.append(random)

        print(randomList)




            # # Random folding of protein
            # protein = randomFold(protein, -10)
            # print(protein)
            # print(protein.stability)
            # visualizer.plotProtein(protein)

            #
            # # Create a visual of the final fold
            # if j == (times - 1):
            #     print(protein.stability)
            #     print(protein)
            #     visualizer.plotProtein(protein)
            #     # visualizer.plotStability(stabilityList)

        # plt.hist(finalStability)
        # plt.xlabel("Stabiliteit")
        # plt.ylabel("Aantal vouwingen")
        # plt.show()
