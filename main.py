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

    if sys.argv[1] == "constructive":
        constructive(proteinString)

    else:
        plane = input("Do you want a 2D or 3D protein? ")

        if sys.argv[1] == "simulated":

            D = int(input("Enter the D: "))
            runnings = int(input("Enter the amount of runnings: "))
            iterations = int(input("Enter the amount of iterations: "))

            if runnings != 1:

                stabilityFile = open("results/stabilitySimulatedAnnealing.txt", "w")
                proteinFile = open("results/proteinSimulatedAnnealing.txt", "w")
                header = [ "Simulated Annealing ", str(proteinString), " ", "D: ", str(D), " ", str(runnings), " x ", str(iterations), '\n']
                stabilityFile.writelines(header)
                proteinFile.writelines(header)

                for i in range(runnings):
                    protein = randomFold(proteinString, plane)
                    protein = simulatedAnnealing(protein, D, iterations, True)
                    stabilityLine = [str(protein.stability), '\n']
                    stabilityFile.writelines(stabilityLine)
                    proteinLine = [str(protein), '\n']
                    proteinFile.writelines(proteinLine)

                stabilityFile.close()
                proteinFile.close()

            else:
                protein = randomFold(proteinString, plane)
                simulatedAnnealing(protein, D, iterations, False)

        elif sys.argv[1] == "hillclimber":

            runnings = int(input("Enter the amount of runnings: "))
            iterations = int(input("Enter the amount of iterations: "))

            if runnings != 1:
                stabilityFile = open("results/stabilityHillClimber.txt", "w")
                proteinFile = open("results/proteinHillClimber.txt", "w")
                header = [ "Hill climber ", str(proteinString), " ", str(runnings), " x ", str(iterations), '\n']
                stabilityFile.writelines(header)
                proteinFile.writelines(header)

                for i in range(runnings):
                    protein = randomFold(proteinString, plane)
                    protein = hillClimber(protein, iterations, True)

                    stabilityLine = [str(protein.stability), '\n']
                    stabilityFile.writelines(stabilityLine)
                    proteinLine = [str(protein), '\n']
                    proteinFile.writelines(proteinLine)

                stabilityFile.close()
                proteinFile.close()

            else:
                protein = randomFold(proteinString, plane)
                hillClimber(protein, iterations, False)

        elif sys.argv[1] == "random":

            iterations = int(input("Enter the amount of iterations: "))

            if iterations != 1:

                stabilityFile = open("results/stabilityRandom.txt", "w")
                proteinFile = open("results/proteinRandom.txt", "w")

                header = [ "Random ", str(proteinString), " ", str(iterations), '\n']
                stabilityFile.writelines(header)
                proteinFile.writelines(header)

                for i in range(iterations):
                    protein = randomFold(proteinString, plane)

                    stabilityLine = [str(protein.stability), '\n']
                    stabilityFile.writelines(stabilityLine)
                    proteinLine = [str(protein), '\n']
                    proteinFile.writelines(proteinLine)

                stabilityFile.close()
                proteinFile.close()

            else:
                protein = randomFold(proteinString, plane)
                visualizer.plotProtein(protein)
