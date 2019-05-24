"""
This script allows the folding of an object of the protein class using
four different algorithms: constructive (depth first), random folding,
hill climber and simulated annealing.

Usage:
python main.py algorithm protein

For more details see README on https://github.com/meikesara/Heuristieken/

Meike Kortleve, Nicole Jansen
"""

import os
import sys
import matplotlib.pyplot as plt
from classes.amino import Amino
from classes.protein import Protein
from algorithms.constructive import constructive
from algorithms.hillclimber import hillClimber
from algorithms.randomfold import randomFold
from algorithms.simulatedannealing import simulatedAnnealing
from helper.estimatestability import estimateStability
import helper.visualizer as visualizer


def checkInput():
    """
    This function checks the input into the main
    """

    # Check if the amount of arguments is 2
    if len(sys.argv) != 3:
        print("A proteinstring and method is needed")
        exit(1)

    methods = {"constructive", "simulated", "hillclimber", "random"}

    # Check if method is valid
    if sys.argv[1].lower() not in methods:
        print("The method should be constructive, simulated, hillclimber" +
                "or random.")
        exit(3)

    # Check if the second argument only contains H, P or C's
    for i in sys.argv[2]:
        if i not in {"H", "P", "C", "h", "p", "c"}:
            print("Protein should only contain H, P or C")
            exit(2)
    return sys.argv[1], sys.argv[2]


if __name__ == "__main__":

    # Check the input and save the protein string and method
    method, proteinString = checkInput()

    if method == "constructive":
        constructive(proteinString)

    else:
        bestStability = []

        # Ask user for the plane
        planes = {"2D", "3D"}
        plane = input("Do you want a 2D or 3D protein (2D/3D)? ").upper()

        # Check input
        while plane not in planes:
            plane = input("Do you want a 2D or 3D protein (2D/3D)? ").upper()

        # estimateStability(proteinString, plane)

        if method == "simulated" or method == "hillclimber":

            # Ask user for the amount of runnings and iterations and check
            while True:
                runnings = input("Enter the amount of runnings: ")
                val = runnings.isdigit()
                if val:
                    runnings = int(runnings)
                    break

            while True:
                iterations = input("Enter the amount of iterations: ")
                val = iterations.isdigit()
                if val:
                    iterations = int(iterations)
                    break

            if method == "simulated":

                # Ask user for D and check
                while True:
                    D = input("Enter D (positive integer): ")
                    val = D.isdigit()
                    if val:
                        D = int(D)
                        break

                if runnings != 1:
                    print("Results will be saved in files.")

                    stabilityFileName = """results/stabilitySimulated
                                        Annealing.txt"""
                    proteinFileName = "results/proteinSimulatedAnnealing.txt"

                    # Make a folder for results if this does not exist
                    os.makedirs(os.path.dirname(stabilityFileName),
                                exist_ok = True)
                    stabilityFile = open(stabilityFileName, "w")
                    proteinFile = open(proteinFileName, "w")

                    # Save header to the files
                    header = ["Simulated Annealing ", str(proteinString),
                              " D: ", str(D), " ", str(runnings), " x ",
                              str(iterations), "\n"]
                    stabilityFile.writelines(header)
                    proteinFile.writelines(header)

                    for i in range(runnings):
                        protein = randomFold(proteinString, plane)
                        protein = simulatedAnnealing(protein, D, iterations,
                                                     True)

                        # Update the files with the stability and protein
                        stabilityLine = [str(protein.stability), "\n"]
                        stabilityFile.writelines(stabilityLine)
                        proteinLine = [str(protein), "\n"]
                        proteinFile.writelines(proteinLine)

                        bestStability.append(protein.stability)

                    stabilityFile.close()
                    proteinFile.close()

                    print("Best stability: ", min(bestStability))

                else:
                    protein = randomFold(proteinString, plane)
                    simulatedAnnealing(protein, D, iterations, False)

            elif method == "hillclimber":

                if runnings != 1:
                    print("Results will be saved in files.")

                    stabilityFileName = "results/stabilityHillClimber.txt"
                    proteinFileName = "results/proteinHillClimber.txt"

                    # Make a folder for results if this does not exist
                    os.makedirs(os.path.dirname(stabilityFileName),
                                exist_ok = True)
                    stabilityFile = open(stabilityFileName, "w")
                    proteinFile = open(proteinFileName, "w")

                    # Write the header to the files
                    header = ["Hill climber ", str(proteinString), " ",
                              str(runnings), " x ", str(iterations), "\n"]
                    stabilityFile.writelines(header)
                    proteinFile.writelines(header)

                    for i in range(runnings):
                        protein = randomFold(proteinString, plane)
                        protein = hillClimber(protein, iterations, True)

                        # Update the files with the stability and protein
                        stabilityLine = [str(protein.stability), "\n"]
                        stabilityFile.writelines(stabilityLine)
                        proteinLine = [str(protein), "\n"]
                        proteinFile.writelines(proteinLine)
                        bestStability.append(protein.stability)

                    stabilityFile.close()
                    proteinFile.close()
                    print("Best stability: ", min(bestStability))

                else:
                    protein = randomFold(proteinString, plane)
                    hillClimber(protein, iterations, False)

        elif method == "random":

            # Ask user for the amount of running
            while True:
                runnings = input("Enter the amount of runnings: ")
                val = runnings.isdigit()
                if val:
                    runnings = int(runnings)
                    break

            if runnings != 1:
                print("Results will be saved in files.")

                stabilityFileName = "results/stabilityRandom.txt"
                proteinFileName = "results/proteinRandom.txt"

                # Make a folder for results if this does not exist
                os.makedirs(os.path.dirname(stabilityFileName), exist_ok=True)
                stabilityFile = open(stabilityFileName, "w")
                proteinFile = open(proteinFileName, "w")

                # Write header to the files
                header = ["Random ", str(proteinString), " ", str(runnings), "\n"]
                stabilityFile.writelines(header)
                proteinFile.writelines(header)

                for i in range(runnings):
                    protein = randomFold(proteinString, plane)

                    stabilityLine = [str(protein.stability), "\n"]
                    stabilityFile.writelines(stabilityLine)

                    proteinLine = [str(protein), "\n"]
                    proteinFile.writelines(proteinLine)

                    bestStability.append(protein.stability)

                stabilityFile.close()
                proteinFile.close()
                print("Best stability: ", min(bestStability))

            else:
                protein = randomFold(proteinString, plane)
                visualizer.plotProtein(protein)

        estimateStability(proteinString, plane)
