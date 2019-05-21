"""
Main file

Meike Kortleve, Nicole Jansen
"""

from amino import Amino
from protein import Protein
import sys
import matplotlib.pyplot as plt
import visualizer
from hillclimber import hillClimber
from randomfold import randomFold

def checkInput():
    """
    This function checks the input into the main
    """

    # Check if the amount of arguments is 2
    if len(sys.argv) != 2:
        print("A proteinstring is needed")
        exit(1)

    # Check if the second argument only contains H, P or C's
    for i in sys.argv[1]:
        if i not in {"H", "P", "C", "h", "p", "c"}:
            print("Protein should only contain H, P or C")
            exit(2)
    return sys.argv[1]


if __name__ == "__main__":

    # Check the input and save the protein string
    proteinString = checkInput()

    finalStability = []
    times = 1
    for j in range(times):

        # Create the protein
        protein = Protein(proteinString, "3D")

        # Fold the protein once
        validFolding = protein.createAminoList()

        while not validFolding:
            newProtein = Protein(proteinString, protein.plane)

            # Fold the protein again
            validFolding = newProtein.createAminoList()
            protein = newProtein

        # Random folding of protein
        protein = randomFold(protein, -49)
        # print(protein)
        # print(protein.stability)
        # visualizer.plotProtein(protein)

        # print(protein.stability)
        # visualizer.plotProtein(protein)
        # # Hill climber
        # protein, stabilityList = hillClimber(protein, 1000, True)

        # Create a visual of the final fold
        print(protein.stability)
        finalStability.append(protein.stability)

        # Create a visual of the final fold
        if j == (times - 1):
            print(protein.stability)
            print(protein)
            visualizer.plotProtein(protein)
            # visualizer.plotStability(stabilityList)

    # plt.hist(finalStability)
    # plt.xlabel("Stabiliteit")
    # plt.ylabel("Aantal vouwingen")
    # plt.show()
