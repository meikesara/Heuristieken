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

    finalStability = []
    sys.setrecursionlimit(1500)

    for j in range(1):

        # Check the input and save the protein string
        proteinString = checkInput()

        # stabilityList = []

        # Create the protein
        protein = Protein(proteinString)

        # Fold the protein once
        validFolding = protein.createAminoList()

        while not validFolding:
            newProtein = Protein(proteinString)

            # Fold the protein again
            validFolding = newProtein.createAminoList(proteinString)
            protein = newProtein

        # Initialise the stability
        stability = protein.stability

        # stabilityList.append(stability)

        # # Loop while the stability is bigger than the minStability
        # while stability > minStability:
        # # for i in range(minStability):
        #
        #     newProtein = Protein(proteinString)
        #
        #     # Fold the protein again
        #     newProtein.createAminoList()
        #
        #     # If the stability of the newly folded protein is higher replace the stability and fold
        #     if newProtein.stability < stability:
        #         protein = newProtein
        #         stability = newProtein.stability
        #         print(stability)
        # print(stability)
        # # print(protein)
        print(protein.stability)
        # visualizer.plotProtein(protein)

        # Hill climber (deze loop zou ook nog in de functie zelf kunnen (of als recursief met extra argument als counter))
        # for i in range(1000):
        #     stabilityList.append(protein.stability)
        #     protein = protein.hillClimber()
        #     # if (i % 10) == 0:
        #     #     print(protein.stability)
        protein, stabilityList = hillClimber(protein, 10000, True)

        # # Create a visual of the final fold
        # # print(protein.stability)
        # finalStability.append(protein.stability)

        # Create a visual of the final fold
        print(protein.stability)
        print(protein)
        visualizer.plotProtein(protein)
        visualizer.plotStability(stabilityList)

    # # plt.hist(finalStability)
    # plt.xlabel("Stabiliteit")
    # plt.ylabel("Aantal vouwingen")
    # plt.show()
