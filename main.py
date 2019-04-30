"""
Main file

Meike, Nicole
"""

from amino import Amino
from protein import Protein
import sys


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
        if i != "H" and i != "P" and i != "C" and i != "h" and i != "p" and i != "c":
            print("Protein should only contain H, P or C")
            exit(2)
    return sys.argv[1]


if __name__ == "__main__":

    # Check the input and save the protein string
    proteinString = checkInput()

    # Initialise the lowest stability
    minStability = input('Stability: ')
    minStability = int(minStability)

    # Create the protein
    protein = Protein(proteinString)

    # Fold the protein once
    protein.createAminoList()

    # Initialise the stability
    stability = protein.stability
    # print(stability)

    # Loop while the stability is bigger than the minStability
    while stability > minStability:
    # for i in range(minStability):

        newProtein = Protein(proteinString)

        # Fold the protein again
        newProtein.createAminoList()

        # If the stability of the newly folded protein is higher replace the stability and fold
        if newProtein.stability < stability:
            protein = newProtein
            stability = newProtein.stability
            print(stability)
    print(stability)
    print(protein)
    protein.createPlot()

    # Hill climber (deze loop zou ook nog in de functie zelf kunnen (of als recursief met extra argument als counter))
    for i in range(10):
        # print(i)
        protein = protein.hillClimber()
        #print(protein.stability)
        print(protein.stability)
        print(protein)
        protein.createPlot()

    # Create a visual of the final fold
    # print(protein.stability)
    # print(protein)
    # protein.createPlot()
