"""
Main file

Meike, Janneke, Nicole
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

    # Create the protein
    protein = Protein(proteinString)

    # Fold the protein once
    protein.createAminoList()

    # Initialise the stability
    stability = protein.stability
    print(stability)

    # Initialise the lowest stability
    minStability = input('Stability: ')

    # Loop while the stability is bigger than the minStability
    while stability > minStability:

        newProtein = Protein(proteinString)

        # Fold the protein again
        newProtein.createAminoList()

        # If the stability of the newly folded protein is higher replace the stability and fold
        if newProtein.stability < stability:
            protein = newProtein
            stability = newProtein.stability
            print(stability)

    #print(stability)
    print(protein)
    # Create a visual of the final fold
    protein.createPlot()
