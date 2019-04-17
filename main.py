"""
Main file

Meike, Janneke, Nicole
"""

from amino import Amino
from protein import Protein
import sys

### Willen we dit nog bewaren?
class Fold(object):
    """
    """

    def __init__(self, proteinString):
        """
        """
        self.initProtein = Protein(proteinString)
        self.initProtein.createAminoList()
        #self.protein = self.initProtein #Check of deze niet in elkaar updaten


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

    # Fold the protein once
    folding = Fold(proteinString)

    # Initialise the stability
    stability = folding.initProtein.stability
    print(stability)

    # Initialise the lowest stability
    minStability = -8

    # Loop while the stability is bigger than the minStability
    while stability > minStability:

        # Fold the protein again
        newfold = Fold(proteinString)

        # If the stability of the newly folded protein is higher replace the stability and fold
        if newfold.initProtein.stability < stability:
            folding = newfold
            stability = folding.initProtein.stability
            print(stability)

    #print(stability)
    print(folding.initProtein)
    # Create a visual of the final fold
    folding.initProtein.createPlot()
