"""
Main file

Meike, Janneke, Nicole
"""

from amino import Amino
from protein import Protein
import sys

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
    if len(sys.argv) != 2:
        print("A proteinstring is needed")
        exit(1)

    for i in sys.argv[1]:
        if i != "H" and i != "P" and i != "h" and i != "p":
            print("Protein should only contain H and P")
            exit(2)

    return sys.argv[1]


if __name__ == "__main__":
    proteinString = checkInput()
    #stabilityList = []
    for i in range(10):
        #print("i = ", i)
        folding = Fold(proteinString)
        print(folding.initProtein.stability)
    #   stabilityList.append(folding.initProtein.stability)

    #print(stabilityList)
