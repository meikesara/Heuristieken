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

    def __init__(self):
        """
        """
        self.protein = Protein(sys.argv[1]).createAminoList()





if __name__ == "__main__":
    proteinString = sys.argv[1]
