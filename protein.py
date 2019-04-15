"""
Class for protein.

Meike, Janneke, Nicole
"""

from amino import Amino

class Protein(object):
    """
    """

    def __init__(self, proteinString):
        """
        """
        self.proteinString = proteinString
        self.aminoList = []
        self.stability = []


    def createAminoList(self):
        """
        """
        for i in range(len(self.proteinString)):
            self.aminoList.append(Amino(i, self.proteinString[i]))
            if i == 0 or i == 1:
                self.aminoList[i].addCoordinate([0, i])

                # coordinate should be (0,0)
