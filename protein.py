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
            if i == 0:
                self.aminoList[0].addCoordinate([0,0])
                # coordinate should be (0,0)