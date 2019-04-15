"""
Class for protein.

Meike, Janneke, Nicole
"""

from amino import Amino
import random

class Protein(object):
    """
    """

    def __init__(self, proteinString):
        """
        """
        self.proteinString = proteinString
        self.aminoList = []
        self.stability = []
        self.occupied = []


    def createAminoList(self):
        """
        """
        for i in range(len(self.proteinString)):
            self.aminoList.append(Amino(i, self.proteinString[i]))
            if i == 0 or i == 1:
                self.aminoList[i].addCoordinate([0, i])
                self.occupied.append([0, i])
            else:
                prevCo = self.aminoList[(i - 1)].getCoordinates()
                posCo = []
                posCo.append([(prevCo[0] - 1), prevCo[1]])
                posCo.append([(prevCo[0] + 1), prevCo[1]])
                posCo.append([prevCo[0], (prevCo[1] - 1)])
                posCo.append([prevCo[0], (prevCo[1] + 1)])

                toRemove = []
                for j in posCo:
                    if j in self.occupied:
                        toRemove.append(j)
                for k in toRemove:
                    if k in self.occupied:
                        posCo.remove(k)

                coordinate = random.choice(posCo)
                self.aminoList[i].addCoordinate(coordinate)
                self.occupied.append(coordinate)

        print(self.occupied)


                # coordinate should be (0,0)
