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
        self.stability = 0
        self.occupied = []


    def getSurroundCo(self, prevCo):
        posCo = []
        posCo.append([(prevCo[0] - 1), prevCo[1]])
        posCo.append([(prevCo[0] + 1), prevCo[1]])
        posCo.append([prevCo[0], (prevCo[1] - 1)])
        posCo.append([prevCo[0], (prevCo[1] + 1)])
        return posCo


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
                posCo = self.getSurroundCo(prevCo)

                toRemove = []
                for j in posCo:
                    if j in self.occupied:
                        toRemove.append(j)
                for k in toRemove:
                    if k in self.occupied:
                        posCo.remove(k)

                if not posCo:
                    print("jeeeej!!!")
                    break
                coordinate = random.choice(posCo)
                self.aminoList[i].addCoordinate(coordinate)
                self.occupied.append(coordinate)

                aroundCo = self.getSurroundCo(coordinate)
                if self.aminoList[i].type == "H":
                    for l in aroundCo:
                        if l in self.occupied:
                            if self.aminoList[self.occupied.index(l)].type == "H":
                                if (self.occupied.index(l) + 1) != i:
                                    self.stability -= 1


        print( self.stability) #"stability = ",
        #print("occupied final = ", self.occupied)
