"""
Class for protein.

Meike, Janneke, Nicole
"""

from amino import Amino
import random
import matplotlib.pyplot as plt

class Protein(object):
    """
    """

    # Initialise variables
    def __init__(self, proteinString):
        """
        """

        self.proteinString = proteinString
        self.aminoList = []
        self.stability = 0
        self.occupied = []

    def __str__(self):
        """
        This method prints the type of amino acids in the protein and their coordinates
        """

        output = ""

        # Loop over the aminoacids in the aminoList
        for amino in self.aminoList:
            # Add the type to the output
            output += amino.type

            # Add the coordinates to the output
            output += str(amino.coordinate) + " "
        return output


    def getSurroundCo(self, prevCo, occupied):
        """
        This method gets the 4 surrounding coordinates
        occupied is true if you want to know which surrounding coordinates are occupied
        occupies is false if you want to know which surrounding coordinates are not occupied
        """

        posCo = []
        coordinates = [[(prevCo[0] - 1), prevCo[1]], [(prevCo[0] + 1), prevCo[1]], [prevCo[0], (prevCo[1] - 1)], [prevCo[0], (prevCo[1] + 1)]]

        for coordinate in coordinates:
            if occupied is True:
                if coordinate in self.occupied:
                    posCo.append(coordinate)
            else:
                if coordinate not in self.occupied:
                    posCo.append(coordinate)

        # posCo.append([(prevCo[0] - 1), prevCo[1]])
        # posCo.append([(prevCo[0] + 1), prevCo[1]])
        # posCo.append([prevCo[0], (prevCo[1] - 1)])
        # posCo.append([prevCo[0], (prevCo[1] + 1)])
        return posCo


    def createAminoList(self):
        """
        This method folds the protein
        """

        # Loop over the letters in the proteinString
        for id in range(len(self.proteinString)):

            # Add amino acid to the aminoList
            self.aminoList.append(Amino(id, self.proteinString[id].upper()))

            # Place the first and second amino-acid
            # The coordinates of the first amino-acid are (0,0)
            # The coordinates of the second amino-acis are (0,1)
            if id == 0 or id == 1:
                self.aminoList[id].addCoordinate([0, id])
                self.occupied.append([0, id])

            # The remaining amino-acids are randomly placed
            else:
                # Get the coordinates of the previous amino-acid
                prevCo = self.aminoList[(id - 1)].coordinate

                # Get the surrounding coordinates
                posCo = self.getSurroundCo(prevCo, False)

                toRemove = []
                # for j in posCo:
                #     if j in self.occupied:
                #         toRemove.append(j)
                # for k in toRemove:
                #     if k in self.occupied:
                #         posCo.remove(k)

                if not posCo:
                    self.stability = 0
                    break
                coordinate = random.choice(posCo)
                self.aminoList[id].addCoordinate(coordinate)
                self.occupied.append(coordinate)

                aroundCo = self.getSurroundCo(coordinate, True)
                typeCo = self.aminoList[id].type
                if typeCo == "H" or typeCo == "C":
                    for l in aroundCo:
                        if l in self.occupied:
                            nextCo = self.aminoList[self.occupied.index(l)].type
                            if nextCo == "H" or nextCo == "C":
                                if (self.occupied.index(l) + 1) != id:
                                    if typeCo == "C" and nextCo == "C":
                                        self.stability -= 5
                                    else:
                                        self.stability -= 1


        #print( self.stability) #"stability = ",
        #print("occupied final = ", self.occupied)

    def createPlot(self):
        """
        This method creates a visual representation of a folded protein
        """

        colorDict = {"P": 'go', "H": 'ro', "C": 'bo'}
        for i in range(len(self.aminoList)):
            amino = self.aminoList[i]
            theseCo = amino.coordinate
            if i != 0:
                prevCo = self.aminoList[(i - 1)].coordinate
                plt.plot([prevCo[0], theseCo[0]], [prevCo[1], theseCo[1]], '-k')
            plt.plot([theseCo[0]], [theseCo[1]], colorDict[amino.type])
        plt.title("P = groen; H = rood; C = blauw")
        plt.show()
