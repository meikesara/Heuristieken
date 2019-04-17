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

        # Create a list of all the possible coordinates
        coordinates = [[(prevCo[0] - 1), prevCo[1]], [(prevCo[0] + 1), prevCo[1]], [prevCo[0], (prevCo[1] - 1)], [prevCo[0], (prevCo[1] + 1)]]

        # Loop over all the coordinates
        for coordinate in coordinates:

            # Only add the coordinates to the possible coordinates list if they are occupied
            if occupied is True:
                if coordinate in self.occupied:
                    posCo.append(coordinate)

            # Only add the coordinates to the possible coordinates list if they are not occupied
            else:
                if coordinate not in self.occupied:
                    posCo.append(coordinate)
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

                # Get the surrounding coordinates that are not occupied
                posCo = self.getSurroundCo(prevCo, False)

                # If there are no surrounding coordinates available break from the loop
                if not posCo:
                    self.stability = 0
                    break

                # Randomly choose one of the possible coordinates
                coordinate = random.choice(posCo)

                # Place the amino-acid on that coordinate
                self.aminoList[id].addCoordinate(coordinate)

                # Add the coordinate to the list of occupied coordinates
                self.occupied.append(coordinate)

                # Get the coordinates of the surrounding amino acids that are occupied
                aroundCos = self.getSurroundCo(coordinate, True)

                # Get the type of the current amino-acid
                typeCo = self.aminoList[id].type

                # Check if the type is H or C
                if typeCo == "H" or typeCo == "C":

                    # Loop over the occupied coordinates around the current amino-acid
                    for aroundCo in aroundCos:

                        # Get the type of the amino-acid
                        nextCo = self.aminoList[self.occupied.index(aroundCo)].type

                        if nextCo == "H" or nextCo == "C":
                            # Wat doet deze regel?
                            if (self.occupied.index(aroundCo) + 1) != id:
                                # If both amino-acids are of type C subtract 5 from the stability
                                if typeCo == "C" and nextCo == "C":
                                    self.stability -= 5
                                # Else subtract 1
                                else:
                                    self.stability -= 1


        #print( self.stability) #"stability = ",
        #print("occupied final = ", self.occupied)

    def createPlot(self):
        """
        This method creates a visual representation of a folded protein
        """

        colorDict = {"P": 'go', "H": 'ro', "C": 'bo'}

        # Loop over the aminoList
        for i in range(len(self.aminoList)):
            amino = self.aminoList[i]

            # Get the coordinates of the previous amino-acid
            theseCo = amino.coordinate

            if i != 0:
                # Get the coordinates of the previous amino-acid
                prevCo = self.aminoList[(i - 1)].coordinate

                # Place a line from between the amino-acids
                plt.plot([prevCo[0], theseCo[0]], [prevCo[1], theseCo[1]], '-k', zorder=-1)

            # Place a dot for the amino-acid
            plt.plot([theseCo[0]], [theseCo[1]], colorDict[amino.type])

        plt.title("P = groen; H = rood; C = blauw")
        plt.axis('off')
        plt.show()
