"""
Class for protein.

Meike, Nicole
"""

from amino import Amino
import random
import matplotlib.pyplot as plt
import copy


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

    def getDiagonalCo(self, currentAmino):

        posDia = []

        currentCo = currentAmino.coordinate

        index = self.aminoList.index(currentAmino)

        if (index + 1) == len(self.aminoList):
            previousAmino = self.aminoList[index - 1]
            previousCo = previousAmino.coordinate

            x = abs(currentCo[0] - previousCo[0])
            y = abs(currentCo[1] - previousCo[1])

            if x == 1:
                diagonals = [[previousCo[0], (previousCo[1] - 1)], [previousCo[0], (previousCo[1] + 1)]]
            else:
                diagonals = [[(previousCo[0] + 1), previousCo[1]], [(previousCo[0] - 1), previousCo[1]]]

            random.shuffle(diagonals)

            for diagonal in diagonals:
                if diagonal not in self.occupied:
                    return [diagonal]

        elif index == 0:
            nextAmino = self.aminoList[index + 1]
            nextCo = nextAmino.coordinate

            x = abs(currentCo[0] - nextCo[0])
            y = abs(currentCo[1] - nextCo[1])

            if x == 1:
                diagonals = [[nextCo[0], (nextCo[1] - 1)], [nextCo[0], (nextCo[1] + 1)]]
            else:
                diagonals = [[(nextCo[0] + 1), nextCo[1]], [(nextCo[0] - 1), nextCo[1]]]

            random.shuffle(diagonals)

            for diagonal in diagonals:
                if diagonal not in self.occupied:
                    return [diagonal]

        else:
            nextAmino = self.aminoList[index + 1]
            nextCo = nextAmino.coordinate

            previousAmino = self.aminoList[index - 1]
            previousCo = previousAmino.coordinate

            diagonals = [[(currentCo[0] + 1), (currentCo[1] + 1)], [(currentCo[0] + 1), (currentCo[1] - 1)], [(currentCo[0] - 1), (currentCo[1] + 1)], [(currentCo[0] - 1), (currentCo[1] - 1)]]
            random.shuffle(diagonals)

            x = abs(currentCo[0] - nextCo[0])
            y = abs(currentCo[1] - nextCo[1])

            for diagonal in diagonals:
                if diagonal not in self.occupied:
                    surroundCo = self.getSurroundCo(diagonal, True)
                    if nextCo in surroundCo:
                        if x == 1:
                            CCo = [currentCo[0], diagonal[1]]
                        else:
                            CCo = [diagonal[0], currentCo[1]]
                        if (CCo not in self.occupied) or (CCo == previousCo):
                            return [diagonal, CCo]


    def getSurroundCo(self, prevCo, occupied):
        """
        This method gets the 4 surrounding coordinates
        occupied is true if you want to know which surrounding coordinates are occupied
        occupied is false if you want to know which surrounding coordinates are not occupied
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
                            # Check if amino is not connected in protein to amino
                            if (self.occupied.index(aroundCo) + 1) != id:
                                # If both amino-acids are of type C subtract 5 from the stability
                                if typeCo == "C" and nextCo == "C":
                                    self.stability -= 5
                                # Else subtract 1
                                else:
                                    self.stability -= 1


        #print( self.stability) #"stability = ",
        #print("occupied final = ", self.occupied)


    def hillClimber(self):
        """
        """

        # Choose random amino to move
        amino = random.choice(self.aminoList)
        print("id = ", amino.id)
        coordinates = self.getDiagonalCo(amino)

        # Make sure amino can be moved
        while not coordinates:
            amino = random.choice(self.aminoList)
            coordinates = self.getDiagonalCo(amino)

        #print(amino.coordinate)
        print("coordinates = ", coordinates)

        # Create copy of self (original protein)
        newProtein = copy.deepcopy(self)

        #check stabiliteit voor verplaatsen
        if amino.type in ["H", "C"]:
            surCo = self.getSurroundCo(amino.coordinate, True)
            #regel 188 protein.py #aannemend dat hier een functie voor is (check of eromheen aminozuren zitten die stabiliteit verlagen)

        chosenCo = coordinates[0] #aannemend dat output getDiagonalCo: [diaCo (L), CCo] (of [diaCo] als animo = eerste of laatste)
        newProtein.aminoList[amino.id].coordinate = chosenCo
        newProtein.occupied[amino.id] = chosenCo

        print(newProtein)
        newProtein.createPlot()
        # TODO: update stability met (nog te maken) functie nu amino verplaatst is

        # if random chosen amino to move is not first or last, make sure other
        # aminoacids are replaced to keep validity of the protein
        if len(coordinates) == 2:
            print(newProtein.aminoList[(amino.id - 1)].coordinate)
            if newProtein.aminoList[(amino.id - 1)].coordinate != coordinates[1]:
                # TODO: update stability met (nog te maken) functie
                newProtein.aminoList[(amino.id - 1)].coordinate = coordinates[1]
                newProtein.occupied[(amino.id - 1)] = coordinates[1]
                # TODO: update stability nu amino verplaatst
                print(newProtein)
                newProtein.createPlot()
            print(newProtein.aminoList[(amino.id - 1)].coordinate)
            surCoPrev = newProtein.getSurroundCo(newProtein.aminoList[(amino.id - 1)].coordinate, True)
            print(surCoPrev)

            if newProtein.aminoList[(amino.id - 2)].coordinate not in surCoPrev:
                newProtein.moveAminos(self, (amino.id - 2))

        if newProtein.stability < self.stability:
            self = newProtein


    def moveAminos(self, oldProtein, idToMove):
        """
        Method for moving aminoacids to create a valid protein.
        oldProtein is the protein from which the new protein (neighbourhood) is
        created.
        idToMove is the id of the aminoacid that needs to be moved.
        """

        surCoPrev = self.getSurroundCo(self.aminoList[(idToMove + 1)].coordinate, True)
        if (self.aminoList[idToMove].coordinate in surCoPrev) or (idToMove < 0):
            return

        # TODO: update stability met (nog te maken) functie
        self.aminoList[idToMove].coordinate = oldProtein.aminoList[(idToMove + 2)].coordinate
        self.occupied[idToMove] = oldProtein.aminoList[(idToMove + 2)]
        # TODO: update stability nu amino verplaatst

        print(self)
        self.createPlot()

        self.moveAminos(oldProtein, (idToMove - 1))


    def createPlot(self):
        """
        This method creates a visual representation of a folded protein
        """

        colorDict = {"P": 'go', "H": 'ro', "C": 'bo'}

        #fig = plt.figure()

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
        plt.axis('equal')
        plt.axis('off')
        plt.show()
