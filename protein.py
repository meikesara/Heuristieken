"""
Class for protein.

Meike Kortleve, Nicole Jansen
"""

from amino import Amino
import random
import matplotlib.pyplot as plt
import copy


class Protein(object):
    """
    Protein class containing necessary attributes and methods for setting up
    and updating protein.
    """

    # Initialise variables
    def __init__(self, proteinString, plane):

        """
        This method initialises the variables of the Protein object
        """

        self.proteinString = proteinString
        self.aminoList = []
        self.stability = 0
        self.occupied = []
        self.proteinLength = len(proteinString)
        self.plane = plane


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
        """
        This method returns a diagonal of the current amino acid, that is not occupied
        and the C if the current amino acid is not the last or the first
        """

        # Get the coordinates and position of the current amino acid
        currentCo = currentAmino.coordinate
        index = self.aminoList.index(currentAmino)

        # If the amino acid is the last get the coordinates of the previous amino acid
        if (index + 1) == len(self.aminoList):
            previousAmino = self.aminoList[index - 1]
            otherCo = previousAmino.coordinate

        # If the amino acid is the first get the coordinates of the next amino acid
        elif index == 0:
            nextAmino = self.aminoList[index + 1]
            otherCo = nextAmino.coordinate

        # If the amino acid is in the middle get the coordinates of the next and previous amino acid
        else:
            nextAmino = self.aminoList[index + 1]
            otherCo = nextAmino.coordinate

            previousAmino = self.aminoList[index - 1]
            previousCo = previousAmino.coordinate

        # Calculate the absolute difference between the x and y coordinates
        x = abs(currentCo[0] - otherCo[0])
        y = abs(currentCo[1] - otherCo[1])
        # print('x,y:', x,y)
        if self.plane == "3D":
            z = abs(currentCo[2] - otherCo[2])
            # print("z:", z)

        # Check if the amino acid is the first or last in the protein
        # NOTE: zou je ipv len(self.aminoList) ook self.proteinLength kunnen gebruiken?
        if (index + 1) == len(self.aminoList) or index == 0:
            # TODO: Voor eerste en laatste nog 3D doen
            # Create a list of diagonals
            if self.plane == "2D":
                if x == 1:
                    diagonals = [[otherCo[0], (otherCo[1] - 1)], [otherCo[0], (otherCo[1] + 1)]]
                else:
                    diagonals = [[(otherCo[0] + 1), otherCo[1]], [(otherCo[0] - 1), otherCo[1]]]
            else:
                if x == 1:
                    diagonals = [[otherCo[0], (otherCo[1] - 1), otherCo[2]],
                                 [otherCo[0], (otherCo[1] + 1), otherCo[2]],
                                 [otherCo[0], otherCo[1], (otherCo[2] - 1)],
                                 [otherCo[0], otherCo[1], (otherCo[2] + 1)]]
                elif y == 1:
                    diagonals = [[(otherCo[0] - 1), otherCo[1], otherCo[2]],
                                 [(otherCo[0] + 1), otherCo[1], otherCo[2]],
                                 [otherCo[0], otherCo[1], (otherCo[2] - 1)],
                                 [otherCo[0], otherCo[1], (otherCo[2] + 1)]]
                else:
                    diagonals = [[(otherCo[0] - 1), otherCo[1], otherCo[2]],
                                 [(otherCo[0] + 1), otherCo[1], otherCo[2]],
                                 [otherCo[0], (otherCo[1] - 1), otherCo[2]],
                                 [otherCo[0], (otherCo[1] + 1), otherCo[2]]]




            # Randomly shuffle the diagonals
            random.shuffle(diagonals)

            # Return the first diagonal that is not occupied
            for diagonal in diagonals:
                if diagonal not in self.occupied:
                    return [diagonal]

        else:
            # Create a list of diagonals of the currentCoordinates
            diagonals = []
            if self.plane == "2D":
                # Only xy-plane
                planes = [[0, 1]]
            else:
                planes = [[0, 1], [0, 2], [1, 2]]

            for plane in planes:
                # TODO: zorgen dat deze code écht beter/netter wordt
                if self.plane == "2D":
                    check = ["0", "1"]
                else:
                    check = ["0", "1", "2"]
                check[plane[0]] = (currentCo[plane[0]] + 1)
                check[plane[1]] = (currentCo[plane[1]] + 1)
                if self.plane == "3D":
                    for i in check:
                        if type(i) == str:
                            check[int(i)] = currentCo[int(i)]
                diagonals.append(check)

                if self.plane == "2D":
                    check = ["0", "1"]
                else:
                    check = ["0", "1", "2"]
                check[plane[0]] = (currentCo[plane[0]] + 1)
                check[plane[1]] = (currentCo[plane[1]] - 1)
                if self.plane == "3D":
                    for i in check:
                        if type(i) == str:
                            check[int(i)] = currentCo[int(i)]
                diagonals.append(check)

                if self.plane == "2D":
                    check = ["0", "1"]
                else:
                    check = ["0", "1", "2"]
                check[plane[0]] = (currentCo[plane[0]] - 1)
                check[plane[1]] = (currentCo[plane[1]] + 1)
                if self.plane == "3D":
                    for i in check:
                        if type(i) == str:
                            check[int(i)] = currentCo[int(i)]
                diagonals.append(check)

                if self.plane == "2D":
                    check = ["0", "1"]
                else:
                    check = ["0", "1", "2"]
                check[plane[0]] = (currentCo[plane[0]] - 1)
                check[plane[1]] = (currentCo[plane[1]] - 1)
                if self.plane == "3D":
                    for i in check:
                        if type(i) == str:
                            check[int(i)] = currentCo[int(i)]
                diagonals.append(check)

                # diagonals.append([(currentCo[plane[0]] + 1), (currentCo[plane[1]] + 1)])
                # diagonals.append([(currentCo[plane[0]] + 1), (currentCo[plane[1]] - 1)])
                # diagonals.append([(currentCo[plane[0]] - 1), (currentCo[plane[1]] + 1)])
                # diagonals.append([(currentCo[plane[0]] - 1), (currentCo[plane[1]] - 1)])
            # diagonals = [[(currentCo[0] + 1), (currentCo[1] + 1)], [(currentCo[0] + 1), (currentCo[1] - 1)],
            #              [(currentCo[0] - 1), (currentCo[1] + 1)], [(currentCo[0] - 1), (currentCo[1] - 1)]]
            # di = [d for dc in dia for d in dc]
            # print(diagonals)
            # Randomly shuffle the diagonals
            random.shuffle(diagonals)

            # Return the first available diagonal and C
            for diagonal in diagonals:

                # Check if the diagonal is not occupied
                if diagonal not in self.occupied:

                    # Get the occupied surrounding coordinates of the diagonal
                    surroundCo = self.getSurroundCo(diagonal, True)

                    # Check if the coordinates of the next aminoacid ar in the surrounding cordinates of the diagonal
                    if otherCo in surroundCo:

                        # Get the coordinates of the C
                        # TODO: zorgen dat dit nog goed gaat worden voor 3D
                        if self.plane == "2D":
                            if x == 1:
                                CCo = [currentCo[0], diagonal[1]]
                            else:
                                CCo = [diagonal[0], currentCo[1]]
                        else:
                            difference = [(currentCo[i] - diagonal[i]) for i in range(len(currentCo))]
                            sameCoIndex = difference.index(0)
                            CCo = [0,0,0]
                            CCo[sameCoIndex] = currentCo[sameCoIndex]

                            xyz = [0,1,2]
                            xyz.remove(xyz[sameCoIndex])
                            CCo[xyz[0]] = currentCo[xyz[0]]
                            CCo[xyz[1]] = diagonal[xyz[1]]
                            # xyzDict = {0:x, 1:y, 2:z}
                            # if xyzDict[xyz[0]] == 1:
                            #      CCo[xyz[0]] = currentCo[xyz[0]]
                            #      CCo[xyz[1]] = diagonal[xyz[1]]
                            # else:
                            #     CCo[]
                            # print("CCo = ", CCo)

                        # Check if C is not occupied or if they are the coordinates of the previous amino acid
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
        if self.plane == "3D":
            [co.append(prevCo[2]) for co in coordinates]
            coordinates.append([prevCo[0], prevCo[1], (prevCo[2] - 1)])
            coordinates.append([prevCo[0], prevCo[1], (prevCo[2] + 1)])

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
        This method folds the protein randomly
        """

        # TODO: misschien kunnen we ook voorkomen dat een eiwit niet goed vouwt

        # Loop over the letters in the proteinString
        for id in range(self.proteinLength):

            # Add amino acid to the aminoList
            self.aminoList.append(Amino(id, self.proteinString[id].upper()))

            # Place the first and second amino-acid
            # The coordinates of first amino-acid are (0,0) (or (0,0,0) for 3D)
            # The coordinates of second amino-acid are (0,1) (or (0,1,0) for 3D)
            if id == 0 or id == 1:
                thisCoordinate = [0, id]
                if self.plane == "3D":
                    thisCoordinate.append(0)
                self.aminoList[id].addCoordinate(thisCoordinate)
                self.occupied.append(thisCoordinate)

            # The remaining amino-acids are randomly placed
            else:
                # Get the coordinates of the previous amino-acid
                prevCo = self.aminoList[(id - 1)].coordinate

                # Get the surrounding coordinates that are not occupied
                posCo = self.getSurroundCo(prevCo, False)

                # If there are no surrounding coordinates available break from the loop
                if not posCo:
                    self.stability = 0
                    return False

                # Randomly choose one of the possible coordinates
                coordinate = random.choice(posCo)

                # Place the amino-acid on that coordinate
                self.aminoList[id].addCoordinate(coordinate)

                # Add the coordinate to the list of occupied coordinates
                self.occupied.append(coordinate)

                # Update the stability
                self.stabilityUpdate(self.aminoList[id])

        return True


    def stabilityUpdate(self, amino, replace=False):
        """
        Method for updating the stability of the protein.
        If replace is set to True, stability will be increased (worsened)
        """

        # Get id, coordinate and type of the current amino-acid
        id = amino.id
        coordinate = amino.coordinate
        typeCo = amino.type

        # Only need to update stability if amino is H or C
        if typeCo in {"H", "C"}:

            # Get surrounding coordinates of given amino acid that are occupied
            aroundCos = self.getSurroundCo(coordinate, True)

            # For each amino next to given amino, check if they create a bond
            for aroundCo in aroundCos:

                # Get the id and type of the neighboring amino acid
                idAround = self.occupied.index(aroundCo)
                aroundType = self.aminoList[idAround].type

                # Bond can only be created when both amino acids are H and/or C
                if aroundType in {"H", "C"}:

                    # Check if amino is not connected in protein to given amino
                    if idAround != (id + 1) and idAround != (id - 1):

                        # Stronger bond created when both aminos are type C
                        if typeCo == "C" and aroundType == "C":
                            if replace:
                                self.stability += 5
                            else:
                                self.stability -= 5

                        # Weaker bond created when at least one amino is type H
                        else:
                            if replace:
                                self.stability += 1
                            else:
                                self.stability -= 1


    def pullMove(self):
        """
        # TODO: comment
        """

        # Choose random amino to move
        amino = random.choice(self.aminoList)

        # Get coordinates of where to move chosen and previous amino to
        coordinates = self.getDiagonalCo(amino)

        # Make sure amino can be moved
        while not coordinates:
            amino = random.choice(self.aminoList)
            coordinates = self.getDiagonalCo(amino)

        # Create copy of protein (original protein)
        newProtein = copy.deepcopy(self)

        newProtein.stabilityUpdate(amino, True)

        chosenCo = coordinates[0]
        newProtein.aminoList[amino.id].coordinate = chosenCo
        newProtein.occupied[amino.id] = chosenCo

        newProtein.stabilityUpdate(newProtein.aminoList[amino.id])

        # if random chosen amino to move is not first or last, make sure other
        # amino acids are replaced to keep validity of the protein
        if len(coordinates) == 2:

            previousAmino = self.aminoList[(amino.id - 1)]
            previousAminoCo = previousAmino.coordinate
            previousAminoId = previousAmino.id


            if previousAminoCo != coordinates[1]:

                newProtein.stabilityUpdate(previousAmino, True)

                newProtein.aminoList[(amino.id - 1)].coordinate = coordinates[1]
                newProtein.occupied[(amino.id - 1)] = coordinates[1]

                newProtein.stabilityUpdate(newProtein.aminoList[(amino.id - 1)])

            surCoPrev = newProtein.getSurroundCo(newProtein.aminoList[(amino.id - 1)].coordinate, True)

            # NOTE: vraag me af of deze if-statement (nog) nodig is,
            # aangezien hetzelfde gecontroleerd wordt in functie
            if newProtein.aminoList[(amino.id - 2)].coordinate not in surCoPrev:
                newProtein.moveAminos(self, (amino.id - 2))

        return newProtein


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

        self.stabilityUpdate(self.aminoList[idToMove], True)

        self.aminoList[idToMove].coordinate = oldProtein.aminoList[(idToMove + 2)].coordinate
        self.occupied[idToMove] = oldProtein.aminoList[(idToMove + 2)].coordinate

        self.stabilityUpdate(self.aminoList[idToMove])

        self.moveAminos(oldProtein, (idToMove - 1))
