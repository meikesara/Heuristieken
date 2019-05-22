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

        self.proteinString = proteinString.upper()
        self.aminoList = []
        self.stability = 0
        self.occupied = []
        self.proteinLength = len(proteinString)
        self.plane = plane.upper()


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
        if (index + 1) == self.proteinLength:
            previousAmino = self.aminoList[index - 1]
            otherCo = previousAmino.coordinate

        # If the amino acid is the first get the coordinates of the next amino acid
        elif index == 0:
            nextAmino = self.aminoList[index + 1]
            otherCo = nextAmino.coordinate

        # If the amino acid is in the middle get the coordinates of the next and previous amino acid
        else:
            nOrP = [index + 1, index - 1]

            # TODO: dit werkend maken
            random.shuffle(nOrP)

            nextAmino = self.aminoList[nOrP[0]]
            otherCo = nextAmino.coordinate

            previousAmino = self.aminoList[nOrP[1]]
            previousCo = previousAmino.coordinate


        # Calculate the absolute difference between the x, y (and z) coordinates
        lengthCo = len(currentCo)
        xyzDif = [abs(currentCo[i] - otherCo[i]) for i in range(lengthCo)]

        # Check if the amino acid is the first or last in the protein
        if index in {0, (self.proteinLength - 1)}:
            # Create a list of diagonals
            if self.plane == "2D":
                if xyzDif[0] == 1:
                    diagonals = [[otherCo[0], (otherCo[1] - 1)], [otherCo[0], (otherCo[1] + 1)]]
                else:
                    diagonals = [[(otherCo[0] + 1), otherCo[1]], [(otherCo[0] - 1), otherCo[1]]]
            else:
                # NOTE: of dit zoals regel 191(?)
                if xyzDif[0] == 1:
                    diagonals = [[otherCo[0], (otherCo[1] - 1), otherCo[2]],
                                 [otherCo[0], (otherCo[1] + 1), otherCo[2]],
                                 [otherCo[0], otherCo[1], (otherCo[2] - 1)],
                                 [otherCo[0], otherCo[1], (otherCo[2] + 1)]]
                elif xyzDif[1] == 1:
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

            # Planes in which the coordinate could lay
            if self.plane == "2D":
                # Only xy-plane
                planes = [[0, 1]]
            else:
                # xy-, xz-, and yz-planes
                planes = [[0, 1], [0, 2], [1, 2]]

            # Get all possible diagonals
            for plane in planes:
                # Places of the diagonal for the given plane
                options = [[(currentCo[plane[0]] + 1), (currentCo[plane[1]] + 1)],
                           [(currentCo[plane[0]] + 1), (currentCo[plane[1]] - 1)],
                           [(currentCo[plane[0]] - 1), (currentCo[plane[1]] + 1)],
                           [(currentCo[plane[0]] - 1), (currentCo[plane[1]] - 1)]]

                # Actually create diagonal and append to list with diagonals
                for option in options:
                    diagonal = copy.copy(currentCo)
                    diagonal[plane[0]] = option[0]
                    diagonal[plane[1]] = option[1]
                    diagonals.append(diagonal)

            # Randomly shuffle the diagonals
            random.shuffle(diagonals)

            # Return the first available diagonal and C
            for diagonal in diagonals:

                # Check if the diagonal is not occupied
                if diagonal not in self.occupied:

                    # Get the occupied surrounding coordinates of the diagonal
                    surroundCo = self.getSurroundCo(diagonal, True)

                    # Check if the coordinates of the next amino acid ar in the surrounding cordinates of the diagonal
                    if otherCo in surroundCo:
                        xyz = list(range(lengthCo))
                        difference = [(currentCo[i] - diagonal[i]) for i in range(lengthCo)]

                        # Initialize list for coordinates of the C
                        CCo = [None] * lengthCo
                        if self.plane == "3D":
                            sameCoIndex = difference.index(0)
                            CCo[sameCoIndex] = currentCo[sameCoIndex]
                            xyz.remove(xyz[sameCoIndex])

                        # Get the coordinates of the C
                        if xyzDif[xyz[0]] == 1:
                            CCo[xyz[0]] = currentCo[xyz[0]]
                            CCo[xyz[1]] = diagonal[xyz[1]]
                        else: #(elif xyzDif[xyz[1]] == 1)
                            CCo[xyz[0]] = diagonal[xyz[0]]
                            CCo[xyz[1]] = currentCo[xyz[1]]

                        # Check if C is not occupied or if they are the coordinates of the previous amino acid
                        if (CCo not in self.occupied) or (CCo == previousCo):
                            return [diagonal, CCo, nOrP[1]]


    def getSurroundCo(self, prevCo, occupied):
        """
        This method gets the 4 surrounding coordinates
        occupied is true if you want to know which surrounding coordinates are occupied
        occupied is false if you want to know which surrounding coordinates are not occupied
        """

        posCo = []

        # Create a list of all the possible coordinates
        coordinates = [[(prevCo[0] - 1), prevCo[1]], [(prevCo[0] + 1), prevCo[1]],
                       [prevCo[0], (prevCo[1] - 1)], [prevCo[0], (prevCo[1] + 1)]]
        if self.plane == "3D":
            [co.append(prevCo[2]) for co in coordinates]
            coordinates.append([prevCo[0], prevCo[1], (prevCo[2] - 1)])
            coordinates.append([prevCo[0], prevCo[1], (prevCo[2] + 1)])

        for coordinate in coordinates:
            # Only add coordinates to possible coordinates list depending on occupied
            # (whether the coordinates should be occupied or not)
            if (coordinate in self.occupied) is occupied:
                posCo.append(coordinate)

        return posCo


    def createAminoList(self):
        """
        This method folds the protein randomly
        Returns True if whole protein could be placed/ folded, False if not
        """
        # TODO: misschien kunnen we ook voorkomen dat een eiwit niet goed vouwt

        self.aminoList = []
        self.occupied = []

        # Loop over the letters in the proteinString
        for id in range(self.proteinLength):

            # Add amino acid to the aminoList
            self.aminoList.append(Amino(id, self.proteinString[id]))

            # Place the first and second amino-acid; the coordinates of the first
            # and second coordinate are (0,0) and (0,1) resp. (or (0,0,0), (0,1,0))
            if id in {0, 1}:
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
                    if idAround not in {(id + 1), (id - 1)}:
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
        Performs pull moves on the protein.
        Returns newly folded protein.

        Based on: Lesh, N., Mitzenmacher, M., & Whitesides, S. (2003).
        A complete and effective move set for simple protein folding. In
        Proceedings of the 7th annual international conference on research
        in computational molecular biology (RECOMB) (pp. 188–195). New York:
        ACM Press.
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
        if len(coordinates) == 3:

            previousAmino = self.aminoList[coordinates[2]]
            previousAminoCo = previousAmino.coordinate
            previousAminoId = previousAmino.id


            if previousAminoCo != coordinates[1]:

                newProtein.stabilityUpdate(previousAmino, True)

                newProtein.aminoList[coordinates[2]].coordinate = coordinates[1]
                newProtein.occupied[coordinates[2]] = coordinates[1]

                newProtein.stabilityUpdate(newProtein.aminoList[coordinates[2]])

            surCoPrev = newProtein.getSurroundCo(newProtein.aminoList[coordinates[2]].coordinate, True)

            dif = coordinates[2] - amino.id
            index = amino.id + (2*dif)
            newProtein.moveAminos(self, index, dif)

        return newProtein


    def moveAminos(self, oldProtein, idToMove, wayToMove):
        """
        Method for moving aminoacids to create a valid protein.
        oldProtein is the protein from which the new protein (neighbor) is created.
        idToMove is the id of the aminoacid that needs to be moved.
        wayToMove -- 1 if to end of protein; -1 if to beginning of protein
        """

        surCoPrev = self.getSurroundCo(self.aminoList[(idToMove - wayToMove)].coordinate, True)
        if ((wayToMove == -1 and idToMove < 0) or
            (wayToMove == 1 and (idToMove > (oldProtein.proteinLength - 1)))):
            return
        if (self.aminoList[idToMove].coordinate in surCoPrev):
            return

        self.stabilityUpdate(self.aminoList[idToMove], True)

        index = idToMove - 2*wayToMove
        self.aminoList[idToMove].coordinate = oldProtein.aminoList[index].coordinate
        self.occupied[idToMove] = oldProtein.aminoList[index].coordinate

        self.stabilityUpdate(self.aminoList[idToMove])

        self.moveAminos(oldProtein, (idToMove + wayToMove), wayToMove)
