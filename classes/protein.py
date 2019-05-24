"""
Class for protein.

This algorithm uses objects of the class Amino.

Meike Kortleve, Nicole Jansen
"""

import copy
import random
import matplotlib.pyplot as plt
from classes.amino import Amino


class Protein(object):
    """
    Protein class containing necessary attributes and methods for setting up
    and updating a protein.
    """

    def __init__(self, proteinString, plane):
        """
        Initializes the variables of the Protein object.

        Arguments:
        proteinString -- string that represents the protein
        plane -- either "2D" or "3D", determining in which dimension protein
                 should be placed/folded
        """

        self.proteinString = proteinString.upper()
        self.aminoList = []
        self.stability = 0
        self.occupied = []
        self.proteinLength = len(proteinString)
        self.plane = plane.upper()


    def __str__(self):
        """
        Prints the type of the amino acids in the protein and their coordinates.
        """

        output = ""
        for amino in self.aminoList:
            output += amino.type
            output += str(amino.coordinate) + " "

        return output


    def initializeProtein(self):
        """
        Initializes the protein by placing it in a straight line parallel to the
        y-axis.
        """

        self.aminoList = []
        self.occupied = []

        for id in range(self.proteinLength):
            self.aminoList.append(Amino(id, self.proteinString[id]))
            thisCoordinate = [0] * int(self.plane[0])
            thisCoordinate[1] = id
            self.aminoList[id].coordinate = thisCoordinate
            self.occupied.append(thisCoordinate)


    def stabilityUpdate(self, amino, replace=False):
        """
        Updates the stability of the protein.

        Argument:
        amino -- amino acid of class Amino that is just (re)placed

        Keyword argument:
        replace -- when True stability will be increased (worsened), when False
                   stability will be decreased (improved) (default False)
        """

        id = amino.id
        coordinate = amino.coordinate
        typeCo = amino.type

        # Bonds only created between amino acids of type H and C
        if typeCo in {"H", "C"}:
            aroundCos = self.getSurroundCo(coordinate, True)

            # For each amino next to given amino, check if they create a bond
            for aroundCo in aroundCos:
                idAround = self.occupied.index(aroundCo)
                aroundType = self.aminoList[idAround].type

                if aroundType in {"H", "C"}:
                    # Check if amino is not connected in protein to given amino
                    if idAround not in {(id + 1), (id - 1)}:
                        if typeCo == "C" and aroundType == "C":
                            if replace:
                                self.stability += 5
                            else:
                                self.stability -= 5
                        else:
                            if replace:
                                self.stability += 1
                            else:
                                self.stability -= 1


    def getSurroundCo(self, givenCo, occupied):
        """
        Returns a list of the surrounding coordinates that are either occupied
        or not.
        A surrounding coordinate is defined as having a difference of one in
        either the x-, y-, or z-coordinate compared to the coordinates of
        givenCo.

        Arguments:
        givenCo -- list of the coordinates (of an amino acid) of which to find
                   to surrounding coordinates
        occupied -- True if surrounding coordinates should be occupied, False
                    if the surrounding coordinates should not be occupied
        """

        possibleCos = []

        # Create a list of all the possible coordinates
        coordinates = [[(givenCo[0] - 1), givenCo[1]],
                       [(givenCo[0] + 1), givenCo[1]],
                       [givenCo[0], (givenCo[1] - 1)],
                       [givenCo[0], (givenCo[1] + 1)]]

        if self.plane == "3D":
            [co.append(givenCo[2]) for co in coordinates]
            coordinates.append([givenCo[0], givenCo[1], (givenCo[2] - 1)])
            coordinates.append([givenCo[0], givenCo[1], (givenCo[2] + 1)])

        for coordinate in coordinates:
            if (coordinate in self.occupied) is occupied:
                possibleCos.append(coordinate)

        return possibleCos


    def getDiagonalCo(self, currentAmino):
        """
        Determines coordinates diagonally to given amino acid.
        Returns a list of only the coordinates an unoccupied diagonal when given
        amino acids is the first or last in the protein; or, when the given
        amino acid is not first or last in protein, a list with an unoccupied
        diagonal, coordinates of the C (as defined by Lesh et al., 2003 (see
        pullMove)), and the index of the amino acid that should also be moved
        (either the next or previous amino acid); or an empty list when no empty
        diagonal could be found.

        Argument:
        currentAmino -- amino acid of class Amino for which a diagonal will be
                        determined
        """

        currentCo = currentAmino.coordinate
        index = self.aminoList.index(currentAmino)

        if (index + 1) == self.proteinLength:
            previousAmino = self.aminoList[index - 1]
            otherCo = previousAmino.coordinate
        elif index == 0:
            nextAmino = self.aminoList[index + 1]
            otherCo = nextAmino.coordinate
        else:
            nextOrPrev = [index + 1, index - 1]
            random.shuffle(nextOrPrev)

            nextAmino = self.aminoList[nextOrPrev[0]]
            otherCo = nextAmino.coordinate

            previousAmino = self.aminoList[nextOrPrev[1]]
            previousCo = previousAmino.coordinate

        # Calculate the absolute difference between the x, y (and z) coordinates
        lengthCo = len(currentCo)
        xyzDif = [abs(currentCo[i] - otherCo[i]) for i in range(lengthCo)]

        if index in {0, (self.proteinLength - 1)}:
            # Create list of diagonals
            if self.plane == "2D":
                if xyzDif[0] == 1:
                    diagonals = [[otherCo[0], (otherCo[1] - 1)],
                                 [otherCo[0], (otherCo[1] + 1)]]
                else:
                    diagonals = [[(otherCo[0] + 1), otherCo[1]],
                                 [(otherCo[0] - 1), otherCo[1]]]
            else:
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

            random.shuffle(diagonals)

            # Return the first diagonal that is not occupied
            for diagonal in diagonals:
                if diagonal not in self.occupied:
                    return [diagonal]

        else:
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
                opts = [[(currentCo[plane[0]] + 1), (currentCo[plane[1]] + 1)],
                        [(currentCo[plane[0]] + 1), (currentCo[plane[1]] - 1)],
                        [(currentCo[plane[0]] - 1), (currentCo[plane[1]] + 1)],
                        [(currentCo[plane[0]] - 1), (currentCo[plane[1]] - 1)]]

                for option in opts:
                    diagonal = copy.copy(currentCo)
                    diagonal[plane[0]] = option[0]
                    diagonal[plane[1]] = option[1]
                    diagonals.append(diagonal)

            random.shuffle(diagonals)

            for diagonal in diagonals:
                if diagonal not in self.occupied:
                    surroundCo = self.getSurroundCo(diagonal, occupied=True)

                    if otherCo in surroundCo:
                        xyz = list(range(lengthCo))
                        difference = [(currentCo[i] - diagonal[i])
                                      for i in range(lengthCo)]

                        CCo = [None] * lengthCo
                        if self.plane == "3D":
                            sameCoIndex = difference.index(0)
                            CCo[sameCoIndex] = currentCo[sameCoIndex]
                            xyz.remove(xyz[sameCoIndex])

                        # Get the coordinates of the C
                        if xyzDif[xyz[0]] == 1:
                            CCo[xyz[0]] = currentCo[xyz[0]]
                            CCo[xyz[1]] = diagonal[xyz[1]]
                        else:
                            CCo[xyz[0]] = diagonal[xyz[0]]
                            CCo[xyz[1]] = currentCo[xyz[1]]

                        if (CCo not in self.occupied) or (CCo == previousCo):
                            return [diagonal, CCo, nextOrPrev[1]]

        return []


    def pullMove(self):
        """
        Performs pull moves on the protein and makes sure that if an amino acid
        that is not the first or the last is chosen to move, that the other
        amino acids will be move as well in order to keep the validity of the
        folding of the protein.
        Returns newly folded protein.

        Based on: Lesh, N., Mitzenmacher, M., & Whitesides, S. (2003).
        A complete and effective move set for simple protein folding. In
        Proceedings of the 7th annual international conference on research
        in computational molecular biology (RECOMB) (pp. 188–195). New York:
        ACM Press.
        """

        amino = random.choice(self.aminoList)
        coordinates = self.getDiagonalCo(amino)

        # Make sure amino can be moved
        while not coordinates:
            amino = random.choice(self.aminoList)
            coordinates = self.getDiagonalCo(amino)

        newProtein = copy.deepcopy(self)
        newProtein.stabilityUpdate(amino, True)

        chosenCo = coordinates[0]
        newProtein.aminoList[amino.id].coordinate = chosenCo
        newProtein.occupied[amino.id] = chosenCo

        newProtein.stabilityUpdate(newProtein.aminoList[amino.id])

        # Make sure validity of folding is preserved
        if len(coordinates) == 3:
            previousAmino = self.aminoList[coordinates[2]]
            previousAminoCo = previousAmino.coordinate
            previousAminoId = previousAmino.id

            if previousAminoCo != coordinates[1]:
                newProtein.stabilityUpdate(previousAmino, True)

                newProtein.aminoList[coordinates[2]].coordinate = coordinates[1]
                newProtein.occupied[coordinates[2]] = coordinates[1]

                newProtein.stabilityUpdate(newProtein.aminoList[coordinates[2]])

            coPrev = newProtein.aminoList[coordinates[2]].coordinate
            surCoPrev = newProtein.getSurroundCo(coPrev, True)

            wayToMove = coordinates[2] - amino.id
            index = amino.id + (2*wayToMove)
            newProtein.moveAminos(self, index, wayToMove)

        return newProtein


    def moveAminos(self, oldProtein, idToMove, wayToMove):
        """
        Moves amino acids to create a valid protein after moving a chosen amino
        acid diagonally (see pullMove; this function is created to only be
        used in combination with the pullMove function).

        Arguments:
        oldProtein -- protein from which the new protein (neighbor) is created
        idToMove -- id of the amino acid that needs to be moved
        wayToMove -- 1 if to end of protein; -1 if to beginning of protein
        """

        prevId = (idToMove - wayToMove)
        surCoPrev = self.getSurroundCo(self.aminoList[prevId].coordinate, True)
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
