from protein import Protein
from amino import Amino
import copy
import visualizer
import matplotlib.pyplot as plt
import visualizer
import time
import cProfile
import re

start_time = time.time()
minStability=[0]
lengthProtein = 0
# xTotal = []
# yTotal = []

def constructive(proteinString):
    """
    """

    # Initialize protein
    protein = Protein(proteinString)
    for i in range(2):
        protein.aminoList.append(Amino(i, protein.proteinString[i].upper()))
        protein.aminoList[i].addCoordinate([0, i])
        protein.occupied.append([0, i])
        # xTotal.append(0)
        # yTotal.append(i)

    createFolded(protein, 2)
    print(proteinString)
    print(minStability[0])


def createFolded(protein, idToMove):

    # NOTE: dit hoeft maar een keer.

    if idToMove > (lengthProtein - 1):
        return

    # NOTE: dit verwijderd de verdere eiwitten
    # TODO: if-statement toevoegen zodat dit niet altijd gebeurt, maar ik kan nu niet nadenken
    # if idToMove < len(protein.aminoList) - 1:
    del protein.aminoList[idToMove + 1:]
    del protein.occupied[idToMove + 1 :]
    # NOTE: dit is hetzelfde als de vorige 2 regels code
    # protein.aminoList = protein.aminoList[0:idToMove + 1]
    # protein.occupied = protein.occupied[0:idToMove + 1]

    prevCo = protein.aminoList[(idToMove - 1)].coordinate
    possibleCos = protein.getSurroundCo(prevCo, False)

    # # NOTE: volgens mij wordt hij hier alleen maar slomer van
    # if not possibleCos:
    #     return

    # NOTE: in een keer hele proteinstring naar upper
    try:
        protein.aminoList[idToMove] = Amino(idToMove, protein.proteinString[idToMove].upper())
    except:
        protein.aminoList.append(Amino(idToMove, protein.proteinString[idToMove].upper()))

    # NOTE: dit zorgt ervoor dat een aantal spiegelbeelden niet wordt gemaakt --> sneller
    # NOTE: dit hoeft niet elke keer hetest worden als het niet meer recht is.
    xTotal = sum([item[0] for item in protein.occupied])
    if xTotal == 0 and len(possibleCos) == 3:
        possibleCos.pop(1)

    # if sum(xTotal) == 0 and len(possibleCos) == 3:
    #     possibleCos.pop(1)



    for possibleCo in possibleCos:

        protein.aminoList[idToMove].addCoordinate(possibleCo)
        try:
            protein.occupied[idToMove] = possibleCo
            # xTotal[idToMove] = possibleCo[0]
            # yTotal[idToMove] = possibleCo[1]
        except:
            protein.occupied.append(possibleCo)
            # xTotal.append(possibleCo[0])
            # yTotal.append(possibleCo[1])

        # NOTE: dit kan nog op een slimmere manier
        xTotal = sum([item[0] for item in protein.occupied])
        yTotal = sum([item[1] for item in protein.occupied])
        # total =  abs(sum(xTotal)) + abs(sum(yTotal))
        total =  abs(xTotal) + abs(yTotal)

        if idToMove == (lengthProtein - 1) and total != lengthProtein * 2:
            getStability(protein)
            # visualizer.plotProtein(protein)

            if protein.stability < minStability[0]:
                minStability[0] = protein.stability
        createFolded(protein, (idToMove + 1))


def getStability(protein):
    # Reset stability voor nieuwe vouwing
    protein.stability = 0

    for amino in protein.aminoList:
        typeCo = amino.type
        if typeCo in {"H", "C"}:
            # Get surrounding coordinates of given amino acid that are occupied
            aroundCos = protein.getSurroundCo(amino.coordinate, True)

            # For each amino next to given amino, check if they create a bond
            for aroundCo in aroundCos:
                # Get id and type of amino acid next to given amino
                idAround = protein.occupied.index(aroundCo)
                aroundType = protein.aminoList[idAround].type

                if aroundType in {"H", "C"}:
                    # Check if amino is not connected in protein to amino
                    id = amino.id
                    if idAround != (id + 1) and idAround != (id - 1):
                        # Bond only needs to be plotted once
                        if id > idAround:
                            # Stronger bond created when both aminos are type C
                            if typeCo == "C" and aroundType == "C":
                                protein.stability -= 5

                            # Weaker bond created when at least one amino is type H
                            else:
                                protein.stability -= 1


if __name__ == "__main__":
    lengthProtein = len("PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP")
    # constructive("HPHPPHHPHPPHPHHPPHPH")
    constructive("PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP")
    # constructive("HHPHHHPH")

    # constructive("HHPHHHPHPHHHPH")

    # constructive("HHPHHHPHPHHHPH")
    # constructive("HHPHHHPHPHHHPH")
    print("--- %s seconds ---" % (time.time() - start_time))
    # cProfile.run('re.compile("constructive("HHPHHHPHPHHHPH")")')
