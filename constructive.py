from protein import Protein
from amino import Amino
import copy
import visualizer
import matplotlib.pyplot as plt
import visualizer
import time

start_time = time.time()

bestFolded=[]
minStability=[0]


def constructive(proteinString): #, protein, validFolds=[]
    """
    """

    # Initialize protein
    protein = Protein(proteinString)
    for i in range(2):
        protein.aminoList.append(Amino(i, protein.proteinString[i].upper()))
        protein.aminoList[i].addCoordinate([0, i])
        protein.occupied.append([0, i])

    createFolded(protein, 2)
    print(proteinString)
    # print(totalFolds)
    # print(bestFolded[0].stability)
    # for bestFold in bestFolded:
    #     print(bestFold)
    #     visualizer.plotProtein(bestFold)
    print(minStability[0])


def createFolded(protein, idToMove):

    lengthProtein = len(protein.proteinString)
    if idToMove > (lengthProtein - 1):
        return

    # NOTE: dit verwijderd de verdere eiwitten
    # TODO: if-statement toevoegen zodat dit niet altijd gebeurt, maar ik kan nu niet nadenken
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

    try:
        protein.aminoList[idToMove] = Amino(idToMove, protein.proteinString[idToMove].upper())
    except:
        protein.aminoList.append(Amino(idToMove, protein.proteinString[idToMove].upper()))

    # NOTE: dit zorgt ervoor dat een aantal spiegelbeelden niet wordt gemaakt --> sneller
    xTotal = sum([item[0] for item in protein.occupied])
    if xTotal == 0 and len(possibleCos) == 3:
        possibleCos.pop(1)


    for possibleCo in possibleCos:
        # Reset stability voor nieuwe vouwing
        protein.stability = 0

        protein.aminoList[idToMove].addCoordinate(possibleCo)
        try:
            protein.occupied[idToMove] = possibleCo
        except:
            protein.occupied.append(possibleCo)

        xTotal = sum([item[0] for item in protein.occupied])
        yTotal = sum([item[1] for item in protein.occupied])
        total =  abs(xTotal) + abs(yTotal)

        if idToMove == (lengthProtein - 1) and total != lengthProtein * 2:
            getStability(protein)
            # visualizer.plotProtein(protein)

            # # Determine current best stability
            # if not bestFolded:
            #     minStability = 0
            # else:
            #     minStability = bestFolded[0].stability

            # if protein.stability < minStability:
            #     # Clear list with best folds, before adding new best
            #     bestFolded.clear()
            #     bestFolded.append(copy.deepcopy(protein))
            #     minStability = protein.stability
            #     print(minStability)
            # if protein.stability == minStability:
            #     bestFolded.append(copy.deepcopy(protein))

            if protein.stability < minStability[0]:
                minStability[0] = protein.stability
        createFolded(protein, (idToMove + 1))


def getStability(protein):
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
                    if idAround != (amino.id + 1) and idAround != (amino.id - 1):
                        # Bond only needs to be plotted once
                        if amino.id > idAround:
                            # Stronger bond created when both aminos are type C
                            if typeCo == "C" and aroundType == "C":
                                protein.stability -= 5

                            # Weaker bond created when at least one amino is type H
                            else:
                                protein.stability -= 1


if __name__ == "__main__":
    constructive("HHPHHHPHPHHHPHHHHPP")
    # constructive("HHPHHHPHPHHHPH")
    # constructive("HHPHHHPH")
    # constructive("HHHHP")
    print("--- %s seconds ---" % (time.time() - start_time))
