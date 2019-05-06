from protein import Protein
from amino import Amino
import copy
import visualizer
import matplotlib.pyplot as plt

bestFolded=[]
minStability=[0]


def constructive(proteinString): #, protein, validFolds=[]
    """
    """

    # Initialize protein
    protein = Protein(proteinString)
    protein.aminoList.append(Amino(0, protein.proteinString[0].upper()))
    protein.aminoList[0].addCoordinate([0, 0])
    protein.occupied.append([0, 0])
    protein.aminoList.append(Amino(1, protein.proteinString[1].upper()))
    protein.aminoList[1].addCoordinate([0, 1])
    protein.occupied.append([0, 1])

    createFolded(protein, 2)
    print(proteinString)
    # print(bestFolded[0].stability)
    # for bestFold in bestFolded:
    #     print(bestFold)
    #     visualizer.plotProtein(bestFold)
    print(minStability[0])

def createFolded(protein, idToMove):

    lengthProtein = len(protein.proteinString)
    if idToMove > (lengthProtein - 1):
        return

    prevCo = protein.aminoList[(idToMove - 1)].coordinate
    possibleCos = protein.getSurroundCo(prevCo, False)
    try:
        protein.aminoList[idToMove] = Amino(idToMove, protein.proteinString[idToMove].upper())
    except:
        protein.aminoList.append(Amino(idToMove, protein.proteinString[idToMove].upper()))

    # Heel misschien maakt dit het (net ietsje) beter
    if not possibleCos:
        return
    print("possibleCos = ", possibleCos)
    print("idToMove = ", idToMove)
    for possibleCo in possibleCos:
        print("possibleCo = ", possibleCo)
        print(protein)
        # Reset stability voor nieuwe vouwing
        protein.stability = 0
        protein.aminoList[idToMove].addCoordinate(possibleCo)
        try:
            protein.occupied[idToMove] = possibleCo
        except:
            protein.occupied.append(possibleCo)
        print(protein)
        print()
        if idToMove == (lengthProtein - 1):
            getStability(protein)

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
        visualizer.plotProtein(protein)
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
    constructive("HHPP")

    """
    HHPHHHPH
    HHPHHHPHPHHHPH
    HPHPPHHPHPPHPHHPPHPH
    PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP

    HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH
    PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP
    CPPCHPPCHPPCPPHHHHHHCCPCHPPCPCHPPHPC
    HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH
    HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH
    """
