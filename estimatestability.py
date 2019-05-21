"""
Function for estimating stability based on the protein as a string.

Nicole Jansen, Meike Kortleve
"""

from math import floor

def estimateStability(proteinString, plane="2D"):
    """
    Returns estimate of stability for protein string.

    Argument:
    proteinString -- string of the protein

    Keyword argument:
    plane -- plane in which protein is placed (either 2D or 3D) (default 2D)
    """

    lengthProtein = len(proteinString)

    allDict = {"evenList": [], "oddList": [], "evenListC": [], "oddListC": []}
    counter = dict()
    counterCC = dict()
    stability = 0

    # Determine on which indices are H and C's
    for aminoIndex in range(lengthProtein):
        amino = proteinString[aminoIndex]

        if amino in {"H", "C"}:
            # Determine if index is even or odd
            useListName = "evenList" if (aminoIndex % 2) == 0 else "oddList"

            allDict[useListName].append(aminoIndex)

            # Extra list for indices of the C's
            if amino == "C":
                allDict[(useListName + "C")].append(aminoIndex)

    allList = allDict["evenList"] + allDict["oddList"]

    # Determine stability for each amino
    for aminoIndex in allList:
        amino = proteinString[aminoIndex]

        # Stability bond only when one amino acid on even and other on odd index
        useListName = "oddList" if (aminoIndex % 2) == 0 else "evenList"
        useList = allDict[useListName]
        useListC = allDict[(useListName + "C")]

        # Need to be at least two aminos inbetween to be able to create bonds
        all = [i for i in useList if (i >= (aminoIndex + 3) or
                                      i <= (aminoIndex - 3))]
        lengthAll = len(all)

        # NOTE: deze regel heb ik nu 2-3 keer precies hetzelfde (functie?)
        # Determine maximum amount of bonds that can be created
        bonds = 3 if aminoIndex in {0, (lengthProtein - 1)} else 2
        if plane == "3D":
            bonds += 2

        stability, counter, counterCC = stabilityUpdate(all, amino, aminoIndex,
                                counter, counterCC, bonds, stability, useListC)

    stability = removeBonds(counter, stability, lengthProtein, 1)
    stability = removeBonds(counterCC, stability, lengthProtein, 5)

    # Stability is doubled because each bond is counted twice
    stability /= 2
    stability = floor(stability)

    print("proteinString ({!s}): {!r}".format(plane, proteinString))
    print(f"stability:  {stability}")


def stabilityUpdate(allList, aminoType, aminoIndex, counterDict, counterDictCC,
                    bonds, stability, useListC):
    """
    Updates the estimate of the stability.
    Returns updated stability and updated dictionaries with counters of bonds.

    Arguments:
    allList -- list with aminos with which the amino acid could make connections
    aminoType -- type of the amino acid (P, H, or C)
    aminoIndex -- place of the amino acid in the protein
    counterDict -- dictionary to count amount of HH/CH-bonds when maximum of
                   bonds can only be achieved by using all possible connections
    counterDictCC -- dictionary to count amount of CC-bonds when maximum of
                     bonds can only be achieved by using all possible connections
    bonds -- maximum amount of bonds that could be created for given amino acid
    stability -- current stability estimate
    useListC -- list of coordinates of C's with which bonds could be created
    """

    if aminoType == "C":
        strengthBond = 5

        # First use list with indices of C's, then with indices of H's
        useFirst = [i for i in useListC if (i >= (aminoIndex + 3) or
                                            i <= (aminoIndex - 3))]
        useSecond = allList

        lengthFirst = len(useFirst)
        lengthSecond = len(useSecond) - lengthFirst

        # First update counters for CC-bonds, then for CH-bonds
        counterDictFirst = counterDictCC
        counterDictSecond = counterDict
    else:
        strengthBond = 1

        # Only bonds of strength 1
        useFirst = allList
        lengthFirst = len(useFirst)
        counterDictFirst = counterDict

    if lengthFirst > bonds:
        # More aminos to create connections with than possible bonds, so use max bonds
        stability -= (bonds * strengthBond)
    else:
        # Use all possible connections
        stability -= (lengthFirst * strengthBond)

        # Count used connections
        for i in useFirst:
            try:
                counterDictFirst[i] += 1
            except:
                counterDictFirst[i] = 1

            # For C-aminos, first checked CC bonds, now check CH-bonds
            if aminoType == "C":
                if lengthSecond >= (bonds - lengthFirst):
                    stability -= (bonds - lengthFirst)
                else:
                    stability -= lengthSecond
                    for i in useSecond:
                        if i not in useFirst:
                            try:
                                counterDictSecond[i] += 1
                            except:
                                counterDictSecond[i] = 1

    if aminoType == "C":
        return stability, counterDictSecond, counterDictFirst
    else:
        return stability, counterDictFirst, counterDictCC


def removeBonds(counterDict, stability, lengthProtein, strengthBond, plane="2D"):
    """
    Remove excessive amount of bonds from stability.
    Returns the updated stability.
    """

    for aminoIndex in counterDict:
        bonds = 3 if aminoIndex in {0, (lengthProtein - 1)} else 2
        if plane == "3D":
            bonds += 2

        impossibleBonds = (counterDict[aminoIndex] - bonds)
        if impossibleBonds > 0:
            stability += (impossibleBonds * strengthBond)

    return stability



if __name__ == "__main__":
    proteinString = "HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH"
    # "HHPHHHPH"
    # "HHPHHHPHPHHHPH"
    # "HPHPPHHPHPPHPHHPPHPH"
    # "PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP"
    # "HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH"

    # "PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP"
    # "CPPCHPPCHPPCPPHHHHHHCCPCHPPCPCHPPHPC"
    # "HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH"
    # "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH"

    estimateStability(proteinString, plane="3D")
