"""
Function for estimating stability based on the protein as a string.

Nicole Jansen, Meike Kortleve
"""

from math import floor

def estimateStability(proteinString):
    """
    Returns estimate of stability for protein string
    """

    lengthProtein = len(proteinString)

    # Initialize variables
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
                # TODO: determine which will be faster
                allDict[(useListName + "C")].append(aminoIndex)
                # allDict["{:}{:}".format(useListName, "C")].append(aminoIndex)

    allList = allDict["evenList"] + allDict["oddList"]

    # Determine stability for each amino
    for aminoIndex in allList:
        amino = proteinString[aminoIndex]

        # If index amino is even only bonds with aminos at odd indices can be made and vice versa
        useListName = "oddList" if (aminoIndex % 2) == 0 else "evenList"
        useList = allDict[useListName]
        useListC = allDict[(useListName + "C")]

        # There should be at least two aminos between aminos to be able to create bonds
        all = [i for i in useList if i >= (aminoIndex + 3) or i <= (aminoIndex - 3)]
        lengthAll = len(all)

        # NOTE: deze regel heb ik nu 2-3 keer precies hetzelfde (functie?)
        # Determine maximum amount of bonds that can be created
        bonds = 3 if aminoIndex in {0, (lengthProtein - 1)} else 2

        # Update stability by appropriate amount
        stability, counter, counterCC = stabilityUpdate(all, amino, aminoIndex,
                                counter, counterCC, bonds, stability, useListC)

    # Remove impossible amount of bonds from stability
    stability = removeBonds(counter, stability, lengthProtein, 1)
    stability = removeBonds(counterCC, stability, lengthProtein, 5)

    # Stability is doubled because each bond is counted twice
    stability /= 2
    stability = floor(stability)

    # Print the proteinString and the stability estimate
    print("proteinString: {!r}".format(proteinString))
    print(f"stability:  {stability}")


def stabilityUpdate(allList, aminoType, aminoIndex, counterDict, counterDictCC, bonds, stability, useListC):
    """
    Returns updated stability and updated dictionaries with counters of type of
    bonds.
    allList is list with aminos with which the amino of the given aminoType at
    the given aminoIndex could make connections.
    counterDict and counterDictCC are dictionaries to count the amount of HH/CH
    and CC-bonds, respectively, of each amino if # TODO: hier uitleggen over wanneer wel niet (kom nu even niet op de woorden)
    bonds is the maximum amount of bonds that could be created for the given
    amino (at aminoIndex)
    stability is the current stability and will be updated
    useListC is list C aminos to determine which one can create connections
    """

    # Initialize variables
    if aminoType == "C":
        strengthBond = 5

        # First use list with indices of C's, then with indices of H's
        useFirst = [i for i in useListC if i >= (aminoIndex + 3) or i <= (aminoIndex - 3)]
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

    # Return updated stability, dictionary with count of HH- and CH-bonds and
    # dictionary with count of CC-bonds
    if aminoType == "C":
        return stability, counterDictSecond, counterDictFirst
    else:
        return stability, counterDictFirst, counterDictCC


def removeBonds(counterDict, stability, lengthProtein, strengthBond):
    """
    Remove excessive amount of bonds from stability.
    Returns the updated stability
    """

    for aminoIndex in counterDict:
        bonds = 3 if aminoIndex in {0, (lengthProtein - 1)} else 2

        impossibleBonds = (counterDict[aminoIndex] - bonds)
        if impossibleBonds > 0:
            stability += (impossibleBonds * strengthBond)

    return stability



if __name__ == "__main__":
    proteinString = "HHPHHHPHPHHHPH"
    # "HHPHHHPH"
    # "HHPHHHPHPHHHPH"
    # "HPHPPHHPHPPHPHHPPHPH"
    # "PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP"
    # "HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH"

    # "PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP"
    # "CPPCHPPCHPPCPPHHHHHHCCPCHPPCPCHPPHPC"
    # "HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH"
    # "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH"

    estimateStability(proteinString)
