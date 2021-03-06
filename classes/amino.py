"""
Class for amino acids.

This class is used as a building block of objects of the Protein class.

Meike Kortleve, Nicole Jansen
"""


class Amino(object):
    """
    Amino class containing necessary attributes and methods for creating an
    amino acid.
    """

    def __init__(self, id, type):
        """
        Initializes an amino acid.

        Arguments:
        id -- integer, position of the amino acid in a protein
        type -- string, should be either H/P/C
        """

        self.id = id
        self.type = type
        self.coordinate = []


    def __str__(self):
        """
        Returns string with id, type, and coordinate of amino acid.s
        """

        return ("id:" + str(self.id) + ", type:" + str(self.type) +
                ", coordinate:" + str(self.coordinate))
