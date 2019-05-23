"""
Class for amino acids.

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
        id -- position of the amino acid in a protein.
        type -- integer, should be either H/P/C.
        """
        self.id = id
        self.type = type
        self.coordinate = []


    def __str__(self):
        """
        Returns string with id, type, and coordinate of amino acid.
        """

        return ("id:" + str(self.id) + ", type:" + str(self.type) +
                ", coordinate:" + str(self.coordinate))
