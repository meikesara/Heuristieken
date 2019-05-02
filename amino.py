"""
Class for Amino-acids.

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
        """
        self.id = id
        self.type = type
        self.coordinate = []


    def __str__(self):
        """
        Return string with id, type, and coordinate of amino acid.
        """

        return ("id:" + str(self.id) + ", type:" + str(self.type) +
                ", coordinate:" + str(self.coordinate))


    def addCoordinate(self, coordinate):
        """
        Add coordinates to the amino-acid
        """
        # TODO: Check of we deze method überhaupt wel nodig hebben
        self.coordinate = coordinate
