"""
Class for Amino-acids.

Meike, Nicole
"""


class Amino(object):
    """
    """

    def __init__(self, id, type):
        """
        Initialise variables
        """
        self.id = id
        self.type = type
        self.coordinate = []


    def __str__(self):
        """
        """

        return ("id:" + str(self.id) + ", type:" + str(self.type) + ", coordinate:" + str(self.coordinate))


    def addCoordinate(self, coordinate):
        """
        Add coordinates to the amino-acid
        """
        self.coordinate = coordinate
