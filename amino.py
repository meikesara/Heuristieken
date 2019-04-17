"""
Class for Amino-acids.

Meike, Janneke, Nicole
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


    def addCoordinate(self, coordinate):
        """
        Add coordinates to the amino-acid
        """
        self.coordinate = coordinate
