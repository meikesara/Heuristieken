"""
Class for Amino-acids.

Meike, Janneke, Nicole
"""


class Amino(object):
    """
    """

    def __init__(self, id, type):
        """
        """
        self.id = id
        self.type = type
        self.coordinate = []


    def addCoordinate(self, coordinate):
        """
        """
        self.coordinate = coordinate
