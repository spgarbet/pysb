from math import *

# Methods from self exporter to import????

class Shape:
    def __init__(self, dimension, volume=None, surface=None):
        self.dimension = dimension
        self.volume    = volume
        self.surface   = surface

# A generic area object
class Area(Shape):
    def __init__(self, area):
        Shape.__init__(self, 2, surface=area)

class Volume(Shape):
    def __init__(self, volume):
        Shape.__init__(self, 3, volume=volume)

class Sphere(Shape):
    def __init__(self,radius,dimension):
        if(not isinstance(radius,int)): raise TypeError("Spherical radius must be an int")
        n  = dimension * 1.0
        Cn = pow(pi, n/2.0) / gamma(n/2.0 + 1.0)
        Shape.__init__(self, dimension, volume=Cn * pow(radius, n), surface=n * Cn * pow(radius, n-1))
        self.radius = radius

class Hypercube(Shape):
    def __init__(self,side,dimension):
        if(not isinstance(side,int)): raise TypeError("Hypercube side length must be an int")
        n  = dimension * 1.0
        Shape.__init__(self, dimension, volume=pow(side, n), surface=2 * n * pow(side, n-1))
        self.side = side

class Geometry(object):
    def __init__(self, shape, location):
        if(not isinstance(location,list)): raise TypeError("location must be a list of parameters")
        for l in location:
            if(not isinstance(l,int)): raise TypeError("each location parameter must be an int")
        self.shape = shape
        self.location = location

class SphericalSurface(Geometry):
    def __init__(self, radius, location):
        Geometry.__init__(self, Sphere(radius, len(location)), location)
        self.shape.dimension = self.shape.dimension - 1

class SphericalSpace(Geometry):
    def __init__(self, radius, location):
        Geometry.__init__(self, Sphere(radius, len(location)), location)

class SquareSpace(Geometry):
    def __init__(self, side, location):
        Geometry.__init__(self, Hypercube(side, len(location)), location)

