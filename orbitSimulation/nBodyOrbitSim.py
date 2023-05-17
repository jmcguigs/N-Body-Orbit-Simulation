import ctypes
import numpy as np
from matplotlib import pyplot as plt


# c++ function calls

lib = ctypes.CDLL("orbitSimulationDLL.dll")

# a pointer to a Body/Simulation object
cppClassHandle = ctypes.POINTER(ctypes.c_int)

# functions to create simulated bodies of various types
# createStaticBody takes mass and radius, returns pointer to BodyStatic object.
lib.createStaticBody.argtypes = [ctypes.c_double, ctypes.c_double]
lib.createStaticBody.restype = cppClassHandle

# createSolvedBody takes mass and radius, returns pointer to BodySolved object.
lib.createSolvedBody.argtypes = [ctypes.c_double, ctypes.c_double]
lib.createSolvedBody.restype = cppClassHandle

# createSpacecraft takes radius, returns pointer to Spacecraft object.
lib.createSpacecraft.argtype = ctypes.c_double
lib.createSpacecraft.restype = cppClassHandle

# create a keplerian body
lib.createKeplerBody.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, cppClassHandle]
lib.createKeplerBody.restype = cppClassHandle

# functions to interact with simulated boedies
# set a body's state vector using keplerian orbital elements and a parent body handle
lib.setKeplerOrbit.argtypes = [cppClassHandle, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, cppClassHandle]
lib.setKeplerOrbit.restype = ctypes.c_void_p

# get a component of a target body's state vector at an index
lib.getStateVectorComponent.argtypes = [cppClassHandle, ctypes.c_int]
lib.getStateVectorComponent.restype = ctypes.c_double

# get a component of a target body's state vector at an index
lib.setStateVectorComponent.argtypes = [cppClassHandle, ctypes.c_int, ctypes.c_double]
lib.setStateVectorComponent.restype = ctypes.c_void_p

# get the mass of a gravitating body
lib.getMass.argtype = cppClassHandle
lib.getMass.restype = ctypes.c_double

# get the mass of a gravitating body
lib.getTimeSinceEpoch.argtype = cppClassHandle
lib.getTimeSinceEpoch.restype = ctypes.c_double

# create a simulation object (pointer to c++ simulation object)
lib.createSimulation.argtype = cppClassHandle
lib.createSimulation.restype = cppClassHandle

# add a body to a simulation object
lib.addToSimulation.argtypes = [cppClassHandle, cppClassHandle]
lib.addToSimulation.restype = ctypes.c_void_p

# run a simulation forward a set number of time steps
lib.runSimulation.argtypes = [cppClassHandle, ctypes.c_int]
lib.runSimulation.restype = ctypes.c_void_p

# set a simulation's time resolution
lib.setTimeResolution.argtypes = [cppClassHandle, ctypes.c_double]
lib.setTimeResolution.restype = ctypes.c_void_p


# FUNCTIONS
def deg2rad(x: float) -> float:
    return np.pi * x / 180

# return mean anomaly from orbital period in seconds
def meanAnomalyFromPeriod(orbitalPeriod: float) -> float:
    return 2 * np.pi / orbitalPeriod


# python wrapper c++ Body class
class Body:
    def __init__(self, mass: float, radius: float) -> None:
        self.pointerToClass: cppClassHandle = None

    # set the position components of the body's state vector
    def setPosition(self, x: float, y: float, z: float) -> None:
        lib.setStateVectorComponent(self.pointerToClass, 0, x)
        lib.setStateVectorComponent(self.pointerToClass, 3, y)
        lib.setStateVectorComponent(self.pointerToClass, 6, z)

     # set the velocity components of the body's state vector
    def setVelocity(self, x: float, y: float, z: float) -> None:
        lib.setStateVectorComponent(self.pointerToClass, 1, x)
        lib.setStateVectorComponent(self.pointerToClass, 4, y)
        lib.setStateVectorComponent(self.pointerToClass, 7, z)

     # set the acceleration components of the body's state vector
    def setPosition(self, x: float, y: float, z: float) -> None:
        lib.setStateVectorComponent(self.pointerToClass, 2, x)
        lib.setStateVectorComponent(self.pointerToClass, 5, y)
        lib.setStateVectorComponent(self.pointerToClass, 8, z)

    # return a numpy ndarray of a body's x, y, z position
    def getPosition(self) -> np.ndarray:
        return np.array([lib.getStateVectorComponent(self.pointerToClass, 0),
                         lib.getStateVectorComponent(self.pointerToClass, 3),
                         lib.getStateVectorComponent(self.pointerToClass, 6)])

    # return a numpy ndarray of a body's x, y, z velocity
    def getVelocity(self) -> np.ndarray:
        return np.array([lib.getStateVectorComponent(self.pointerToClass, 1),
                         lib.getStateVectorComponent(self.pointerToClass, 4),
                         lib.getStateVectorComponent(self.pointerToClass, 7)])

    # return a numpy ndarray of a body's x, y, z acceleration
    def getAcceleration(self) -> np.ndarray:
        return np.array([lib.getStateVectorComponent(self.pointerToClass, 2),
                         lib.getStateVectorComponent(self.pointerToClass, 5),
                         lib.getStateVectorComponent(self.pointerToClass, 8)])

    # set the body's state vector using Keplerian orbital elements
    def setKeplerOrbit(self, lonAscendingNode: float, semiMajorAxis: float, eccentricity: float, argOfPeriapsis: float, inclination: float, meanAnomaly: float, time: float, parent: 'Body') -> None:
        lib.setKeplerOrbit(self.pointerToClass, lonAscendingNode, semiMajorAxis, eccentricity, argOfPeriapsis, inclination, meanAnomaly, time, parent.pointerToClass)


# subclass of Body that follows a rigid Keplerian orbit around a parent Body object
class BodyKeplerian(Body):
    def __init__(self, mass: float, radius: float, lonAscendingNode: float, semiMajorAxis: float, eccentricity: float, argOfPeriapsis: float, inclination: float, meanAnomaly: float, parent: Body) -> None:
        super().__init__(mass, radius)
        self.pointerToClass = lib.createKeplerBody(mass, radius, lonAscendingNode, semiMajorAxis, eccentricity, argOfPeriapsis, inclination, meanAnomaly, parent.pointerToClass)


# subclass of Body that updates its' state vector using gravitational forces from all other bodies in the simulation
class BodySolved(Body):
    def __init__(self, mass: float, radius: float) -> None:
        super().__init__(mass, radius)
        self.pointerToClass = lib.createSolvedBody(mass, radius)

# essentially a BodySolved object that has negligible mass and is not factored into BodySolved state updates.
class Spacecraft(Body):
    def __init__(self, radius: float) -> None:
        super().__init__(0.0, radius)
        self.pointerToClass = lib.createSpacecraft(radius)

# subclass of Body that experiences no motion- still contributes gravitational forces to BodySolved objects.
class BodyStatic(Body):
    def __init__(self, mass: float, radius: float) -> None:
        super().__init__(mass, radius)
        self.pointerToClass = lib.createStaticBody(mass, radius)


# warpper for c++ simulation class- propagates gravitating bodies forward in time
class Simulation:
    def __init__(self) -> None:
        self.pointerToClass: cppClassHandle = lib.createSimulation()

    # set the simulation's time resolution in seconds
    def setTimeResolution(self, timeResolution: float) -> None:
        lib.setTimeResolution(self.pointerToClass, timeResolution)

    # add a body to the simulation
    def addBody(self, bodyToAdd: Body) -> None:
        lib.addToSimulation(self.pointerToClass, bodyToAdd.pointerToClass)

    # step the simulation forward numberOfSteps time steps
    def run(self, numberOfSteps: int) -> None:
        lib.runSimulation(self.pointerToClass, numberOfSteps)
