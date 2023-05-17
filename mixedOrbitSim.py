import ctypes
import numpy as np


lib = ctypes.CDLL("C:\\Users\\jemcg\\source\\repos\\orbitSimulation\\Debug\\orbitSimulationDLL.dll")

simulatedBodyHandle = ctypes.POINTER(ctypes.c_int)

# createStaticBody takes mass and radius, returns pointer to BodyStatic object.
lib.createStaticBody.argtypes = [ctypes.c_double, ctypes.c_double]
lib.createStaticBody.restypes = simulatedBodyHandle

# get a component of a target body's state vector at an index
lib.getStateVectorComponent.argtypes = [simulatedBodyHandle, ctypes.c_int]
lib.getStateVectorComponent.restypes = ctypes.c_double


jupiterMass: float = 1.8982e27
jupiterRadius: float = 69911e3

jupiter = lib.createStaticBody(jupiterMass, jupiterRadius)

jupiterPosX = lib.getStateVectorComponent(jupiter, 0)
print(jupiterPosX)