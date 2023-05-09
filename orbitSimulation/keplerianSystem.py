from orbitSimulation import Body, StaticBody, KeplerianBody, Spacecraft, Simulation
from matplotlib import pyplot as plt
import numpy as np


jovianSystem: Simulation = Simulation(10.0)

# center of simulation, modeled as static
jupiter: Body = StaticBody()
jupiter.mass = 1.8982e27
jovianSystem.addBody(jupiter)

# Io
io: Body = KeplerianBody()
io.argOfPeriapsis