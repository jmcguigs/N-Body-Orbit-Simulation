from orbitSimulation import Body, StaticBody, KeplerianBody, Spacecraft, Simulation, MonteCarloOrbitTrace, deg2rad, meanAnomalyFromPeriod
from matplotlib import pyplot as plt
import numpy as np

# set up a monte carlo orbit trace simulation with 10s time resolution
jovianSystem: Simulation = Simulation(10.0)

# center of simulation, modeled as static
jupiter: Body = StaticBody()
jupiter.mass = 1.8982e27
jupiter.radius = 69911e3
jovianSystem.addBody(jupiter)

# Io
io: Body = KeplerianBody()
io.mass = 8.931938e23
io.radius = 1821e6
io.argOfPeriapsis = deg2rad(49.1)
io.meanAnomaly = deg2rad(330.9)
io.inclination = deg2rad(0.0)
io.lonAscendingNode = deg2rad(0.0)
io.semiMajorAxis = 421800e3
io.eccentricity = 0.004
io.parentBody = jupiter
io.updateState(np.zeros([9, 9]), 0)   # initialize the state vector (empty 9x9 transition matrix since it's not needed)
jovianSystem.addBody(io)

# Europa
europa: Body = KeplerianBody()
europa.mass = 4.799844e22
europa.radius = 1561e3
europa.argOfPeriapsis = deg2rad(45.0)
europa.meanAnomaly = deg2rad(345.4)
europa.inclination = deg2rad(0.5)
europa.lonAscendingNode = deg2rad(184.0)
europa.semiMajorAxis = 671100e3
europa.eccentricity = 0.009
europa.parentBody = jupiter
europa.updateState(np.zeros([9, 9]), 0)   # initialize the state vector (empty 9x9 transition matrix since it's not needed)
jovianSystem.addBody(europa)

# Ganymede
ganymede: Body = KeplerianBody()
ganymede.mass = 1.4819e23
ganymede.radius = 2634e3
ganymede.argOfPeriapsis = deg2rad(198.3)
ganymede.meanAnomaly = deg2rad(324.8)
ganymede.inclination = deg2rad(0.2)
ganymede.lonAscendingNode = deg2rad(58.5)
ganymede.semiMajorAxis = 1070400e3
ganymede.eccentricity = 0.001
ganymede.parentBody = jupiter
ganymede.updateState(np.zeros([9, 9]), 0)   # initialize the state vector (empty 9x9 transition matrix since it's not needed)
jovianSystem.addBody(ganymede)

# Callisto
callisto: Body = KeplerianBody()
callisto.mass = 1.075938e23
callisto.radius = 2410e3
callisto.argOfPeriapsis = deg2rad(43.8)
callisto.meanAnomaly = deg2rad(87.4)
callisto.inclination = deg2rad(0.3)
callisto.lonAscendingNode = deg2rad(309.1)
callisto.semiMajorAxis = 1882700e3
callisto.eccentricity = 0.007
callisto.parentBody = jupiter
callisto.updateState(np.zeros([9, 9]), 0)   # initialize the state vector (empty 9x9 transition matrix since it's not needed)
jovianSystem.addBody(callisto)

# spacecraft
satellite: Body = Spacecraft()
satellite.mass = 100
a: float = 800000e3
e: float = 0.6
w: float = deg2rad(0)
M: float = meanAnomalyFromPeriod(90*60)
i: float = deg2rad(10)
O: float = deg2rad(90)

satellite.setStateFromKepler(O, a, e, w, i, M, 0, jupiter)
jovianSystem.addBody(satellite)

# simulate the system
iterations: int = 1000000
reduction: int = 50    # plot every 1 in N points
redCount: int = int(iterations / reduction)

jXY: np.ndarray = np.zeros([2, redCount])
ioXY: np.ndarray = np.zeros([2, redCount])
euXY: np.ndarray = np.zeros([2, redCount])
gaXY: np.ndarray = np.zeros([2, redCount])
caXY: np.ndarray = np.zeros([2, redCount])
sXY: np.ndarray = np.zeros([2, redCount])
ind: int = 0

for i in range(iterations):
    jovianSystem.step()
    if i % reduction == 0:
        jXY[:, ind] = jupiter.stateVector[[0, 3]]
        ioXY[:, ind] = io.stateVector[[0, 3]]
        euXY[:, ind] = europa.stateVector[[0, 3]]
        gaXY[:, ind] = ganymede.stateVector[[0, 3]]
        caXY[:, ind] = callisto.stateVector[[0, 3]]
        sXY[:, ind] = satellite.stateVector[[0, 3]]
        ind += 1

plt.plot(ioXY[0, :], ioXY[1, :], ':g', alpha=0.2)
plt.plot(euXY[0, :], euXY[1, :], ':g', alpha=0.2)
plt.plot(gaXY[0, :], gaXY[1, :], ':g', alpha=0.2)
plt.plot(caXY[0, :], caXY[1, :], ':g', alpha=0.2)
plt.scatter(jXY[0, :], jXY[1, :], marker='o', color='r', alpha=0.5)
plt.plot(sXY[0, :], sXY[1, :], '-b', alpha=0.5)
plt.show()
