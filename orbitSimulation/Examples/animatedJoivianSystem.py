import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from nBodyOrbitSim import deg2rad, meanAnomalyFromPeriod, BodyKeplerian, BodySolved, BodyStatic, Spacecraft, Simulation


jovianSystem: Simulation = Simulation()
jovianSystem.setTimeResolution(36.0)  # 36-second time resolution in simulation

# model Jupiter as a static body in the center of the simulation
jupiter_mass: float = 1.8982e27
jupiter_radius: float = 69911e3
jupiter: BodyStatic = BodyStatic(jupiter_mass, jupiter_radius)
jovianSystem.addBody(jupiter)

# model Io as a Keplerian body
io_mass: float = 8.931938e23
io_radius: float = 1821e6
io_argOfPeriapsis: float = deg2rad(49.1)
io_meanAnomaly: float = deg2rad(330.9)
io_inclination: float = deg2rad(0.0)
io_lonAscendingNode: float = deg2rad(0.0)
io_semiMajorAxis: float = 421800e3
io_eccentricity: float = 0.004
io: BodyKeplerian = BodyKeplerian(io_mass, io_radius, io_lonAscendingNode, io_semiMajorAxis, io_eccentricity, io_argOfPeriapsis, io_inclination, io_meanAnomaly, jupiter)
jovianSystem.addBody(io)

# model Europa as a Keplerian body
europa_mass: float = 4.799844e22
europa_radius: float = 1561e3
europa_argOfPeriapsis: float = deg2rad(45.0)
europa_meanAnomaly: float = deg2rad(345.4)
europa_inclination: float = deg2rad(0.5)
europa_lonAscendingNode: float = deg2rad(184.0)
europa_semiMajorAxis: float = 671100e3
europa_eccentricity: float = 0.009
#europa = BodySolved(europa_mass, europa_radius)
#europa.setKeplerOrbit(europa_lonAscendingNode, europa_semiMajorAxis, europa_eccentricity, europa_argOfPeriapsis, europa_inclination, europa_meanAnomaly, 0.0, jupiter)
europa = BodyKeplerian(europa_mass, europa_radius, europa_lonAscendingNode, europa_semiMajorAxis, europa_eccentricity, europa_argOfPeriapsis, europa_inclination, europa_meanAnomaly, jupiter)
jovianSystem.addBody(europa)

# model Ganymede as a Keplerian body
ganymede_mass: float = 1.4819e23
ganymede_radius: float = 2634e3
ganymede_argOfPeriapsis: float = deg2rad(198.3)
ganymede_meanAnomaly: float = deg2rad(324.8)
ganymede_inclination: float = deg2rad(0.2)
ganymede_lonAscendingNode: float = deg2rad(58.5)
ganymede_semiMajorAxis: float = 1070400e3
ganymede_eccentricity: float = 0.001
ganymede: BodyKeplerian = BodyKeplerian(ganymede_mass, ganymede_radius, ganymede_lonAscendingNode, ganymede_semiMajorAxis, ganymede_eccentricity, ganymede_argOfPeriapsis, ganymede_inclination, ganymede_meanAnomaly, jupiter)
jovianSystem.addBody(ganymede)

# model Callisto as a Keplerian body
callisto_mass: float = 1.075938e23
callisto_radius: float = 2410e3
callisto_argOfPeriapsis: float = deg2rad(43.8)
callisto_meanAnomaly: float = deg2rad(87.4)
callisto_inclination: float = deg2rad(0.3)
callisto_lonAscendingNode: float = deg2rad(309.1)
callisto_semiMajorAxis: float = 1882700e3
callisto_eccentricity: float = 0.007
callisto: BodyKeplerian = BodyKeplerian(callisto_mass, callisto_radius, callisto_lonAscendingNode, callisto_semiMajorAxis, callisto_eccentricity, callisto_argOfPeriapsis, callisto_inclination, callisto_meanAnomaly, jupiter)
jovianSystem.addBody(callisto)


# spacecraft
satellite_radius:float = 2  # keep-out (collision) boundary
satellite: Spacecraft = Spacecraft(satellite_radius)
a: float = 6000e3
e: float = 0.1
w: float = deg2rad(45)
M: float = deg2rad(90)
i: float = deg2rad(10)
O: float = deg2rad(310)
satellite.setKeplerOrbit(O, a, e, w, i, M, 0.0, io)
jovianSystem.addBody(satellite)

iterations: int = 10000

jPosition: np.ndarray = np.array([[0], [0], [0]])

ioPositions: np.ndarray = np.zeros((iterations, 3))
europaPositions: np.ndarray = np.zeros((iterations, 3))
ganymedePositions: np.ndarray = np.zeros((iterations, 3))
callistoPositions: np.ndarray = np.zeros((iterations, 3))
satellitePositions: np.ndarray = np.zeros((iterations, 3))

# precompute positions
for i in range(iterations):
    ioPositions[i, :] = io.getPosition()
    europaPositions[i, :] = europa.getPosition()
    ganymedePositions[i, :] = ganymede.getPosition()
    callistoPositions[i, :] = callisto.getPosition()
    satellitePositions[i, :] = satellite.getPosition()

    # step forward 100 time steps - 3600s (1hr) total between graphics updates
    jovianSystem.run(100)


def update(i: int) -> None:
    ax.clear()

    # plot orbit traces
    ax.plot(ioPositions[:, 0], ioPositions[:, 1], ioPositions[:, 2], alpha = 0.5, linestyle = '-', label = 'Io')
    ax.plot(europaPositions[:, 0], europaPositions[:, 1], europaPositions[:, 2], alpha = 0.5, linestyle = '-', label = 'Europa')
    ax.plot(ganymedePositions[:, 0], ganymedePositions[:, 1], ganymedePositions[:, 2], alpha = 0.5, linestyle = '-',label = 'Ganymede')
    ax.plot(callistoPositions[:, 0], callistoPositions[:, 1], callistoPositions[:, 2], alpha = 0.5, linestyle = '-',label = 'Callisto')
    ax.plot(satellitePositions[:, 0], satellitePositions[:, 1], satellitePositions[:, 2], alpha = 0.5, linestyle = '-',label = 'Satellite')

    # plot current body positions
    ax.scatter(jPosition[0], jPosition[1], jPosition[2], label = 'Jupiter')
    ax.scatter(ioPositions[i, 0], ioPositions[i, 1], ioPositions[i, 2], alpha = 0.8, label = 'Io')
    ax.scatter(europaPositions[i, 0], europaPositions[i, 1], europaPositions[i, 2], alpha = 0.8, label = 'Europa')
    ax.scatter(ganymedePositions[i, 0], ganymedePositions[i, 1], ganymedePositions[i, 2], alpha = 0.8, label = 'Ganymede')
    ax.scatter(callistoPositions[i, 0], callistoPositions[i, 1], callistoPositions[i, 2], alpha = 0.8, label = 'Callisto')
    ax.scatter(satellitePositions[i, 0], satellitePositions[i, 1], satellitePositions[i, 2], alpha = 0.8, label = 'Satellite')

    # set equal axis scales
    ax.set_xlim([-2e9, 2e9])
    ax.set_ylim([-2e9, 2e9])
    ax.set_zlim([-2e9, 2e9])

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

ani = FuncAnimation(fig, update, frames=range(1000))
plt.show()
