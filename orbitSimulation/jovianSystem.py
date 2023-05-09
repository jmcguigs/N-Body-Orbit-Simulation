from orbitSimulation import Body, Spacecraft, Simulation
from matplotlib import pyplot as plt
import numpy as np

# simulation with 10-second time steps
jovianSystem: Simulation = Simulation(10)

# Jupiter - center of the simulation
jupiter: Body = Body()
jupiter.mass = 1.8982e27
jovianSystem.addBody(jupiter)

# Io
io: Body = Body()
io.setPosition(421700e3, 0, 0)
io.mass = 8.931938e23
io.setVelocity(0, 17.334e3, 0)
jovianSystem.addBody(io)

# Europa
europa: Body = Body()
europa.setPosition(670900e3, 0, 0)
europa.setVelocity(0, 13743.36, 0)
europa.mass = 4.799844e22
jovianSystem.addBody(europa)

# Ganymede
ganymede: Body = Body()
ganymede.setPosition(1070000e3, 0, 0)
ganymede.setVelocity(0, 10.88e3, 0)
ganymede.mass = 1.4819e23
jovianSystem.addBody(ganymede)

# Callisto
callisto: Body = Body()
callisto.setPosition(1885000e3, 0, 0)
callisto.setVelocity(0, 8.204e3, 0)
callisto.mass = 1.075938e23
jovianSystem.addBody(callisto)

# spacecraft
satellite: Spacecraft = Spacecraft()
satellite.setStateFromKepler(0, 8000e3, 0.2, np.pi, 0, 2*np.pi/90000, 0, io)
satellite.mass = 100
jovianSystem.addBody(satellite)

iterations: int = 10000
reduction: int = 100
redCount: int = int(iterations / reduction)

jXY: np.ndarray = np.zeros([2, redCount])
ioXY: np.ndarray = np.zeros([2, redCount])
sXY: np.ndarray = np.zeros([2, redCount])
ind: int = 0
for i in range(iterations):
    jovianSystem.step()
    if i % reduction == 0:
        ioXY[:, ind] = io.stateVector[[0, 3]]
        jXY[:, ind] = jupiter.stateVector[[0, 3]]
        sXY[:, ind] = satellite.stateVector[[0, 3]]
        ind += 1

plt.plot(ioXY[0, :], ioXY[1, :], ':g')
plt.scatter(jXY[0, :], jXY[1, :], marker='o')
plt.plot(sXY[0, :], sXY[1, :], '-b')
plt.show()

