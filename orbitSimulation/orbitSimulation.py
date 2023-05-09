import numpy as np

# global variables
gravitationalConstant: float = 6.6743e-11

# return mean anomaly from orbital period in seconds
def meanAnomalyFromPeriod(orbitalPeriod: float) -> float:
    return 2 * np.pi / orbitalPeriod


# a gravitating body in 3D space
class Body:
    def __init__(self) -> None:
        self.mass: float = 1
        self.type: string = "solved"

        # column vector of x, dx, ddx, y, dy, ddy, z, dz, ddz
        self.stateVector: np.ndarray = np.transpose(np.zeros(9))

    # set the positon component of the state vector
    def setPosition(self, x: float, y: float, z: float) -> None:
        self.stateVector[[0, 3, 6]] = [x, y, z]

    # set the velocity component of the state vector
    def setVelocity(self, x: float, y: float, z: float) -> None:
        self.stateVector[[1, 4, 7]] = [x, y, z]

    # set the acceleration component of the state vector
    def setAcceleration(self, x: float, y: float, z: float) -> None:
        self.stateVector[[2, 5, 8]] = [x, y, z]

    # get the gravitational force in 3 dimensions between bodies (force direction referenced to self)
    def getForceFrom(self, other: 'Body') -> np.ndarray:
        global gravitationalConstant

        # magnitude of position difference vector
        distance: float = np.sqrt(np.sum(np.square(other.stateVector[[0, 3, 6]] - self.stateVector[[0, 3, 6]])))
        if distance == 0:
            return np.zeros(3)

        # unit vector pointing from self to other
        direction: np.ndarray = (other.stateVector[[0, 3, 6]] - self.stateVector[[0, 3, 6]]) / distance

        return direction * (gravitationalConstant * self.mass * other.mass / distance**2)

    # get forces to/from all other bodies in list
    def getTotalForceFromBodies(self, others: list) -> np.ndarray:
        totalForce: np.ndarray = np.zeros(3)
        for otherBody in others:
            totalForce += self.getForceFrom(otherBody)

        return totalForce

    def getTotalAcceleration(self, others: list) -> None:
        self.stateVector[[2, 5, 8]] = self.getTotalForceFromBodies(others) / self.mass

    # update the body's state vector using the state transition matrix
    def updateState(self, stateTransitionMatrix: np.ndarray, time: float) -> None:
        self.stateVector = stateTransitionMatrix @ self.stateVector

    # set a body's state vector using a kepler orbit
    def setStateFromKepler(self, lonAscendingNode: float, semiMajorAxis: float, eccentricity: float, argOfPeriapsis: float,
                           inclination: float, meanAnomaly: float, time: float, parentBody: 'Body') -> None:

        mu: float = 6.67430e-11 * parentBody.mass
        anomaly: float = meanAnomaly + time * np.sqrt(mu / semiMajorAxis**3)
        M: float = 2 * np.arctan(np.tan(anomaly)) + np.pi   # limit to 0-2pi

        # initial value of eccentric anomaly
        E = M + eccentricity*(np.sin(M + eccentricity) + np.sin(M))
        correction: float = 1.0
        tolerance: float = 1e-10

        # solve for eccentric anomaly with Newton-Raphson
        while abs(correction) > tolerance:
            F = E - eccentricity*np.sin(E) - M
            dF = 1 - eccentricity*np.cos(E)
            if abs(dF) < tolerance:
                break

            correction = F/dF
            E = E - correction

        trueAnomaly: float = 2 * np.arctan2(np.sqrt(1 + eccentricity) * np.sin(E/2),
                                     np.sqrt(1 - eccentricity) * np.cos(E/2))

        radius: float = semiMajorAxis * (1 - eccentricity * np.cos(E))

        # position vector x/y/z
        o: np.ndarray = np.array([radius * np.cos(trueAnomaly), radius * np.sin(trueAnomaly), 0])
        rx: float = (o[0] * (np.cos(argOfPeriapsis) * np.cos(lonAscendingNode) - np.sin(argOfPeriapsis) * np.cos(inclination) * np.sin(lonAscendingNode))
                   - o[1] * (np.sin(argOfPeriapsis) * np.cos(lonAscendingNode) + np.cos(argOfPeriapsis)  *np.cos(inclination) * np.sin(lonAscendingNode)))

        ry: float = (o[0] * (np.cos(argOfPeriapsis) * np.sin(lonAscendingNode) + np.sin(argOfPeriapsis) * np.cos(inclination) * np.cos(lonAscendingNode))
                   + o[1] * (np.cos(argOfPeriapsis) * np.cos(inclination) * np.cos(lonAscendingNode) - np.sin(argOfPeriapsis) * np.sin(lonAscendingNode)))

        rz: float = o[0] * np.sin(argOfPeriapsis) * np.sin(inclination) + o[1] * np.cos(argOfPeriapsis) * np.sin(inclination)


        # velocity vector x/y/z
        do: np.ndarray = (np.sqrt(mu * semiMajorAxis) / radius) * np.array([-np.sin(E), np.sqrt(1 - eccentricity**2) * np.cos(E), 0])

        drx: float = (do[0] * (np.cos(argOfPeriapsis) * np.cos(lonAscendingNode) - np.sin(argOfPeriapsis) * np.cos(inclination) * np.sin(lonAscendingNode))
                      - do[1] * (np.sin(argOfPeriapsis) * np.cos(lonAscendingNode) + np.cos(argOfPeriapsis) * np.cos(inclination) * np.sin(lonAscendingNode)))

        dry: float = (do[0] * (np.cos(argOfPeriapsis) * np.sin(lonAscendingNode) + np.sin(argOfPeriapsis) * np.cos(inclination) * np.cos(lonAscendingNode))
                      + do[1] * (np.cos(argOfPeriapsis) * np.cos(inclination) * np.cos(lonAscendingNode) - np.sin(argOfPeriapsis) * np.sin(lonAscendingNode)))

        drz: float = do[0] * np.sin(argOfPeriapsis) * np.sin(inclination) + do[1] * np.cos(argOfPeriapsis) * np.sin(inclination)


        # center about parent body
        rx += parentBody.stateVector[0]
        ry += parentBody.stateVector[3]
        rz += parentBody.stateVector[6]

        drx += parentBody.stateVector[1]
        dry += parentBody.stateVector[4]
        drz += parentBody.stateVector[7]

        # update the body's position
        self.stateVector[[0, 3, 6]] = [rx, ry, rz]

        # update the body's velocity
        self.stateVector[[1, 4, 7]] = [drx, dry, drz]


# a body that follows a rigid Keplerian orbit- helpful to speed up simulations
class KeplerianBody(Body):
    def __init__(self) -> None:
        super().__init__()
        self.type = "keplerian"
        self.parentBody: Body = None
        self.lonAscendingNode: float = 0
        self.semiMajorAxis: float = 0
        self.eccentricity: float = 0
        self.argOfPeriapsis: float = 0
        self.inclination: float = 0
        self.meanAnomaly: float = 0

    def updateState(self, stateTransitionMatrix: np.ndarray, time: float) -> None:
        super().setStateFromKepler(self.lonAscendingNode, self.semiMajorAxis, self.eccentricity, self.argOfPeriapsis, self.inclination,
                                   self.meanAnomaly, time, self.parentBody)


# a stationary body that does not change with time- helpful for improved sim speed / central body in simulations
class StaticBody(Body):
    def __init__(self) -> None:
        super().__init__()
        self.type = "static"

    def updateState(self, stateTransitionMatrix: np.ndarray, time: float) -> None:
        return

# a body capable of maneuvering
class Spacecraft(Body):
    def __init__(self) -> None:
        super().__init__()
        self.type = "solved"
        self.pointing: np.ndarray = np.array([0, 0, 1])
        self.thrust: float = 0

    # set the spacecraft's pointing angle
    def setPointing(rollAngle: float, yawAngle: float) -> None:
        return

    # start burn of constant thrust (in Newtons) at current pointing angle
    def startBurn(self, thrust: float) -> None:
        self.thrust = thrust

    def stopBurn(self) -> None:
        self.thrust = 0

    # set the acceleration portion of the state vector using thrust and gravitational forces
    def getTotalAcceleration(self, others: list) -> None:
        gravAcceleration: np.ndarray = super().getTotalForceFromBodies(others) / self.mass
        self.stateVector[[2, 5, 8]] = gravAcceleration + (self.thrust * self.pointing / self.mass)


# N-body orbit simulation
class Simulation:
    def __init__(self, timeStep: float) -> None:
        # get the state transition matrix used to propagate state vectors
        self.stateTransitionMatrix: np.ndarray = self.stateTransition3DCA(timeStep)
        self.bodies: list = []
        self.timeStep: float = timeStep
        self.time: float = 0
        self.assumeSpacecraftNegligible: bool = True

    def addBody(self, body: Body) -> None:
        self.bodies.append(body)

    # step the simulation forward in time
    def step(self) -> None:
        bodiesAtPreviousStep: list = self.bodies
        time += self.timeStep

        # update accelerations in state vectors
        for body in self.bodies:
            # only compute forces/accelerations for solved bodies
            if body.type != "keplerian":
                body.getTotalAcceleration(bodiesAtPreviousStep)

            body.updateState(self.stateTransitionMatrix, time)
            

    # function to generate the state transition matrix for a 3D constant acceleration motion model
    def stateTransition3DCA(self, dt: float) -> np.ndarray:
        # 1D constant-acceleration state transition matrix
        ca1d: np.ndarray = np.array([[1, dt, dt**2/2],
                                        [0, 1, dt],
                                        [0, 0, 1]])

        # 3x3 zero matrix
        o3: np.ndarray = np.zeros([3, 3])

        # rows of 3D state transition matrix
        r1: np.ndarray = np.hstack((np.hstack((ca1d, o3)), o3))
        r2: np.ndarray = np.roll(r1, 3, axis=1)
        r3: np.ndarray = np.roll(r2, 3, axis=1)

        # full 3D state transition matrix
        ca3d: np.ndarray = np.vstack((np.vstack((r1, r2)), r3))

        return ca3d
