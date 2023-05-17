#include "pch.h"
#include "orbitSimulation.h"
#include <iostream>

using namespace orbitSimulation;


// ----------------- gravitating body parent class --------------------------
orbitSimulation::Body::Body(void)
{
	this->stateVector = std::vector<double>(9, 0);
}

orbitSimulation::Body::~Body(void)
{
}

// set a Body's position in its state vector
void orbitSimulation::Body::setPosition(double x, double y, double z)
{
	this->stateVector[0] = x;
	this->stateVector[3] = y;
	this->stateVector[6] = z;
	return;
}

// base class state update- should always be overridden by a derived class
void orbitSimulation::Body::updateState(std::vector<Body*> others, double dt)
{
	std::cout << "Base class 'Body' updateState method called." << std::endl;
}

// set a Body's velocity in its state vector
void orbitSimulation::Body::setVelocity(double x, double y, double z)
{
	this->stateVector[1] = x;
	this->stateVector[4] = y;
	this->stateVector[7] = z;
	return;
}

// set a Body's acceleration in its state vector
void orbitSimulation::Body::setAcceleration(double x, double y, double z)
{
	this->stateVector[2] = x;
	this->stateVector[5] = y;
	this->stateVector[8] = z;
	return;
}

// set a Body's state vector using Keplerian orbital elements referenced to a parent body
void orbitSimulation::Body::setStateFromKepler(double lonAscendingNode, double semiMajorAxis, double eccentricity, double argOfPeriapsis,
	double inclination, double meanAnomaly, double time, Body *parentBody)
{
	double mu = gravitationalConstant * (*parentBody).mass;
	double anomaly = meanAnomaly + time * sqrt(mu / pow(semiMajorAxis, 3));
	double M = 2 * atan(tan(anomaly)) + pi;		// limit to 0-2pi

	// initial value for eccentric anomaly
	double E = M + eccentricity * (sin(M + eccentricity) + sin(M));
	double correction = 1.0;
	double tolerance = 1e-10;

	// solve for eccentric anomaly with Newton-Raphson
	while (abs(correction) > tolerance)
	{
		double F = E - eccentricity * sin(E) - M;
		double dF = 1 - eccentricity * cos(E);
		if (abs(dF) < tolerance)
		{
			break;
		}
		correction = F / dF;
		E = E - correction;
	}

	double trueAnomaly = 2 * atan2(sqrt(1 + eccentricity) * sin(E / 2), sqrt(1 - eccentricity) * cos(E / 2));
	double radius = semiMajorAxis * (1 - eccentricity * cos(E));

	// used for position vector x/y/z
	double ox = radius * cos(trueAnomaly);
	double oy = radius * sin(trueAnomaly);

	// we'll use the sin/cos of several parameters a lot, so compute them once here (trig functions are slow- repeat as few times as possible)
	double sinw = sin(argOfPeriapsis);	// "w" for lowercase omega
	double cosw = cos(argOfPeriapsis);

	double sinOmega = sin(lonAscendingNode);
	double cosOmega = cos(lonAscendingNode);

	double sini = sin(inclination);
	double cosi = cos(inclination);

	// calculate the position in 3 axes relative to the parent body
	double rx = ox * (cosw * cosOmega - sinw * cosi * sinOmega) - oy * (sinw * cosOmega + cosw * cosi * sinOmega);

	double ry = ox * (cosw * sinOmega + sinw * cosi * cosOmega) + oy * (cosw * cosi * cosOmega - sinw * sinOmega);
	double rz = ox * sinw * sini + oy * cosw * sini;

	// used for velocity vector x/y/z
	double dox = -(sqrt(mu * semiMajorAxis) / radius) * sin(E);
	double doy = (sqrt(mu * semiMajorAxis) / radius) * sqrt(1 - pow(eccentricity, 2)) * cos(E);

	// calculate the velocity in 3 axes relative to the parent body
	double drx = dox * (cosw * cosOmega - sinw * cosi * sinOmega) - doy * (sinw * cosOmega + cosw * cosi * sinOmega);
	double dry = dox * (cosw * sinOmega + sinw * cosi * cosOmega) + doy * (cosw * cosi * cosOmega - sinw * sinOmega);
	double drz = dox * sinw * sini + doy * cosw * sini;

	// center about the parent body
	rx += parentBody->stateVector[0];
	ry += parentBody->stateVector[3];
	rz += parentBody->stateVector[6];

	drx += parentBody->stateVector[1];
	dry += parentBody->stateVector[4];
	drz += parentBody->stateVector[7];

	// update the state vector (avoiding method calls for speed)
	this->stateVector[0] = rx;
	this->stateVector[3] = ry;
	this->stateVector[6] = rz;
	this->stateVector[1] = drx;
	this->stateVector[4] = dry;
	this->stateVector[7] = drz;

	return;
}

// ----------------- full-solve body derived class --------------------------

// solved body subclass
orbitSimulation::BodySolved::BodySolved(void)
{
	this->stateVector = std::vector<double>(9, 0);
}

// updates a Body's state vector using gravitational accelerations from all other bodies
void orbitSimulation::BodySolved::updateState(std::vector<Body*> others, double dt)
{
	// iterate over all other bodies in simulation

	// reset accelerations
	this->stateVector[2] = 0;
	this->stateVector[5] = 0;
	this->stateVector[8] = 0;

	for (auto other : others)
	{
		// don't calculate force between self and self
		if (other == this)
		{
			continue;
		}

		else if (other->mass == 0)
		{
			continue;
		}

		double dx = other->stateVector[0] - this->stateVector[0];
		double dy = other->stateVector[3] - this->stateVector[3];
		double dz = other->stateVector[6] - this->stateVector[6];

		double distance = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));

		// direction vector of force
		std::vector<double> direction = { dx / distance, dy / distance, dz / distance };

		// compute gravitational acceleration magnitude
		double Gm2divR = gravitationalConstant * other->mass / pow(distance, 2);	// compute this once instead of 3 times

		// update the state vector to include the accelerations calculated above
		this->stateVector[2] += Gm2divR * direction[0];
		this->stateVector[5] += Gm2divR * direction[1];
		this->stateVector[8] += Gm2divR * direction[2];
	}

	// perform 3D constant-acceleration motion using the total acceleration from all bodies
	double dtSquared = pow(dt, 2) / 2;	// compute this once (pow is slow!)

	// position along one axis: x1 = x0 + vx * dt + ax * (dt^2 / 2)
	this->stateVector[0] += this->stateVector[1] * dt + this->stateVector[2] * dtSquared;
	// valocity along one axis: v1 = v0 + ax * dt
	this->stateVector[1] += this->stateVector[2] * dt;

	// repeat for other two axes
	this->stateVector[3] += this->stateVector[4] * dt + this->stateVector[5] * dtSquared;
	this->stateVector[4] += this->stateVector[5] * dt;
	this->stateVector[6] += this->stateVector[7] * dt + this->stateVector[8] * dtSquared;
	this->stateVector[7] += this->stateVector[8] * dt;

	return;
}

// ----------------- keplerian body derived class --------------------------
orbitSimulation::BodyKeplerian::BodyKeplerian(double lonAscendingNode, double semiMajorAxis, double eccentricity, double argOfPeriapsis, double inclination, double meanAnomaly, double time, Body* parentBody)
{
	this->stateVector = std::vector<double>(9, 0);
	this->lonAscendingNode = lonAscendingNode;
	this->semiMajorAxis = semiMajorAxis;
	this->eccentricity = eccentricity;
	this->argOfPeriapsis = argOfPeriapsis;
	this->inclination = inclination;
	this->meanAnomaly = meanAnomaly;
	this->timeSinceEpoch = time;
	this->parentBody = parentBody;
	this->setStateFromKepler(this->lonAscendingNode, this->semiMajorAxis, this->eccentricity, this->argOfPeriapsis, this->inclination, this->meanAnomaly, 0.0, this->parentBody);
}

// keplerian body subclass
void orbitSimulation::BodyKeplerian::updateState(std::vector<Body*> others, double dt)
{
	// update the time since epoch
	this->timeSinceEpoch = this->timeSinceEpoch + dt;

	// use the parent class setStateFromKepler method to update the state vector at the current time
	this->setStateFromKepler(this->lonAscendingNode, this->semiMajorAxis, this->eccentricity, this->argOfPeriapsis,
		this->inclination, this->meanAnomaly, this->timeSinceEpoch, this->parentBody);

	return;
}


// ----------------- static body derived class --------------------------
orbitSimulation::BodyStatic::BodyStatic(void)
{
	this->stateVector = std::vector<double> (9, 0);
}

// static body subclass update state method
void orbitSimulation::BodyStatic::updateState(std::vector<Body*> others, double dt)
{
	// static, so the state never updates
	return;
}

// ----------------- spacecraft derived class --------------------------
orbitSimulation::Spacecraft::Spacecraft(void)
{
	this->stateVector = std::vector<double>(9, 0);
	this->mass = 0;
}

// ----------------- simulation class --------------------------
// simulation class to propagate multiple bodies simultaneously
orbitSimulation::Simulation::Simulation()
{
	this->bodies.clear();
	this->negligibleBodies.clear();
}

// step the simulation forward a given number of time steps
void orbitSimulation::Simulation::step(int numberOfSteps)
{
	int steps = 0;
	while (steps < numberOfSteps)
	{
		// update states of bodies using all others (excluding negligible-mass bodies)
		for (Body *currentBody : this->bodies)
		{
			currentBody->updateState(this->bodies, this->timeResolution);
		}

		// update states of negligible-mass bodies
		for (Body* currentBody : this->negligibleBodies)
		{
			currentBody->updateState(this->bodies, this->timeResolution);
		}

		steps++;
	}
	return;
}


// External interface to Python ctypes
extern "C"
{
	// return a pointer to a newly-created keplerian body
	__declspec(dllexport) Body* createKeplerBody(double mass, double radius, double lonAscendingNode, double semiMajorAxis, double eccentricity, double argOfPeriapsis,
		double inclination, double meanAnomaly, Body* pointerToParent)
	{
		Body* newBody = new BodyKeplerian(lonAscendingNode, semiMajorAxis, eccentricity, argOfPeriapsis, inclination, meanAnomaly, 0.0, pointerToParent);
		newBody->mass = mass;
		newBody->radius = radius;
		return newBody;
	}

	// return a pointer to a newly-created static
	__declspec(dllexport) Body* createStaticBody(double mass, double radius)
	{
		Body* newBody = new BodyStatic;
		newBody->mass = mass;
		newBody->radius = radius;
		return newBody;
	}

	// return a pointer to a newly-created solved body
	__declspec(dllexport) Body* createSolvedBody(double mass, double radius)
	{
		Body* newBody = new BodySolved;
		newBody->mass = mass;
		newBody->radius = radius;
		return newBody;
	}

	// return a pointer to a newly-created solved body
	__declspec(dllexport) Body* createSpacecraft(double radius)
	{
		Body* newBody = new Spacecraft;
		newBody->radius = radius;
		return newBody;
	}

	// set a body's state using keplerian orbital elements
	__declspec(dllexport) void setKeplerOrbit(Body *pointerToBody, double lonAscendingNode, double semiMajorAxis, double eccentricity, double argOfPeriapsis,
		double inclination, double meanAnomaly, double time, Body *pointerToParent)
	{
		(*pointerToBody).setStateFromKepler(lonAscendingNode, semiMajorAxis, eccentricity, argOfPeriapsis,
			inclination, meanAnomaly, time, pointerToParent);
	}

	// create a simulation
	__declspec(dllexport) Simulation* createSimulation(void)
	{
		Simulation* sim = new Simulation;
		return sim;
	}

	// add a body to a simulation
	__declspec(dllexport) void addToSimulation(Simulation *pointerToSimulation, Body *pointerToBody)
	{
		if (pointerToBody->mass <= 0)
		{
			pointerToSimulation->negligibleBodies.push_back(pointerToBody);
		}

		else 
		{
			pointerToSimulation->bodies.push_back(pointerToBody);
		}
		std::cout << "Successfully added 1 body to simulation." << std::endl;
	}

	// step the simulation forward a specified number of steps
	__declspec(dllexport) void runSimulation(Simulation* pointerToSimulation, int numberOfSteps)
	{
		pointerToSimulation->step(numberOfSteps);
	}

	// add a body to a simulation
	__declspec(dllexport) double getStateVectorComponent(Body* pointerToBody, int stateVectorIndex)
	{
		return pointerToBody->stateVector.at(stateVectorIndex);
	}

	// set a body's state vector component
	__declspec(dllexport) void setStateVectorComponent(Body* pointerToBody, int stateVectorIndex, double value)
	{
		pointerToBody->stateVector[stateVectorIndex] = value;
	}

	// return the mass of a gravitating body- mostly for debug purposes
	__declspec(dllexport) double getMass(Body* pointerToBody)
	{
		return pointerToBody->mass;
	}

	// set the time resolution of a simulation
	__declspec(dllexport) void setTimeResolution(Simulation* pointerToSimulation, double timeResolution)
	{
		pointerToSimulation->timeResolution = timeResolution;
	}

	// set the time resolution of a simulation
	__declspec(dllexport) double getTimeSinceEpoch(BodyKeplerian *pointerToBody)
	{
		return pointerToBody->timeSinceEpoch;
	}
};
