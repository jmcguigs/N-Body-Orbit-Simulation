#pragma once
#include <math.h>
#include <vector>
#include <algorithm>

constexpr double pi = 3.14159265;
constexpr double gravitationalConstant = 6.6743e-11;

// header file for orbitSimulation classes
namespace orbitSimulation
{
	// functions / utilities
	double meanAnomalyFromPeriod(double orbitalPeriod)
	{
		return 2 * pi / orbitalPeriod;
	}

	// convert an angle in degrees to one in radians
	double deg2rad(double angleDegrees)
	{
		return pi * angleDegrees / 180;
	}

	// classes

	// gravitating body parent class
	class Body
	{
	public:
		Body(void);
		~Body(void);

		// empty state vector (3D constant-acceleration)
		std::vector<double> stateVector{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		double mass;
		double radius;

		void setPosition(double x, double y, double z);
		virtual void updateState(std::vector<Body*> others, double dt);
		void setVelocity(double x, double y, double z);
		void setAcceleration(double x, double y, double z);
		//bool checkCollision(Body other);
		void setStateFromKepler(double lonAscendingNode, double semiMajorAxis, double eccentricity, double argOfPeriapsis,
			double inclination, double meanAnomaly, double time, Body *parentBody);
	};

	// full-solve gravitating body
	class BodySolved : public Body
	{
	public:
		BodySolved(void);
		void updateState(std::vector<Body*> others, double dt) override;
	};

	// body that follows a rigid keplerian orbit centered on a parent body
	class BodyKeplerian : public Body
	{
	public:
		BodyKeplerian(double lonAscendingNode, double semiMajorAxis, double eccentricity, double argOfPeriapsis,
			double inclination, double meanAnomaly, double time, Body* parentBody);
		// orbital parameters: 
		double timeSinceEpoch;
		double lonAscendingNode;
		double semiMajorAxis;
		double eccentricity;
		double argOfPeriapsis;
		double inclination;
		double meanAnomaly;
		Body *parentBody;

		void updateState(std::vector<Body*> others, double dt) override;
	};

	// a static body that does not move over time
	class BodyStatic : public Body
	{
	public:
		BodyStatic(void);
		void updateState(std::vector<Body*> others, double dt) override;
	};

	// a full-solve body that is treated as having negligible mass relative to other bodies
	class Spacecraft : public BodySolved
	{
	public:
		Spacecraft(void);
	};

	// simulation class to handle propagating multiple bodies
	class Simulation
	{
	public:
		std::vector<Body*> bodies;

		// bodies with negligible mass (not factored into body state updates)
		std::vector<Body*> negligibleBodies;
		int time = 0;
		double timeResolution;

		Simulation();
		void step(int numberOfSteps);
	};
}