# N-Body-Orbit-Simulation
### Classes &amp; functions for simulating n-body systems

This package is used for simulating how systems with multiple interacting gravitating bodies evolve over time. Since this is a fairly computationally-intensive task,
the heavy lifting is done with a compiled DLL from the c++ code in orbitSimulation.cpp. To allow quick development and experimentation, these classes and their methods
have Python wrappers in nBodyOrbitSim.py.

### Examples
There are two example files included in this repository, both of which simulate the Jovian system with a spacecraft interacting with Jupiter & its 4 largest moons.
