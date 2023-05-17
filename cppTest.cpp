#include "pch.h"
#include "cppHeader.h"

// generate a random number on a standard normal distribution
void RandomGenerator::getSample(void)
{
	std::random_device rd{};
	std::mt19937 gen{ rd() };
	std::normal_distribution<> d{ 0, 1 };
	this->previousValue = d(gen);
}