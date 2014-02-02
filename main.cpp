/*
 * main.cpp
 *
 *  Created on: 30 Jan 2014
 *      Author: genway
 */

#include <cmath>
#include <iostream>
#include <fstream>

#include "defs.hpp"
#include "PdeGrid.hpp"

#define DELTA 0.3
#define RADIUS 15.0
#define ALPHA 5.0


double constantTerm(double X,double Y){
	return (- DELTA * Y / RADIUS * sin(X) ) + ALPHA;
}

double xDerivative(double X, double Y){
	return (DELTA * Y / RADIUS * cos(X) );
}

double yDerivative(double X, double Y){
	return (2 * DELTA * RADIUS * sin(X) ) + ALPHA * (Y);
}

double xxDerivative(double X, double Y){
	return 0.1;
}

double yyDerivative(double X, double Y){
	return 0.5;
}



int main(int argc, char *argv[]){

	double (*functions[5]) (double,double);

	functions[0]=&constantTerm;
	functions[1]=&xDerivative;
	functions[2]=&yDerivative;
	functions[3]=&xxDerivative;
	functions[4]=&yyDerivative;

	std::cout << (*functions[0])(2.0,5.0) << std::endl;
	std::cout << functions[1](2.0,5.0) << std::endl;
	std::cout << functions[4](2.0,5.0) << std::endl;

	PdeGrid * ring = new PdeGrid(functions,0,2*3.14159,-3.0,3.0,50,50,0.002,true);

	std::cout << ring->func << std::endl;
	std::cout << std::endl;

	for (int i = 0; i<10000; ++i){

		ring->iterateTimeSweep();

	}

	std::cout << ring->func << std::endl;
	std::cout << std::endl;

	std::ofstream file("steadystate.dat");
	file << ring->func;

	return 0;
}


