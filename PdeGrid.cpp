/*
 * PdeGrid.cpp
 *
 *  Created on: 29 Jan 2014
 *      Author: genway
 */

#include "PdeGrid.hpp"

PdeGrid::PdeGrid(double (*f[5])(double,double), double Xmin, double Xmax, double Ymin,
		double Ymax, unsigned Xsteps, unsigned Ysteps, double timestep, bool Xpbc) {

	X0=Xmin;
	Y0=Ymin;
	X1=Xmax;
	Y1=Ymax;
	iSize=Xsteps;
	jSize=Ysteps;
	hi = (X1-X0) / Xsteps;
	hj = (Y1-Y0) / Ysteps;
	pbcX = Xpbc;

	func = DoubleMatrix::Zero(Xsteps,Ysteps)/(Ysteps*Ysteps);

	func_temp = DoubleMatrix::Zero(Xsteps,Ysteps)/(Ysteps*Ysteps);

	if(pbcX){
	func.block(0,1,Xsteps,Ysteps-2) =
			DoubleMatrix::Ones(Xsteps,Ysteps-2)/((Xsteps)*(Ysteps-2));
	//func(Xsteps/2,Ysteps/2)=1.0;
	}
	//func(Xsteps/2,Ysteps/2)=1.0;

	timeStep = timestep;
	for(unsigned i=0;i<5;++i) coeff[i]=f[i];

}

PdeGrid::~PdeGrid() {

}



/* Sweep across grid leaving the edges untouched */
void PdeGrid::iterateTimeSweep(){

	for (unsigned i=1; i<iSize-1; ++i){
		for (unsigned j=1; j<jSize-1; ++j){
			iterate_at_site(i,j);
		}
	}

	if(pbcX) {

		for (unsigned j=1; j<jSize-1; ++j){
			iterate_pbc_X(j);
		}

	}

	func = func_temp;

}





void PdeGrid::iterate_at_site(unsigned i, unsigned j){

	double newFuncVal = coeff[0](X(i),Y(j))*func(i,j) +
			coeff[1](X(i),Y(j))*(func(i+1,j)-func(i-1,j))/(2*hi) +
			coeff[2](X(i),Y(j))*(func(i,j+1)-func(i,j-1))/(2*hj) +
			coeff[3](X(i),Y(j))*(func(i+1,j)-2*func(i,j)+func(i-1,j))/(hi*hi)+
			coeff[4](X(i),Y(j))*(func(i,j+1)-2*func(i,j)+func(i,j-1))/(hj*hj);

	func_temp(i,j) = func(i,j) + newFuncVal*timeStep;
}


void PdeGrid::iterate_pbc_X(unsigned j){

	//Update points at iSize-1
	double newFuncVal1 = coeff[0](X(iSize-1),Y(j))*func(iSize-1,j) +
			coeff[1](X(iSize-1),Y(j))*(func(0,j)-func(iSize-2,j))/(2*hi) +
			coeff[2](X(iSize-1),Y(j))*(func(iSize-1,j+1)-func(iSize-1,j-1))/(2*hj) +
			coeff[3](X(iSize-1),Y(j))*(func(0,j)-2*func(iSize-1,j)+func(iSize-2,j))/(hi*hi)+
			coeff[4](X(iSize-1),Y(j))*(func(iSize-1,j+1)-2*func(iSize-1,j)+func(iSize-1,j-1))/(hj*hj);

	func_temp(iSize-1,j) = func(iSize-1,j) + newFuncVal1*timeStep;

	//Update points at 0
	double newFuncVal2 = coeff[0](X(0),Y(j))*func(0,j) +
			coeff[1](X(0),Y(j))*(func(1,j)-func(iSize-1,j))/(2*hi) +
			coeff[2](X(0),Y(j))*(func(0,j+1)-func(0,j-1))/(2*hj) +
			coeff[3](X(0),Y(j))*(func(1,j)-2*func(0,j)+func(iSize-1,j))/(hi*hi)+
			coeff[4](X(0),Y(j))*(func(0,j+1)-2*func(0,j)+func(0,j-1))/(hj*hj);

	func_temp(0,j) = func(0,j) + newFuncVal2*timeStep;

}





double PdeGrid::X(unsigned i){

	return X0 + i*hi;

}


double PdeGrid::Y(unsigned j){

	return Y0 + j*hi;

}
