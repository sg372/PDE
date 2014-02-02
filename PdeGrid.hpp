/*
 * PdeGrid.hpp
 *
 *  Created on: 29 Jan 2014
 *      Author: genway
 */

#ifndef PDEGRID_HPP_
#define PDEGRID_HPP_

#include "defs.hpp"

class PdeGrid {


public:
	PdeGrid(double (*[])(double,double), double, double, double, double,
			unsigned, unsigned, double, bool);
	virtual ~PdeGrid();

	DoubleMatrix func, func_temp;
	double X0, Y0;
	double X1, Y1;
	double hi, hj;

	double timeStep;

	unsigned iSize, jSize;

	double (*coeff[5])(double, double);

	bool pbcX;

	void iterateTimeSweep();

	void iterate_at_site(unsigned, unsigned);

	void iterate_pbc_X(unsigned);

	double X(unsigned);

	double Y(unsigned);

};

#endif /* PDEGRID_HPP_ */
