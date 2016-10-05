/*
 * data.cpp
 *
 *  Created on: Oct 30, 2015
 *      Author: starside
 */
#include "data.h"
#include <cstring>

//Function Definitions
void printjVector(const jVector *vec){
	for(int i=0; i < 3; i++) {
		std::cout << vec->vector[i] << ",";
	}
	std::cout << "\n";
}

//prints a column major order (Fortran Style) matrix
void printjMatrix(const jMatrix *mat){
	int index;
	for(int row = 0; row < 3; row++){
		for(int col = 0; col < 3; col++){
			index = col*3 + row;
			std::cout << mat->matrix[index] << ",";
		}
		std::cout << "\n";
	}
}

void vectorAddPd(double *res, const double *a, const double *b, const double bsign) {
	res[0] = a[0] + bsign*b[0];
	res[1] = a[1] + bsign*b[1];
	res[2] = a[2] + bsign*b[2];
}

/*The return value is a location in memory.  I plan on rotating many
vectors, so not using the stack seems beneficial*/
void transformVectorP(jVector* res, const jMatrix *mat, const jVector *vec){
	//The LDA paramater is the stride.  In other words, how many
	//items per column.  It is used to map to linear memory
	cblas_dgemv(CblasColMajor, CblasNoTrans, 3, 3, 1.0, mat->matrix, 3,
		vec->vector, 1, 0, res->vector, 1);
}

/*The return value is a location in memory.  Perform matrix multiplication*/
void transformMatrixP(jMatrix *res, const jMatrix *mLeft, const jMatrix *mRight) {
	cblas_dgemm(	CblasColMajor, CblasNoTrans,
			CblasNoTrans, 3, 3,
			3, 1.0, mLeft->matrix,
			3, mRight->matrix, 3,
			0, res->matrix, 3);
}

//I am using Jose and Saletan p.528 for this definition
jMatrix eulerAngleZRotCW(const double theta){
	jMatrix result;
	/*
	cos t -sin t 0
	sin t  cos t 0
	0      0     1
	*/
	result.matrix[0] = cos(theta);
	result.matrix[1] = sin(theta);
	result.matrix[2] = 0;
	result.matrix[3] = -sin(theta);
	result.matrix[4] = cos(theta);
	result.matrix[5] = 0;
	result.matrix[6] = 0; result.matrix[7] = 0; result.matrix[8] = 1.0;
	return result;
}

//I am using Jose and Saletan p.528 for this definition
jMatrix eulerAngleXRotCW(const double theta){
	jMatrix result;
	/*
	1      0     0
	0    cos t  -sin t
	0    sin t   cos t
	*/
	result.matrix[0] = 1;result.matrix[1] = 0;result.matrix[2] = 0;
	result.matrix[3] = 0;
	result.matrix[4] = cos(theta);
	result.matrix[5] = sin(theta);
	result.matrix[6] = 0;
	result.matrix[7] = -1.0*sin(theta);
	result.matrix[8] = cos(theta);
	return result;
}

/*This is a body->lab frame transform
Pretty wasteful with the stack, but who cares.  Seems
a little safe than dynamic allocation
x_lab = M x_body.  To go the other direction, swap
psi and phi, then invert all signs.  Needs testing
*/
jMatrix eulerAngleMatrix(const double phi, const double theta, const double psi){
	jMatrix mPhi, mTheta, mPsi,mEuler,mTemp;
	mPhi   = eulerAngleZRotCW(phi);
	mTheta = eulerAngleXRotCW(theta);
	mPsi   = eulerAngleZRotCW(psi);
	transformMatrixP(&mTemp,  &mTheta, &mPsi);
	transformMatrixP(&mEuler, &mPhi,   &mTemp);
	return mEuler;
}


int testSomeRotations(){
	jVector test = { {0,0,1} };
	jVector res = { {6,6,6} };
	jMatrix resM;
	jMatrix rot1 = eulerAngleXRotCW(45*M_PI/180.0);
	jMatrix rot2 = eulerAngleZRotCW(45*M_PI/180.0);
	transformVectorP(&res, &rot1, &test);

	transformMatrixP(&resM,&rot2,&rot1);
	//printjVector(&res);
	std::cout << "----M1----\n";
	printjMatrix(&rot1);
	std::cout << "----M2----\n";
	printjMatrix(&rot2);
	std::cout << "---M1*M2---\n";
	transformMatrixP(&resM,&rot1,&rot2);
	printjMatrix(&resM);
	std::cout << "---M2*M1---\n";
	transformMatrixP(&resM,&rot2,&rot1);
	printjMatrix(&resM);
	return 0;
}


