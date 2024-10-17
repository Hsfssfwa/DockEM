/*************************************************************************
Improvements:
Inline Functions should place in head file. 2007.5.11.
BUG: PointsOnLine();
c3[0]=c2[0]-(c1[0]-c2[0])*dis32/dis12;=>c3[0]=c2[0]-(c1[0]-c2[0])*dis32/dis12;
Date 2008.7.11.
Author: Cao Yang
Date:2008.7.11.

Improvements:
Add others Functions
Author: Biao zhang
Data:2018.3.5 
*************************************************************************/

#ifndef Geometryh
#define Geometryh

// C++ headers
#include <vector>
#include <cmath>
#include <iostream>
#include <cstdlib> 
using namespace std;


const float PI=3.141592653589793;
const float Extra=1.0e-4;
const float UpMax=1.0e+10;
typedef vector<vector<float> > matrix;


void ShowMyvector(const vector<float> &cc);

//inline void 
//Functions using vectors
void crossproduct(const vector<float> &c1, const vector<float> &c2, vector<float> &cc);

float innerproduct(const vector<float> &c1, const vector<float> &c2);

float innerproduct_univ(const vector<float> &c1, const vector<float> &c2);

void vectorsum(const vector<float> &c1, const vector<float> &c2, vector<float> &cc);

void subtract_univ(const vector<float> &c1, const vector<float> &c2, vector<float> &cc);

bool norm(const vector<float> &c, vector<float> &cc);

bool norm_univ(const vector<float> &c, vector<float> &cc);

bool norm_warnless(vector<float> &c);

bool norm_univ_warnless(vector<float> &c);

float vectorlength(const vector<float> &c);

float vectorlength_univ(const vector<float> &c);

float vectorlength2(const vector<float> &c);

float vectorlength2_univ(const vector<float> &c);

void multi(float coefficient, vector<float> &c, vector<float> &cc);

void multi_univ(float coefficient, vector<float> &c, vector<float> &cc);

void multi(float coefficient, vector<float> &cc);

void multi_univ(float coefficient, vector<float> &cc);

float deg2rad(float deg);

float rad2deg(float rad);

float Points2Distance2(const vector<float> &c1, const vector<float> &c2);

float Points2Distance(const vector<float> &c1, const vector<float> &c2);

//Return the angle of c1-c2-c3. Unit: radian
//Angle <c1c2c3 
float Points2Angle(const vector<float> &c1, const vector<float> &c2, const vector<float> &c3);
float Points2Dihedral(const vector<float> &c1, const vector<float> &c2, const vector<float> &c3, const vector<float> &c4);

//////////////////////////////////Matrix Tool/////////////////////////////////

//Method to initialize a matrix
void SetMatrix(matrix &sm, int m, int n);

void MatrixTimesMatrix(const matrix &a, const matrix &v, matrix &x, int m, int n, int l);

bool TransVectorTimesVector(const vector<float> &trans, const vector<float> &vctor, matrix &mtx);

bool MatrixTimesTransVector(const matrix &mtx, const vector<float> &tvt, vector<float> &vct);

void RealTimesMatrix(float at, const matrix &mx, matrix &mc);

bool MatrixAddMatrix(const matrix &ma, const matrix &mb, matrix &mc);

bool norm2(const vector<float> &c, vector<float> &cc);

bool norm2_warnless(vector<float> &c);

/************************Rotation Matrix****************************
* Give the Rotation Axis and the Rotation Angle
* Generate the Rotation Matrix
* Author: C.Y.
* Date 2006.11.30.
* BUG: romtx[0][2]=0 should be romtx[3][2]=0
*************************End***************************************/
bool RotationMatrixA(const vector<float> &axis, float angle, matrix &romtx);

/************************Rotation Matrix****************************
* Give the Rotation Axis and the Rotation Angle
* Generate the Rotation Matrix
* Author: C.Y.
* Date 2006.12.6.
*************************End***************************************/
bool RotationMatrixB(const vector<float> &axis, float angle, matrix &romtx);

/***************************CoordinateRotation********************
* Give the Rotation Axis and the Rotation Angle for a PointA
* Generate the coordinate of PointB after operation
* Author: C.Y.
* Date 2006.12.2.  
*****************************END**********************************/
bool CoordinateRotation(vector<float> &pointA, const vector<float> &axis, float angle, vector<float> &pointB);

/***************************CoordinateRotation********************
* For a PointA, Given the Coordinates of 2 points(axisA, axisB) 
* in line of Rotation Axis and the Rotation Angle 
* Generate the coordinate of PointB after operation
* Author: C.Y.
* Date 2006.12.2.  
*
*****************************END**********************************/
bool CoordinateRotation(const vector<float> &pointA, const vector<float> &axisA, const vector<float> &axisB, float angle, vector<float> &pointB);

/******************************************************************************
功能：用于大量空间点坐标的旋转变换
参数：axisA-axisB两个点构成空间旋转轴，angle是旋转角度，
     pointB是原来的空间点坐标，运算完毕后更新为新的坐标。
作者: CaoYang
时间：2008.6.14.
******************************************************************************/
bool GroupRotation(const vector<float> &axisA, const vector<float> &axisB, float angle, vector<vector<float> > &pointB);
//对指定部分进行旋转
bool GroupRotation(const vector<float> &axisA, const vector<float> &axisB, float angle, vector<vector<float> > &pointB, const vector<short>& index);


/******************************************************************************
功能：用于大量空间点坐标的平移变换
参数：trans构成空间平移量，pointB是原来的空间点坐标，运算完毕后更新为新的坐标。
时间：2008.6.14.
******************************************************************************/
bool GroupTranslation(const vector<float> &trans, vector<vector<float> > &pointB);

/*****************************************************************************
功能：使用欧拉角对三个坐标轴（x, y, z）进行旋转操作。
参数：phi,theta, psi分别为Z轴，新X轴，  
时间：2007.7.7.
******************************************************************************/
/*
bool EulerAngleApproach(const float phi, const float theta, const float psi, vector<vector<float>> &pointB);
*/


#endif
