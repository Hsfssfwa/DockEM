#ifndef Rotbuilder_h
#define Rotbuilder_h

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string.h>
#include <cstring>
#include <iomanip>
#include <vector>
#include <math.h>
#include <time.h>
#include "GeometryTools.h"
#include "readpdb.h"

using namespace std;

//#define PI 3.1415926
struct Rotm
{
	vector<vector<float> > Rotcd;
	vector<string > Rotnm;
	vector<float> Ocd;
};

vector<float> calculateCoordinates(vector<float>& refA,vector<float>& refB,vector<float>& refC,float L,float ang,float di);

Rotm makeGly(vector<float> N, vector<float> CA, vector<float> C);
Rotm makeAla(vector<float> N, vector<float> CA, vector<float> C);
Rotm makeSer(vector<float> N, vector<float> CA, vector<float> C);
Rotm makeCys(vector<float> N, vector<float> CA, vector<float> C);
Rotm makeVal(vector<float> N, vector<float> CA, vector<float> C);
Rotm makeIle(vector<float> N, vector<float> CA, vector<float> C);
Rotm makeLeu(vector<float> N, vector<float> CA, vector<float> C);  
Rotm makeThr(vector<float> N, vector<float> CA, vector<float> C);
Rotm makeArg(vector<float> N, vector<float> CA, vector<float> C);
Rotm makeLys(vector<float> N, vector<float> CA, vector<float> C);
Rotm makeAsp(vector<float> N, vector<float> CA, vector<float> C);
Rotm makeAsn(vector<float> N, vector<float> CA, vector<float> C);
Rotm makeGlu(vector<float> N, vector<float> CA, vector<float> C);
Rotm makeGln(vector<float> N, vector<float> CA, vector<float> C);
Rotm makeMet(vector<float> N, vector<float> CA, vector<float> C);
Rotm makeHis(vector<float> N, vector<float> CA, vector<float> C);
Rotm makePro(vector<float> N, vector<float> CA, vector<float> C);
Rotm makePhe(vector<float> N, vector<float> CA, vector<float> C);
Rotm makeTyr(vector<float> N, vector<float> CA, vector<float> C);
Rotm makeTrp(vector<float> N, vector<float> CA, vector<float> C);

Rotm CalculateRotm(string AA,vector<float> N,vector<float> CA,vector<float> C);
#endif
