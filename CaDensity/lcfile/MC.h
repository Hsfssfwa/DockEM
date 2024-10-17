#ifndef INCLUDED_MC_H
#define INCLUDED_MC_H

#include <iostream>
#include <sstream>
#include <vector>
#include <cmath>
#include <time.h>
#include <stdlib.h>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "readpdb.h"
#include "GeometryTools.h"

using namespace std;

#define MX_PI 3.1415926

vector<float> Generate_random_point(float A,float B,float r);
float randomAB(float A,float B);
float Distance_point(vector<float> A,vector<float> B);
void MC_sampley();

#endif