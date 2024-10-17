#ifndef randomx_h
#define randomx_h

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <cstring>
#include <iomanip>
#include <vector>
#include <math.h>
#include <time.h>
#include <random>

using namespace std;


// std::uniform_real_distribution<double> fxdist;
// std::uniform_int_distribution<int> ixdist;  
//double randf0and1(){ return fdist(mtRNG); }
double randf0and1();
int randIntCustom(int a, int b);


#endif