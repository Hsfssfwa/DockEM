#ifndef vdwRadiiHead
#define vdwRadiiHead
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <cstring>
#include <cstdlib>
#include "readpdb.h"

using namespace std;

float GetVdwRadius(const string& atm);
float GetVdwEg(string& R1,string& R2,float Dij);

#endif