#include "GetVdwRadius.h"

using namespace std;

float GetVdwRadius(const string& atm)
{
    if(atm[0] =='H')  return 1.00; 
    if(atm[0] =='C')  return 1.80;
    if(atm[0] =='N')  return 1.65;
    if(atm[0] =='O')  return 1.40;
    if(atm[0] =='S')  return 1.85;
    if(atm[0] =='P')  return 1.90;
    else return 1.70; 
}

float GetVdwEg(string& R1,string& R2,float Dij)
{
	float A = GetVdwRadius(R1)/Dij;
	float B = GetVdwRadius(R2)/Dij;
	float E = 4.0*pow(A,12.0)-4.0*pow(B,6.0);
	return E;
}