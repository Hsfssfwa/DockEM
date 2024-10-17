#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <string.h>
#include <random>

#include "mol.h"
#include "flex_dock.h"
#include "CaDensityx.h"
#if !defined(REMC_H)
    #define REMC_H 1

//bool REMC(AMBER_TYPER & amber, vector <DOCKMol> & mol,Energy_Score & c_nrg, string filevdw, string fileflex,string fileflex_drive_tbl,string grid_file,int cutoff,vector<vector <float> > &x,vector<vector <float> > &y, vector<vector <float> > &z,vector<float> &,float,float,string,mt19937 &);
bool sortmol(vector<DOCKMol> & vectormols, vector<int> & vectorindex);
bool sortmol_em(vector<DOCKMol>& vectormols, vector<int>& vectorindex);
bool sortmol_CC(vector<float>& vectormols, vector<int>& vectorindex);

vector<int> randomgeneratenum(int num,int maxvalue);

bool REMC(AMBER_TYPER & amber, vector <DOCKMol> & REmolvec, DOCKMol & receptor, Energy_Score & c_nrg,string filevdw, string fileflex,string fileflex_drive_tbl,string grid_file,int cutoff,vector<vector <float> > &x,vector<vector <float> > &y, vector<vector <float> > &z,vector<float> &accenergy,float Tmin, float Tmax,string outfile,mt19937 &gen, bool flexible_flag,vector<vector<int>> biolip_matrix,int swapnumber,int MC_steps,vector<float>sphx,vector<float>sphy,vector<float> sphz,float bindsite_weight,vector<int> protein_atom_index,vector<vector <float> > ave, vector<vector <float> > std, ElectronDensity theDensityMap, vector <DOCKMol> &best_model, vector<float> &best_E, vector<poseCoord> &binding_site, vector<float>& centers);

using namespace std;
#endif
