#ifndef readpdb_h
#define readpdb_h

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <cstring>
#include <iomanip>
#include <vector>
#include <math.h>
#include <time.h>
#include "GeometryTools.h"

using namespace std;

//#define PI 3.1415926


struct Atom
{
	int serial;
	string chara;
	string recna;
	string atona;  //atom name
	float tempfac;
	vector<float> xyzVect; //  Orthogonal coordinates for X,Y,Z in Angstroms.
};

struct Residue
{
//	bool het; // ture - HEATM, false -ATOM
	int resseq;  // residue sequence number
	string icode;  // inseriton code
	string resname; // residue name
	vector<Atom> heavyatoms; // nheavyatom
	vector<Atom> atoms; // list of atoms
};

struct Chain
{
	string chaid;  // chain ID
//	string sequence;  // sequence converted from CA coordinate
//	string ss; //secondary structure
	vector<Residue> residues; //list of residues
};

struct Model
{
	vector<Chain> chains; // list of chains
};

struct coordtype {
	vector<string> costr;   // before "ATOM"
	vector<string> recna;  // record name
	vector<int> serial;    // atom serial number
	vector<string> atona;  // atom name
	vector<string> chara;    // alternate location indicaator
	vector<string> resname; //residue name
	vector<string> chaid;  // chain ID
	vector<int> resseq;   // residue sequence number
	vector<string> icode; // code for insetion of residue
	vector<float> cx;   // Orthogonal coordinates for X in Angstroms.
	vector<float> cy;   //Orthogonal coordinates for Y in Angstroms.
	vector<float> cz;    // Orthogonal coordinates for X in Angstroms.
	vector<float> occ;   // occupancy
	vector<float> tempfac; //Temperature  factor
	vector<string> ele; //Element symbol, right-justified.
	vector<string> chag; //Charge  on the atom.
//	stirng endstr; // follow the "ATOM"
} ;

struct RT
{
	vector<vector<float> > Rotation;
	vector<float> Translation;
};

struct poseCoord{
	vector<float> x_;
	float B_;   // temperature factor (Default =0.0 in PDB file). this parmeters is corellation to B factor
	std::string elt_;
};

struct Rot
{
	vector<float> ca;
	vector<vector<float> > res_atm; // all sidechain atom in residue
	vector<string> resatmnm;
	vector<float> cn;	
	vector<float> cc;
	vector<float> co;

	vector<float> cb;
	vector<float> cg;
	vector<float> cg1;
	vector<float> cg2;
	vector<float> og;
	vector<float> sg;
	vector<float> cd;
	vector<float> cd1;
	vector<float> cd2;
	vector<float> sd;
	vector<float> nd1;
	vector<float> nd2;
	vector<float> od2;
	vector<float> ce;
	vector<float> oe2;
	vector<float> ne;
	vector<float> ne2;
	vector<float> nz;
	vector<float> cz;
};

typedef vector<poseCoord> poseCoords;
typedef vector<vector<float> > FMatrix;

string AAtoS(string animoacide);
string StoAA(string animoacide);
void readpdbstructure(const char *filename,Model &proteins);
void readpdbstructurex(const char *filename,Model &proteins);
void writePDBStructure(const char *filename,Model &protein);
void writePDBStructure(const char *filename,Model &protein,string outname);
void readcoordinatefrompdb(const char *filename,coordtype &coord);
void writePDB(const char *filename,coordtype &coordx);
void rotatepdb(coordtype &m);
poseCoords rotatePDB(poseCoords &pose);
float RandomDoubleX(double numx,double numy);
//int RandomDoubleX(int numx,int numy);
void readPDBcoords(std::string filename, poseCoords &atmlist);
std::string name2elt( std::string line );
vector<string> string_splitx(string const & in, char splitchar);
void writePDBcoords(std::string filename, poseCoords &atmlist);

bool GroupRotationx(const vector<float> &axisA, const vector<float> &axisB, float angle, vector<Rot> &pointB);
bool GroupRotationy(const vector<float> &axisA, const vector<float> &axisB, float angle, vector<Rot> &pointB);
// rotate the amount of point beginning from index 
bool GroupRotationx(const vector<float> &axisA, const vector<float> &axisB, float angle, vector<Rot> &pointB, int index);
// Translation
bool GroupTranslationx(const vector<float> &trans, vector<Rot> &pointB);
// translate the amount of point beginning from index 
bool GroupTranslationx(const vector<float> &trans, vector<Rot> &pointB, int index);
bool GroupRotationt(const vector<float> &axisA, const vector<float> &axisB, float angle, vector<float> &pointB);
bool GroupRotationid(const vector<float> &axisA, const vector<float> &axisB, float angle, vector<vector<float> > &pointB, int index);
bool GroupRotationp(const vector<float> &axisA, const vector<float> &axisB, float angle, vector<poseCoord> &pointB);
std::string name2eltx( std::string atmid);

#endif
