#ifndef vdwRadiiHead
#define vdwRadiiHead
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <string.h>
#include <vector>
#include <cstring>
#include <cstdlib>
//#include "readpdb.h"
#include "randomx.h"
#include "GeometryTools.h"

using namespace std;

#define nosegunit 9
#define nosegunit2 3
#define nosegdh 30
#define lennca 1.460f
#define delnca 0.004f*20
#define lencac 1.525f
#define delcac 0.004f*20
#define lencn 1.338f
#define delcn 0.005f*15
#define lennn 3.131f
#define lencc 3.231f
#define lencaca 3.813f
#define delcaca 0.019f*10
#define lencan1 2.441f
#define delcan1 0.036f*4
#define lencca1 2.446f
#define delcca1 0.036f*4
#define lennc 2.460f
#define angcacn 116.617f
#define angcnca 121.614f
#define delnc 0.012f*10
#define lenco 1.2314f
#define lencacb 1.533f
#define lenccb 2.506f
#define lennh 1.0000f //0.9666f
#define lencaha 1.0505f
#define degrad 180.0/PI
#define raddeg PI/180.0
#define tornccaha 240.4438f
#define lencaha 1.0505f
#define angccaha 108.7950f
#define angncac 111.008f
#define topno 200

//#define PI 3.1415926
//typedef vector<vector<float> > matrix;

static double casta[][3]={
//i,i+2, angle, torsion
{5.460699,1.598969,50.143396},//alpha
{6.640796,2.151540,190.00000},//beta
{6.196809,1.882359,250.00000}//coil
};
static double cbsta[][3]={
1.52369,1.92448,122.35124, 1.52962,1.91875,122.28332, 1.53149,1.92096,122.53073, 1.53132,1.92149,122.55859,
1.53219,1.91936,122.36077, 1.51371,1.90607,121.58025, 1.53172,1.92135,122.58755, 1.54507,1.92240,122.99168,
1.53180,1.92265,122.48313, 1.53108,1.92040,122.28572, 1.53078,1.91922,122.34940, 1.53148,1.92241,122.84907,
1.52996,1.94084,115.54662, 1.53057,1.92128,122.53531, 1.53085,1.92046,122.42439, 1.52991,1.91734,122.39611,
1.54070,1.91400,122.79225, 1.54500,1.92132,123.02119, 1.53172,1.92002,122.56818, 1.53251,1.91842,122.36359
};//n c ca cb  len ang tor

static double sglatavg2[][3]={//exclude cb
0.00000, 1.00000,  180.00, 2.35390, 2.14468,  130.58, 2.41571, 2.13243,  135.58, 3.15619, 2.13210,  136.99, 2.93226, 2.16722,  141.46,
0.00000, 1.00000,    0.00, 2.80663, 2.20671,  139.13, 2.51831, 2.13046,  131.80, 3.63956, 2.13159,  138.74, 2.84465, 2.18538,  144.33,
3.11849, 2.14252,  140.96, 2.44354, 2.17746,  137.26, 2.23972, 2.20653,   57.35, 3.14531, 2.15497,  138.77, 4.01554, 2.10163,  135.50,
2.05520, 1.86421,  116.85, 2.07510, 2.02899,  124.14, 2.13220, 2.02264,  134.07, 3.12109, 2.06139,  143.86, 3.11752, 2.15751,  142.69,
};
static double stddelta=3.9;
static double stdangle[][4]={
{111.0,2.7*stddelta,111.0-2.7*stddelta,111.0+2.7*stddelta},//ncac
{117.2,2.2*stddelta,117.2-2.2*stddelta,117.2+2.2*stddelta},//cacn
{121.7,2.5*stddelta,121.7-2.5*stddelta,121.7+2.5*stddelta},//cnca

{113.1,2.5*stddelta,113.1-2.5*stddelta,113.1+2.5*stddelta},//ncac gly
{116.2,2.0*stddelta,116.2-2.0*stddelta,116.2+2.0*stddelta},//cacn gly 
{122.3,2.1*stddelta,122.3-2.1*stddelta,122.3+2.1*stddelta},//cnca gly  

{112.1,2.6*stddelta,112.1-2.6*stddelta,112.1+2.6*stddelta},//ncac pro  
{117.1,2.8*stddelta,117.1-2.8*stddelta,117.1+2.8*stddelta},//cacn pro 
{119.3,1.5*stddelta,119.3-1.5*stddelta,119.3+1.5*stddelta},//cnca pro

{120.1,2.1*stddelta,120.1-2.1*stddelta,120.1+2.1*stddelta},//caco  
{120.6,1.8*stddelta,120.6-1.8*stddelta,120.6+1.8*stddelta},//caco gly  
{120.2,2.4*stddelta,120.2-2.4*stddelta,120.2+2.4*stddelta},//caco pro  

{122.7,1.6*stddelta,122.7-1.6*stddelta,122.7+1.6*stddelta},//ocn  
{123.2,1.7*stddelta,123.2-1.7*stddelta,123.2+1.7*stddelta},//ocn gly  
{121.1,1.9*stddelta,121.1-1.9*stddelta,121.1+1.9*stddelta},//ocn pro  
};
static double stdlength[][4]={
	{1.459,0.020*stddelta,1.459-0.020*stddelta,1.459+0.020*stddelta},//nca  
	{1.525,0.026*stddelta,1.525-0.026*stddelta,1.525+0.026*stddelta},//cac  
	{1.336,0.023*stddelta,1.336-0.023*stddelta,1.336+0.023*stddelta},//cn  
	{3.813,0.080*stddelta,3.813-0.080*stddelta,3.813+0.080*stddelta},//caca

	{1.456,0.015*stddelta,1.456-0.015*stddelta,1.456+0.015*stddelta},//nca gly  
	{1.514,0.016*stddelta,1.514-0.016*stddelta,1.514+0.016*stddelta},//cac gly  
	{1.326,0.018*stddelta,1.326-0.018*stddelta,1.326+0.018*stddelta},//cn gly  

	{1.468,0.017*stddelta,1.468-0.017*stddelta,1.468+0.017*stddelta},//nca pro  
	{1.524,0.020*stddelta,1.524-0.020*stddelta,1.524+0.020*stddelta},//cac pro  
	{1.338,0.019*stddelta,1.338-0.019*stddelta,1.338+0.019*stddelta},//cn pro  

	{1.229,0.019*stddelta,1.229-0.019*stddelta,1.229+0.019*stddelta},//co 
	{1.232,0.016*stddelta,1.232-0.016*stddelta,1.232+0.016*stddelta},//co gly  
	{1.228,0.020*stddelta,1.228-0.020*stddelta,1.228+0.020*stddelta},//co pro  
};

static double angleaa[30]={0.01371,0.05788,0.115,0.17441,0.25209,0.3115,0.38842,0.45163,0.52094,0.5773,
0.62909,0.68469,0.73343,0.77684,0.80807,0.8431,0.86443,0.89337,0.9086,0.92764,
0.93906,0.95582,0.9642,0.97258,0.97791,0.98857,0.99238,0.99466,0.99999,1.0};
static int torideal[][2]={
148,159, 148,158, 148,159, 148,159, 149,156, 
148,159, 148,158, 148,157, 149,158, 148,158, 
148,158, 148,159, 150,161, 148,159, 148,159, 
148,158, 148,157, 148,157, 149,157, 149,156, 

104, 75, 119, 65, 135, 53, 113, 71, 117, 71, 
 90, 90, 121, 66, 119, 64, 113, 71, 126, 62, 
115, 70, 125, 64, 149, 70, 116, 69, 118, 67, 
110, 76, 119, 64, 119, 64, 120, 68, 114, 76, 

148, 72, 148, 67, 136,  0, 150,164, 148, 67, 
 41,  0, 148, 68, 145, 64, 147,170, 143, 74, 
148, 70,  25, 21, 149, 73, 135,  0, 149, 68, 
144, 78, 142, 81, 146, 66, 144, 71, 148, 65
};//h20 e20 c20

static double subratio[][21]={//20-i,3.0	
0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.002,0.011,0.040,0.118,0.288,0.585,1.000,
0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.002,0.007,0.027,0.081,0.197,0.400,0.683,1.000,
0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.001,0.006,0.021,0.063,0.154,0.312,0.533,0.779,1.000,
0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.001,0.005,0.018,0.054,0.133,0.269,0.460,0.673,0.863,1.000,
0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.001,0.005,0.017,0.050,0.123,0.250,0.426,0.624,0.801,0.927,1.000,
0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.001,0.005,0.016,0.049,0.119,0.242,0.413,0.604,0.775,0.897,0.968,1.000,
0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.001,0.004,0.016,0.048,0.118,0.239,0.408,0.597,0.766,0.887,0.956,0.988,1.000,
0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.001,0.004,0.016,0.048,0.117,0.238,0.406,0.595,0.763,0.884,0.953,0.985,0.997,1.000,
0.000,0.000,0.000,0.000,0.000,0.000,0.001,0.004,0.016,0.048,0.117,0.238,0.406,0.594,0.762,0.883,0.952,0.984,0.996,0.999,1.000,
0.000,0.000,0.000,0.000,0.000,0.001,0.004,0.016,0.048,0.117,0.238,0.406,0.594,0.762,0.883,0.952,0.984,0.996,0.999,1.000,1.000,
0.000,0.000,0.000,0.000,0.001,0.004,0.016,0.048,0.117,0.238,0.406,0.594,0.762,0.883,0.952,0.984,0.996,0.999,1.000,1.000,1.000,
0.000,0.000,0.000,0.001,0.004,0.016,0.048,0.117,0.238,0.406,0.594,0.762,0.883,0.952,0.984,0.996,0.999,1.000,1.000,1.000,1.000,
0.000,0.000,0.001,0.004,0.016,0.048,0.117,0.238,0.406,0.594,0.762,0.883,0.952,0.984,0.996,0.999,1.000,1.000,1.000,1.000,1.000,
0.000,0.001,0.004,0.016,0.048,0.117,0.238,0.406,0.594,0.762,0.883,0.952,0.984,0.996,0.999,1.000,1.000,1.000,1.000,1.000,1.000,
0.000,0.003,0.015,0.047,0.116,0.237,0.405,0.594,0.762,0.883,0.952,0.984,0.996,0.999,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
0.000,0.012,0.044,0.113,0.234,0.403,0.592,0.761,0.882,0.952,0.984,0.996,0.999,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
0.000,0.032,0.103,0.225,0.396,0.587,0.758,0.881,0.951,0.984,0.995,0.999,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
0.000,0.073,0.199,0.376,0.574,0.750,0.877,0.950,0.983,0.995,0.999,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
0.000,0.137,0.327,0.540,0.731,0.867,0.946,0.982,0.995,0.999,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
0.000,0.221,0.467,0.688,0.846,0.937,0.979,0.994,0.999,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
0.000,0.317,0.600,0.803,0.919,0.973,0.993,0.998,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
};

static	double vdwd[8][8]={
3.600,2.300,3.700,2.900,3.500,3.200,3.500,1.000,//ca
3.700,2.500,3.500,2.500,3.500,2.100,3.600,1.000,//n
2.300,1.200,2.700,2.700,2.300,1.800,2.200,1.000,//c
2.500,2.100,2.400,2.300,2.600,1.700,2.200,1.000,//o
3.500,2.300,3.500,2.800,3.300,2.300,3.400,1.000,//cb
3.500,2.200,3.000,1.800,3.500,2.200,3.500,1.000,//h
3.500,2.200,3.300,2.800,3.300,2.100,3.500,1.000,//ha
1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000};//sg ori 1.5
static	double vdwds[8][8]={
12.96000, 5.29000,13.69000, 8.41000,12.25000,10.24000,12.25000, 1.00000,
13.69000, 6.25000,12.25000, 6.25000,12.25000, 4.41000,12.96000, 1.00000,
 5.29000, 1.44000, 7.29000, 7.29000, 5.29000, 3.24000, 4.84000, 1.00000,
 6.25000, 4.41000, 5.76000, 5.29000, 6.76000, 2.89000, 4.84000, 1.00000,
12.25000, 5.29000,12.25000, 7.84000,10.89000, 5.29000,11.56000, 1.00000,
12.25000, 4.84000, 9.00000, 3.24000,12.25000, 4.84000,12.25000, 1.00000,
12.25000, 4.84000,10.89000, 7.84000,10.89000, 4.41000,12.25000, 1.00000,
 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000};
static	double vdwds2[8][8]={
12.96000, 5.29000,13.69000, 8.41000,12.25000,10.24000,12.25000, 0.64000,
13.69000, 6.25000,12.25000, 6.25000,12.25000, 4.41000,12.96000, 0.64000,
 5.29000, 1.44000, 7.29000, 7.29000, 5.29000, 3.24000, 4.84000, 0.64000,
 6.25000, 4.41000, 5.76000, 5.29000, 6.76000, 2.89000, 4.84000, 0.64000,
12.25000, 5.29000,12.25000, 7.84000,10.89000, 5.29000,11.56000, 0.64000,
12.25000, 4.84000, 9.00000, 3.24000,12.25000, 4.84000,12.25000, 0.64000,
12.25000, 4.84000,10.89000, 7.84000,10.89000, 4.41000,12.25000, 0.64000,
 0.64000, 0.64000, 0.64000, 0.64000, 0.64000, 0.64000, 0.64000, 0.64000};

static char aad1[]= {
'A','C','D','E','F',
'G','H','I','K','L',
'M','N','P','Q','R',
'S','T','V','W','Y',
'J','B','Z','X','O','U'};


static char aad3[][4] = {
"ALA","CYS","ASP","GLU","PHE",
"GLY","HIS","ILE","LYS","LEU",
"MET","ASN","PRO","GLN","ARG",
"SER","THR","VAL","TRP","TYR",
"XLE","ASX","GLX","UNK","XAA","SEC"};

static char dnarnaatoms[8][24][4] ={
	{"P","OP1","OP2","O5'","C5'","C4'","O4'","C3'","O3'","C2'","C1'","N9","C8","N7","C5","C6","N6","N1","C2","N3","C4"},//a21tot
	{"P","OP1","OP2","O5'","C5'","C4'","O4'","C3'","O3'","C2'","C1'","N9","C8","N7","C5","C6","O6","N1","C2","N2","N3","C4"},//g22toc
	{"P","OP1","OP2","O5'","C5'","C4'","O4'","C3'","O3'","C2'","C1'","N1","C2","O2","N3","C4","N4","C5","C6"},//c19
	{"P","OP1","OP2","O5'","C5'","C4'","O4'","C3'","O3'","C2'","C1'","N1","C2","O2","N3","C4","O4","C5","C7","C6"},//t20
	{"P","OP1","OP2","O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'","N9","C8","N7","C5","C6","N6","N1","C2","N3","C4"},//a22tou
	{"P","OP1","OP2","O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'","N9","C8","N7","C5","C6","O6","N1","C2","N2","N3","C4"},//g23toc
	{"P","OP1","OP2","O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'","N1","C2","O2","N3","C4","N4","C5","C6"},//c20
	{"P","OP1","OP2","O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'","N1","C2","O2","N3","C4","O4","C5","C6"},//u20
};
static char dnarnares[8][4] ={
	" DA"," DG"," DC"," DT",
	"  A","  G","  C","  U"};

struct vwdt
{
	float radius;
	float well_depth;
};

struct point3d
{
	float x,y,z;
};

struct point3s
{
	float x,y,z;
};

struct boneinfo
{
	int indn,indca,indc,indo,indcb,indf,indh,indha;//pro no h gly has 1ha no cb
	int indall[84];
	int resind;
	char resid;//ACD
	int istart,iend;
	unsigned char sst;//0 coil 1 helix 2 sheet 3 turn
};
 struct lableproblematic
{
	int nn[60];//wrong pos
	double fn[60];//totval
	bool *flap[60];//pos concern
	int *indn[60];//ene int
	double *indf[60];//ene dou
};
 struct point3f
{
	///////////////////////////////////////////////////////aaa,ssm,ss2 not exist in fragments
	int resind;//from fragments
	float x,y,z;
	float len[3],ang[3],tor[3];//n ca c  0psi-1 n-1 ca-1 c-1 n; 1ome ca-1 c-1 n ca; 2phi c-1 n ca c
	float vpos,tpos,prob;//v means co-pair in the left; t means co-pair in the right
	int indl,indr,tpl,tpr;
	float leng,angl,phi;
	point3s ptn,ptc,pto,ptb,pth,ptha,ptsg,ptg;
	char name[8];//from fragments
	char residueid;//ACD from fragments
	char stype;//from fragments
	char ss2;//psipred in the fragment
	char ssm;//mycalculation
	char aaa;//realseq ACD
	unsigned char iaa;//realACD index
};
 struct residueatoms
{
	int scpr,mcpr;
	int indl,indr,tpl,tpr;
	int aopt;//0 move 1 use initial then to move 2 fixed 3 move only close to 0
	int indrot;//index for rotamer
	int residueid;//sequence
	int bintor[7];//phi psi [0,71] 
	point3d pc[6],pn[6];
	float dof[10],tor[3],dih[3],len[3],ang[3];//[iphi psii] [i-1psi iome iphi]
	float vpos,tpos,prob;//v means co-pair in the left; t means co-pair in the right
	float sspro[4];//c h e solve
	point3s ptv[15];
	point3s pth[16];
	point3s pbb[2],pvb[2],ptsg,ptg;//all; only ncacocb	
	char aaa,resind;//acd 1-20
	char stype;//from reference
	char ss2;//psipred
	char ssm;//mycalculation
//	char ti;//aass
	//0 rama 1 rotout 1-5 2 bondlengangle 6-13 hd[3] clash 14-16
	bool cw[20],hd[20];//wrong or not,need to process if hd is true	
};
typedef struct pairaa
{
	bool ispair;
	double dist;
	double dstd;
	double pval;
	int pbin,tnum;
	double ratio[20];
}pairaa;
typedef struct sssegment
{
	int init;
	int term;
	char ss;
}sssegment;
typedef struct segse
{
	int indss2;
	int init;
	int term;
}segse;
typedef struct alphahelix
{
	segse seg;
	double dist1,dist2;
	point3d cen,dir;
}alphahelix;
typedef struct betastrand
{
	segse seg;
	double dist1,dist2;
	point3d cen,dir,pla;
}betastrand;

float GetVdwRadius(const string& atm);
float GetVdwEg(string& R1,string& R2,float Dij);
float GetVdwEg(float& R1,float& R2,float Dij);
float GetVdwEgC(string& R1,string& R2,float Dij);
float GetVdwEgCT(string& R1,string& R2,float Dij);
float GetVdwEgCx(string& R1,string& R2,float Dij);
float GetVdwEgCG(string& R1,string& R2,float Dij);
float GetVdwEgCH(string& R1,string& R2,float Dij);
vwdt GetVdwRadiusC(const string& atm);
vwdt GetVdwRadiusx(const string& atm);
void Spicker(vector<vector<point3f>> vect_decstr,vector<vector<point3f>> &out_model);
void Spickerx(vector<vector<point3f>> vect_decstr,vector<vector<point3f>> &out_model);
int energyclash(vector<point3f> &decstr,int numseq,float &fene);
double energyexcludedvolume(vector<point3f> decstr,int numseq);
int sec_str(float dis13, float dis14, float dis15, float dis24, float dis25, float dis35);
string sec_strx(float dis13, float dis14, float dis15, float dis24, float dis25, float dis35);
void make_secx(vector<vector<float> > x, int len, vector<string> &sec);
void make_sec(vector<vector<float> > x, int len, vector<int> &sec);
int getmovetype(vector<float> tmov,int totmov,float trandnum);

void v2rot(point3d p1,double rot[]);
point3d rotv(point3d p1,double rot[]);
point3d rana(double theta);
point3d setv(float tx,float ty,float tz);
point3d minu(point3d p1,point3d p2);
float norm(point3d p1);
point3d prod(point3d p1,point3d p2);
point3d unit(point3d p1);
double phi(double xi,double yi,double zi,double xj,double yj,double zj,double xk,
           double yk,double zk,double xl,double yl,double zl);
float phi(float xi,float yi,float zi,float xj,float yj,float zj,float xk,
           float yk,float zk,float xl,float yl,float zl);
float squgaussian(float x,float sigma,float miu);
float expgaussian2(float x,float sigma,float miu);
float expgaussian(float x,float sigma,float miu);
float expgaussianq(float xx,float sigsig);
float angv(point3d p1,point3d p2);
float energyhbondcanc(vector<point3f> &decstr,int numseq);
void calcssennhoc(vector<boneinfo> &bb, int numbb, vector<point3f> &decstr);
void str2tor(vector<point3f> &decstr,int seqnum,int type);
char numtoss(int num);
void getLinLen(int Nlen, vector<point3f> xyzCa, vector<int> &a);
bool GroupRotationpid(const vector<float> &axisA, const vector<float> &axisB, float angle, vector<point3f> &pointB,int index0,int indexn);
void randomwalk(vector<vector<float> > &LinxyzCa, int istart, int iend);
float *get_bond(float al);
vector<float> get_bondx(float al);
int mcheck(int i0, vector<vector<float>> LinxyzCa, int amcheck_dis);
void singlemoveLMPf(vector<point3f> &tmstr,int k,int mtype);
bool mcfragsweepLMP2(vector<point3f> &tmstr,int numseq,int poss,int pose);
void singlemoveLMPb(vector<point3f> &tmstr,int k,int mtype);
point3d scalx(point3d p1,double f);
int fenergybondangle(vector<residueatoms> &decstr,int numseq,double fene,bool flagacc);
int aminoid(vector<char> &aminoname);
int aminoid(char aminoname);
int Eenergybondangle(vector<point3f> &decstr,int numseq,float &fene);
int energybondlength(vector<point3f> &decstr,int numseq,float &fene);
float energypartialrmsd(point3f *decstr,int numseq);
double energyhbondnhoc2(vector<point3f> &decstr,int numseq);
double energyhbondnhoc3(vector<point3f> &decstr,int numseq);
double energyhbondnhoc4(vector<point3f> &decstr,int numseq);
double energyhbondnhoc2(vector<point3f> &decstr,int numseq,float &hbalpha,float &hbbeta);
bool tor2pos22(float xi,float yi,float zi,float xj,float yj,float zj,float xk,
           float yk,float zk,float tang,float tleng,float tinner,float *xl,float *yl,float *zl);
int maxnormal(point3d p1);
void  q2rot(double *q,double rmax[]);
point3d mmat(double mat[3][3],point3d p1);
point3d mmat(double mat[9],point3d p1);
bool tor2stroh(vector<point3f> &decstr,int seqnum);
bool movementhel(vector<point3f> &tmstr,int numseq);
//bool genesse(vector<point3f> &decstr,int numseq,sssegment *sse,int &numsse);
bool genesse(vector<point3f> &decstr,int numseq,vector<sssegment> &sse,int &numsse);
bool genessex(vector<point3f> &decstr,int numseq,vector<sssegment> &sse,int &numsse,int &numshelix);
void calcabind(vector<int> &alphasind, vector<int> &alphaind, vector<int> &betaind,int numsse,vector<sssegment> &sse);
void str2torp(vector<point3f> &decstr,int seqnum,int istart,int iend);
double energyhelixpacking(vector<point3f> &decstr,int numseq,vector<vector<float> > hhda,vector<vector<vector<float>>> haad,float &enehpk);
int linecross(point3d ps1,point3d pe1,point3d ps2,point3d pe2,point3d *pc1,point3d *pc2,double  *dist);
double sdet(double a[],int n);
point3d addv(point3d p1,point3d p2);
double footpoint(point3d ps1,point3d pe1, point3d tp,point3d *tfp,double *tdist);
void getahelix(vector<point3f> &decstr,int numseq,vector<alphahelix> &ahelix);
vector<vector<float>> loadhelixhelix(const char *filename);
vector<vector<vector<float>>> loadhelixpairsg(const char *filename);
bool mcmovementLMP(vector<point3f> &tmstr,int seqlen,int nums,int nume);
point3d ranv(double fac);
bool mcfragsweepaaa(vector<point3f> &tmstr,int numseq,int numshelix);
bool mcfragsweepCCD6(vector<point3f> &decstr,int numseq,int lps,int lpe,point3d pt1,int ind1,point3d pt2,int ind2);
bool mcfragsweepCCD4(vector<point3f> &tmstr2,int numseq,int lps,int lpe,point3d pt1,int ind1,point3d pt2,int ind2);
bool tor2str(vector<point3f> &decstr,int seqnum,int type);
bool tor2strp(vector<point3f> &decstr,int seqnum,int istart);
bool tor2strp(vector<point3f> &decstr,int seqnum,int istart,int iend);
bool itor2strp(vector<point3f> &decstr,int seqnum,int iss,int istart);
bool itor2str(vector<point3f> &decstr,int seqnum);
bool tor2strsg2(vector<point3f> &decstr,int seqnum,vector<vector<vector<vector<float>>>> &sgposdat);
bool loadsgpos2(char *filename,int ndim,vector<vector<vector<vector<float>>>> &sgpos2);
int energyrama(vector<point3f> &decstr,int numseq,float &fene,vector<vector<vector<float>>> &rlogduke,vector<vector<vector<int>>> &ramaduke);
bool loadca2ncbins(char *filename,vector<vector<vector<vector<float> > > > &cancbinsx);
bool loadca2ncbins2(char *filename,vector<vector<vector<vector<float> > > > &cancbinsx);
void ca2nc(vector<point3f> &decstr,int numseq,vector<vector<vector<vector<float> > > > &cancbins);
double energytorsion2(vector<point3f> &decstr,int numseq,vector<vector<vector<double>>> phipsidisss);
bool loadphipsidisss(char *filename,vector<vector<vector<double>>> &phipsidisss);
bool loadramadukebn(char *datafile,vector<vector<vector<float>>> &rlogdukex,vector<vector<vector<int>>> &ramadukex);
int getmovetype(double tmov[],int totmov,double trandnum);
bool mcfragsweepome(vector<point3f> &decstr,int numseq,float *phipsiprob[8][20]);
bool mcfragsweepphi(vector<point3f> &decstr,int numseq,float *phipsiprob[8][20]);
bool mcfragsweeppsi(vector<point3f> &decstr,int numseq,float *phipsiprob[8][20]);
int findpos2(double *p,int is,int ie,double pt);
int findpos2(float *p,int is,int ie,float pt);
int findpos2(vector<float> p,int is,int ie,float pt);
bool loadphipsiprob(char *filename,double *phipsiprob[8][20]);
bool loadphipsiprob(char *filename,vector<vector<vector<float>>> &phipsiprob);
bool mcfragsweeplen(vector<point3f> &decstr,int numseq);
bool mcfragsweepang(vector<point3f> &decstr,int numseq);

#endif