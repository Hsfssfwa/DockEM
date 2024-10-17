// C++ headers
#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <cstring>
#include <iomanip>


// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArray4D.hh>
#include <fourier/FFT.hh>
#include "SplineInterp.h"
#include "xray_scattering.h"

// My headers
#include "readpdb.h"
#include "GeometryTools.h"

using namespace std;


// x mod y, returns z in [0,y-1]
int pos_mod(int x,int y);
float pos_mod(float x,float y);
double pos_mod(double x,double y);
float square_len(vector<float> vec_s);

static void swap4_aligned(void *v, long ndata);
vector<string> string_split(string const & in, char splitchar);

int findSampling(float MINSMP, int NMUL);
template<class S, class T>
void resample(ObjexxFCL::FArray3D< S > const &density, ObjexxFCL::FArray3D< T > &newDensity,vector<int> newDims);
void spline_coeffs(ObjexxFCL::FArray3D< double > const & data , ObjexxFCL::FArray3D< double > & coeffs);
void spline_coeffs(ObjexxFCL::FArray3D< float > const & data , ObjexxFCL::FArray3D< double > & coeffs);


class ElectronDensity
{
public:
	/// null constructor
	ElectronDensity();

	bool readMRCandResize(
		std::string mapfile,
		float reso,
		float gridSpacing
	);

	bool readMRCandResize(
		std::istream & mapin,
		std::string mapfile,
		float reso,
		float gridSpacing
	);

		/// @brief constructor from an FArray3D (debugging only!)
	template<class T>
	ElectronDensity( ObjexxFCL::FArray3D< T > const &map,
		float apix = 1.0,
		vector<float> new_origin=vector<float>(3,0.0),
		bool fftshift=false) {
		init();

		isLoaded = true;
		origin = new_origin;
		grid = vector<int>(3,map.u1());
		grid[1] = map.u2();
		grid[2] = map.u3();
		cellDimensions = vector<float>(3,apix*map.u1());
		cellDimensions[1] = apix*map.u2();
		cellDimensions[2] = apix*map.u3();
		cellAngles = vector<float>(3,90.0);
		density.dimension( map.u1(),map.u2(),map.u3() );
		//if ( fftshift ) origin -= grid/2;
		if ( fftshift ) 
		{
			origin[0] -= grid[0]/2;
			origin[1] -= grid[1]/2;
			origin[2] -= grid[2]/2;
		}
		for ( int i=1; i<=(int)map.u1(); ++i ) {
			int tmpi = pos_mod( i-(map.u1()/2)-1 , (int)map.u1())+1;
			int fi = (int)(fftshift ? tmpi : i);
			for ( int j=1; j<=(int)map.u2(); ++j ) {
				int tmpj = pos_mod( j-(map.u2()/2)-1 , (int)map.u2())+1;
				int fj = (int)(fftshift ? tmpj : j);
				for ( int k=1; k<=(int)map.u3(); ++k ) {
					int tmpk = pos_mod( k-(map.u3()/2)-1 , (int)map.u3())+1;
					int fk = (int)(fftshift ? tmpk : k);
					density(fi,fj,fk) = (float)map(i,j,k);
				}
			}
		}
	}

	/// @brief initialize vars from command line options
	void init();

	template<class Q>
	void idx2cart( vector<Q> const & idxX , vector<float> &cartX ) const {
		vector<float> fracX(3,(idxX[0]  + origin[0] -1 ) / grid[0]);
		fracX[1]= (idxX[1]  + origin[1] -1 ) / grid[1] ;
		fracX[2]= (idxX[2]  + origin[2] -1 ) / grid[2] ;
//		cartX = f2c*fracX;
		MatrixTimesTransVector(f2c,fracX,cartX);
	}

	/// @brief (debugging) Write MRC mapfile
	bool writeMRC(std::string mapfilestem);

	///@brief set voxel spacing of the map
	void set_voxel_spacing(vector<float> apix );

	///@brief set voxel spacing of the map
	void
	set_voxel_spacing( float apix ) {
		set_voxel_spacing( vector<float>(3,apix) );
	}	

	///@brief resize the map via FFT resampling
	void resize( float approxGridSpacing );

	///@brief  access raw density data
	inline ObjexxFCL::FArray3D< float > const & get_data() const { return density; };	

	/// @brief get Rho Calc
	void calcRhoC(
		poseCoords const &pose,
		float highres_limit,
		ObjexxFCL::FArray3D< float > &rhoC,
		ObjexxFCL::FArray3D< float > &mask,
		float forceB = -1.0,
		float B_upper_limit = 600,
		float force_mask = -1.0 );	

	void calcRhoCy(
		Model const &pose,
		float highres_limit,
		ObjexxFCL::FArray3D< float > &rhoC,
		ObjexxFCL::FArray3D< float > &mask,
		float fixed_mask_B  = -1,
		float B_upper_limit = 600,
		float force_mask =-1 );

	void calcRhoCx(
		Model const &pose,
		float highres_limit,
		ObjexxFCL::FArray3D< float > &rhoC,
		ObjexxFCL::FArray3D< float > &mask,
		float fixed_mask_B = -1 ,
		float B_upper_limit = 600 ,
		float force_mask =-1);	

	float S2(int h, int k, int l) {
		return ( h*h*RcellDimensions[0]*RcellDimensions[0]
			+ k*k*RcellDimensions[1]*RcellDimensions[1]
			+ l*l*RcellDimensions[2]*RcellDimensions[2]
			+ 2*h*k*RcellDimensions[0]*RcellDimensions[1]*cosRcellAngles[2]
			+ 2*h*l*RcellDimensions[0]*RcellDimensions[2]*cosRcellAngles[1]
			+ 2*k*l*RcellDimensions[1]*RcellDimensions[2]*cosRcellAngles[0] );
	}

	/// @brief Real-space correlation
	float getRSCC( ObjexxFCL::FArray3D< float > const &density2,  ObjexxFCL::FArray3D< float > const &mask);
	float getRSCCX( ObjexxFCL::FArray3D< float > const &density2,  ObjexxFCL::FArray3D< float > const &mask); 

	/// @brief Real-space correlation
	float getRSCC( ObjexxFCL::FArray3D< float > &rhoC ) {
		ObjexxFCL::FArray3D< float > nullmask;
		return getRSCC( rhoC, nullmask );
	}		


public:
	// Controllable parameters
	float preso, pATOM_MASK, pCA_MASK, pforce_apix_on_map_load_, pSC_scaling_;
	int pWINDOW_;
	bool pscore_window_context_, premap_symm_;
	int pnkbins_;
	bool pde_edensity;




private:
	bool isLoaded;

	bool de_edenisty;

	// the density data array and spline coeffs
	ObjexxFCL::FArray3D< float > density;
	ObjexxFCL::FArray3D< double > coeffs_density_;

		// fft of density -- only used in FSC calc
	ObjexxFCL::FArray3D< std::complex<double> > Fdensity;

	// symmetric transforms (in frac. coords)
	// this is only used to expand the density data outside the ASU
	//    and is unrelated Rosetta's symmetric modelling
	vector< RT > symmOps;
	vector< FMatrix > symmRotOps;

	// min multiples in each dim
	vector<int> MINMULT;

	// map info
	float minimumB;    // minimum B factor allowed by map
	float effectiveB;  // B factor blurring based on map resolution
	vector< int > grid;	
	vector<float> origin;
	bool use_altorigin;   // which field to write origin to ... only affects map outputting!

	// unit cell, reciprocal unit cell parameters, volume
	vector<float> cellDimensions, cellAngles;
	vector<float> RcellDimensions, cosRcellAngles;
	float V, RV;

	// Controllable parameters
	float reso, ATOM_MASK, CA_MASK, force_apix_on_map_load_, SC_scaling_;
	float ATOM_MASK_PADDING;
	int WINDOW_;
	bool score_window_context_, remap_symm_;


	// (fast scoring) precomputed rhocrhoo, d_rhocrhoo
	ObjexxFCL::FArray4D< double > fastdens_score;
	int nkbins_;	
	vector< int > fastgrid;           // grid & origin		

	///////////////////
	/// CRYSTAL INFO
	///////////////////
	// TO DO --- put all this in a self-contained class
	// converting fractional to cartesian coords
	//numeric::xyzMatrix<core::Real> f2c, c2f;
	vector<vector<float> > f2c,c2f;



private:
	// helper functions for symmetry
	void initializeSymmOps( vector< string > const & symList );
	void computeCrystParams();
	void expandToUnitCell();
	/// @brief The function is called everytime the density changes
	void
	density_change_trigger();

	// helper functions for map statistics
	void computeGradients();
	// volume of 1 voxel
	inline float voxel_volume() {
		return float(float(V) / float((grid[0]*grid[1]*grid[2])));
	};
};


