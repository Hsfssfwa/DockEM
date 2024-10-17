#include "CaDensity.h"

using namespace std;

inline float d2r(float d) {return (d*PI/180.0);}
inline double d2r(double d) {return (d*PI/180.0);}
inline float square(float x) {return (x*x);}
inline double square(double x) {return (x*x);}

bool factorsLTE19(int X) {
	while ( X != 1 && X%2 == 0 ) X /= 2;
	while ( X != 1 && X%3 == 0 ) X /= 3;
	while ( X != 1 && X%5 == 0 ) X /= 5;
	while ( X != 1 && X%7 == 0 ) X /= 7;
	while ( X != 1 && X%11 == 0 ) X /= 11;
	while ( X != 1 && X%13 == 0 ) X /= 13;
	while ( X != 1 && X%17 == 0 ) X /= 17;
	while ( X != 1 && X%19 == 0 ) X /= 19;

	return (X == 1);
}

// x mod y, returns z in [0,y-1]
int pos_mod(int x,int y) {
	int r=x%y; if ( r<0 ) r+=y;
	return r;
}
float pos_mod(float x,float y) {
	float r=std::fmod(x,y); if ( r<0 ) r+=y;
	return r;
}
double pos_mod(double x,double y) {
	double r=std::fmod(x,y); if ( r<0 ) r+=y;
	return r;
}

static void swap4_aligned(void *v, long ndata) {
	int *data = (int *) v;
	long i;
	for ( i=0; i<ndata; i++ ) {
		int *N;
		N = data + i;
		*N=(((*N>>24)&0xff) | ((*N&0xff)<<24) | ((*N>>8)&0xff00) | ((*N&0xff00)<<8));
	}
}

int findSampling(float MINSMP, int NMUL) {
	if ( MINSMP <= 0 ) return NMUL;

	// multiple of nmul nearest minsmp
	int N = (int) floor( MINSMP/NMUL + 0.5 ) * NMUL;

	// increment until no factors >= 19
	while ( !factorsLTE19(N) )
			N += NMUL;

	return N;
}

float square_len(vector<float> vec_s)
{
	float squ_v=0.0;
	for(int i=0;i<vec_s.size();i++)
	{
		squ_v +=vec_s[i]*vec_s[i];
	}
	return squ_v;
}
/*
bool pose_has_nonzero_Bs( Model const & pose ) {  //////
	for(int i=0;i<pose.chains.size();i++)
	{
		Chain chainxx=pose.chains[i];
		for(int j=0;j<chainxx.residues.size();j++)
		{
			Residue residuexx=chainxx.residues[j];
			for(int k=0;k<residuexx.atoms.size();k++)
			{
				Atom atomxx=residuexx.atoms[k];
				float B = atomxx.tempfac;
				if( B>0 ) return true;
			}
		}
	}
	return false;
}  */

bool pose_has_nonzero_Bs( Model const & pose ) {
	for(int i=0;i<pose.chains.size();i++)
	{
		for(int j=0;j<pose.chains[i].residues.size();j++)
		{
			for(int t=0; t<pose.chains[i].residues[j].atoms.size();t++)
			{
				float tf=pose.chains[i].residues[j].atoms[t].tempfac;
				if(tf>0) return true;
			}
		}
	}
//	for ( int i=0 ; i<(int)pose.size(); ++i ) {
//		if ( pose[i].B_ > 0 ) {
//			return true;
//		}
//	}
	return false;
}

bool pose_has_nonzero_Bs( poseCoords const & pose ) {
	for ( int i=0 ; i<(int)pose.size(); ++i ) {
		if ( pose[i].B_ > 0 ) {
			return true;
		}
	}
	return false;
}


void
ElectronDensity::init() {
	isLoaded = false;

	grid = vector< int >(3,0);
	origin = vector< float >(3,0);
	cellDimensions = vector< float >(3,1.0);
	cellAngles = vector< float >(3,90.0);
	use_altorigin =  false;

	// command line overrides defaults
	reso = preso;
	ATOM_MASK = pATOM_MASK;
	ATOM_MASK_PADDING = 2;
	CA_MASK = pCA_MASK;
	WINDOW_ = pWINDOW_;
	score_window_context_ = pscore_window_context_;
	remap_symm_ = premap_symm_;
	force_apix_on_map_load_ = pforce_apix_on_map_load_;
	nkbins_ = pnkbins_;

	// use-specified B factor (may be overridden)
	effectiveB = 0;
}

/// null constructor
ElectronDensity::ElectronDensity() 
{
	init();
}

// ElectronDensity::readMRC(std::string mapfile)
//      read an MRC/CCP4 density map
bool
ElectronDensity::readMRCandResize(
	std::string mapfile,
	float reso,
	float gridSpacing
) {
	std::ifstream mapin(mapfile.c_str() , std::ios::binary | std::ios::in );
	bool isLoaded(readMRCandResize(mapin, mapfile, reso, gridSpacing));
	mapin.close();

	return isLoaded;
}

float ElectronDensity::getRSCC( ObjexxFCL::FArray3D< float > const &density2,  ObjexxFCL::FArray3D< float > const &mask) {
//	runtime_assert( density.u1()==density2.u1() && density.u2()==density2.u2() && density.u3()==density2.u3() );

	float sumC_i=0, sumO_i=0, sumCO_i=0, vol_i=0, CC_i=0;
	float sumO2_i=0.0, sumC2_i=0.0, varC_i=0, varO_i=0;
	float clc_x, obs_x, mask_x;
	for ( int x=0; x<density.u1()*density.u2()*density.u3(); ++x ) {
		clc_x = density2[x];
		obs_x = density[x];
		mask_x = (mask.u1() == 0) ? 1.0 : mask[x];

		sumCO_i += mask_x*clc_x*obs_x; sumO_i  += mask_x*obs_x;
		sumO2_i += mask_x*obs_x*obs_x; sumC_i  += mask_x*clc_x;
		sumC2_i += mask_x*clc_x*clc_x;
		vol_i   += mask_x;
	}
	varC_i = (sumC2_i - sumC_i*sumC_i / vol_i );
	varO_i = (sumO2_i - sumO_i*sumO_i / vol_i ) ;
	if ( varC_i == 0 || varO_i == 0 ) {
		CC_i = 0;
	} else {
		CC_i = (sumCO_i - sumC_i*sumO_i/ vol_i) / sqrt( varC_i * varO_i );
	}

	return CC_i;
}

void ElectronDensity::calcRhoC(
	poseCoords const &pose,
	float highres_limit,
	ObjexxFCL::FArray3D< float > &rhoC,
	ObjexxFCL::FArray3D< float > &mask,
	float fixed_mask_B /* = -1 */,
	float B_upper_limit /* = 600 */,
	float force_mask /*=-1*/ ) {

	// get rho_c
	rhoC.dimension(density.u1() , density.u2() , density.u3());
	mask.dimension(density.u1() , density.u2() , density.u3());
	rhoC = 0.0;
	mask = 0.0;

	bool use_Bs = pose_has_nonzero_Bs( pose );
	if ( !use_Bs ) {
		cout << "Input pose has no nonzero B factors ... setting to baseline" << std::endl;
	}

	vector<float> per_atom_ks(pose.size(), 0.0);
	vector<float> per_atom_Cs(pose.size(), 0.0);
	vector<float> ATOM_MASK_SQS(pose.size(), 0.0);
	poseCoords pose_grid = pose;
	for ( int i=0 ; i<(int)pose.size(); ++i ) {
		std::string elt_i = pose[i].elt_;
		OneGaussianScattering sig_j = get_A( elt_i );
		per_atom_ks[i] = sig_j.k( pose[i].B_, B_upper_limit );
		if ( use_Bs ) {
			per_atom_ks[i] = std::min ( per_atom_ks[i], (float) (4*M_PI*M_PI/minimumB) );
		} else {
			per_atom_ks[i] = std::min ( per_atom_ks[i], (float) (4*M_PI*M_PI/effectiveB) );
		}
		per_atom_Cs[i] = sig_j.C( per_atom_ks[i] );

		float C = per_atom_Cs[i], k = per_atom_ks[i];
		if ( C < 1e-6 ) continue;
		ATOM_MASK_SQS[i] = (1.0/k) * ( std::log( C ) - std::log(1e-4) );  // 1e-4 is the density value where the mask goes to 0

		// for b factor minimization, we don't want the mask to change.
		// so we fix the mask at a high B value
		if ( fixed_mask_B>0 ) {
			float mK = sig_j.k( fixed_mask_B, B_upper_limit );
			float mC = sig_j.C( mK );
			ATOM_MASK_SQS[i] = (1.0/mK) * (std::log( mC ) - std::log(1e-4));
		}

		// force mask
		if ( force_mask>0 ) {
			ATOM_MASK_SQS[i] = force_mask*force_mask;
		}

		// TO DO: OPTIONALLY control periodic/nonperiodic boundaries
		//atm_idx[0] = pos_mod (fracX[0]*grid[0] - origin[0] + 1 , (double)grid[0]);
		//atm_idx[1] = pos_mod (fracX[1]*grid[1] - origin[1] + 1 , (double)grid[1]);
		//atm_idx[2] = pos_mod (fracX[2]*grid[2] - origin[2] + 1 , (double)grid[2]);
		vector<float> cartX = pose[i].x_; // - getTransform();
	//	vector<float> fracX = c2f*cartX;
		vector<float> fracX;
		MatrixTimesTransVector(c2f,cartX,fracX);
		vector<float> atm_idx(3,fracX[0]*grid[0] - origin[0] + 1);
		atm_idx[1] = fracX[1]*grid[1] - origin[1] + 1;
		atm_idx[2] = fracX[2]*grid[2] - origin[2] + 1;
		pose_grid[i].x_ = atm_idx;
	}

	// PASS 1: rho_calc
	for ( int i=1 ; i<=(int)pose.size(); ++i ) {
		std::cout << "\rrho_calc " << i-1 << " of " << pose.size() << std::flush;

		float C = per_atom_Cs[i-1], k = per_atom_ks[i-1];
		if ( C < 1e-6 ) continue;

		// dist cutoff for mask
		float ATOM_MASK_SQ = ATOM_MASK_SQS[i-1];

		// dist cutoff for density
		float ATOM_DENS_SQ = (1.0/k) * (std::log( C ) - std::log(1e-4));

		vector<float> atm_j(3,0.0), del_ij(3,0.0), atm_idx;
		atm_idx =  pose_grid[i-1].x_;

		if ( atm_idx[0]<1.0 || atm_idx[0]>grid[0] ) continue;
		if ( atm_idx[1]<1.0 || atm_idx[1]>grid[1] ) continue;
		if ( atm_idx[2]<1.0 || atm_idx[2]>grid[2] ) continue;

		for ( int z=1; z<=density.u3(); ++z ) {
			atm_j[2] = z;
			del_ij[2] = (atm_idx[2] - atm_j[2]) / grid[2];
			del_ij[0] = del_ij[1] = 0.0;
			vector<float> frac_tmpz;
		//	if ( (f2c*del_ij).length_squared() > (ATOM_MASK_SQ) ) continue;
			MatrixTimesTransVector(f2c,del_ij,frac_tmpz);
			if(square_len(frac_tmpz)> (ATOM_MASK_SQ) ) continue;
			for ( int y=1; y<=density.u2(); ++y ) {
				atm_j[1] = y;
				del_ij[1] = (atm_idx[1] - atm_j[1]) / grid[1] ;
				del_ij[0] = 0.0;
			//	if ( (f2c*del_ij).length_squared() > (ATOM_MASK_SQ) ) continue;
				vector<float> frac_tmpy;
				MatrixTimesTransVector(f2c,del_ij,frac_tmpy);
				if(square_len(frac_tmpy)> (ATOM_MASK_SQ) ) continue;
				for ( int x=1; x<=density.u1(); ++x ) {
					atm_j[0] = x;
					del_ij[0] = (atm_idx[0] - atm_j[0]) / grid[0];
					vector<float> cart_del_ij;// = (f2c*del_ij);  // cartesian offset from (x,y,z) to atom_i
					MatrixTimesTransVector(f2c,del_ij,cart_del_ij);
				//	float d2 = (cart_del_ij).length_squared();
					float d2 = square_len(cart_del_ij);

					if ( d2 <= (ATOM_MASK_SQ) ) {
						mask(x,y,z) = 1.0; // problem?
						if ( d2 <= (ATOM_DENS_SQ) ) {
							float atm = C*exp(-k*d2);
							rhoC(x,y,z) += atm;
						}
					}
				}
			}
		}
	}
	std::cout << std::endl;

	if ( highres_limit == 0 ) return;

	// bandlimit mask at 'radius'
	ObjexxFCL::FArray3D< std::complex<double> > Fmask;
	fourier::fft3(mask, Fmask);
	//int H,K,L;
	for ( int z=1; z<=(int)grid[2]; ++z ) {
		int H = (z < (int)grid[2]/2) ? z-1 : z-grid[2] - 1;
		for ( int y=1; y<=(int)grid[1]; ++y ) {
			int K = (y < (int)grid[1]/2) ? y-1 : y-grid[1]-1;
			for ( int x=1; x<=(int)grid[0]; ++x ) {
				int L = (x < (int)grid[0]/2) ? x-1 : x-grid[0]-1;
				float S2c =  S2(H,K,L);

				// exp fade
				float scale = exp(-S2c*(highres_limit*highres_limit));
				Fmask(x,y,z) *= scale;
			}
		}
	}
	fourier::ifft3(Fmask, mask);
}

void ElectronDensity::calcRhoCx(
	Model const &pose,
	float highres_limit,
	ObjexxFCL::FArray3D< float > &rhoC,
	ObjexxFCL::FArray3D< float > &mask,
	float fixed_mask_B /* = -1 */,
	float B_upper_limit /* = 600 */,
	float force_mask /*=-1*/ )
	{
	// get rho_c
	rhoC.dimension(density.u1() , density.u2() , density.u3());
	mask.dimension(density.u1() , density.u2() , density.u3());
	rhoC = 0.0;
	mask = 0.0;
	bool use_Bs = pose_has_nonzero_Bs( pose );
	if ( !use_Bs ) {
		cout << "Input pose has no nonzero B factors ... setting to baseline" << std::endl;
	}

	vector<float> per_atom_ks;
	vector<float> per_atom_Cs;
	vector<float> ATOM_MASK_SQS;
	vector<vector<float> > pose_grid;
	int tmp_p=0;
	for(int i=0;i<pose.chains.size();i++)
	{
		Chain chain_x=pose.chains[i];
		for(int j=0;j<chain_x.residues.size();j++)
		{
			Residue residue_x = chain_x.residues[j];
			for(int t=0; t<residue_x.atoms.size();t++)
			{
				Atom atom_x = residue_x.atoms[t];
				std::string elt_i = name2eltx(atom_x.atona);
				OneGaussianScattering sig_j = get_A( elt_i );
				per_atom_ks.push_back(sig_j.k( atom_x.tempfac, B_upper_limit ));
				if ( use_Bs ) {
					per_atom_ks[tmp_p] = std::min ( per_atom_ks[tmp_p], (float) (4*M_PI*M_PI/minimumB ));
				} else {
					per_atom_ks[tmp_p] = std::min ( per_atom_ks[tmp_p], (float) (4*M_PI*M_PI/effectiveB ));
				}
				per_atom_Cs.push_back(sig_j.C( per_atom_ks[tmp_p] ));

				float C = per_atom_Cs[tmp_p], k = per_atom_ks[tmp_p];
				if ( C < 1e-6 ) continue;
				ATOM_MASK_SQS.push_back((1.0/k) * ( std::log( C ) - std::log(1e-4) ));  // 1e-4 is the density value where the mask goes to 0

				// for b factor minimization, we don't want the mask to change.
				// so we fix the mask at a high B value
				if ( fixed_mask_B>0 ) {
					float mK = sig_j.k( fixed_mask_B, B_upper_limit );
					float mC = sig_j.C( mK );
					ATOM_MASK_SQS[tmp_p] = (1.0/mK) * (std::log( mC ) - std::log(1e-4));
				}

				// force mask
				if ( force_mask>0 ) {
					ATOM_MASK_SQS[tmp_p] = force_mask*force_mask;
				}

				// TO DO: OPTIONALLY control periodic/nonperiodic boundaries
				//atm_idx[0] = pos_mod (fracX[0]*grid[0] - origin[0] + 1 , (double)grid[0]);
				//atm_idx[1] = pos_mod (fracX[1]*grid[1] - origin[1] + 1 , (double)grid[1]);
				//atm_idx[2] = pos_mod (fracX[2]*grid[2] - origin[2] + 1 , (double)grid[2]);
				vector<float> cartX = atom_x.xyzVect; // - getTransform();
			//	vector<float> fracX = c2f*cartX;
				vector<float> fracX;
				MatrixTimesTransVector(c2f,cartX,fracX);
				vector<float> atm_idx(3,fracX[0]*grid[0] - origin[0] + 1);
				atm_idx[1] = fracX[1]*grid[1] - origin[1] + 1;
				atm_idx[2] = fracX[2]*grid[2] - origin[2] + 1;
				pose_grid.push_back( atm_idx );				
				tmp_p = tmp_p+1;
			}
		}
	}

	// PASS 1: rho_calc
	int tmp_pp=0;
	for(int i=0;i<pose.chains.size();i++)
	{
		Chain chain_x=pose.chains[i];
		for(int j=0;j<chain_x.residues.size();j++)
		{
			Residue residue_x = chain_x.residues[j];
			for(int t=0; t<residue_x.atoms.size();t++)
			{

				Atom atom_x = residue_x.atoms[t];
				std::cout << "\rrho_calc " << tmp_pp << " of " << std::flush;	
				float C = per_atom_Cs[tmp_pp], k = per_atom_ks[tmp_pp];
				if ( C < 1e-6 ) continue;	
				
				// dist cutoff for mask
				float ATOM_MASK_SQ = ATOM_MASK_SQS[tmp_pp];

				// dist cutoff for density
				float ATOM_DENS_SQ = (1.0/k) * (std::log( C ) - std::log(1e-4));

				vector<float> atm_j(3,0.0), del_ij(3,0.0), atm_idx;
				atm_idx =  pose_grid[tmp_pp];

				if ( atm_idx[0]<1.0 || atm_idx[0]>grid[0] ) continue;
				if ( atm_idx[1]<1.0 || atm_idx[1]>grid[1] ) continue;
				if ( atm_idx[2]<1.0 || atm_idx[2]>grid[2] ) continue;

				for(int z=1;z<=density.u3();z++)
				{
					atm_j[2] = z;
					del_ij[2] =(atm_idx[2]-atm_j[2])/grid[2];
					vector<float> frac_tmpz;
					MatrixTimesTransVector(f2c,del_ij,frac_tmpz);
					if(square_len(frac_tmpz)> (ATOM_MASK_SQ) ) continue;					
					for(int y=1;y<=density.u2();y++)
					{
						atm_j[1] = y;
						del_ij[1] = (atm_idx[1] - atm_j[1])/grid[1];
						vector<float> frac_tmpy;
						MatrixTimesTransVector(f2c,del_ij,frac_tmpy);
						if(square_len(frac_tmpy)> (ATOM_MASK_SQ) ) continue;						
						for(int x=1;x<=density.u1();x++)
						{
							atm_j[0] = x;
							del_ij[0] = (atm_idx[0] - atm_j[0])/grid[0];
							vector<float> cart_del_ij;
							MatrixTimesTransVector(f2c,del_ij,cart_del_ij);
							float d2 = square_len(cart_del_ij);
						
							if ( d2 <= (ATOM_MASK_SQ) ) {
								mask(x,y,z) = 1.0; // problem?
								if ( d2 <= (ATOM_DENS_SQ) ) {
									float atm = C*exp(-k*d2);
									rhoC(x,y,z) += atm;
								}
							}
						}

					}
				}				
				tmp_pp = tmp_pp +1;
			}
		}
	}
	cout<<endl;

	if ( highres_limit == 0 ) return;

	// bandlimit mask at 'radius'
	ObjexxFCL::FArray3D< std::complex<double> > Fmask;
	fourier::fft3(mask, Fmask);
	//int H,K,L;
	for ( int z=1; z<=(int)grid[2]; ++z ) {
		int H = (z < (int)grid[2]/2) ? z-1 : z-grid[2] - 1;
		for ( int y=1; y<=(int)grid[1]; ++y ) {
			int K = (y < (int)grid[1]/2) ? y-1 : y-grid[1]-1;
			for ( int x=1; x<=(int)grid[0]; ++x ) {
				int L = (x < (int)grid[0]/2) ? x-1 : x-grid[0]-1;
				float S2c =  S2(H,K,L);

				// exp fade
				float scale = exp(-S2c*(highres_limit*highres_limit));
				Fmask(x,y,z) *= scale;
			}
		}
	}
	fourier::ifft3(Fmask, mask);
 }

void ElectronDensity::calcRhoCy(
	Model const &pose,
	float highres_limit,
	ObjexxFCL::FArray3D< float > &rhoC,
	ObjexxFCL::FArray3D< float > &mask,
	float fixed_mask_B /* = -1 */,
	float B_upper_limit /* = 600 */,
	float force_mask /*=-1*/ )
	{
	// get rho_c
	rhoC.dimension(density.u1() , density.u2() , density.u3());
	mask.dimension(density.u1() , density.u2() , density.u3());
	rhoC = 0.0;
	mask = 0.0;
	bool use_Bs = pose_has_nonzero_Bs( pose );
	if ( !use_Bs ) {
		cout << "Input pose has no nonzero B factors ... setting to baseline" << std::endl;
	}

	// PASS 1: rho_calc
	for(int i=0;i<pose.chains.size();i++)
	{
		Chain chain_x=pose.chains[i];
		for(int j=0;j<chain_x.residues.size();j++)
		{
			Residue residue_x = chain_x.residues[j];
			for(int t=0; t<residue_x.atoms.size();t++)
			{
				Atom atom_x = residue_x.atoms[t];
				vector<float> cartX = atom_x.xyzVect;
				vector<float> fracX ;
				MatrixTimesTransVector(c2f,cartX,fracX);
				vector<float> atm_idx(3,fracX[0]*grid[0] - origin[0] + 1);
				atm_idx[1] = fracX[1]*grid[1] - origin[1] + 1;
				atm_idx[2] = fracX[2]*grid[2] - origin[2] + 1;
				vector<float> atm_j(3,0.0), del_ij(3,0.0);

				if ( atm_idx[0]<1.0 || atm_idx[0]>grid[0] ) continue;
				if ( atm_idx[1]<1.0 || atm_idx[1]>grid[1] ) continue;
				if ( atm_idx[2]<1.0 || atm_idx[2]>grid[2] ) continue;

				std::string elt_i = name2eltx(atom_x.atona);
				OneGaussianScattering sig_j = get_A( elt_i );
				float per_atom_ks = sig_j.k( atom_x.tempfac, B_upper_limit );
				if ( use_Bs ) {
					per_atom_ks = std::min ( per_atom_ks, (float) (4*M_PI*M_PI/minimumB ));
				} else {
					per_atom_ks = std::min ( per_atom_ks, (float) (4*M_PI*M_PI/effectiveB ));
				}				
				float C= sig_j.C( per_atom_ks);
				float k = per_atom_ks;

				float ATOM_MASK_SQ=0.0;
				ATOM_MASK_SQ = (1.0/k) * ( std::log( C ) - std::log(1e-4) );
				// dist cutoff for density
				float ATOM_DENS_SQ = (1.0/k) * (std::log( C ) - std::log(1e-4));
				// for b factor minimization, we don't want the mask to change.
				// so we fix the mask at a high B value				
				if ( fixed_mask_B>0 ) {
					float mK = sig_j.k( fixed_mask_B, B_upper_limit );
					float mC = sig_j.C( mK );
					ATOM_MASK_SQ = (1.0/mK) * (std::log( mC ) - std::log(1e-4));
				}
				// force mask
				if ( force_mask>0 ) {
					ATOM_MASK_SQ = force_mask*force_mask;
				}

				for(int z=1;z<=density.u3();z++)
				{
					atm_j[2] = z;
					del_ij[2] =(atm_idx[2]-atm_j[2])/grid[2];
					vector<float> frac_tmpz;
					MatrixTimesTransVector(f2c,del_ij,frac_tmpz);
					if(square_len(frac_tmpz)> (ATOM_MASK_SQ) ) continue;					
					for(int y=1;y<=density.u2();y++)
					{
						atm_j[1] = y;
						del_ij[1] = (atm_idx[1] - atm_j[1])/grid[1];
						vector<float> frac_tmpy;
						MatrixTimesTransVector(f2c,del_ij,frac_tmpy);
						if(square_len(frac_tmpy)> (ATOM_MASK_SQ) ) continue;						
						for(int x=1;x<=density.u1();x++)
						{
							atm_j[0] = x;
							del_ij[0] = (atm_idx[0] - atm_j[0])/grid[0];
							vector<float> cart_del_ij;
							MatrixTimesTransVector(f2c,del_ij,cart_del_ij);
							float d2 = square_len(cart_del_ij);
						
							if ( d2 <= (ATOM_MASK_SQ) ) {
								mask(x,y,z) = 1.0; // problem?
								if ( d2 <= (ATOM_DENS_SQ) ) {
									float atm = C*exp(-k*d2);
									rhoC(x,y,z) += atm;
								}
							}
						}

					}
				}				

			}
		}
	}
	cout<<endl;

	if ( highres_limit == 0 ) return;

	// bandlimit mask at 'radius'
	ObjexxFCL::FArray3D< std::complex<double> > Fmask;
	fourier::fft3(mask, Fmask);
	//int H,K,L;
	for ( int z=1; z<=(int)grid[2]; ++z ) {
		int H = (z < (int)grid[2]/2) ? z-1 : z-grid[2] - 1;
		for ( int y=1; y<=(int)grid[1]; ++y ) {
			int K = (y < (int)grid[1]/2) ? y-1 : y-grid[1]-1;
			for ( int x=1; x<=(int)grid[0]; ++x ) {
				int L = (x < (int)grid[0]/2) ? x-1 : x-grid[0]-1;
				float S2c =  S2(H,K,L);

				// exp fade
				float scale = exp(-S2c*(highres_limit*highres_limit));
				Fmask(x,y,z) *= scale;
			}
		}
	}
	fourier::ifft3(Fmask, mask);
 }


// ElectronDensity::readMRC(std::istream mapfile)
//      read an MRC/CCP4 density map
bool
ElectronDensity::readMRCandResize(
	std::istream & mapin,
	std::string mapfile,
	float reso,
	float gridSpacing
) {
	char mapString[4], symData[81];

	int  crs2xyz[3], extent[3], mode, symBytes, grid[3], origin[3];
	int  xyz2crs[3], vol_xsize, vol_ysize, vol_zsize;
	int xIndex, yIndex, zIndex, vol_xySize, coord[3];
	unsigned long long dataOffset, filesize;
	float *rowdata;

	cellDimensions=vector<float>(3,0.0);
	cellAngles=vector<float>(3,0.0);

	bool swap=false;

	if ( !mapin ) {
		cout << "Error opening MRC map " << mapfile << ".  Not loading map." << std::endl;
		return false;
	}

	if (    !mapin.read(reinterpret_cast <char*> (extent), 3*sizeof(int))
			|| !mapin.read(reinterpret_cast <char*> (&mode), 1*sizeof(int))
			|| !mapin.read(reinterpret_cast <char*> (&origin[0]), 3*sizeof(int))
			|| !mapin.read(reinterpret_cast <char*> (&grid[0]), 3*sizeof(int))
			|| !mapin.read(reinterpret_cast <char*> (&cellDimensions[0]), 3*sizeof(float))
			|| !mapin.read(reinterpret_cast <char*> (&cellAngles[0]), 3*sizeof(float))
			|| !mapin.read(reinterpret_cast <char*> (crs2xyz), 3*sizeof(int)) )  {
		cout << "Improperly formatted line in MRC map.  Not loading map." << std::endl;
		return false;
	}

	// Check the number of bytes used for storing symmetry operators
	mapin.seekg(92, std::ios::beg);
	if ( !mapin.read(reinterpret_cast <char*> (&symBytes), 1*sizeof(int)) ) {
		cout << "Failed reading symmetry bytes record.  Not loading map." << "\n";
		return false;
	}

	// alt: MRC files have floating-point origin info at byte 196
	// read this and try to figure out if it is used
	float altorigin[3];
	mapin.seekg(196, std::ios::beg);
	if ( !mapin.read(reinterpret_cast <char*> (altorigin), 3*sizeof(float)) ) {
		cout << "Improperly formatted line in MRC map.  Not loading map." << std::endl;
		return false;
	}

	// Check for the string "MAP" at byte 208, indicating a CCP4 file.
	mapin.seekg(208, std::ios::beg);
	mapString[3] = '\0';
	if ( !mapin.read(mapString, 3) || (std::string(mapString) != "MAP") ) {
		cout << "'MAP' string missing, not a valid MRC map.  Not loading map." << std::endl;
		return false;
	}
	cout<<"model: "<<mode<<endl;
	// Check the file endianness
	if ( mode != 2 ) {
		swap4_aligned(&mode, 1);
		if ( mode != 2 ) {
			cout << "Non-real (32-bit float) data types are unsupported.  Not loading map." << std::endl;
			return false;
		} else {
			swap = true; // enable byte swapping
		}
	}

	// Swap all the information obtained from the header
	if ( swap ) {
		swap4_aligned(extent, 3);
		swap4_aligned(&origin[0], 3);
		swap4_aligned(&altorigin[0], 3);
		swap4_aligned(&grid[0], 3);
		swap4_aligned(&cellDimensions[0], 3);
		swap4_aligned(&cellAngles[0], 3);
		swap4_aligned(crs2xyz, 3);
		swap4_aligned(&symBytes, 1);
	}

	if ( reso != 0 ) cout << " Setting resolution to " << reso << "A" << std::endl;
	else cout << " Setting resolution to AUTO" << std::endl;
	cout << "          atom mask to " << ATOM_MASK << "A" << std::endl;
	cout << "            CA mask to " << CA_MASK << "A" << std::endl;
	cout << " Read density map'" << mapfile << "'" << std::endl;
	cout << "     extent: " << extent[0] << " x " << extent[1] << " x " << extent[2] << std::endl;
	cout << "     origin: " << origin[0] << " x " << origin[1] << " x " << origin[2] << std::endl;
	cout << "  altorigin: " << altorigin[0] << " x " << altorigin[1] << " x " << altorigin[2] << std::endl;
	cout << "       grid: " << grid[0] << " x " << grid[1] << " x " << grid[2] << std::endl;
	cout << "    celldim: " << cellDimensions[0] << " x " << cellDimensions[1] << " x " << cellDimensions[2] << std::endl;
	cout << " cellangles: " << cellAngles[0] << " x " << cellAngles[1] << " x " << cellAngles[2] << std::endl;
	cout << " crs2xyz: " << crs2xyz[0] << " x " << crs2xyz[1] << " x " << crs2xyz[2] << std::endl;

	// Check the dataOffset: this fixes the problem caused by files claiming
	// to have symmetry records when they do not.
	mapin.seekg(0, std::ios::end);
	filesize = mapin.tellg();
	dataOffset = filesize - 4L *
		(static_cast<long long>(extent[0]) * static_cast<long long>(extent[1]) * static_cast<long long>(extent[2]));
	if ( dataOffset != static_cast<unsigned long long>(1024 + symBytes) ) {
		if ( dataOffset == static_cast<unsigned long long>(1024) ) {
			// Bogus symmetry record information
			cout << "File contains bogus symmetry record.  Continuing." << std::endl;
			symBytes = 0;
		} else if ( dataOffset < static_cast<unsigned long long>(1024) ) {
			cout << "File appears truncated and doesn't match header.  Not loading map." << std::endl;
			return false;
		} else if ( (dataOffset > static_cast<unsigned long long>(1024)) && (dataOffset < (1024*1024)) ) {
			// Fix for loading SPIDER files which are larger than usual
			// In this specific case, we must absolutely trust the symBytes record
			dataOffset = 1024 + symBytes;
			cout << "Warning,File is larger than expected and doesn't match header.  Reading anyway." <<
				std::endl;
		} else {
			cout << "File is MUCH larger than expected and doesn't match header.  Not loading map." <<
				std::endl;
			cout << dataOffset  << std::endl;
			return false;
		}
	}

	// Read symmetry records -- organized as 80-byte lines of text.
//	utility::vector1< std::string > symList;
	vector<string> symList;
	symData[80]='\0';
	if ( symBytes != 0 ) {
		cout << "Symmetry records found:" << std::endl;
		mapin.seekg(1024, std::ios::beg);
		for ( int i = 0; i < symBytes/80; i++ ) {
			mapin.read(symData, 80);
			symList.push_back(symData);
			cout << symData << std::endl;
		}
	} else {
		// no symm info; assume P 1
		symList.push_back("X,  Y,  Z");
	}
	initializeSymmOps( symList );

	// check extent and grid interval counts
	if ( grid[0] == 0 && extent[0] > 0 ) {
		grid[0] = extent[0] - 1;
		cout << "Warning,Fixed X interval count.  Continuing." << std::endl;
	}
	if ( grid[1] == 0 && extent[1] > 0 ) {
		grid[1] = extent[1] - 1;
		cout << "Warning, Fixed Y interval count.  Continuing." << std::endl;
	}
	if ( grid[2] == 0 && extent[2] > 0 ) {
		grid[2] = extent[2] - 1;
		cout << "Warning, Fixed Z interval count.  Continuing." << std::endl;
	}

	// Mapping between CCP4 column, row, section and Cartesian x, y, z.
	if ( crs2xyz[0] == 0 && crs2xyz[1] == 0 && crs2xyz[2] == 0 ) {
		cout << "Warning, All crs2xyz records are zero." << std::endl;
		cout << "Warning, Setting crs2xyz to 1, 2, 3 and continuing." << std::endl;
		crs2xyz[0] = 1; crs2xyz[1] = 2; crs2xyz[2] = 3;
	}

	xyz2crs[crs2xyz[0]-1] = 0; xyz2crs[crs2xyz[1]-1] = 1; xyz2crs[crs2xyz[2]-1] = 2;
	xIndex = xyz2crs[0]; yIndex = xyz2crs[1]; zIndex = xyz2crs[2];

	vol_xsize = extent[xIndex];
	vol_ysize = extent[yIndex];
	vol_zsize = extent[zIndex];
	vol_xySize = vol_xsize * vol_ysize;

	// coord = <col, row, sec>
	// extent = <colSize, rowSize, secSize>
	rowdata = new float[extent[0]];
	mapin.seekg(dataOffset, std::ios::beg);

	// 'alloc' changes ordering of "extent"
	density.dimension( vol_xsize,vol_ysize,vol_zsize );

	for ( coord[2] = 1; coord[2] <= extent[2]; coord[2]++ ) {
		for ( coord[1] = 1; coord[1] <= extent[1]; coord[1]++ ) {
			// Read an entire row of data from the file, then write it into the
			// datablock with the correct slice ordering.
			if ( mapin.eof() ) {
				cout << "Error, Unexpected end-of-file. Not loading map." << std::endl;
				delete [] rowdata;
				return false;
			}
			if ( mapin.fail() ) {
				cout << "Error, Problem reading the file. Not loading map." << std::endl;
				delete [] rowdata;
				return false;
			}
			if ( !mapin.read( reinterpret_cast< char* >(rowdata), sizeof(float)*extent[0]) ) {
				cout << "Error reading data row. Not loading map." << std::endl;
				delete [] rowdata;
				return false;
			}


			for ( coord[0] = 1; coord[0] <= extent[0]; coord[0]++ ) {
				density( coord[xyz2crs[0]], coord[xyz2crs[1]], coord[xyz2crs[2]]) = rowdata[coord[0]-1];
//				cout<< rowdata[coord[0]-1];
			}
		}
	}

//	cout<<"rowdata: "<<rowdata[32]<<endl;

	if ( swap == 1 ) {
		swap4_aligned( &density[0], vol_xySize * vol_zsize);
	}
	delete [] rowdata;


	this->origin=vector<float>(3,0.0);
	this->origin[0] =(float) origin[xyz2crs[0]];
	this->origin[1] =(float) origin[xyz2crs[1]];
	this->origin[2] =(float) origin[xyz2crs[2]];
//	this->origin.push_back((float)origin[xyz2crs[0]]); //error
//	this->origin.push_back((float)origin[xyz2crs[1]]);
//	this->origin.push_back((float)origin[xyz2crs[2]]);
//	cout<<"this grid: "<<this->origin[0]<<" "<<this->origin[1]<<" "<<this->origin[2]<<endl;

	// grid doesnt seemed to get remapped in ccp4 maps
	this->grid=vector<int>(3,1);
	this->grid[0] = grid[0];
	this->grid[1] = grid[1];
	this->grid[2] = grid[2];
//	this->grid.push_back(grid[0]); // error
//	this->grid.push_back(grid[1]);
//	this->grid.push_back(grid[2]);

	///////////////////////////////////
	/// POST PROCESSING
	// expand to unit cell
	this->computeCrystParams();
	cout << " voxel vol.: " << voxel_volume() << std::endl;
	// mrc format maps occasionally specify a real-valued origin in a different spot in the header
	if (  altorigin[0]!=0 &&  altorigin[1]!=0 &&  altorigin[2]!=0 &&
			( altorigin[0] > -10000 && altorigin[0] < 10000) &&
			( altorigin[1] > -10000 && altorigin[1] < 10000) &&
			( altorigin[2] > -10000 && altorigin[2] < 10000)
			) {
		this->origin[0] = altorigin[xyz2crs[0]];
		this->origin[1] = altorigin[xyz2crs[1]];
		this->origin[2] = altorigin[xyz2crs[2]];
//		numeric::xyzVector<core::Real> fracX = c2f*(this->origin);
		vector<float> fracX;
		MatrixTimesTransVector(c2f,(this->origin),fracX);
//		this->origin = numeric::xyzVector<core::Real>( fracX[0]*grid[0] , fracX[1]*grid[1] , fracX[2]*grid[2] );
		this->origin[0] = fracX[0]*grid[0];
		this->origin[1] = fracX[1]*grid[1];
		this->origin[2] = fracX[2]*grid[2];

		cout << "Using ALTERNATE origin\n";
		cout << "     origin =" << this->origin[0] << " x " << this->origin[1] << " x " << this->origin[2] << std::endl;

		use_altorigin = true;
	} else {
		use_altorigin = false;
	}

	this->expandToUnitCell();

	// optionally, if force_apix_on_map_load_ is specified
	// override input apix with supplied value
	if ( force_apix_on_map_load_ > 0 ) {
		set_voxel_spacing( force_apix_on_map_load_ );
	}

	// resample the map
	if ( gridSpacing > 0 ) {
		this->resize( gridSpacing );
	}

	// Resolution logic is a bit convoluted
	//   We define a "minimum" resolution and an "effective" resolution
	//   The minimum resolution is the absolute minimum resolution component we can represent (grid_sampling/2)
	//   Effective resolution is the "actual" resolution of the data, and comes from -map_reso (if provided)
	//      or is guessed at grid_sampling
	float max_del_grid = std::max( cellDimensions[0]/((double)this->grid[0]) , cellDimensions[1]/((double)this->grid[1]) );
	max_del_grid = std::max( max_del_grid ,(float) (cellDimensions[2]/((double)this->grid[2])) );
	minimumB = 16*max_del_grid*max_del_grid;

	if ( reso == 0 ) {
		max_del_grid *= 1.5;
	} else {
		max_del_grid = std::max( max_del_grid, reso/2 );
	}

	cout << "Effective resolution = " << 2*max_del_grid << std::endl;
	effectiveB = 16*max_del_grid*max_del_grid;

	// mask width logic is also convoluted
	// make sure mask extends >= 3 carbon stdevs at effectiveB
	float mask_min = 3.0 * sqrt( effectiveB / (2*M_PI*M_PI) );
//	if ( basic::options::option[ basic::options::OptionKeys::edensity::atom_mask_min ].user() ) {
//		mask_min = basic::options::option[ basic::options::OptionKeys::edensity::atom_mask_min ];
//	}

	if ( ATOM_MASK < mask_min ) {
		ATOM_MASK = mask_min;
	}
	mask_min = 3.0 * sqrt(2.0) * std::max( (float)(2.4+1.6*reso) , reso ) / M_PI;
	if ( CA_MASK < mask_min ) {
		CA_MASK = mask_min;
	}

	// set map resolution
	this->reso = 2*max_del_grid;

	// update/clear derived data
	density_change_trigger();

//	ElectronDensity(density,1.0, vector<float>(3,0.0), false ).writeMRC( "rho_before_expand1.mrc");
//	test_density=density;

	// we're done!
	isLoaded = true;
	return isLoaded;
}

void ElectronDensity::density_change_trigger() {
	fastdens_score.clear();
	fastgrid = vector<int>(3,0);

	Fdensity.clear();
	coeffs_density_.clear();

	computeGradients();  // visualization only
}

void spline_coeffs(
	ObjexxFCL::FArray3D< double > const & data ,
	ObjexxFCL::FArray3D< double > & coeffs) {
	int dims[3] = { data.u3(), data.u2(), data.u1() };
	coeffs = data;
	SplineInterp::compute_coefficients3( const_cast<double*>(&coeffs[0]) , dims );  // external code wants nonconst even though array is unchanged
}

void spline_coeffs(
	ObjexxFCL::FArray3D< float > const & data ,
	ObjexxFCL::FArray3D< double > & coeffs) {
	int N = data.u3()*data.u2()*data.u1();
	ObjexxFCL::FArray3D< double > data_d(data.u1(),data.u2(),data.u3()) ;
	for ( int i=0; i<N; ++i ) {
		data_d[i] = (double)data[i];
	}
	spline_coeffs( data_d, coeffs );
}

/////////////////////////////////////
//  compute gradients of density --- used to calculate surface normals in visualizer
void ElectronDensity::computeGradients() {
#ifdef GL_GRAPHICS
	vector<int> dims(3,density.u1());
	dims[1] = density.u2();
	dims[2] = density.u3();

	// Allocate arrays
	ObjexxFCL::FArray3D< float > grad_x, grad_y, grad_z;
	grad_x.dimension( dims[0] , dims[1] , dims[2] ); grad_x = 0;
	grad_y.dimension( dims[0] , dims[1] , dims[2] ); grad_y = 0;
	grad_z.dimension( dims[0] , dims[1] , dims[2] ); grad_z = 0;

	for (int x=2; x<dims[0]; ++x) {
		for (int y=2; y<dims[1]; ++y) {
			for (int z=2; z<dims[2]; ++z) {
				grad_x(x,y,z) = 0.5*(density(x+1,y,z) - density(x-1,y,z));
				grad_x(x,y,z) = 0.125*(density(x+1,y+1,z) - density(x-1,y+1,z));
				grad_x(x,y,z) = 0.125*(density(x+1,y-1,z) - density(x-1,y-1,z));
				grad_x(x,y,z) = 0.125*(density(x+1,y,z+1) - density(x-1,y,z+1));
				grad_x(x,y,z) = 0.125*(density(x+1,y,z-1) - density(x-1,y,z-1));

				grad_y(x,y,z) = 0.5*(density(x,y+1,z) - density(x,y-1,z));
				grad_y(x,y,z) = 0.125*(density(x+1,y+1,z) - density(x+1,y-1,z));
				grad_y(x,y,z) = 0.125*(density(x-1,y+1,z) - density(x-1,y-1,z));
				grad_y(x,y,z) = 0.125*(density(x,y+1,z+1) - density(x,y-1,z+1));
				grad_y(x,y,z) = 0.125*(density(x,y+1,z-1) - density(x,y-1,z-1));

				grad_z(x,y,z) = 0.5*(density(x,y,z+1) - density(x,y,z-1));
				grad_z(x,y,z) = 0.125*(density(x+1,y,z+1) - density(x+1,y,z-1));
				grad_z(x,y,z) = 0.125*(density(x-1,y,z+1) - density(x-1,y,z-1));
				grad_z(x,y,z) = 0.125*(density(x,y+1,z+1) - density(x,y+1,z-1));
				grad_z(x,y,z) = 0.125*(density(x,y-1,z+1) - density(x,y-1,z-1));
			}
		}
	}
	// spline coeff of dCOdx's
	spline_coeffs( grad_x , coeff_grad_x );
	spline_coeffs( grad_y , coeff_grad_y );
	spline_coeffs( grad_z , coeff_grad_z );

	cout << "Finished computing gradient maps" << std::endl;
#endif
}

template<class S, class T>
void resample(ObjexxFCL::FArray3D< S > const &density, ObjexxFCL::FArray3D< T > &newDensity,vector<int> newDims)
{
	if ( density.u1() == newDims[0] && density.u2() == newDims[1] && density.u3() == newDims[2] ) {
		newDensity = density;
		return;
	}

	newDensity.dimension( newDims[0], newDims[1], newDims[2] );

	// convert map to complex<double>
	ObjexxFCL::FArray3D< std::complex<double> > Foldmap, Fnewmap;
	Fnewmap.dimension( newDims[0], newDims[1], newDims[2] );

	// fft
	ObjexxFCL::FArray3D< std::complex<double> > cplx_density;
	cplx_density.dimension( density.u1() , density.u2() , density.u3() );
	for ( int i=0; i<density.u1()*density.u2()*density.u3(); ++i ) cplx_density[i] = (double)density[i];
	fourier::fft3(cplx_density, Foldmap);	 

	// reshape (handles both shrinking and growing in each dimension)
	for ( int i=0; i<Fnewmap.u1()*Fnewmap.u2()*Fnewmap.u3(); ++i ) Fnewmap[i] = std::complex<double>(0.0,0.0);

	int tmp_t=(int)std::min(Foldmap.u1(), Fnewmap.u1())/2;
	vector<int> nyq(3,tmp_t);
	nyq[1]=(int)std::min(Foldmap.u2(), Fnewmap.u2())/2;
	nyq[2]=(int)std::min(Foldmap.u3(), Fnewmap.u3())/2;

	int tmp_k=std::max(Foldmap.u1() - (std::min(Foldmap.u1(),Fnewmap.u1())-nyq[0]) + 1 , nyq[0]+1);
	vector<int> nyqplus1_old(3,tmp_k);
	nyqplus1_old[1]=std::max(Foldmap.u2() - (std::min(Foldmap.u2(),Fnewmap.u2())-nyq[1]) + 1 , nyq[1]+1);
	nyqplus1_old[2]=std::max(Foldmap.u3() - (std::min(Foldmap.u3(),Fnewmap.u3())-nyq[2]) + 1 , nyq[2]+1);

	int tmp_u = std::max(Fnewmap.u1() - (std::min(Foldmap.u1(),Fnewmap.u1())-nyq[0]) + 1 , nyq[0]+1);
	vector<int> nyqplus1_new(3,tmp_u);
	nyqplus1_new[1] = std::max(Fnewmap.u2() - (std::min(Foldmap.u2(),Fnewmap.u2())-nyq[1]) + 1 , nyq[1]+1);
	nyqplus1_new[2] = std::max(Fnewmap.u3() - (std::min(Foldmap.u3(),Fnewmap.u3())-nyq[2]) + 1 , nyq[2]+1);

	for ( int i=1; i<=Fnewmap.u1(); i++ ) {
		for ( int j=1; j<=Fnewmap.u2(); j++ ) {
			for ( int k=1; k<=Fnewmap.u3(); k++ ) {
				if ( i-1<=nyq[0] ) {
					if ( j-1<=nyq[1] ) {
						if ( k-1<=nyq[2] ) {
							Fnewmap(i,j,k) = Foldmap(i, j, k);
						} else if ( k-1>=nyqplus1_new[2] ) {
							Fnewmap(i,j,k) = Foldmap(i, j, k-nyqplus1_new[2]+nyqplus1_old[2]);
						}

					} else if ( j-1>=nyqplus1_new[1] ) {
						if ( k-1<=nyq[2] ) {
							Fnewmap(i,j,k) = Foldmap(i, j-nyqplus1_new[1]+nyqplus1_old[1], k);
						} else if ( k-1>=nyqplus1_new[2] ) {
							Fnewmap(i,j,k) = Foldmap(i, j-nyqplus1_new[1]+nyqplus1_old[1], k-nyqplus1_new[2]+nyqplus1_old[2]);
						}
					}
				} else if ( i-1>=nyqplus1_new[0] ) {
					if ( j-1<=nyq[1] ) {
						if ( k-1<=nyq[2] ) {
							Fnewmap(i,j,k) = Foldmap(i-nyqplus1_new[0]+nyqplus1_old[0],j,k);
						} else if ( k-1>=nyqplus1_new[2] ) {
							Fnewmap(i,j,k) = Foldmap(i-nyqplus1_new[0]+nyqplus1_old[0],j, k-nyqplus1_new[2]+nyqplus1_old[2]);
						}

					} else if ( j-1>=nyqplus1_new[1] ) {
						if ( k-1<=nyq[2] ) {
							Fnewmap(i,j,k) = Foldmap(i-nyqplus1_new[0]+nyqplus1_old[0],j-nyqplus1_new[1]+nyqplus1_old[1],k);
						} else if ( k-1>=nyqplus1_new[2] ) {
							Fnewmap(i,j,k) = Foldmap(i-nyqplus1_new[0]+nyqplus1_old[0],
								j-nyqplus1_new[1]+nyqplus1_old[1],
								k-nyqplus1_new[2]+nyqplus1_old[2]);
						}
					}
				}
			}
		}
	}

	// ifft
	fourier::ifft3(Fnewmap, newDensity);		
}

/////////////////////////////////////
// resize a map (using FFT-interpolation)
void ElectronDensity::resize( float approxGridSpacing ) {
	// potentially expand map to cover entire unit cell
	if ( grid[0] != density.u1() || grid[1] != density.u2() || grid[2] != density.u3() ) {
		cout << "Error, resize() not supported for maps not covering the entire unit cell."<< std::endl;
		cout << "   " << grid[0] << " != " << density.u1()
			<< " || " << grid[1] << " != " << density.u2()
			<< " || " << grid[2] << " != " << density.u3() << std::endl;
//		utility_exit();
		return;
	}

	// compute new dimensions & origin
	vector<int> newDims(3,0),  newGrid(3,0);
//	numeric::xyzVector<int> newDims,  newGrid;
//	numeric::xyzVector<double> newOri;
	vector<float> newOri(3,0.0);

	//fpd since we're doing a bunch of FFTs now resize this to something with no large prime factors
	newDims[0] = findSampling( float(cellDimensions[0] / approxGridSpacing ), MINMULT[0] );
	newDims[1] = findSampling( float(cellDimensions[1] / approxGridSpacing ), MINMULT[1] );
	newDims[2] = findSampling( float(cellDimensions[2] / approxGridSpacing ), MINMULT[2] );

	newOri[0] = newDims[0]*origin[0] / ((float)grid[0]);
	newOri[1] = newDims[1]*origin[1] / ((float)grid[1]);
	newOri[2] = newDims[2]*origin[2] / ((float)grid[2]);
	newGrid = newDims;

	ObjexxFCL::FArray3D< double > newDensity;
	ObjexxFCL::FArray3D< double > xdensity(density.u1(),density.u2(),density.u3());
/*	for(int i=1;i<=density.u1();i++)
	{
		for(int j=1;j<=density.u2();j++)
		{
			for(int k=1;k<=density.u3();k++)
			{
				xdensity = double (density(i,j,k));
			}
		}
	}*/

	resample( xdensity, newDensity, newDims );
/*	for(int i=1;i<=density.u1();i++)
	{
		for(int j=1;j<=density.u2();j++)
		{
			for(int k=1;k<=density.u3();k++)
			{
				density = float (xdensity(i,j,k));
			}
		}
	} */	
	cout << "Resizing " << density.u1() << "x" << density.u2() << "x" << density.u3() << " to "
		<< newDensity.u1() << "x" << newDensity.u2() << "x" << newDensity.u3() << std::endl;

	// update density
	density.dimension( newDims[0], newDims[1], newDims[2] );
	for ( int i=0; i< newDims[0]*newDims[1]*newDims[2] ; ++i ) {
		density[i] = (float)newDensity[i];
	}
	this->grid = newGrid;
	this->origin = newOri;

	cout << " new extent: " << density.u1() << " x " << density.u2() << " x " << density.u3() << std::endl;
	cout << " new origin: " << origin[0] << " x " << origin[1] << " x " << origin[2] << std::endl;
	cout << "   new grid: " << grid[0] << " x " << grid[1] << " x " << grid[2] << std::endl;
}

// expand density to cover complete unit cell
// maintain origin
void ElectronDensity::expandToUnitCell() {
//	numeric::xyzVector< int > extent( density.u1(), density.u2(), density.u3() );
	vector<int> extent(3,density.u1());
	extent[1] = density.u2();
	extent[2] = density.u3();

	// if it already covers unit cell do nothing
	if ( grid[0] == extent[0] && grid[1] == extent[1] && grid[2] == extent[2] ) {
		return;
	}

	ObjexxFCL::FArray3D< float > newDensity( grid[0],grid[1],grid[2], 0.0 );

	// copy the block
	int limX=std::min(extent[0],grid[0]),
		limY=std::min(extent[1],grid[1]),
		limZ=std::min(extent[2],grid[2]);
	for ( int x=1; x<=limX; ++x ) {
		for ( int y=1; y<=limY; ++y ) {
			for ( int z=1; z<=limZ; ++z ) {
				newDensity( x,y,z ) = density( x,y,z );
			}
		}
	}

	// apply symmetry
	// why backwards? it is a mystery
	for ( int x=grid[0]; x>=1; --x ) {
		for ( int y=grid[1]; y>=1; --y ) {
			for ( int z=grid[2]; z>=1; --z ) {
				if ( x <= limX && y <= limY && z <= limZ ) {
					continue;
				}

				vector<float> fracX=vector<float>(3,0.0);
				fracX[0] = ( (float)x + origin[0] - 1 ) / grid[0];
				fracX[1] = ( (float)y + origin[1] - 1 ) / grid[1];
				fracX[2] = ( (float)z + origin[2] - 1 ) / grid[2];

				for ( int symm_idx=0; symm_idx<(int)symmOps.size(); symm_idx++ ) {
					//numeric::xyzVector<core::Real> SfracX =
					vector<float> tmp_fracX;
					MatrixTimesTransVector(symmOps[symm_idx].Rotation,fracX,tmp_fracX);
					vector<float> SfracX(3,0.0);
					vectorsum(tmp_fracX,symmOps[symm_idx].Translation,SfracX);
					//vector<float> SfracX =
						 //tmp_fracX +  symmOps[symm_idx].Translation;
						//symmOps[symm_idx].Rotation * fracX +  symmOps[symm_idx].Translation;

						//symmOps[symm_idx].get_rotation() * fracX +  symmOps[symm_idx].get_translation();

					// indices of symm copy
					int Sx = pos_mod((int)floor(SfracX[0]*grid[0]+0.5 - origin[0]) , grid[0]) + 1;
					int Sy = pos_mod((int)floor(SfracX[1]*grid[1]+0.5 - origin[1]) , grid[1]) + 1 ;
					int Sz = pos_mod((int)floor(SfracX[2]*grid[2]+0.5 - origin[2]) , grid[2]) + 1 ;

					if ( Sx <= limX && Sy <= limY && Sz <= limZ ) {
						newDensity( x,y,z ) = density(Sx,Sy,Sz);
					}
				}
			}
		}
	}

	de_edenisty=pde_edensity;
	if ( de_edenisty ) {
		ElectronDensity(density,1.0, vector<float>(3,0.0), false ).writeMRC( "rho_before_expand.mrc");
		ElectronDensity(newDensity,1.0, vector<float>(3,0.0), false ).writeMRC( "rho_after_expand.mrc");
	}
//	if ( basic::options::option[ basic::options::OptionKeys::edensity::debug ]() ) {
//		ElectronDensity(density,1.0, vector<float>(3,0.0), false ).writeMRC( "rho_before_expand.mrc");
//		ElectronDensity(newDensity,1.0, vector<float>(3,0), false ).writeMRC( "rho_after_expand.mrc");
//	}

	// new map!
	density = newDensity;
}

bool ElectronDensity::writeMRC(std::string mapfilename) {
	std::fstream outx( (mapfilename).c_str() , std::ios::binary | std::ios::out );

	float buff_f;
	int buff_i;
	float buff_vf[3];
	int buff_vi[3];
	int symBytes = 0;

	if ( !outx ) {
		cout << "Error opening MRC map for writing." << std::endl;
		return false;
	}

	// extent
	buff_vi[0] = density.u1(); buff_vi[1] = density.u2(); buff_vi[2] = density.u3();
	outx.write(reinterpret_cast <char*>(buff_vi), sizeof(int)*3);

	// mode
	buff_i = 2;
	outx.write(reinterpret_cast <char*>(&buff_i), sizeof(int));

	// origin
	int ori_int[3];
	ori_int[0] = use_altorigin? 0 : (int)std::floor( origin[0] );
	ori_int[1] = use_altorigin? 0 : (int)std::floor( origin[1] );
	ori_int[2] = use_altorigin? 0 : (int)std::floor( origin[2] );
	outx.write(reinterpret_cast <char*>(ori_int), sizeof(int)*3);

	// grid
	outx.write(reinterpret_cast <char*>(&grid[0]), sizeof(int)*3);

	// cell params
	outx.write(reinterpret_cast <char*>(&cellDimensions), sizeof(float)*3);
	outx.write(reinterpret_cast <char*>(&cellAngles), sizeof(float)*3);

	// crs2xyz
	buff_vi[0] = 1; buff_vi[1] = 2; buff_vi[2] = 3;
	outx.write(reinterpret_cast <char*>(buff_vi), sizeof(int)*3);

	// min, max, mean dens
	buff_vf[0] = -100.0; buff_vf[1] = 100.0; buff_vf[2] = 0.0;
	outx.write(reinterpret_cast <char*>(buff_vf), sizeof(float)*3);

	// 4 bytes junk
	buff_i = 0;
	outx.write(reinterpret_cast <char*>(&buff_i), sizeof(int));

	// symmops (to do!)
	outx.write(reinterpret_cast <char*>(&symBytes), sizeof(int));

	// 104 bytes junk
	buff_i = 0;
	for ( int i=0; i<25; ++i ) {
		outx.write(reinterpret_cast <char*>(&buff_i), sizeof(int));
	}

	// alt origin (MRC)
	float ori_float[3]={0,0,0};
//	numeric::xyzVector<core::Real> origin_realspace(0,0,0);
	vector<float> origin_realspace(3,0);
	if ( use_altorigin ) {
		idx2cart (vector<float>(3,1), origin_realspace );
		ori_float[0] = (float)( origin_realspace[0] );
		ori_float[1] = (float)( origin_realspace[1] );
		ori_float[2] = (float)( origin_realspace[2] );
	}
	outx.write(reinterpret_cast <char*>(ori_float), sizeof(float)*3);

	// Write "MAP" at byte 208, indicating a CCP4 file.
	char buff_s[80]; strcpy(buff_s, "MAP DD  ");
	outx.write(reinterpret_cast <char*>(buff_s), 8);

	// fill remainder of head with junk
	int nJunkWords = (1024 - 216) /4;
	buff_i = 0;
	for ( int i=0; i<nJunkWords; ++i ) {
		outx.write(reinterpret_cast <char*>(&buff_i), sizeof(int));
	}

	// data
	int coord[3];
	for ( coord[2] = 1; coord[2] <= density.u3(); coord[2]++ ) {
		for ( coord[1] = 1; coord[1] <= density.u2(); coord[1]++ ) {
			for ( coord[0] = 1; coord[0] <= density.u1(); coord[0]++ ) {
				buff_f = (float) density(coord[0],coord[1],coord[2]);
				outx.write(reinterpret_cast <char*>(&buff_f), sizeof(float));
			}
		}
	}

	return true;
}

void
ElectronDensity::set_voxel_spacing( vector<float> apix ) {
	// 1) scale origin
	//origin[0] *= (apix[0] * this->grid[0]) / cellDimensions[0];
	//origin[1] *= (apix[1] * this->grid[1]) / cellDimensions[1];
	//origin[2] *= (apix[2] * this->grid[2]) / cellDimensions[2];

	// 2) scale cell dimensions
	cellDimensions[0] = apix[0] * this->grid[0];
	cellDimensions[1] = apix[1] * this->grid[1];
	cellDimensions[2] = apix[2] * this->grid[2];

	// 3) update f2c,c2f
	this->computeCrystParams();

	cout << "Forcing apix to " << apix[0] << "," << apix[1] << "," << apix[2] << std::endl;
	//TR << "new origin:" << this->origin[0] << "," << this->origin[1] << "," << this->origin[2] << std::endl;

	//fpd no need for full density change trigger, just invalidate a few things
	// Fdensity and coeffs_density_ are still valid
	fastdens_score.clear();
	fastgrid = vector<int>(3,0);
}


void ElectronDensity::initializeSymmOps( vector< string > const & symList ) {
//	using core::kinematics::RT;
	RT RT_sym;

	RT_sym.Rotation = vector<vector<float> >(3,vector<float>(3,0.0));
	RT_sym.Rotation[0][0]=1.0;
	RT_sym.Rotation[1][1]=1.0;
	RT_sym.Rotation[2][2]=1.0;
	RT_sym.Translation = vector<float>(3,0.0);	

	symmOps.clear();

	if ( symList.size() == 0 ) { // no symminfo in header, assume P 1
		symmOps.push_back( RT_sym );
	}

	MINMULT=vector<int>(3,2);
	MINMULT[0] = MINMULT[1] = MINMULT[2] = 2;
	for ( int i=0; i<(int)symList.size(); ++i ) {
		std::string line = symList[i];
//		utility::vector1< std::string > rows = utility::string_split(line, ',');
		vector< string > rows = string_splitx(line, ',');
		if ( rows.size() != 3 ) {
			cout << "Error, invalid symmop in map file" << std::endl;
			cout << line << std::endl;
			cout << "Error, Setting symmetry to P1 and continuing!" << line << std::endl;

			// should we throw an exception here????  nah, just set symm to P1 and continue
			symmOps.clear();
			symmOps.push_back( RT_sym );

			return;
		}

		// _REALLY_ simple parser
		//int k;
		for ( int j=0; j<3; ++j ) {
			int k = rows[j].find('/');
			if ( k != (int)std::string::npos ) {
				// numerator
				int startDenom = k+1;
				float denom = std::atof( &rows[j][startDenom]);

				// make sure this shift corresponds to a point in the map
				int oldMinMult = MINMULT[j-1];
				while ( std::fmod(MINMULT[j-1] , denom) > 1e-6 ) MINMULT[j-1] += oldMinMult;
			}
		}
	}

	return;
}

/////////////////////////////////////
// void ElectronDensity::computeCrystParams()
void ElectronDensity::computeCrystParams() {
	// recompute reciprocal cell
	// f2c, c2f
	float ca = cos(d2r(cellAngles[0])), cb = cos(d2r(cellAngles[1])), cg = cos(d2r(cellAngles[2]));
	float sa = sin(d2r(cellAngles[0])), sb = sin(d2r(cellAngles[1])), sg = sin(d2r(cellAngles[2]));

	// conversion from fractional cell coords to cartesian coords
	this->f2c=vector<vector<float> >(3, vector<float>(3,0.0));
	this->f2c[0][0] = cellDimensions[0];
	this->f2c[0][1] = cellDimensions[1] * cg;
	this->f2c[0][2] = cellDimensions[2] * cb;
	this->f2c[1][1] = cellDimensions[1] * sg;
	this->f2c[1][2] = cellDimensions[2] * (ca - cb*cg) / sg;
	this->f2c[2][2] = cellDimensions[2] * sb * sqrt(1.0 - square((cb*cg - ca)/(sb*sg)));

	float D = f2c[0][0] * ( (f2c[1][1]*f2c[2][2]) - (f2c[2][1]*f2c[1][2]) ) -
		 	 f2c[0][1] * ( (f2c[1][0]*f2c[2][2]) - (f2c[2][0]*f2c[1][2]) ) +
		 	 f2c[0][2] * ( (f2c[1][0]*f2c[2][1]) - (f2c[2][0]*f2c[1][1]) );
	if(D == 0 ) {
		cout<< "Warning, Invalid crystal cell dimensions. " << endl;
		return;
	}
		
	// c2f is inverse of f2c
	this->c2f=vector<vector<float> >(3, vector<float>(3,0.0));
	this->c2f[0][0] = (f2c[1][1] * f2c[2][2]-f2c[1][2] * f2c[2][1])/D;
	this->c2f[0][1] = -(f2c[0][1] * f2c[2][2]-f2c[0][2] * f2c[2][1])/D;
	this->c2f[0][2] = (f2c[0][1] * f2c[1][2]-f2c[0][2] * f2c[1][1])/D;
	this->c2f[1][0] = -(f2c[1][0] * f2c[2][2]-f2c[2][0] * f2c[1][2])/D;
	this->c2f[1][1] = (f2c[0][0] * f2c[2][2]-f2c[0][2] * f2c[2][0])/D;
	this->c2f[1][2] = -(f2c[0][0] * f2c[1][2]-f2c[0][2] * f2c[1][0])/D;
	this->c2f[2][0] = (f2c[1][0] * f2c[2][1]-f2c[2][0] * f2c[1][1])/D;
	this->c2f[2][1] = -(f2c[0][0] * f2c[2][1]-f2c[0][1] * f2c[2][0])/D;
	this->c2f[2][2] = (f2c[0][0] * f2c[1][1]-f2c[0][1] * f2c[1][0])/D;
	this->V = cellDimensions[0] * cellDimensions[1] * cellDimensions[2]* sqrt(1-square(ca)-square(cb)-square(cg)+2*ca*cb*cg);

	// reciprocal space cell dimensions
	RcellDimensions = vector<float>(3,0.0);
	this->RcellDimensions[0] = cellDimensions[1]*cellDimensions[2]*sa/V;
	this->RcellDimensions[1] = cellDimensions[0]*cellDimensions[2]*sb/V;
	this->RcellDimensions[2] = cellDimensions[0]*cellDimensions[1]*sg/V;
	cosRcellAngles = vector<float>(3,0.0);
	this->cosRcellAngles[0] = cos(  asin( std::min( std::max( V/(cellDimensions[0]*cellDimensions[1]*cellDimensions[2]*sb*sg) , (float)-1.0) , (float)1.0) )  );
	this->cosRcellAngles[1] = cos(  asin( std::min( std::max( V/(cellDimensions[0]*cellDimensions[1]*cellDimensions[2]*sa*sg) , (float)-1.0) , (float)1.0) )  );
	this->cosRcellAngles[2] = cos(  asin( std::min( std::max( V/(cellDimensions[0]*cellDimensions[1]*cellDimensions[2]*sa*sb) , (float)-1.0) , (float)1.0) )  );
	this->RV = 1.0/V;
}