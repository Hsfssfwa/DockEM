#include "CaDensityx.h"

//using namespace std;
const float PI = 3.141592653589793;

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

vector<string> string_splitx(string const & in, char splitchar)
{
    vector<string> parts;
    int i=0,j=0;
    while( j!= std::string::npos)
    {
        j=in.find(splitchar,i);
        std::string const part = in.substr(i,j-i);
        parts.push_back( part);
        i=j+1;
    }
    return parts;
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

/*
void idxoffset2cart( vector<float> & idxX , vector<float> &cartX ) {
	vector<float> fracX(3,0.0);
	fracX[0] = ( (float) idxX[0] ) / grid[0];
	fracX[1] = ( (float) idxX[1] ) / grid[1];
	fracX[2] = ( (float) idxX[2] ) / grid[2];
	numeric::xyzVector<core::Real> fracX(
		( (core::Real) idxX[0] ) / grid[0],
		( (core::Real) idxX[1] ) / grid[1],
		( (core::Real) idxX[2] ) / grid[2] );
	cartX = f2c*fracX;
} */

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


//bool pose_has_nonzero_Bs( vector<vector<molecule> > const &moleculex ) {
//	for(int i=0;i<moleculex.size();i++)
//	{
//	    for (int j=0;j<moleculex[i].size();j++)
//	    {
//	        for (int k = 0; k < moleculex[i][j].atoms.size(); i++)
//	        {
//				float tf=moleculex[i][j].atoms[k].temfactor;
//				if(tf>0) return true;           
//	        }
//	    }		
//	}
//
//	return false;
//}
//bool pose_has_nonzero_Bs( vector<molecule> const &moleculex ) {
//    for (int i=0;i<moleculex.size();i++)
//    {
//        for (int j = 0; j < moleculex[i].atoms.size(); j++)
//        {
//			float tf=moleculex[i].atoms[i].temfactor;
//			if(tf>0) return true;           
//        }
//    }
//	return false;
//}

//bool pose_has_nonzero_Bs( molecule const & moleculex ) {
//	for ( int i=0 ; i<(int) moleculex.atoms.size(); ++i ) {
//		if ( moleculex.atoms[i].temfactor > 0 ) {
//			return true;
//		}
//	}
//	return false;
//}


//void ElectronDensity::init() {
//	isLoaded = false;
//
//	grid = vector< int >(3,0);
//	origin = vector< float >(3,0);
//	cellDimensions = vector< float >(3,1.0);
//	cellAngles = vector< float >(3,90.0);
//	use_altorigin =  false;
//
//	// command line overrides defaults
//	reso = preso;
//	ATOM_MASK = pATOM_MASK;
//	ATOM_MASK_PADDING = 2;
//	CA_MASK = pCA_MASK;
//	WINDOW_ = pWINDOW_;
//	score_window_context_ = pscore_window_context_;
//	remap_symm_ = premap_symm_;
//	force_apix_on_map_load_ = pforce_apix_on_map_load_;
//	nkbins_ = pnkbins_;
//
////	cout<<"ATOM_MASK: "<<ATOM_MASK<<endl;
////	cout<<"111311"<<endl;
//	// use-specified B factor (may be overridden)
//	effectiveB = 0;
//}

/// null constructor
ElectronDensity::ElectronDensity() 
{
	init();
}

void ElectronDensity::init() {
	// 进行初始化工作的具体代码
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


 float ElectronDensity::maxNominalRes() {
	float S = (1/sqrt(3.0)) * sqrt(S2( density.u1()/2, density.u2()/2, density.u3()/2 ));
	return 1.0/S;
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

//float ElectronDensity::getRSCCX( ObjexxFCL::FArray3D< float > const &density2,  ObjexxFCL::FArray3D< float > const &mask) {
////	runtime_assert( density.u1()==density2.u1() && density.u2()==density2.u2() && density.u3()==density2.u3() );
//
//	float sumC_i=0, sumO_i=0, sumCO_i=0, vol_i=0, CC_i=0;
//	float sumO2_i=0.0, sumC2_i=0.0, varC_i=0, varO_i=0;
//	float clc_x, obs_x, mask_x;
//	for ( int x=0; x<density.u1()*density.u2()*density.u3(); ++x ) {
//		clc_x = density2[x];
//		obs_x = density[x];
//		mask_x = (mask.u1() == 0) ? 1.0 : mask[x];
//
//		sumCO_i += mask_x*clc_x*mask_x*obs_x; sumO_i  += mask_x*obs_x;
//		sumO2_i += mask_x*obs_x*mask_x*obs_x; sumC_i  += mask_x*clc_x;
//		sumC2_i += mask_x*clc_x*mask_x*clc_x;
//		vol_i   += mask_x;
//	}
//	sumO_i = sumO_i/vol_i;
//	sumC_i = sumC_i/vol_i;
//	sumCO_i = sumCO_i/vol_i;
//	sumC2_i = sumC2_i/vol_i;
//	sumO2_i = sumO2_i/vol_i;
//	varC_i = (sumC2_i - sumC_i*sumC_i);
//	varO_i = (sumO2_i - sumO_i*sumO_i) ;
//	if ( varC_i == 0 || varO_i == 0 ) {
//		CC_i = 0;
//	} else {
//		CC_i = (sumCO_i - sumC_i*sumO_i) / sqrt( varC_i * varO_i );
//	}
//
//	return CC_i;
//}

//void ElectronDensity::calcRhoC(
//	molecule const &moleculex,
//	float highres_limit,
//	ObjexxFCL::FArray3D< float > &rhoC,
//	ObjexxFCL::FArray3D< float > &mask,
//	float fixed_mask_B /* = -1 */,
//	float B_upper_limit /* = 600 */,
//	float force_mask /*=-1*/ ) {
//
//	// get rho_c
//	rhoC.dimension(density.u1() , density.u2() , density.u3());
//	mask.dimension(density.u1() , density.u2() , density.u3());
//	rhoC = 0.0;
//	mask = 0.0;
//
//	bool use_Bs = pose_has_nonzero_Bs( moleculex );
//	if ( !use_Bs ) {
////		cout << "Input pose has no nonzero B factors ... setting to baseline" << std::endl;
//	}
//
//	vector<float> per_atom_ks(moleculex.atoms.size(), 0.0);
//	vector<float> per_atom_Cs(moleculex.atoms.size(), 0.0);
//	vector<float> ATOM_MASK_SQS(moleculex.atoms.size(), 0.0);
//	vector<vector<float>> pose_grid ;//= moleculex;
//	for ( int i=0 ; i<(int)moleculex.atoms.size(); ++i ) {
//		string elt_i = moleculex.atoms[i].atomname_left;
////		cout<<"elt_i: "<<elt_i<<endl;
//		OneGaussianScattering sig_j = get_A( elt_i );
//		per_atom_ks[i] = sig_j.k( moleculex.atoms[i].temfactor, B_upper_limit );
//		if ( use_Bs ) {
//			per_atom_ks[i] = std::min ( per_atom_ks[i], (float) (4*M_PI*M_PI/minimumB) );
//		} else {
//			per_atom_ks[i] = std::min ( per_atom_ks[i], (float) (4*M_PI*M_PI/effectiveB) );
//		}
//		per_atom_Cs[i] = sig_j.C( per_atom_ks[i] );
//
//		float C = per_atom_Cs[i], k = per_atom_ks[i];
//		if ( C < 1e-6 ) continue;
//		ATOM_MASK_SQS[i] = (1.0/k) * ( std::log( C ) - std::log(1e-4) );  // 1e-4 is the density value where the mask goes to 0
//
//		// for b factor minimization, we don't want the mask to change.
//		// so we fix the mask at a high B value
//		if ( fixed_mask_B>0 ) {
//			float mK = sig_j.k( fixed_mask_B, B_upper_limit );
//			float mC = sig_j.C( mK );
//			ATOM_MASK_SQS[i] = (1.0/mK) * (std::log( mC ) - std::log(1e-4));
//		}
//
//		// force mask
//		if ( force_mask>0 ) {
//			ATOM_MASK_SQS[i] = force_mask*force_mask;
//		}
//
//		// TO DO: OPTIONALLY control periodic/nonperiodic boundaries
//		//atm_idx[0] = pos_mod (fracX[0]*grid[0] - origin[0] + 1 , (double)grid[0]);
//		//atm_idx[1] = pos_mod (fracX[1]*grid[1] - origin[1] + 1 , (double)grid[1]);
//		//atm_idx[2] = pos_mod (fracX[2]*grid[2] - origin[2] + 1 , (double)grid[2]);
//		vector<float> cartX(3,0.0); // - getTransform();
//		cartX[0]=moleculex.atoms[i].X;
//		cartX[1]=moleculex.atoms[i].Y;
//		cartX[2]=moleculex.atoms[i].Z;
//	//	vector<float> fracX = c2f*cartX;
//		vector<float> fracX;
//		MatrixTimesTransVector(c2f,cartX,fracX);
//		vector<float> atm_idx(3,fracX[0]*grid[0] - origin[0] + 1);
//		atm_idx[1] = fracX[1]*grid[1] - origin[1] + 1;
//		atm_idx[2] = fracX[2]*grid[2] - origin[2] + 1;
//		pose_grid.push_back(atm_idx) ;
//	}
//
//	// PASS 1: rho_calc
//	for ( int i=1 ; i<=(int)moleculex.atoms.size(); ++i ) {
////		std::cout << "\rrho_calc " << i-1 << " of " << pose.size() << std::flush;
//
//		float C = per_atom_Cs[i-1], k = per_atom_ks[i-1];
//		if ( C < 1e-6 ) continue;
//
//		// dist cutoff for mask
//		float ATOM_MASK_SQ = ATOM_MASK_SQS[i-1];
//
//		// dist cutoff for density
//		float ATOM_DENS_SQ = (1.0/k) * (std::log( C ) - std::log(1e-4));
//
//		vector<float> atm_j(3,0.0), del_ij(3,0.0), atm_idx;
//		atm_idx =  pose_grid[i-1];
//
//		if ( atm_idx[0]<1.0 || atm_idx[0]>grid[0] ) continue;
//		if ( atm_idx[1]<1.0 || atm_idx[1]>grid[1] ) continue;
//		if ( atm_idx[2]<1.0 || atm_idx[2]>grid[2] ) continue;
//
//		for ( int z=1; z<=density.u3(); ++z ) {
//			atm_j[2] = z;
//			del_ij[2] = (atm_idx[2] - atm_j[2]) / grid[2];
//			del_ij[0] = del_ij[1] = 0.0;
//			vector<float> frac_tmpz;
//		//	if ( (f2c*del_ij).length_squared() > (ATOM_MASK_SQ) ) continue;
//			MatrixTimesTransVector(f2c,del_ij,frac_tmpz);
//			if(square_len(frac_tmpz)> (ATOM_MASK_SQ) ) continue;
//			for ( int y=1; y<=density.u2(); ++y ) {
//				atm_j[1] = y;
//				del_ij[1] = (atm_idx[1] - atm_j[1]) / grid[1] ;
//				del_ij[0] = 0.0;
//			//	if ( (f2c*del_ij).length_squared() > (ATOM_MASK_SQ) ) continue;
//				vector<float> frac_tmpy;
//				MatrixTimesTransVector(f2c,del_ij,frac_tmpy);
//				if(square_len(frac_tmpy)> (ATOM_MASK_SQ) ) continue;
//				for ( int x=1; x<=density.u1(); ++x ) {
//					atm_j[0] = x;
//					del_ij[0] = (atm_idx[0] - atm_j[0]) / grid[0];
//					vector<float> cart_del_ij;// = (f2c*del_ij);  // cartesian offset from (x,y,z) to atom_i
//					MatrixTimesTransVector(f2c,del_ij,cart_del_ij);
//				//	float d2 = (cart_del_ij).length_squared();
//					float d2 = square_len(cart_del_ij);
//
//					if ( d2 <= (ATOM_MASK_SQ) ) {
//						mask(x,y,z) = 1.0; // problem?
//						if ( d2 <= (ATOM_DENS_SQ) ) {
//							float atm = C*exp(-k*d2);
//							rhoC(x,y,z) += atm;
//						}
//					}
//				}
//			}
//		}
//	}
////	std::cout << std::endl;
//
//	if ( highres_limit == 0 ) return;
//
//	// bandlimit mask at 'radius'
//	ObjexxFCL::FArray3D< std::complex<double> > Fmask;
//	fourier::fft3(mask, Fmask);
//	//int H,K,L;
//	for ( int z=1; z<=(int)grid[2]; ++z ) {
//		int H = (z < (int)grid[2]/2) ? z-1 : z-grid[2] - 1;
//		for ( int y=1; y<=(int)grid[1]; ++y ) {
//			int K = (y < (int)grid[1]/2) ? y-1 : y-grid[1]-1;
//			for ( int x=1; x<=(int)grid[0]; ++x ) {
//				int L = (x < (int)grid[0]/2) ? x-1 : x-grid[0]-1;
//				float S2c =  S2(H,K,L);
//
//				// exp fade
//				float scale = exp(-S2c*(highres_limit*highres_limit));
//				Fmask(x,y,z) *= scale;
//			}
//		}
//	}
//	fourier::ifft3(Fmask, mask);
//}

//float ElectronDensity::matchposex(vector<molecule> PB)
//{
//    vector<float> del_ijx(3,0.0);
//    vector<float> atm_jx(3,0.0); 
//    string elt_i;   
//    int pnum = PB.size(); 
//    vector<float> atm_idx;
//    ObjexxFCL::FArray3D< float > rhoC;
//    ObjexxFCL::FArray3D< float > inv_rho_mask0;
//    rhoC.dimension(density.u1() , density.u2() , density.u3());
//    inv_rho_mask0.dimension(density.u1() , density.u2() , density.u3());
//    for ( int t=0; t<density.u1()*density.u2()*density.u3(); ++t ) {
//        rhoC[t]=0.0;
//        inv_rho_mask0[t]=1.0;
//    }    
////    cout<<"ATOM_MASK_PADDING,CA_MASK: "<<ATOM_MASK_PADDING<<" "<<CA_MASK<<endl;
//    for(int i=0; i<pnum;i++)
//    {
//    	for(int j=0;j<PB[i].atoms.size();j++)
//    	{
////        cout<<"YYYY: "<<pointsBx[i].elt_<<endl;;
////        if(pointsBx[i].elt_ == "CA")
////        {
//	        atm_idx = vector<float>(3,0.0);
//	        vector<float> cartX1;
//	        vector<float> fracX1;
//	        elt_i = PB[i].atoms[j].atomname_left;
//	        elt_i = elt_i[0];
//	        OneGaussianScattering sig_j = get_A( elt_i );
//	        float k = sig_j.k( effectiveB );
//	        float C = sig_j.C( k );
//	        if ( C < 1e-6 ) continue;  
//
//	        cartX1[0] = PB[i].atoms[j].X; 
//	        cartX1[1] = PB[i].atoms[j].Y; 
//	        cartX1[2] = PB[i].atoms[j].Z;         
//	        cartX1[0] = cartX1[0] - origin[0];
//	        cartX1[1] = cartX1[1] - origin[1];
//	        cartX1[2] = cartX1[2] - origin[2];           
//	        MatrixTimesTransVector(c2f,cartX1,fracX1); 
//	        // the location of atom in grid ?
//	        atm_idx[0] = pos_mod (double(fracX1[0]*grid[0]  + 1) , (double)grid[0]);
//	        atm_idx[1] = pos_mod (double(fracX1[1]*grid[1]  + 1) , (double)grid[1]);
//	        atm_idx[2] = pos_mod (double(fracX1[2]*grid[2]  + 1) , (double)grid[2]);   
//	//        atm_idx[tmp_i][0] = fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1 ;
//	//        atm_idx[tmp_i][1] = fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1 ;
//	//        atm_idx[tmp_i][2] = fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1 ;
//	//        cout<<"RRRR: "<<atm_idx[tmp_i][0]<<" "<<atm_idx[tmp_i][1]<<" "<<atm_idx[tmp_i][2]<<endl; 
//	        for(int z=1;z<=density.u3();z++)
//	        {
//	            atm_jx[2] = z;
//	            del_ijx[2] =(atm_idx[2]-atm_jx[2])/grid[2];
//	            if ( del_ijx[2] > 0.5 ) del_ijx[2]-=1.0;
//	            if ( del_ijx[2] < -0.5 ) del_ijx[2]+=1.0;                    
//	            del_ijx[0] = del_ijx[1] = 0.0;
//	            vector<float> frac_tmpz;
//	            MatrixTimesTransVector(f2c,del_ijx,frac_tmpz);
//	            if(square_len(frac_tmpz)> (ATOM_MASK_PADDING+CA_MASK)*(ATOM_MASK_PADDING+CA_MASK) ) continue;                    
//	            for(int y=1;y<=density.u2();y++)
//	            {
//	                atm_jx[1] = y;
//	                del_ijx[1] = (atm_idx[1] - atm_jx[1])/grid[1];
//	                // wrap-around??
//	                if ( del_ijx[1] > 0.5 ) del_ijx[1]-=1.0;
//	                if ( del_ijx[1] < -0.5 ) del_ijx[1]+=1.0;                        
//	                del_ijx[0] = 0.0;
//	                vector<float> frac_tmpy;
//	                MatrixTimesTransVector(f2c,del_ijx,frac_tmpy);
//	                if(square_len(frac_tmpy)> (ATOM_MASK_PADDING+CA_MASK)*(ATOM_MASK_PADDING+CA_MASK) ) continue;                        
//	                for(int x=1;x<=density.u1();x++)
//	                {
//	                    atm_jx[0] = x;
//	                    del_ijx[0] = (atm_idx[0] - atm_jx[0])/grid[0];
//	                    // wrap-around??
//	                    if ( del_ijx[0] > 0.5 ) del_ijx[0]-=1.0;
//	                    if ( del_ijx[0] < -0.5 ) del_ijx[0]+=1.0;                            
//	                    vector<float> cart_del_ij2;
//	                    MatrixTimesTransVector(f2c,del_ijx,cart_del_ij2);
//	                    float d2 = square_len(cart_del_ij2);
//	                    if(d2 > (ATOM_MASK_PADDING+CA_MASK)*(ATOM_MASK_PADDING+CA_MASK) ) continue;
//	                
//	                    float atm = C*exp(-k*d2);
//	//                            cout<<"MASK: "<<theDensityMap.ATOM_MASK<<endl;
//	                    float sigmoid_msk = exp( d2 - (ATOM_MASK)*(ATOM_MASK)  );
//	                //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.ATOM_MASK)  );
//	                //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.CA_MASK)  );
//	                    float inv_msk = 1/(1+sigmoid_msk);
//	                    rhoC(x,y,z) += atm;
//	                    inv_rho_mask0(x,y,z) *= (1 - inv_msk);
//
//	//                            if ( d2 <= (theDensityMap.ATOM_MASK_SQ) ) {
//	//                                mask2(x,y,z) = 1.0; // problem?
//	//                                if ( d2 <= (theDensityMap.ATOM_DENS_SQ) ) {
//	//                                    float atm = C*exp(-k*d2);
//	//                                    rhoC(x,y,z) += atm;
//	//                                }
//	//                            }
//	                }
//	            }
//	        }
//        } 
////            cout<<" "<<tmp_i<<" "<<endl;                                       
//    }
//
//    float sumC_i2=0.0, sumO_i2=0.0, sumCO_i2=0.0, vol_i2=0.0, CC_i2=0.0;
//    float sumO2_i2=0.0, sumC2_i2=0.0, varC_i2=0.0, varO_i2=0.0;
//    float clc_x2=0.0, obs_x2=0.0, eps_x2=0.0;
//    for ( int x=0; x<density.u1()*density.u2()*density.u3(); ++x ) {
////                tmp_den = 0.0; 
////                tmp_rhc=0.0;
////              if(theDensityMap.density[x]>(2.0/3.0)*max_den) tmp_den=theDensityMap.density[x];
////                if(rhoC2[x]>(2.0/3.0)*max_rhc) tmp_rhc=rhoC2[x];
////                if(theDensityMap.density[x]<(3.0/5.0)*max_den) continue;
//        // fetch this point
//        clc_x2 = rhoC[x];
////                clc_x2 = tmp_rhc;
//        obs_x2 = density[x];
////                obs_x2 = tmp_den;
//        eps_x2 = 1-inv_rho_mask0[x]; //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
////        eps_x2 = 1.0;
//
//        // SMOOTHED
//        sumCO_i2 += eps_x2*clc_x2*obs_x2;
//        sumO_i2  += eps_x2*obs_x2;
//        sumO2_i2 += eps_x2*obs_x2*obs_x2;
//        sumC_i2  += eps_x2*clc_x2;
//        sumC2_i2 += eps_x2*clc_x2*clc_x2;
//        vol_i2   += eps_x2;
//    }
//    varC_i2 = (sumC2_i2 - sumC_i2*sumC_i2 / vol_i2 );
//    varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2 / vol_i2 ) ;
//    if ( varC_i2 == 0 || varO_i2 == 0 ) {
//        CC_i2 = 0;
//    } else {
//        CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2/ vol_i2) / sqrt( varC_i2 * varO_i2 );
//    }
//
//    return CC_i2;        
//}

void ElectronDensity::cart(vector <float> &car , vector <float>& carf,float max)
{
	vector<float> cartX1(3, 0.0), cartX2(3, 0.0), carf1(3,0.0), carf2(3,0.0);
	vector<float> fracX1, fracX2;
	int x, x1, y, y1, z, z1;

	cartX1[0] = car[0];
	cartX1[1] = car[1];
	cartX1[2] = car[2];
	MatrixTimesTransVector(c2f, cartX1, fracX1);
	// the location of atom in grid ?
	carf1[0] = (double(fracX1[0] * grid[0] - origin[0] + 1));
	carf1[1] = (double(fracX1[1] * grid[1] - origin[1] + 1));
	carf1[2] = (double(fracX1[2] * grid[2] - origin[2] + 1));
	

	//cartX2[0] = car[0] + max * 0.5;
	//cartX2[1] = car[1] + max * 0.5;
	//cartX2[2] = car[2] + max * 0.5;
	cartX2[0] = car[0] + max ;
	cartX2[1] = car[1] + max ;
	cartX2[2] = car[2] + max ;
	MatrixTimesTransVector(c2f, cartX2, fracX2);
	// the location of atom in grid ?
	carf2[0] = (double(fracX2[0] * grid[0] - origin[0] + 1));
	carf2[1] = (double(fracX2[1] * grid[1] - origin[1] + 1));
	carf2[2] = (double(fracX2[2] * grid[2] - origin[2] + 1));

	x = (carf1[0] + carf2[0] - carf1[0]);
	if (((carf1[0] + carf2[0] - carf1[0]) - x) > 0.1) x_max = x + 1;
	else x_max = x;
	x1 = (carf1[0] - (carf2[0] - carf1[0]));
	if (((carf1[0] - (carf2[0] - carf1[0])) - x1) < 0.9) x_min = x1;
	else x_min = x1 + 1;
	y = (carf1[1] + carf2[1] - carf1[1]);
	if (((carf1[1] + carf2[1] - carf1[1]) - y) > 0.1) y_max = y + 1;
	else y_max = y;
	y1 = (carf1[1] - (carf2[1] - carf1[1]));
	if (((carf1[1] - (carf2[1] - carf1[1])) - y1) < 0.9) y_min = y1;
	else y_min = y1 + 1;
	z = (carf1[2] + carf2[2] - carf1[2]);
	if (((carf1[2] + carf2[2] - carf1[2]) - z) > 0.1) z_max = z + 1;
	else z_max = z;
	z1 = (carf1[2] - (carf2[2] - carf1[2]));
	if (((carf1[2] - (carf2[2] - carf1[2])) - z1) < 0.9) z_min = z1;
	else z_min = z1 + 1;

	//x_max = (carf1[0] + carf2[0] - carf1[0]);
	//x_min = (carf1[0] - (carf2[0] - carf1[0]));
	//y_max = (carf1[1] + carf2[1] - carf1[1]);
	//y_min = (carf1[1] - (carf2[1] - carf1[1]));
	//z_max = (carf1[2] + carf2[2] - carf1[2]);
	//z_min = (carf1[2] - (carf2[2] - carf1[2]));
	//z_min = (int)(carf1[2] - (carf2[2] - carf1[2]));
	cout << "DDDDDDDDD" << (carf1[0] - (carf2[0] - carf1[0])) << "," << (carf1[1] - (carf2[1] - carf1[1])) << "," << (carf1[2] - (carf2[2] - carf1[2])) << "," << (carf1[0] + carf2[0] - carf1[0]) << "," << (carf1[1] + carf2[1] - carf1[1]) << "," << (carf1[2] + carf2[2] - carf1[2]) << endl;
	cout << "DDDDDDDDD" << x_min <<"," << y_min << "," << z_min << "," << x_max << "," << y_max << "," << z_max << endl;
	cout << "Carf1:  " << carf1[0] << "  " << carf1[1] << "  " << carf1[2] << "  " << endl;
	cout << "Carf2:  " << carf2[0] << "  " << carf2[1] << "  " << carf2[2] << "  " << endl;
}
float ElectronDensity::clCC(DOCKMol lig)
{
	vector<float> del_ijx(3, 0.0);
	vector<float> atm_jx(3, 0.0);
	string elt_i;
	int pnum = 1;
	vector<float> atm_idx;
	ObjexxFCL::FArray3D< float > rhoC;
	ObjexxFCL::FArray3D< float > inv_rho_mask0;
	int Denx = density.u1();
	int Deny = density.u2();
	int Denz = density.u3();
	int Denxyz = Denx * Deny * Denz;
	rhoC.dimension(Denx, Deny, Denz, 0.0);
	inv_rho_mask0.dimension(Denx, Deny, Denz, 1.0);

	//    cout<<"ATOM_MASK_PADDING,CA_MASK: "<<ATOM_MASK_PADDING<<" "<<CA_MASK<<endl;

	for (int i = 0; i < lig.num_atoms; i++)
	{
		//        cout<<"YYYY: "<<pointsBx[i].elt_<<endl;;
		//        if(pointsBx[i].elt_ == "CA")
		//        {
		elt_i = lig.atom_names[i];
		elt_i = elt_i[0];
		{
			atm_idx = vector<float>(3, 0.0);
			vector<float> cartX1(3, 0.0);
			vector<float> fracX1;
			//elt_i = PB[i].elt_;
			//elt_i = elt_i[0];
			OneGaussianScattering sig_j = get_A(elt_i);
			float k = sig_j.k(effectiveB);
			float C = sig_j.C(k);
			if (C < 1e-6) continue;
			//cartX1 = PB[i].x_;
			cartX1[0] = lig.x[i];
			cartX1[1] = lig.y[i];
			cartX1[2] = lig.z[i];
			MatrixTimesTransVector(c2f, cartX1, fracX1);
			// the location of atom in grid ?
			atm_idx[0] = (double(fracX1[0] * grid[0] - origin[0] + 1));
			atm_idx[1] = (double(fracX1[1] * grid[1] - origin[1] + 1));
			atm_idx[2] = (double(fracX1[2] * grid[2] - origin[2] + 1));
			//        atm_idx[tmp_i][0] = fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1 ;
			//        atm_idx[tmp_i][1] = fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1 ;
			//        atm_idx[tmp_i][2] = fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1 ;
			//        cout<<"RRRR: "<<atm_idx[tmp_i][0]<<" "<<atm_idx[tmp_i][1]<<" "<<atm_idx[tmp_i][2]<<endl; 
			
			for (int z = 1; z <= Denz; z++)
			{
				if (z > z_min && z < z_max)
				{
					atm_jx[2] = z;
					del_ijx[2] = (atm_idx[2] - atm_jx[2]) / grid[2];
					del_ijx[0] = del_ijx[1] = 0.0;
					vector<float> frac_tmpz;
					MatrixTimesTransVector(f2c, del_ijx, frac_tmpz);
					if (square_len(frac_tmpz) > (ATOM_MASK_PADDING + CA_MASK) * (ATOM_MASK_PADDING + CA_MASK)) continue;
					for (int y = 1; y <= Deny; y++)
					{
						if (y > y_min && y < y_max)
						{
							atm_jx[1] = y;
							del_ijx[1] = (atm_idx[1] - atm_jx[1]) / grid[1];
							del_ijx[0] = 0.0;
							vector<float> frac_tmpy;
							MatrixTimesTransVector(f2c, del_ijx, frac_tmpy);
							if (square_len(frac_tmpy) > (ATOM_MASK_PADDING + CA_MASK) * (ATOM_MASK_PADDING + CA_MASK)) continue;
							for (int x = 1; x <= Denx; x++)
							{
								if (x > x_min && x < x_max)
								{
									atm_jx[0] = x;
									del_ijx[0] = (atm_idx[0] - atm_jx[0]) / grid[0];
									vector<float> cart_del_ij2;
									MatrixTimesTransVector(f2c, del_ijx, cart_del_ij2);
									float d2 = square_len(cart_del_ij2);
									if (d2 > (ATOM_MASK_PADDING + CA_MASK) * (ATOM_MASK_PADDING + CA_MASK)) continue;
									float atm = C * exp(-k * d2);
									//                            cout<<"MASK: "<<theDensityMap.ATOM_MASK<<endl;
									float sigmoid_msk = exp(d2 - (ATOM_MASK) * (ATOM_MASK));
									//    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.ATOM_MASK)  );
									//    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.CA_MASK)  );
									float inv_msk = 1 / (1 + sigmoid_msk);
									rhoC(x, y, z) += atm;
									inv_rho_mask0(x, y, z) *= (1 - inv_msk);
								}
								else continue;
								//                            if ( d2 <= (theDensityMap.ATOM_MASK_SQ) ) {
								//                                mask2(x,y,z) = 1.0; // problem?
								//                                if ( d2 <= (theDensityMap.ATOM_DENS_SQ) ) {
								//                                    float atm = C*exp(-k*d2);
								//                                    rhoC(x,y,z) += atm;
								//                                }
								//                            }
							}
						}
						else continue;
					}
				}
				else continue;
			}

			//for (int z = static_cast<int>(z_min); z <= z_max; z++)
			//{
			//	atm_jx[2] = z;
			//	del_ijx[2] = (atm_idx[2] - atm_jx[2]) / grid[2];
			//	del_ijx[0] = del_ijx[1] = 0.0;
			//	vector<float> frac_tmpz;
			//	MatrixTimesTransVector(f2c, del_ijx, frac_tmpz);
			//	if (square_len(frac_tmpz) > (ATOM_MASK_PADDING + CA_MASK) * (ATOM_MASK_PADDING + CA_MASK)) continue;
			//	for (int y = static_cast<int>(y_min); y <= y_max; y++)
			//	{
			//		atm_jx[1] = y;
			//		del_ijx[1] = (atm_idx[1] - atm_jx[1]) / grid[1];
			//		del_ijx[0] = 0.0;
			//		vector<float> frac_tmpy;
			//		MatrixTimesTransVector(f2c, del_ijx, frac_tmpy);
			//		if (square_len(frac_tmpy) > (ATOM_MASK_PADDING + CA_MASK) * (ATOM_MASK_PADDING + CA_MASK)) continue;
			//		for (int x = static_cast<int>(x_min); x <= x_max; x++)
			//		{
			//				atm_jx[0] = x;
			//				del_ijx[0] = (atm_idx[0] - atm_jx[0]) / grid[0];
			//				vector<float> cart_del_ij2;
			//				MatrixTimesTransVector(f2c, del_ijx, cart_del_ij2);
			//				float d2 = square_len(cart_del_ij2);
			//				if (d2 > (ATOM_MASK_PADDING + CA_MASK) * (ATOM_MASK_PADDING + CA_MASK)) continue;
			//				float atm = C * exp(-k * d2);
			//				//                            cout<<"MASK: "<<theDensityMap.ATOM_MASK<<endl;
			//				float sigmoid_msk = exp(d2 - (ATOM_MASK) * (ATOM_MASK));
			//				//    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.ATOM_MASK)  );
			//				//    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.CA_MASK)  );
			//				float inv_msk = 1 / (1 + sigmoid_msk);
			//				rhoC(x, y, z) += atm;
			//				inv_rho_mask0(x, y, z) *= (1 - inv_msk);
			//		}
			//	}

			//}
                                  
		}
	}

	float sumC_i2 = 0.0, sumO_i2 = 0.0, sumCO_i2 = 0.0, vol_i2 = 0.0, CC_i2 = 0.0;
	float sumO2_i2 = 0.0, sumC2_i2 = 0.0, varC_i2 = 0.0, varO_i2 = 0.0;
	float clc_x2 = 0.0, obs_x2 = 0.0, eps_x2 = 0.0;
	//for (int x = 0; x < Denxyz; ++x) {
	//	//                tmp_den = 0.0; 
	//	//                tmp_rhc=0.0;
	//	//              if(theDensityMap.density[x]>(2.0/3.0)*max_den) tmp_den=theDensityMap.density[x];
	//	//                if(rhoC2[x]>(2.0/3.0)*max_rhc) tmp_rhc=rhoC2[x];
	//	//                if(theDensityMap.density[x]<(3.0/5.0)*max_den) continue;
	//			// fetch this point
	//	clc_x2 = rhoC[x];
	//	//                clc_x2 = tmp_rhc;
	//	obs_x2 = density[x];
	//	//                obs_x2 = tmp_den;
	//	//cout << "clc:  " << clc_x2 << "  obs:  " << obs_x2 << endl;

	//	eps_x2 = 1 - inv_rho_mask0[x]; //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
	//	//        eps_x2 = 1.0;

	//			// SMOOTHED
	//	sumCO_i2 += eps_x2 * clc_x2 * obs_x2;
	//	sumO_i2 += eps_x2 * obs_x2;
	//	sumO2_i2 += eps_x2 * obs_x2 * obs_x2;
	//	sumC_i2 += eps_x2 * clc_x2;
	//	sumC2_i2 += eps_x2 * clc_x2 * clc_x2;
	//	vol_i2 += eps_x2;
	//}
	int idx = 0;
	/*for (int z = static_cast<int>(z_min); z <= z_max; ++z)
	{
		for (int y = static_cast<int>(y_min); y <= y_max; ++y)
		{
			for (int x = static_cast<int>(x_min); x <= x_max; ++x)
			{
					clc_x2 = rhoC(x,y,z);
					obs_x2 = density(x, y, z);
					eps_x2 = 1 - inv_rho_mask0(x,y,z);
					sumCO_i2 += eps_x2 * clc_x2 * obs_x2;
					sumO_i2 += eps_x2 * obs_x2;
					sumO2_i2 += eps_x2 * obs_x2 * obs_x2;
					sumC_i2 += eps_x2 * clc_x2;
					sumC2_i2 += eps_x2 * clc_x2 * clc_x2;
					vol_i2 += eps_x2;
					idx++;
			}
		}
	}*/

	for (int z = 1; z <= Denz; ++z)
	{
		if (z > z_min && z < z_max)
		{
			for (int y = 1; y <= Deny; ++y)
			{
				if (y > y_min && y < y_max)
				{
					for (int x = 1; x <= Denx; ++x)
					{
						if (x > x_min && x < x_max)
						{
							clc_x2 = rhoC(x, y, z);
							obs_x2 = density(x, y, z);
							eps_x2 = 1 - inv_rho_mask0(x, y, z);
							sumCO_i2 += eps_x2 * clc_x2 * obs_x2;
							sumO_i2 += eps_x2 * obs_x2;
							sumO2_i2 += eps_x2 * obs_x2 * obs_x2;
							sumC_i2 += eps_x2 * clc_x2;
							sumC2_i2 += eps_x2 * clc_x2 * clc_x2;
							vol_i2 += eps_x2;
							idx++;
						}
						else continue;
					}
				}
				else continue;
			}
		}
		else continue;
	}
	//cout << "IDXIDXIDXIDX" << idx << endl;
	varC_i2 = (sumC2_i2 - sumC_i2 * sumC_i2 / vol_i2);
	varO_i2 = (sumO2_i2 - sumO_i2 * sumO_i2 / vol_i2);
	if (varC_i2 == 0 || varO_i2 == 0 || vol_i2 == 0) {
		CC_i2 = 0;
	}
	else {
		CC_i2 = (sumCO_i2 - sumC_i2 * sumO_i2 / vol_i2) / sqrt(varC_i2 * varO_i2);
	}
	return CC_i2;
}
float ElectronDensity::clCC_1(DOCKMol lig)
{
	vector<float> del_ijx(3, 0.0);
	vector<float> atm_jx(3, 0.0);
	string elt_i;
	int pnum = 1;
	vector<float> atm_idx;
	ObjexxFCL::FArray3D< float > rhoC;
	ObjexxFCL::FArray3D< float > inv_rho_mask0;
	int Denx = density.u1();
	int Deny = density.u2();
	int Denz = density.u3();
	int Denxyz = Denx * Deny * Denz;
	rhoC.dimension(Denx, Deny, Denz, 0.0);
	inv_rho_mask0.dimension(Denx, Deny, Denz, 1.0);

	//    cout<<"ATOM_MASK_PADDING,CA_MASK: "<<ATOM_MASK_PADDING<<" "<<CA_MASK<<endl;

	for (int i = 0; i < lig.num_atoms; i++)
	{
		//        cout<<"YYYY: "<<pointsBx[i].elt_<<endl;;
		//        if(pointsBx[i].elt_ == "CA")
		//        {
		elt_i = lig.atom_names[i];
		elt_i = elt_i[0];
		if(elt_i != "H")
		{
			atm_idx = vector<float>(3, 0.0);
			vector<float> cartX1(3, 0.0);
			vector<float> fracX1;
			//elt_i = PB[i].elt_;
			//elt_i = elt_i[0];
			OneGaussianScattering sig_j = get_A(elt_i);
			float k = sig_j.k(effectiveB);
			float C = sig_j.C(k);
			if (C < 1e-6) continue;

			//cartX1 = PB[i].x_;
			cartX1[0] = lig.x[i];
			cartX1[1] = lig.y[i];
			cartX1[2] = lig.z[i];
			MatrixTimesTransVector(c2f, cartX1, fracX1);
			// the location of atom in grid ?
			atm_idx[0] = (double(fracX1[0] * grid[0] - origin[0] + 1));
			atm_idx[1] = (double(fracX1[1] * grid[1] - origin[1] + 1));
			atm_idx[2] = (double(fracX1[2] * grid[2] - origin[2] + 1));
			//        atm_idx[tmp_i][0] = fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1 ;
			//        atm_idx[tmp_i][1] = fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1 ;
			//        atm_idx[tmp_i][2] = fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1 ;
			//        cout<<"RRRR: "<<atm_idx[tmp_i][0]<<" "<<atm_idx[tmp_i][1]<<" "<<atm_idx[tmp_i][2]<<endl; 
			for (int z = z_min; z <= z_max; z++)
			//for (int z = 1; z <= Denz; z++)
			{
				{
					atm_jx[2] = z;
					del_ijx[2] = (atm_idx[2] - atm_jx[2]) / grid[2];
					del_ijx[0] = del_ijx[1] = 0.0;
					vector<float> frac_tmpz;
					MatrixTimesTransVector(f2c, del_ijx, frac_tmpz);
					if (square_len(frac_tmpz) > (ATOM_MASK_PADDING + CA_MASK) * (ATOM_MASK_PADDING + CA_MASK)) continue;
					for (int y = y_min; y <= y_max; y++)
					//for (int y = 1; y <= Deny; y++)
					{
						{
							atm_jx[1] = y;
							del_ijx[1] = (atm_idx[1] - atm_jx[1]) / grid[1];
							del_ijx[0] = 0.0;
							vector<float> frac_tmpy;
							MatrixTimesTransVector(f2c, del_ijx, frac_tmpy);
							if (square_len(frac_tmpy) > (ATOM_MASK_PADDING + CA_MASK) * (ATOM_MASK_PADDING + CA_MASK)) continue;
							for (int x = x_min; x <= x_max; x++)
							//for (int x = 1; x <= Denx; x++)
							{
								{
									atm_jx[0] = x;
									del_ijx[0] = (atm_idx[0] - atm_jx[0]) / grid[0];
									vector<float> cart_del_ij2;
									MatrixTimesTransVector(f2c, del_ijx, cart_del_ij2);
									float d2 = square_len(cart_del_ij2);
									if (d2 > (ATOM_MASK_PADDING + CA_MASK) * (ATOM_MASK_PADDING + CA_MASK)) continue;

									float atm = C * exp(-k * d2);
									//                            cout<<"MASK: "<<theDensityMap.ATOM_MASK<<endl;
									float sigmoid_msk = exp(d2 - (ATOM_MASK) * (ATOM_MASK));
									//    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.ATOM_MASK)  );
									//    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.CA_MASK)  );
									float inv_msk = 1 / (1 + sigmoid_msk);
									rhoC(x, y, z) += atm;
									inv_rho_mask0(x, y, z) *= (1 - inv_msk);
								}
								//                            if ( d2 <= (theDensityMap.ATOM_MASK_SQ) ) {
								//                                mask2(x,y,z) = 1.0; // problem?
								//                                if ( d2 <= (theDensityMap.ATOM_DENS_SQ) ) {
								//                                    float atm = C*exp(-k*d2);
								//                                    rhoC(x,y,z) += atm;
								//                                }
								//                            }
							}
						}
					}
				}
			}
			//            cout<<" "<<tmp_i<<" "<<endl;                                       
		}
	}

	float sumC_i2 = 0.0, sumO_i2 = 0.0, sumCO_i2 = 0.0, vol_i2 = 0.0, CC_i2 = 0.0;
	float sumO2_i2 = 0.0, sumC2_i2 = 0.0, varC_i2 = 0.0, varO_i2 = 0.0;
	float clc_x2 = 0.0, obs_x2 = 0.0, eps_x2 = 0.0;
	//for (int x = 0; x < Denxyz; ++x) {
	//	//                tmp_den = 0.0; 
	//	//                tmp_rhc=0.0;
	//	//              if(theDensityMap.density[x]>(2.0/3.0)*max_den) tmp_den=theDensityMap.density[x];
	//	//                if(rhoC2[x]>(2.0/3.0)*max_rhc) tmp_rhc=rhoC2[x];
	//	//                if(theDensityMap.density[x]<(3.0/5.0)*max_den) continue;
	//			// fetch this point
	//	clc_x2 = rhoC[x];
	//	//                clc_x2 = tmp_rhc;
	//	obs_x2 = density[x];
	//	//                obs_x2 = tmp_den;
	//	//cout << "clc:  " << clc_x2 << "  obs:  " << obs_x2 << endl;
	/*for (int z = 1; z <= Denz; ++z)
	{
		for (int y = 1; y <= Deny; ++y)
		{
			for (int x = 1; x <= Denx; ++x)
			{
				if (x > x_min && x < x_max && y > y_min && y < y_max && z > z_min && z < z_max)
				{
					clc_x2 = rhoC(x, y, z);
					obs_x2 = density(x, y, z);
					cout << "AAAAAAAAAAAAAAAAAAAAAA";
						cout << "clc:  " << clc_x2 << "obs:  " << obs_x2 << endl;
				}
				else
				{
					clc_x2 = rhoC(x, y, z);
					obs_x2 = density(x, y, z);
					cout << "clc:  " << clc_x2 << "obs:  " << obs_x2 << endl;
				}
			}
		}

	}*/
	//	eps_x2 = 1 - inv_rho_mask0[x]; //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
	//	//        eps_x2 = 1.0;

	//			// SMOOTHED
	//	sumCO_i2 += eps_x2 * clc_x2 * obs_x2;
	//	sumO_i2 += eps_x2 * obs_x2;
	//	sumO2_i2 += eps_x2 * obs_x2 * obs_x2;
	//	sumC_i2 += eps_x2 * clc_x2;
	//	sumC2_i2 += eps_x2 * clc_x2 * clc_x2;
	//	vol_i2 += eps_x2;
	//}
	int idx = 0;
	//for (int z = 1; z <= Denz; ++z)
	for (int z = z_min; z <= z_max; z++)
	{
		{
			//for (int y = 1; y <= Deny; ++y)
			for (int y = y_min; y <=y_max; y++)
			{
				{
					for (int x = x_min; x <= x_max; x++)
					//for (int x = 1; x <= Denx; ++x)
					{
						{
							clc_x2 = rhoC(x, y, z);
							obs_x2 = density(x, y, z);
							eps_x2 = 1 - inv_rho_mask0(x, y, z);
							sumCO_i2 += eps_x2 * clc_x2 * obs_x2;
							sumO_i2 += eps_x2 * obs_x2;
							sumO2_i2 += eps_x2 * obs_x2 * obs_x2;
							sumC_i2 += eps_x2 * clc_x2;
							sumC2_i2 += eps_x2 * clc_x2 * clc_x2;
							vol_i2 += eps_x2;
							idx++;
						}
					}
				}
			}
		}
	}
	//cout << "IDXIDXIDXIDX" << idx << endl;
	varC_i2 = (sumC2_i2 - sumC_i2 * sumC_i2 / vol_i2);
	varO_i2 = (sumO2_i2 - sumO_i2 * sumO_i2 / vol_i2);
	if (varC_i2 == 0 || varO_i2 == 0 || vol_i2 == 0) {
		CC_i2 = 0;
	}
	else {
		CC_i2 = (sumCO_i2 - sumC_i2 * sumO_i2 / vol_i2) / sqrt(varC_i2 * varO_i2);
	}
	return CC_i2;
}
float ElectronDensity::clCC_2(DOCKMol lig, vector<poseCoord>lig1)
{
	vector<float> del_ijx(3, 0.0);
	vector<float> atm_jx(3, 0.0);
	string elt_i;
	int pnum = 1;
	vector<float> atm_idx;
	ObjexxFCL::FArray3D< float > rhoC;
	ObjexxFCL::FArray3D< float > inv_rho_mask0;
	int Denx = density.u1();
	int Deny = density.u2();
	int Denz = density.u3();
	int Denxyz = Denx * Deny * Denz;
	rhoC.dimension(Denx, Deny, Denz, 0.0);
	inv_rho_mask0.dimension(Denx, Deny, Denz, 1.0);

	for (int i = 0; i < lig1.size(); i++)
	{
		elt_i = lig1[i].atom_name;
		elt_i = elt_i[0];
		{
			atm_idx = vector<float>(3, 0.0);
			vector<float> cartX1(3, 0.0);
			vector<float> fracX1;
			OneGaussianScattering sig_j = get_A(elt_i);
			float k = sig_j.k(effectiveB);
			float C = sig_j.C(k);
			if (C < 1e-6) continue;

			//cartX1 = PB[i].x_;
			cartX1[0] = lig1[i].x;
			cartX1[1] = lig1[i].y;
			cartX1[2] = lig1[i].z;
			MatrixTimesTransVector(c2f, cartX1, fracX1);
			// the location of atom in grid ?
			atm_idx[0] = (double(fracX1[0] * grid[0] - origin[0] + 1));
			atm_idx[1] = (double(fracX1[1] * grid[1] - origin[1] + 1));
			atm_idx[2] = (double(fracX1[2] * grid[2] - origin[2] + 1));
			for (int z = z_min; z <= z_max; z++)
			{
				{
					atm_jx[2] = z;
					del_ijx[2] = (atm_idx[2] - atm_jx[2]) / grid[2];
					del_ijx[0] = del_ijx[1] = 0.0;
					vector<float> frac_tmpz;
					MatrixTimesTransVector(f2c, del_ijx, frac_tmpz);
					if (square_len(frac_tmpz) > (ATOM_MASK_PADDING + CA_MASK) * (ATOM_MASK_PADDING + CA_MASK)) continue;

					for (int y = y_min; y <= y_max; y++)
					{
						{
							atm_jx[1] = y;
							del_ijx[1] = (atm_idx[1] - atm_jx[1]) / grid[1];
							del_ijx[0] = 0.0;
							vector<float> frac_tmpy;
							MatrixTimesTransVector(f2c, del_ijx, frac_tmpy);
							if (square_len(frac_tmpy) > (ATOM_MASK_PADDING + CA_MASK) * (ATOM_MASK_PADDING + CA_MASK)) continue;

							for (int x = x_min; x <= x_max; x++)
							{
								{
									atm_jx[0] = x;
									del_ijx[0] = (atm_idx[0] - atm_jx[0]) / grid[0];
									vector<float> cart_del_ij2;
									MatrixTimesTransVector(f2c, del_ijx, cart_del_ij2);
									float d2 = square_len(cart_del_ij2);
									if (d2 > (ATOM_MASK_PADDING + CA_MASK) * (ATOM_MASK_PADDING + CA_MASK)) continue;

									float atm = C * exp(-k * d2);
									//                            cout<<"MASK: "<<theDensityMap.ATOM_MASK<<endl;
									float sigmoid_msk = exp(d2 - (ATOM_MASK) * (ATOM_MASK));
									//    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.ATOM_MASK)  );
									//    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.CA_MASK)  );
									float inv_msk = 1 / (1 + sigmoid_msk);
									rhoC(x, y, z) += atm;
									inv_rho_mask0(x, y, z) *= (1 - inv_msk);
								}
							}
						}
					}
				}
			}
		}
	}
	for (int i = 0; i < lig.num_atoms; i++)
	{
		//        cout<<"YYYY: "<<pointsBx[i].elt_<<endl;;
		//        if(pointsBx[i].elt_ == "CA")
		//        {
		elt_i = lig.atom_names[i];
		elt_i = elt_i[0];
		{
			atm_idx = vector<float>(3, 0.0);
			vector<float> cartX1(3, 0.0);
			vector<float> fracX1;
			//elt_i = PB[i].elt_;
			//elt_i = elt_i[0];
			OneGaussianScattering sig_j = get_A(elt_i);
			float k = sig_j.k(effectiveB);
			float C = sig_j.C(k);
			if (C < 1e-6) continue;

			//cartX1 = PB[i].x_;
			cartX1[0] = lig.x[i];
			cartX1[1] = lig.y[i];
			cartX1[2] = lig.z[i];
			MatrixTimesTransVector(c2f, cartX1, fracX1);
			// the location of atom in grid ?
			atm_idx[0] = (double(fracX1[0] * grid[0] - origin[0] + 1));
			atm_idx[1] = (double(fracX1[1] * grid[1] - origin[1] + 1));
			atm_idx[2] = (double(fracX1[2] * grid[2] - origin[2] + 1));
			for (int z = z_min; z <= z_max; z++)
			{
				{
					atm_jx[2] = z;
					del_ijx[2] = (atm_idx[2] - atm_jx[2]) / grid[2];
					del_ijx[0] = del_ijx[1] = 0.0;
					vector<float> frac_tmpz;
					MatrixTimesTransVector(f2c, del_ijx, frac_tmpz);
					if (square_len(frac_tmpz) > (ATOM_MASK_PADDING + CA_MASK) * (ATOM_MASK_PADDING + CA_MASK)) continue;

					for (int y = y_min; y <= y_max; y++)
					{
						{
							atm_jx[1] = y;
							del_ijx[1] = (atm_idx[1] - atm_jx[1]) / grid[1];
							del_ijx[0] = 0.0;
							vector<float> frac_tmpy;
							MatrixTimesTransVector(f2c, del_ijx, frac_tmpy);
							if (square_len(frac_tmpy) > (ATOM_MASK_PADDING + CA_MASK) * (ATOM_MASK_PADDING + CA_MASK)) continue;

							for (int x = x_min; x <= x_max; x++)
							{
								{
									atm_jx[0] = x;
									del_ijx[0] = (atm_idx[0] - atm_jx[0]) / grid[0];
									vector<float> cart_del_ij2;
									MatrixTimesTransVector(f2c, del_ijx, cart_del_ij2);
									float d2 = square_len(cart_del_ij2);
									if (d2 > (ATOM_MASK_PADDING + CA_MASK) * (ATOM_MASK_PADDING + CA_MASK)) continue;

									float atm = C * exp(-k * d2);
									//                            cout<<"MASK: "<<theDensityMap.ATOM_MASK<<endl;
									float sigmoid_msk = exp(d2 - (ATOM_MASK) * (ATOM_MASK));
									//    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.ATOM_MASK)  );
									//    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.CA_MASK)  );
									float inv_msk = 1 / (1 + sigmoid_msk);
									rhoC(x, y, z) += atm;
									inv_rho_mask0(x, y, z) *= (1 - inv_msk);
								}
							}
						}
					}
				}
			}
		}
	}

	float sumC_i2 = 0.0, sumO_i2 = 0.0, sumCO_i2 = 0.0, vol_i2 = 0.0, CC_i2 = 0.0;
	float sumO2_i2 = 0.0, sumC2_i2 = 0.0, varC_i2 = 0.0, varO_i2 = 0.0;
	float clc_x2 = 0.0, obs_x2 = 0.0, eps_x2 = 0.0;

	int idx = 0;
	for (int z = z_min; z <= z_max; z++)
	{
		{
			for (int y = y_min; y <= y_max; y++)
			{
				{

					for (int x = x_min; x <= x_max; x++)
					{
						{
							clc_x2 = rhoC(x, y, z);
							obs_x2 = density(x, y, z);
							eps_x2 = 1 - inv_rho_mask0(x, y, z);
							sumCO_i2 += eps_x2 * clc_x2 * obs_x2;
							sumO_i2 += eps_x2 * obs_x2;
							sumO2_i2 += eps_x2 * obs_x2 * obs_x2;
							sumC_i2 += eps_x2 * clc_x2;
							sumC2_i2 += eps_x2 * clc_x2 * clc_x2;
							vol_i2 += eps_x2;
							idx++;
						}
					}
				}
			}
		}
	}
	varC_i2 = (sumC2_i2 - sumC_i2 * sumC_i2 / vol_i2);
	varO_i2 = (sumO2_i2 - sumO_i2 * sumO_i2 / vol_i2);
	if (varC_i2 == 0 || varO_i2 == 0 || vol_i2 == 0) {
		CC_i2 = 0;
	}
	else {
		CC_i2 = (sumCO_i2 - sumC_i2 * sumO_i2 / vol_i2) / sqrt(varC_i2 * varO_i2);
	}
	return CC_i2;
}
float ElectronDensity::clCC_3(DOCKMol lig, DOCKMol lig1 )
{
	vector<float> del_ijx(3, 0.0);
	vector<float> atm_jx(3, 0.0);
	string elt_i;
	int pnum = 1;
	vector<float> atm_idx;
	ObjexxFCL::FArray3D< float > rhoC;
	ObjexxFCL::FArray3D< float > inv_rho_mask0;
	int Denx = density.u1();
	int Deny = density.u2();
	int Denz = density.u3();
	int Denxyz = Denx * Deny * Denz;
	rhoC.dimension(Denx, Deny, Denz, 0.0);
	inv_rho_mask0.dimension(Denx, Deny, Denz, 1.0);

	for (int i = 0; i < lig1.num_atoms; i++)
	{
		//        cout<<"YYYY: "<<pointsBx[i].elt_<<endl;;
		//        if(pointsBx[i].elt_ == "CA")
		//        {
		elt_i = lig1.atom_names[i];
		elt_i = elt_i[0];
		{
			atm_idx = vector<float>(3, 0.0);
			vector<float> cartX1(3, 0.0);
			vector<float> fracX1;
			//elt_i = PB[i].elt_;
			//elt_i = elt_i[0];
			OneGaussianScattering sig_j = get_A(elt_i);
			float k = sig_j.k(effectiveB);
			float C = sig_j.C(k);
			if (C < 1e-6) continue;

			//cartX1 = PB[i].x_;
			cartX1[0] = lig1.x[i];
			cartX1[1] = lig1.y[i];
			cartX1[2] = lig1.z[i];
			MatrixTimesTransVector(c2f, cartX1, fracX1);
			// the location of atom in grid ?
			atm_idx[0] = (double(fracX1[0] * grid[0] - origin[0] + 1));
			atm_idx[1] = (double(fracX1[1] * grid[1] - origin[1] + 1));
			atm_idx[2] = (double(fracX1[2] * grid[2] - origin[2] + 1));
			for (int z = z_min; z <= z_max; z++)
			{
				{
					atm_jx[2] = z;
					del_ijx[2] = (atm_idx[2] - atm_jx[2]) / grid[2];
					del_ijx[0] = del_ijx[1] = 0.0;
					vector<float> frac_tmpz;
					MatrixTimesTransVector(f2c, del_ijx, frac_tmpz);
					if (square_len(frac_tmpz) > (ATOM_MASK_PADDING + CA_MASK) * (ATOM_MASK_PADDING + CA_MASK)) continue;

					for (int y = y_min; y <= y_max; y++)
					{
						{
							atm_jx[1] = y;
							del_ijx[1] = (atm_idx[1] - atm_jx[1]) / grid[1];
							del_ijx[0] = 0.0;
							vector<float> frac_tmpy;
							MatrixTimesTransVector(f2c, del_ijx, frac_tmpy);
							if (square_len(frac_tmpy) > (ATOM_MASK_PADDING + CA_MASK) * (ATOM_MASK_PADDING + CA_MASK)) continue;

							for (int x = x_min; x <= x_max; x++)
							{
								{
									atm_jx[0] = x;
									del_ijx[0] = (atm_idx[0] - atm_jx[0]) / grid[0];
									vector<float> cart_del_ij2;
									MatrixTimesTransVector(f2c, del_ijx, cart_del_ij2);
									float d2 = square_len(cart_del_ij2);
									if (d2 > (ATOM_MASK_PADDING + CA_MASK) * (ATOM_MASK_PADDING + CA_MASK)) continue;

									float atm = C * exp(-k * d2);
									//                            cout<<"MASK: "<<theDensityMap.ATOM_MASK<<endl;
									float sigmoid_msk = exp(d2 - (ATOM_MASK) * (ATOM_MASK));
									//    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.ATOM_MASK)  );
									//    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.CA_MASK)  );
									float inv_msk = 1 / (1 + sigmoid_msk);
									rhoC(x, y, z) += atm;
									inv_rho_mask0(x, y, z) *= (1 - inv_msk);
								}
							}
						}
					}
				}
			}
		}
	}
	for (int i = 0; i < lig.num_atoms; i++)
	{
		//        cout<<"YYYY: "<<pointsBx[i].elt_<<endl;;
		//        if(pointsBx[i].elt_ == "CA")
		//        {
		elt_i = lig.atom_names[i];
		elt_i = elt_i[0];
		{
			atm_idx = vector<float>(3, 0.0);
			vector<float> cartX1(3, 0.0);
			vector<float> fracX1;
			//elt_i = PB[i].elt_;
			//elt_i = elt_i[0];
			OneGaussianScattering sig_j = get_A(elt_i);
			float k = sig_j.k(effectiveB);
			float C = sig_j.C(k);
			if (C < 1e-6) continue;

			//cartX1 = PB[i].x_;
			cartX1[0] = lig.x[i];
			cartX1[1] = lig.y[i];
			cartX1[2] = lig.z[i];
			MatrixTimesTransVector(c2f, cartX1, fracX1);
			// the location of atom in grid ?
			atm_idx[0] = (double(fracX1[0] * grid[0] - origin[0] + 1));
			atm_idx[1] = (double(fracX1[1] * grid[1] - origin[1] + 1));
			atm_idx[2] = (double(fracX1[2] * grid[2] - origin[2] + 1));
			for (int z = z_min; z <= z_max; z++)
			{
				{
					atm_jx[2] = z;
					del_ijx[2] = (atm_idx[2] - atm_jx[2]) / grid[2];
					del_ijx[0] = del_ijx[1] = 0.0;
					vector<float> frac_tmpz;
					MatrixTimesTransVector(f2c, del_ijx, frac_tmpz);
					if (square_len(frac_tmpz) > (ATOM_MASK_PADDING + CA_MASK) * (ATOM_MASK_PADDING + CA_MASK)) continue;

					for (int y = y_min; y <= y_max; y++)
					{
						{
							atm_jx[1] = y;
							del_ijx[1] = (atm_idx[1] - atm_jx[1]) / grid[1];
							del_ijx[0] = 0.0;
							vector<float> frac_tmpy;
							MatrixTimesTransVector(f2c, del_ijx, frac_tmpy);
							if (square_len(frac_tmpy) > (ATOM_MASK_PADDING + CA_MASK) * (ATOM_MASK_PADDING + CA_MASK)) continue;

							for (int x = x_min; x <= x_max; x++)
							{
								{
									atm_jx[0] = x;
									del_ijx[0] = (atm_idx[0] - atm_jx[0]) / grid[0];
									vector<float> cart_del_ij2;
									MatrixTimesTransVector(f2c, del_ijx, cart_del_ij2);
									float d2 = square_len(cart_del_ij2);
									if (d2 > (ATOM_MASK_PADDING + CA_MASK) * (ATOM_MASK_PADDING + CA_MASK)) continue;

									float atm = C * exp(-k * d2);
									//                            cout<<"MASK: "<<theDensityMap.ATOM_MASK<<endl;
									float sigmoid_msk = exp(d2 - (ATOM_MASK) * (ATOM_MASK));
									//    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.ATOM_MASK)  );
									//    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.CA_MASK)  );
									float inv_msk = 1 / (1 + sigmoid_msk);
									rhoC(x, y, z) += atm;
									inv_rho_mask0(x, y, z) *= (1 - inv_msk);
								}
							}
						}
					}
				}
			}
		}
	}

	float sumC_i2 = 0.0, sumO_i2 = 0.0, sumCO_i2 = 0.0, vol_i2 = 0.0, CC_i2 = 0.0;
	float sumO2_i2 = 0.0, sumC2_i2 = 0.0, varC_i2 = 0.0, varO_i2 = 0.0;
	float clc_x2 = 0.0, obs_x2 = 0.0, eps_x2 = 0.0;

	int idx = 0;
	for (int z = z_min; z <= z_max; z++)
	{
		{
			for (int y = y_min; y <= y_max; y++)
			{
				{

					for (int x = x_min; x <= x_max; x++)
					{
						{
							clc_x2 = rhoC(x, y, z);
							obs_x2 = density(x, y, z);
							eps_x2 = 1 - inv_rho_mask0(x, y, z);
							sumCO_i2 += eps_x2 * clc_x2 * obs_x2;
							sumO_i2 += eps_x2 * obs_x2;
							sumO2_i2 += eps_x2 * obs_x2 * obs_x2;
							sumC_i2 += eps_x2 * clc_x2;
							sumC2_i2 += eps_x2 * clc_x2 * clc_x2;
							vol_i2 += eps_x2;
							idx++;
						}
					}
				}
			}
		}
	}
	varC_i2 = (sumC2_i2 - sumC_i2 * sumC_i2 / vol_i2);
	varO_i2 = (sumO2_i2 - sumO_i2 * sumO_i2 / vol_i2);
	if (varC_i2 == 0 || varO_i2 == 0 || vol_i2 == 0) {
		CC_i2 = 0;
	}
	else {
		CC_i2 = (sumCO_i2 - sumC_i2 * sumO_i2 / vol_i2) / sqrt(varC_i2 * varO_i2);
	}
	return CC_i2;
}
float ElectronDensity::clCC_4( vector<poseCoord>lig1)
{
	vector<float> del_ijx(3, 0.0);
	vector<float> atm_jx(3, 0.0);
	string elt_i;
	int pnum = 1;
	vector<float> atm_idx;
	ObjexxFCL::FArray3D< float > rhoC;
	ObjexxFCL::FArray3D< float > inv_rho_mask0;
	int Denx = density.u1();
	int Deny = density.u2();
	int Denz = density.u3();
	int Denxyz = Denx * Deny * Denz;
	rhoC.dimension(Denx, Deny, Denz, 0.0);
	inv_rho_mask0.dimension(Denx, Deny, Denz, 1.0);

	for (int i = 0; i < lig1.size(); i++)
	{
		elt_i = lig1[i].atom_name;
		elt_i = elt_i[0];
		{
			atm_idx = vector<float>(3, 0.0);
			vector<float> cartX1(3, 0.0);
			vector<float> fracX1;
			OneGaussianScattering sig_j = get_A(elt_i);
			float k = sig_j.k(effectiveB);
			float C = sig_j.C(k);
			if (C < 1e-6) continue;

			//cartX1 = PB[i].x_;
			cartX1[0] = lig1[i].x;
			cartX1[1] = lig1[i].y;
			cartX1[2] = lig1[i].z;
			MatrixTimesTransVector(c2f, cartX1, fracX1);
			// the location of atom in grid ?
			atm_idx[0] = (double(fracX1[0] * grid[0] - origin[0] + 1));
			atm_idx[1] = (double(fracX1[1] * grid[1] - origin[1] + 1));
			atm_idx[2] = (double(fracX1[2] * grid[2] - origin[2] + 1));
			for (int z = z_min; z <= z_max; z++)
			{
				{
					atm_jx[2] = z;
					del_ijx[2] = (atm_idx[2] - atm_jx[2]) / grid[2];
					del_ijx[0] = del_ijx[1] = 0.0;
					vector<float> frac_tmpz;
					MatrixTimesTransVector(f2c, del_ijx, frac_tmpz);
					if (square_len(frac_tmpz) > (ATOM_MASK_PADDING + CA_MASK) * (ATOM_MASK_PADDING + CA_MASK)) continue;

					for (int y = y_min; y <= y_max; y++)
					{
						{
							atm_jx[1] = y;
							del_ijx[1] = (atm_idx[1] - atm_jx[1]) / grid[1];
							del_ijx[0] = 0.0;
							vector<float> frac_tmpy;
							MatrixTimesTransVector(f2c, del_ijx, frac_tmpy);
							if (square_len(frac_tmpy) > (ATOM_MASK_PADDING + CA_MASK) * (ATOM_MASK_PADDING + CA_MASK)) continue;

							for (int x = x_min; x <= x_max; x++)
							{
								{
									atm_jx[0] = x;
									del_ijx[0] = (atm_idx[0] - atm_jx[0]) / grid[0];
									vector<float> cart_del_ij2;
									MatrixTimesTransVector(f2c, del_ijx, cart_del_ij2);
									float d2 = square_len(cart_del_ij2);
									if (d2 > (ATOM_MASK_PADDING + CA_MASK) * (ATOM_MASK_PADDING + CA_MASK)) continue;

									float atm = C * exp(-k * d2);
									//                            cout<<"MASK: "<<theDensityMap.ATOM_MASK<<endl;
									float sigmoid_msk = exp(d2 - (ATOM_MASK) * (ATOM_MASK));
									//    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.ATOM_MASK)  );
									//    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.CA_MASK)  );
									float inv_msk = 1 / (1 + sigmoid_msk);
									rhoC(x, y, z) += atm;
									inv_rho_mask0(x, y, z) *= (1 - inv_msk);
								}
							}
						}
					}
				}
			}
		}
	}

	float sumC_i2 = 0.0, sumO_i2 = 0.0, sumCO_i2 = 0.0, vol_i2 = 0.0, CC_i2 = 0.0;
	float sumO2_i2 = 0.0, sumC2_i2 = 0.0, varC_i2 = 0.0, varO_i2 = 0.0;
	float clc_x2 = 0.0, obs_x2 = 0.0, eps_x2 = 0.0;

	int idx = 0;
	for (int z = z_min; z <= z_max; z++)
	{
		{
			for (int y = y_min; y <= y_max; y++)
			{
				{

					for (int x = x_min; x <= x_max; x++)
					{
						{
							clc_x2 = rhoC(x, y, z);
							obs_x2 = density(x, y, z);
							eps_x2 = 1 - inv_rho_mask0(x, y, z);
							sumCO_i2 += eps_x2 * clc_x2 * obs_x2;
							sumO_i2 += eps_x2 * obs_x2;
							sumO2_i2 += eps_x2 * obs_x2 * obs_x2;
							sumC_i2 += eps_x2 * clc_x2;
							sumC2_i2 += eps_x2 * clc_x2 * clc_x2;
							vol_i2 += eps_x2;
							idx++;
						}
					}
				}
			}
		}
	}
	varC_i2 = (sumC2_i2 - sumC_i2 * sumC_i2 / vol_i2);
	varO_i2 = (sumO2_i2 - sumO_i2 * sumO_i2 / vol_i2);
	if (varC_i2 == 0 || varO_i2 == 0 || vol_i2 == 0) {
		CC_i2 = 0;
	}
	else {
		CC_i2 = (sumCO_i2 - sumC_i2 * sumO_i2 / vol_i2) / sqrt(varC_i2 * varO_i2);
	}
	return CC_i2;
}

float ElectronDensity::matchcentroidposex(DOCKMol lig)
{
    vector<float> del_ijx(3,0.0);
    vector<float> atm_jx(3,0.0); 
    string elt_i;   
    int pnum = 1;
    vector<float> atm_idx;
    ObjexxFCL::FArray3D< float > rhoC;
    ObjexxFCL::FArray3D< float > inv_rho_mask0;
    rhoC.dimension(density.u1() , density.u2() , density.u3(),0.0);
    inv_rho_mask0.dimension(density.u1() , density.u2() , density.u3(),1.0);
    for ( int t=0; t<density.u1()*density.u2()*density.u3(); ++t ) {
        rhoC[t]=0.0;
        inv_rho_mask0[t]=1.0;
    }    
    vector<float> cartX1(3,0.0);
    vector<float> fracX1;
    vector<vector<vector<int> > > rho_dx_pt(pnum);
    vector<vector<vector<vector<float > > > > rho_dx_mask(pnum),rho_dx_atm(pnum);
//    cout<<"ATOM_MASK_PADDING,CA_MASK: "<<ATOM_MASK_PADDING<<" "<<CA_MASK<<endl;
//    cout<<"PPP1"<<endl;
    for(int i=0; i<pnum;i++)
    {
    	int tmp_j=lig.num_atoms;
    	rho_dx_pt[i] = vector<vector<int>>(tmp_j);
    	rho_dx_mask[i]=vector<vector<vector<float>>>(tmp_j);
    	rho_dx_atm[i]=vector<vector<vector<float>>>(tmp_j);
    	for(int j=0;j<lig.num_atoms; j++)
    	{
//        cout<<"YYYY: "<<pointsBx[i].elt_<<endl;;
//        if(pointsBx[i].elt_ == "CA")
//        {
	        atm_idx = vector<float>(3,0.0);

	        elt_i = lig.atom_names[j];
	        elt_i = elt_i[0];
	        OneGaussianScattering sig_j = get_A( elt_i );
	        float k = sig_j.k( effectiveB );
	        float C = sig_j.C( k );
//	         cout<<"YYYY: "<<C<<endl;;
	        if ( C < 1e-6 ) continue; 

	        cartX1[0] = lig.x[j];
	        cartX1[1] = lig.y[j];
	        cartX1[2] = lig.z[j];
//	        cout<<"CCCC: "<<cartX1[0]<<" "<<cartX1[1]<<" "<<cartX1[2]<<endl;     
//	        cartX1[0] = cartX1[0] - origin[0];
//	        cartX1[1] = cartX1[1] - origin[1];
//	        cartX1[2] = cartX1[2] - origin[2];           
	        MatrixTimesTransVector(c2f,cartX1,fracX1); 
	        // the location of atom in grid ?
	        atm_idx[0] = pos_mod (float(fracX1[0]*grid[0]- origin[0]  + 1.0) , (float)grid[0]);
	        atm_idx[1] = pos_mod (float(fracX1[1]*grid[1]- origin[1]  + 1.0) , (float)grid[1]);
	        atm_idx[2] = pos_mod (float(fracX1[2]*grid[2]- origin[2]  + 1.0) , (float)grid[2]); 
//	        atm_idx[0] =  (double(fracX1[0]*grid[0]- origin[0]  + 1) );
//	        atm_idx[1] =  (double(fracX1[1]*grid[1]- origin[1]  + 1) );
//	        atm_idx[2] =  (double(fracX1[2]*grid[2]- origin[2]  + 1) ); 	         
	//        atm_idx[tmp_i][0] = fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1 ;
	//        atm_idx[tmp_i][1] = fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1 ;
	//        atm_idx[tmp_i][2] = fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1 ;
	//        cout<<"RRRR: "<<atm_idx[tmp_i][0]<<" "<<atm_idx[tmp_i][1]<<" "<<atm_idx[tmp_i][2]<<endl; 
	        for(int z=1;z<=density.u3();z++)
	        {
	            atm_jx[2] = z;
	            del_ijx[2] =float(atm_idx[2]-atm_jx[2])/float(grid[2]);
	            if ( del_ijx[2] > 0.5 ) del_ijx[2]-=1.0;
	            if ( del_ijx[2] < -0.5 ) del_ijx[2]+=1.0;                    
	            del_ijx[0] = del_ijx[1] = 0.0;
	            vector<float> frac_tmpz;
	            MatrixTimesTransVector(f2c,del_ijx,frac_tmpz);
	            if(square_len(frac_tmpz)> (ATOM_MASK_PADDING+ATOM_MASK)*(ATOM_MASK_PADDING+ATOM_MASK) ) continue;
	    //        if(square_len(frac_tmpz)> (ATOM_MASK_PADDING+CA_MASK)*(ATOM_MASK_PADDING+CA_MASK) ) continue;                     
				for(int y=1;y<=density.u2();y++)
	            {
	                atm_jx[1] = y;
	                del_ijx[1] = float(atm_idx[1] - atm_jx[1])/float(grid[1]);
	                // wrap-around??
	                if ( del_ijx[1] > 0.5 ) del_ijx[1]-=1.0;
	                if ( del_ijx[1] < -0.5 ) del_ijx[1]+=1.0;                        
	                del_ijx[0] = 0.0;
	                vector<float> frac_tmpy;
	                MatrixTimesTransVector(f2c,del_ijx,frac_tmpy);
	                if(square_len(frac_tmpy)> (ATOM_MASK_PADDING+ATOM_MASK)*(ATOM_MASK_PADDING+ATOM_MASK) ) continue;  
	            //    if(square_len(frac_tmpy)> (ATOM_MASK_PADDING+CA_MASK)*(ATOM_MASK_PADDING+CA_MASK) ) continue; 
	                for(int x=1;x<=density.u1();x++)
	                {
	                    atm_jx[0] = x;
	                    del_ijx[0] = float(atm_idx[0] - atm_jx[0])/float(grid[0]);
	                    // wrap-around??
	                    if ( del_ijx[0] > 0.5 ) del_ijx[0]-=1.0;
	                    if ( del_ijx[0] < -0.5 ) del_ijx[0]+=1.0;                            
	                    vector<float> cart_del_ij2;
	                    MatrixTimesTransVector(f2c,del_ijx,cart_del_ij2);
	                    float d2 = square_len(cart_del_ij2);
	                    if(d2 > (ATOM_MASK_PADDING+ATOM_MASK)*(ATOM_MASK_PADDING+ATOM_MASK) ) continue;
	            //        if(d2 > (ATOM_MASK_PADDING+CA_MASK)*(ATOM_MASK_PADDING+CA_MASK) ) continue;
	                
	                    float atm = C*exp(-k*d2);
	//                            cout<<"MASK: "<<theDensityMap.ATOM_MASK<<endl;
	                    float sigmoid_msk = exp( d2 - (ATOM_MASK)*(ATOM_MASK)  );
	//                    cout<<"ATOM_MASK: "<<ATOM_MASK<<" CA_MASK "<<CA_MASK<<" ATOM_MASK_PADDING "<<ATOM_MASK_PADDING<<endl;
	                //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.ATOM_MASK)  );
	                //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.CA_MASK)  );
	                    float inv_msk = 1.0/(1.0+sigmoid_msk);
	                    rhoC(x,y,z) += atm;
	                    inv_rho_mask0(x,y,z) *= (1.0 - inv_msk);
	                    int idx = (z-1)*density.u2()*density.u1() + (y-1)*density.u1() + x-1;

	                    float eps_i = 1.0-inv_msk, inv_eps_i;
						if ( eps_i == 0 ) { // divide-by-zero
							inv_eps_i = sigmoid_msk;
						} else {
							inv_eps_i = 1.0/eps_i;
						}
						//cout << "677" << endl;
						rho_dx_pt[i][j].push_back(idx);
						//cout << "679" << endl;
						vector<float> cart_del_ijm(3,0),cart_del_ijn(3,0);
				//		cout << "ddd " <<sigmoid_msk<<" "<<inv_msk<<" "<<inv_eps_i<<" "<<" "<<cart_del_ij2[0]<<" "<<cart_del_ij2[1]<<" "<<cart_del_ij2[2]<< endl;
						for(int pp=0;pp<3;pp++)
						{
							cart_del_ijm[pp]=(-2.0*k*atm)*cart_del_ij2[pp];
							cart_del_ijn[pp]=(-2.0*sigmoid_msk*inv_msk*inv_msk*inv_eps_i)*cart_del_ij2[pp];
//							cout << "ddd " <<cart_del_ijm[0]<<" "<<cart_del_ijm[1]<<" "<<cart_del_ijm[2]<<" "<<" "<<cart_del_ijn[0]<<" "<<cart_del_ijn[1]<<" "<<cart_del_ijn[2]<< endl;
						}
						rho_dx_atm[i][j].push_back ( cart_del_ijm);
						rho_dx_mask[i][j].push_back( cart_del_ijn);	                
	//                            if ( d2 <= (theDensityMap.ATOM_MASK_SQ) ) {
	//                                mask2(x,y,z) = 1.0; // problem?
	//                                if ( d2 <= (theDensityMap.ATOM_DENS_SQ) ) {
	//                                    float atm = C*exp(-k*d2);
	//                                    rhoC(x,y,z) += atm;
	//                                }
	//                            }
	                }
	            }
	        }
        } 
//            cout<<" "<<tmp_i<<" "<<endl;                                       
    }
    float sumC_i2=0.0, sumO_i2=0.0, sumCO_i2=0.0, vol_i2=0.0, CC_i2=0.0;
    float sumO2_i2=0.0, sumC2_i2=0.0, varC_i2=0.0, varO_i2=0.0;
    float clc_x2=0.0, obs_x2=0.0, eps_x2=0.0;
    for ( int x=0; x<density.u1()*density.u2()*density.u3(); ++x ) {
//                tmp_den = 0.0; 
//                tmp_rhc=0.0;
//              if(theDensityMap.density[x]>(2.0/3.0)*max_den) tmp_den=theDensityMap.density[x];
//                if(rhoC2[x]>(2.0/3.0)*max_rhc) tmp_rhc=rhoC2[x];
//                if(theDensityMap.density[x]<(3.0/5.0)*max_den) continue;
        // fetch this point
        clc_x2 = rhoC[x];
//                clc_x2 = tmp_rhc;
        obs_x2 = density[x];
//                obs_x2 = tmp_den;
        eps_x2 = 1-inv_rho_mask0[x]; //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
//        eps_x2 = 1.0;

        // SMOOTHED
        sumCO_i2 += eps_x2*clc_x2*obs_x2;
        sumO_i2  += eps_x2*obs_x2;
        sumO2_i2 += eps_x2*obs_x2*obs_x2;
        sumC_i2  += eps_x2*clc_x2;
        sumC2_i2 += eps_x2*clc_x2*clc_x2;
        vol_i2   += eps_x2;
    }
    varC_i2 = (sumC2_i2 - sumC_i2*sumC_i2 / vol_i2 );
    varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2 / vol_i2 ) ;
    if ( varC_i2 == 0 || varO_i2 == 0 || vol_i2==0) {
//    	cout<<"CCCCC"<<endl;
        CC_i2 = 0;
    } else {
//    	cout<<"sumCO_i2:"<<sumCO_i2<<" sumC_i2:"<<sumC_i2 <<" sumO_i2:"<<sumO_i2<<" vol_i2:"<<vol_i2<<" varC_i2:"<<varC_i2<<" varO_i2:"<<varO_i2<<endl;
        CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2/ vol_i2) / sqrt( varC_i2 * varO_i2 );
    }
//    cout<<"PPP"<<endl;
    /// 4  CALCULATE PER-ATOM DERIVATIVES
    
 	//std::map< core::Size , numeric::xyzMatrix< core::Real > > symmRots;
 	dCCdxs_aacen=vector<vector<vector<float>>>(pnum);
	for ( int i=0 ; i<pnum; i++ ) {
		
		int tmp_xx= lig.num_atoms;
		dCCdxs_aacen[i]=vector<vector<float>> (tmp_xx,vector<float>(3,0.0));
//		vector<vector<float>> tmp_op(tmp_xx);
		for(int j=0;j<lig.num_atoms;j++) {
			
			vector<float> dVdx_ij(3,0), dOdx_ij(3,0), dO2dx_ij(3,0), dCOdx_ij(3,0), dC2dx_ij(3,0);

			vector<int> rho_dx_pt_ij   = rho_dx_pt[i][j];
			vector<vector<float>> rho_dx_mask_ij = rho_dx_mask[i][j];
			vector<vector<float>> rho_dx_atm_ij  = rho_dx_atm[i][j];
			
			int npoints = rho_dx_pt_ij.size();
			for ( int n=0; n<npoints; ++n ) {
				const int x(rho_dx_pt_ij[n]);
//				cout<<"VVV: "<<rho_dx_pt_ij<<endl;
				float clc_x = rhoC[x];
				float obs_x = density[x];
				float inv_eps_x = inv_rho_mask0[x];

				vector<float> del_mask(3,0); 
				for(int tt=0;tt<3;tt++)
				{
					del_mask[tt]=inv_eps_x*rho_dx_mask_ij[n][tt];
				}
				vector<float> del_rhoc = rho_dx_atm_ij[n];

				for(int tt=0;tt<3;tt++)
				{
					dVdx_ij[tt]  += del_mask[tt];
					dOdx_ij[tt]  += obs_x*del_mask[tt];
					dO2dx_ij[tt] += obs_x*obs_x*del_mask[tt];
					dCOdx_ij[tt] += obs_x*del_rhoc[tt];
					dC2dx_ij[tt] += 2.0*del_rhoc[tt]*clc_x;					
				}


			}

			// finally compute dCC/dx_ij
			float f = ( sumCO_i2 - sumC_i2*sumO_i2 / vol_i2 );
			float g = sqrt ( varO_i2 * varC_i2 );

			vector<float> fprime(3,0.0);
			vector<float> gprime(3,0.0);
			for(int kk=0;kk<3;kk++)
			{
				fprime[kk] = dCOdx_ij[kk] - 1/(vol_i2*vol_i2) * ( dOdx_ij[kk]*sumC_i2*vol_i2 - sumO_i2*sumC_i2*dVdx_ij[kk]);
				gprime[kk] = 0.5 * (sqrt(varO_i2)/sqrt(varC_i2) * ( dC2dx_ij[kk] + ( sumC_i2*sumC_i2*dVdx_ij[kk]/(vol_i2*vol_i2) ) )  +
				sqrt(varC_i2)/sqrt(varO_i2) * ( dO2dx_ij[kk] - ( 1/(vol_i2*vol_i2) * ( 2*vol_i2*sumO_i2*dOdx_ij[kk] - sumO_i2*sumO_i2*dVdx_ij[kk] ) ) ) );
			}
//			cout << "mmm" << endl;
			vector<float> tmp_vct(3,0.0);
			for(int kk=0;kk<3;kk++)
			{
				tmp_vct[kk]=(-1.0)*(g*fprime[kk] - f*gprime[kk]) / (g*g);//  梯度为-CC的梯度
//				cout << "ddd " <<fprime[0]<<" "<<fprime[1]<<" "<<fprime[2]<<" "<<" "<<gprime[0]<<" "<<gprime[1]<<" "<<gprime[2]<< endl;
			}
//			tmp_op[j]=tmp_vct;
			dCCdxs_aacen[i][j]=tmp_vct;
//			cout << "ddd" <<tmp_vct[0]<<" "<<tmp_vct[1]<<" "<<tmp_vct[2]<<" "<< endl;
		}
//		dCCdxs_aacen.push_back(tmp_op);
	}
//	cout<<"PPP2"<<endl;
    return CC_i2;        
}
float ElectronDensity::matchcentroidposexz(DOCKMol lig)
{
	vector<float> del_ijx(3, 0.0);
	vector<float> atm_jx(3, 0.0);
	string elt_i;
	int pnum = 1;
	vector<float> atm_idx;
	ObjexxFCL::FArray3D< float > rhoC;
	ObjexxFCL::FArray3D< float > inv_rho_mask0;
	rhoC.dimension(density.u1(), density.u2(), density.u3(), 0.0);
	inv_rho_mask0.dimension(density.u1(), density.u2(), density.u3(), 1.0);
	for (int t = 0; t < density.u1() * density.u2() * density.u3(); ++t) {
		rhoC[t] = 0.0;
		inv_rho_mask0[t] = 1.0;
	}
	vector<float> cartX1(3, 0.0);
	vector<float> fracX1;
	vector<vector<vector<int> > > rho_dx_pt(pnum);
	vector<vector<vector<vector<float > > > > rho_dx_mask(pnum), rho_dx_atm(pnum);
	//    cout<<"ATOM_MASK_PADDING,CA_MASK: "<<ATOM_MASK_PADDING<<" "<<CA_MASK<<endl;
	//    cout<<"PPP1"<<endl;
	for (int i = 0; i < pnum; i++)
	{
		int tmp_j = lig.num_atoms;
		rho_dx_pt[i] = vector<vector<int>>(tmp_j);
		rho_dx_mask[i] = vector<vector<vector<float>>>(tmp_j);
		rho_dx_atm[i] = vector<vector<vector<float>>>(tmp_j);
		for (int j = 0; j < lig.num_atoms; j++)
		{
			//        cout<<"YYYY: "<<pointsBx[i].elt_<<endl;;
			//        if(pointsBx[i].elt_ == "CA")
			//        {
			atm_idx = vector<float>(3, 0.0);

			elt_i = lig.atom_names[j];
			elt_i = elt_i[0];
			OneGaussianScattering sig_j = get_A(elt_i);
			float k = sig_j.k(effectiveB);
			float C = sig_j.C(k);
			//	         cout<<"YYYY: "<<C<<endl;;
			if (C < 1e-6) continue;

			cartX1[0] = lig.x[j];
			cartX1[1] = lig.y[j];
			cartX1[2] = lig.z[j];
			//	        cout<<"CCCC: "<<cartX1[0]<<" "<<cartX1[1]<<" "<<cartX1[2]<<endl;     
			//	        cartX1[0] = cartX1[0] - origin[0];
			//	        cartX1[1] = cartX1[1] - origin[1];
			//	        cartX1[2] = cartX1[2] - origin[2];           
			MatrixTimesTransVector(c2f, cartX1, fracX1);
			// the location of atom in grid ?
			atm_idx[0] = pos_mod(float(fracX1[0] * grid[0] - origin[0] + 1.0), (float)grid[0]);
			atm_idx[1] = pos_mod(float(fracX1[1] * grid[1] - origin[1] + 1.0), (float)grid[1]);
			atm_idx[2] = pos_mod(float(fracX1[2] * grid[2] - origin[2] + 1.0), (float)grid[2]);
			//	        atm_idx[0] =  (double(fracX1[0]*grid[0]- origin[0]  + 1) );
			//	        atm_idx[1] =  (double(fracX1[1]*grid[1]- origin[1]  + 1) );
			//	        atm_idx[2] =  (double(fracX1[2]*grid[2]- origin[2]  + 1) ); 	         
				//        atm_idx[tmp_i][0] = fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1 ;
				//        atm_idx[tmp_i][1] = fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1 ;
				//        atm_idx[tmp_i][2] = fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1 ;
				//        cout<<"RRRR: "<<atm_idx[tmp_i][0]<<" "<<atm_idx[tmp_i][1]<<" "<<atm_idx[tmp_i][2]<<endl; 
			for (int z = z_min; z <= z_max; z++)
			{
				atm_jx[2] = z;
				del_ijx[2] = float(atm_idx[2] - atm_jx[2]) / float(grid[2]);
				if (del_ijx[2] > 0.5) del_ijx[2] -= 1.0;
				if (del_ijx[2] < -0.5) del_ijx[2] += 1.0;
				del_ijx[0] = del_ijx[1] = 0.0;
				vector<float> frac_tmpz;
				MatrixTimesTransVector(f2c, del_ijx, frac_tmpz);
				if (square_len(frac_tmpz) > (ATOM_MASK_PADDING + ATOM_MASK) * (ATOM_MASK_PADDING + ATOM_MASK)) continue;
				//        if(square_len(frac_tmpz)> (ATOM_MASK_PADDING+CA_MASK)*(ATOM_MASK_PADDING+CA_MASK) ) continue;                     
				for (int y = y_min; y <= y_max; y++)
				{
					atm_jx[1] = y;
					del_ijx[1] = float(atm_idx[1] - atm_jx[1]) / float(grid[1]);
					// wrap-around??
					if (del_ijx[1] > 0.5) del_ijx[1] -= 1.0;
					if (del_ijx[1] < -0.5) del_ijx[1] += 1.0;
					del_ijx[0] = 0.0;
					vector<float> frac_tmpy;
					MatrixTimesTransVector(f2c, del_ijx, frac_tmpy);
					if (square_len(frac_tmpy) > (ATOM_MASK_PADDING + ATOM_MASK) * (ATOM_MASK_PADDING + ATOM_MASK)) continue;
					//    if(square_len(frac_tmpy)> (ATOM_MASK_PADDING+CA_MASK)*(ATOM_MASK_PADDING+CA_MASK) ) continue; 
					for (int x = x_min; x <= x_max; x++)
					{
						atm_jx[0] = x;
						del_ijx[0] = float(atm_idx[0] - atm_jx[0]) / float(grid[0]);
						// wrap-around??
						if (del_ijx[0] > 0.5) del_ijx[0] -= 1.0;
						if (del_ijx[0] < -0.5) del_ijx[0] += 1.0;
						vector<float> cart_del_ij2;
						MatrixTimesTransVector(f2c, del_ijx, cart_del_ij2);
						float d2 = square_len(cart_del_ij2);
						if (d2 > (ATOM_MASK_PADDING + ATOM_MASK) * (ATOM_MASK_PADDING + ATOM_MASK)) continue;
						//        if(d2 > (ATOM_MASK_PADDING+CA_MASK)*(ATOM_MASK_PADDING+CA_MASK) ) continue;

						float atm = C * exp(-k * d2);
						//                            cout<<"MASK: "<<theDensityMap.ATOM_MASK<<endl;
						float sigmoid_msk = exp(d2 - (ATOM_MASK) * (ATOM_MASK));
						//                    cout<<"ATOM_MASK: "<<ATOM_MASK<<" CA_MASK "<<CA_MASK<<" ATOM_MASK_PADDING "<<ATOM_MASK_PADDING<<endl;
										//    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.ATOM_MASK)  );
										//    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.CA_MASK)  );
						float inv_msk = 1.0 / (1.0 + sigmoid_msk);
						rhoC(x, y, z) += atm;
						inv_rho_mask0(x, y, z) *= (1.0 - inv_msk);
						int idx = (z - 1) * density.u2() * density.u1() + (y - 1) * density.u1() + x - 1;

						float eps_i = 1.0 - inv_msk, inv_eps_i;
						if (eps_i == 0) { // divide-by-zero
							inv_eps_i = sigmoid_msk;
						}
						else {
							inv_eps_i = 1.0 / eps_i;
						}
						//cout << "677" << endl;
						rho_dx_pt[i][j].push_back(idx);
						//cout << "679" << endl;
						vector<float> cart_del_ijm(3, 0), cart_del_ijn(3, 0);
						//		cout << "ddd " <<sigmoid_msk<<" "<<inv_msk<<" "<<inv_eps_i<<" "<<" "<<cart_del_ij2[0]<<" "<<cart_del_ij2[1]<<" "<<cart_del_ij2[2]<< endl;
						for (int pp = 0; pp < 3; pp++)
						{
							cart_del_ijm[pp] = (-2.0 * k * atm) * cart_del_ij2[pp];
							cart_del_ijn[pp] = (-2.0 * sigmoid_msk * inv_msk * inv_msk * inv_eps_i) * cart_del_ij2[pp];
							//							cout << "ddd " <<cart_del_ijm[0]<<" "<<cart_del_ijm[1]<<" "<<cart_del_ijm[2]<<" "<<" "<<cart_del_ijn[0]<<" "<<cart_del_ijn[1]<<" "<<cart_del_ijn[2]<< endl;
						}
						rho_dx_atm[i][j].push_back(cart_del_ijm);
						rho_dx_mask[i][j].push_back(cart_del_ijn);
						//                            if ( d2 <= (theDensityMap.ATOM_MASK_SQ) ) {
						//                                mask2(x,y,z) = 1.0; // problem?
						//                                if ( d2 <= (theDensityMap.ATOM_DENS_SQ) ) {
						//                                    float atm = C*exp(-k*d2);
						//                                    rhoC(x,y,z) += atm;
						//                                }
						//                            }

					}

				}

			}
		}
		//            cout<<" "<<tmp_i<<" "<<endl;                                       
	}
	float sumC_i2 = 0.0, sumO_i2 = 0.0, sumCO_i2 = 0.0, vol_i2 = 0.0, CC_i2 = 0.0;
	float sumO2_i2 = 0.0, sumC2_i2 = 0.0, varC_i2 = 0.0, varO_i2 = 0.0;
	float clc_x2 = 0.0, obs_x2 = 0.0, eps_x2 = 0.0;
	for (int z = z_min; z <= z_max; z++)
	{
		for (int y = y_min; y <= y_max; y++)
		{
			for (int x = x_min; x <= x_max; x++)
			{
				clc_x2 = rhoC(x, y, z);
				//                clc_x2 = tmp_rhc;
				obs_x2 = density(x, y, z);
				//                obs_x2 = tmp_den;
				eps_x2 = 1 - inv_rho_mask0(x, y, z); //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
				//        eps_x2 = 1.0;

						// SMOOTHED
				sumCO_i2 += eps_x2 * clc_x2 * obs_x2;
				sumO_i2 += eps_x2 * obs_x2;
				sumO2_i2 += eps_x2 * obs_x2 * obs_x2;
				sumC_i2 += eps_x2 * clc_x2;
				sumC2_i2 += eps_x2 * clc_x2 * clc_x2;
				vol_i2 += eps_x2;

			}

		}
		//                tmp_den = 0.0; 
		//                tmp_rhc=0.0;
		//              if(theDensityMap.density[x]>(2.0/3.0)*max_den) tmp_den=theDensityMap.density[x];
		//                if(rhoC2[x]>(2.0/3.0)*max_rhc) tmp_rhc=rhoC2[x];
		//                if(theDensityMap.density[x]<(3.0/5.0)*max_den) continue;
				// fetch this point
		
	}
	varC_i2 = (sumC2_i2 - sumC_i2 * sumC_i2 / vol_i2);
	varO_i2 = (sumO2_i2 - sumO_i2 * sumO_i2 / vol_i2);
	if (varC_i2 == 0 || varO_i2 == 0 || vol_i2 == 0) {
		//    	cout<<"CCCCC"<<endl;
		CC_i2 = 0;
	}
	else {
		//    	cout<<"sumCO_i2:"<<sumCO_i2<<" sumC_i2:"<<sumC_i2 <<" sumO_i2:"<<sumO_i2<<" vol_i2:"<<vol_i2<<" varC_i2:"<<varC_i2<<" varO_i2:"<<varO_i2<<endl;
		CC_i2 = (sumCO_i2 - sumC_i2 * sumO_i2 / vol_i2) / sqrt(varC_i2 * varO_i2);
	}
	//    cout<<"PPP"<<endl;
		/// 4  CALCULATE PER-ATOM DERIVATIVES

		//std::map< core::Size , numeric::xyzMatrix< core::Real > > symmRots;
	dCCdxs_aacen = vector<vector<vector<float>>>(pnum);
	for (int i = 0; i < pnum; i++) {

		int tmp_xx = lig.num_atoms;
		dCCdxs_aacen[i] = vector<vector<float>>(tmp_xx, vector<float>(3, 0.0));
		//		vector<vector<float>> tmp_op(tmp_xx);
		for (int j = 0; j < lig.num_atoms; j++) {

			vector<float> dVdx_ij(3, 0), dOdx_ij(3, 0), dO2dx_ij(3, 0), dCOdx_ij(3, 0), dC2dx_ij(3, 0);

			vector<int> rho_dx_pt_ij = rho_dx_pt[i][j];
			vector<vector<float>> rho_dx_mask_ij = rho_dx_mask[i][j];
			vector<vector<float>> rho_dx_atm_ij = rho_dx_atm[i][j];

			int npoints = rho_dx_pt_ij.size();
			for (int n = 0; n < npoints; ++n) {
				const int x(rho_dx_pt_ij[n]);
				//				cout<<"VVV: "<<rho_dx_pt_ij<<endl;
				float clc_x = rhoC[x];
				float obs_x = density[x];
				float inv_eps_x = inv_rho_mask0[x];

				vector<float> del_mask(3, 0);
				for (int tt = 0; tt < 3; tt++)
				{
					del_mask[tt] = inv_eps_x * rho_dx_mask_ij[n][tt];
				}
				vector<float> del_rhoc = rho_dx_atm_ij[n];

				for (int tt = 0; tt < 3; tt++)
				{
					dVdx_ij[tt] += del_mask[tt];
					dOdx_ij[tt] += obs_x * del_mask[tt];
					dO2dx_ij[tt] += obs_x * obs_x * del_mask[tt];
					dCOdx_ij[tt] += obs_x * del_rhoc[tt];
					dC2dx_ij[tt] += 2.0 * del_rhoc[tt] * clc_x;
				}


			}

			// finally compute dCC/dx_ij
			float f = (sumCO_i2 - sumC_i2 * sumO_i2 / vol_i2);
			float g = sqrt(varO_i2 * varC_i2);

			vector<float> fprime(3, 0.0);
			vector<float> gprime(3, 0.0);
			for (int kk = 0; kk < 3; kk++)
			{
				fprime[kk] = dCOdx_ij[kk] - 1 / (vol_i2 * vol_i2) * (dOdx_ij[kk] * sumC_i2 * vol_i2 - sumO_i2 * sumC_i2 * dVdx_ij[kk]);
				gprime[kk] = 0.5 * (sqrt(varO_i2) / sqrt(varC_i2) * (dC2dx_ij[kk] + (sumC_i2 * sumC_i2 * dVdx_ij[kk] / (vol_i2 * vol_i2))) +
					sqrt(varC_i2) / sqrt(varO_i2) * (dO2dx_ij[kk] - (1 / (vol_i2 * vol_i2) * (2 * vol_i2 * sumO_i2 * dOdx_ij[kk] - sumO_i2 * sumO_i2 * dVdx_ij[kk]))));
			}
			//			cout << "mmm" << endl;
			vector<float> tmp_vct(3, 0.0);
			for (int kk = 0; kk < 3; kk++)
			{
				tmp_vct[kk] = (-1.0) * (g * fprime[kk] - f * gprime[kk]) / (g * g);//  梯度为-CC的梯度
				//				cout << "ddd " <<fprime[0]<<" "<<fprime[1]<<" "<<fprime[2]<<" "<<" "<<gprime[0]<<" "<<gprime[1]<<" "<<gprime[2]<< endl;
			}
			//			tmp_op[j]=tmp_vct;
			dCCdxs_aacen[i][j] = tmp_vct;
			//			cout << "ddd" <<tmp_vct[0]<<" "<<tmp_vct[1]<<" "<<tmp_vct[2]<<" "<< endl;
		}
		//		dCCdxs_aacen.push_back(tmp_op);
	}
	//	cout<<"PPP2"<<endl;
	return CC_i2;
}

//float ElectronDensity::matchpose(DOCKMol& lig)
//{
//    vector<float> del_ijx(3,0.0);
//    vector<float> atm_jx(3,0.0); 
//    string elt_i;   
////    int pnum = PB.size(); 
//    vector<float> atm_idx;
//    ObjexxFCL::FArray3D< float > rhoC;
//    ObjexxFCL::FArray3D< float > inv_rho_mask0;
//    rhoC.dimension(density.u1() , density.u2() , density.u3());
//    inv_rho_mask0.dimension(density.u1() , density.u2() , density.u3());
////    for(int i=0; i<pnum;i++)
//    {
////        cout<<"YYYY: "<<pointsBx[i].elt_<<endl;;
////        if(pointsBx[i].elt_ == "CA")
////        {
////        int nheavyatm= PB[i].size();
//        int tmp_i=0;
//        for(int j=0;j<1;j++)
//        {
//            atm_idx = vector<float>(3,0.0);
//            vector<float> cartX1;
//            vector<float> fracX1;
//            elt_i = PB.atomname_left;
//            elt_i = elt_i[0];
//            OneGaussianScattering sig_j = get_A( elt_i );
//            float k = sig_j.k( effectiveB );
//            float C = sig_j.C( k );
//            if ( C < 1e-6 ) continue;  
//
//            cartX1[0] = PB.X;     
//            cartX1[1] = PB.Y; 
//            cartX1[2] = PB.Z;      
//            MatrixTimesTransVector(c2f,cartX1,fracX1); 
//            // the location of atom in grid ?
////            atm_idx[0] = pos_mod (double(fracX1[0]*grid[0] - origin[0] + 1) , (double)grid[0]);
////            atm_idx[1] = pos_mod (double(fracX1[1]*grid[1] - origin[1] + 1) , (double)grid[1]);
////            atm_idx[2] = pos_mod (double(fracX1[2]*grid[2] - origin[2] + 1) , (double)grid[2]);   
//            atm_idx[0] = fracX1[0]*grid[0] - origin[0] + 1 ;
//            atm_idx[1] = fracX1[1]*grid[1] - origin[1] + 1 ;
//            atm_idx[2] = fracX1[2]*grid[2] - origin[2] + 1 ;
//    //        cout<<"RRRR: "<<atm_idx[tmp_i][0]<<" "<<atm_idx[tmp_i][1]<<" "<<atm_idx[tmp_i][2]<<endl; 
//            for(int z=1;z<=density.u3();z++)
//            {
//                atm_jx[2] = z;
//                del_ijx[2] =(atm_idx[2]-atm_jx[2])/grid[2];
//                if ( del_ijx[2] > 0.5 ) del_ijx[2]-=1.0;
//                if ( del_ijx[2] < -0.5 ) del_ijx[2]+=1.0;                    
//                del_ijx[0] = del_ijx[1] = 0.0;
//                vector<float> frac_tmpz;
//                MatrixTimesTransVector(f2c,del_ijx,frac_tmpz);
//                if(square_len(frac_tmpz)> (ATOM_MASK_PADDING+CA_MASK)*(ATOM_MASK_PADDING+CA_MASK) ) continue;                    
//                for(int y=1;y<=density.u2();y++)
//                {
//                    atm_jx[1] = y;
//                    del_ijx[1] = (atm_idx[1] - atm_jx[1])/grid[1];
//                    // wrap-around??
//                    if ( del_ijx[1] > 0.5 ) del_ijx[1]-=1.0;
//                    if ( del_ijx[1] < -0.5 ) del_ijx[1]+=1.0;                        
//                    del_ijx[0] = 0.0;
//                    vector<float> frac_tmpy;
//                    MatrixTimesTransVector(f2c,del_ijx,frac_tmpy);
//                    if(square_len(frac_tmpy)> (ATOM_MASK_PADDING+CA_MASK)*(ATOM_MASK_PADDING+CA_MASK) ) continue;                        
//                    for(int x=1;x<=density.u1();x++)
//                    {
//                        atm_jx[0] = x;
//                        del_ijx[0] = (atm_idx[0] - atm_jx[0])/grid[0];
//                        // wrap-around??
//                        if ( del_ijx[0] > 0.5 ) del_ijx[0]-=1.0;
//                        if ( del_ijx[0] < -0.5 ) del_ijx[0]+=1.0;                            
//                        vector<float> cart_del_ij2;
//                        MatrixTimesTransVector(f2c,del_ijx,cart_del_ij2);
//                        float d2 = square_len(cart_del_ij2);
//                        if(d2 > (ATOM_MASK_PADDING+CA_MASK)*(ATOM_MASK_PADDING+CA_MASK) ) continue;
//                    
//                        float atm = C*exp(-k*d2);
////                            cout<<"MASK: "<<theDensityMap.ATOM_MASK<<endl;
//                        float sigmoid_msk = exp( d2 - (ATOM_MASK)*(ATOM_MASK)  );
//                    //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.ATOM_MASK)  );
//                    //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.CA_MASK)  );
//                        float inv_msk = 1/(1+sigmoid_msk);
//                        rhoC(x,y,z) += atm;
//                        inv_rho_mask0(x,y,z) *= (1 - inv_msk);
//
////                            if ( d2 <= (theDensityMap.ATOM_MASK_SQ) ) {
////                                mask2(x,y,z) = 1.0; // problem?
////                                if ( d2 <= (theDensityMap.ATOM_DENS_SQ) ) {
////                                    float atm = C*exp(-k*d2);
////                                    rhoC(x,y,z) += atm;
////                                }
////                            }
//                    }
//                }
//            } 
////            cout<<" "<<tmp_i<<" "<<endl;           
//            tmp_i = tmp_i + 1;                                    
//        }
//    }
//
//    float sumC_i2=0.0, sumO_i2=0.0, sumCO_i2=0.0, vol_i2=0.0, CC_i2=0.0;
//    float sumO2_i2=0.0, sumC2_i2=0.0, varC_i2=0.0, varO_i2=0.0;
//    float clc_x2=0.0, obs_x2=0.0, eps_x2=0.0;
//    for ( int x=0; x<density.u1()*density.u2()*density.u3(); ++x ) {
////                tmp_den = 0.0; 
////                tmp_rhc=0.0;
////              if(theDensityMap.density[x]>(2.0/3.0)*max_den) tmp_den=theDensityMap.density[x];
////                if(rhoC2[x]>(2.0/3.0)*max_rhc) tmp_rhc=rhoC2[x];
////                if(theDensityMap.density[x]<(3.0/5.0)*max_den) continue;
//        // fetch this point
//        clc_x2 = rhoC[x];
////                clc_x2 = tmp_rhc;
//        obs_x2 = density[x];
////                obs_x2 = tmp_den;
//        eps_x2 = 1-inv_rho_mask0[x]; //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
////        eps_x2 = 1.0;
//
//        // SMOOTHED
//        sumCO_i2 += eps_x2*clc_x2*obs_x2;
//        sumO_i2  += eps_x2*obs_x2;
//        sumO2_i2 += eps_x2*obs_x2*obs_x2;
//        sumC_i2  += eps_x2*clc_x2;
//        sumC2_i2 += eps_x2*clc_x2*clc_x2;
//        vol_i2   += eps_x2;
//    }
//    varC_i2 = (sumC2_i2 - sumC_i2*sumC_i2 / vol_i2 );
//    varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2 / vol_i2 ) ;
//    if ( varC_i2 == 0 || varO_i2 == 0 ) {
//        CC_i2 = 0;
//    } else {
//        CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2/ vol_i2) / sqrt( varC_i2 * varO_i2 );
//    }
//
//    return CC_i2;        
//}
float ElectronDensity::matchpose(DOCKMol lig)
{
	vector<float> del_ijx(3, 0.0);
	vector<float> atm_jx(3, 0.0);
	string elt_i,elt_i1;
	int pnum = 1;
	vector<float> atm_idx;
	ObjexxFCL::FArray3D< float > rhoC;
	ObjexxFCL::FArray3D< float > inv_rho_mask0;
	int Denx = density.u1();
	int Deny = density.u2();
	int Denz = density.u3();
	int Denxyz = Denx * Deny * Denz;
	rhoC.dimension(Denx, Deny, Denz, 0.0);
	inv_rho_mask0.dimension(Denx, Deny, Denz, 1.0);

	//    cout<<"ATOM_MASK_PADDING,CA_MASK: "<<ATOM_MASK_PADDING<<" "<<CA_MASK<<endl;
	
	for (int i = 0; i < lig.num_atoms; i++)
	{
		//        cout<<"YYYY: "<<pointsBx[i].elt_<<endl;;
		//        if(pointsBx[i].elt_ == "CA")
		//        {
		elt_i = lig.atom_names[i];
		elt_i1 = elt_i[0];
		if(elt_i == "CA")
		{
			atm_idx = vector<float>(3, 0.0);
			vector<float> cartX1(3, 0.0);
			vector<float> fracX1;
			//elt_i = PB[i].elt_;
			//elt_i = elt_i[0];
			OneGaussianScattering sig_j = get_A(elt_i1);
			float k = sig_j.k(effectiveB);
			float C = sig_j.C(k);
			if (C < 1e-6) continue;

			//cartX1 = PB[i].x_;
			cartX1[0] = lig.x[i];
			cartX1[1] = lig.y[i];			
			cartX1[2] = lig.z[i];
			MatrixTimesTransVector(c2f, cartX1, fracX1);
			// the location of atom in grid ?
			atm_idx[0] = (double(fracX1[0] * grid[0] - origin[0] + 1));
			atm_idx[1] = (double(fracX1[1] * grid[1] - origin[1] + 1));
			atm_idx[2] = (double(fracX1[2] * grid[2] - origin[2] + 1));
			//        atm_idx[tmp_i][0] = fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1 ;
			//        atm_idx[tmp_i][1] = fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1 ;
			//        atm_idx[tmp_i][2] = fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1 ;
			//        cout<<"RRRR: "<<atm_idx[tmp_i][0]<<" "<<atm_idx[tmp_i][1]<<" "<<atm_idx[tmp_i][2]<<endl; 
			for (int z = 1; z <= Denz; z++)
			{
				atm_jx[2] = z;
				del_ijx[2] = (atm_idx[2] - atm_jx[2]) / grid[2];
				del_ijx[0] = del_ijx[1] = 0.0;
				vector<float> frac_tmpz;
				MatrixTimesTransVector(f2c, del_ijx, frac_tmpz);
				if (square_len(frac_tmpz) > (ATOM_MASK_PADDING + CA_MASK) * (ATOM_MASK_PADDING + CA_MASK)) continue;
				for (int y = 1; y <= Deny; y++)
				{
					atm_jx[1] = y;
					del_ijx[1] = (atm_idx[1] - atm_jx[1]) / grid[1];
					del_ijx[0] = 0.0;
					vector<float> frac_tmpy;
					MatrixTimesTransVector(f2c, del_ijx, frac_tmpy);
					if (square_len(frac_tmpy) > (ATOM_MASK_PADDING + CA_MASK) * (ATOM_MASK_PADDING + CA_MASK)) continue;
					for (int x = 1; x <= Denx; x++)
					{
						atm_jx[0] = x;
						del_ijx[0] = (atm_idx[0] - atm_jx[0]) / grid[0];
						vector<float> cart_del_ij2;
						MatrixTimesTransVector(f2c, del_ijx, cart_del_ij2);
						float d2 = square_len(cart_del_ij2);
						if (d2 > (ATOM_MASK_PADDING + CA_MASK) * (ATOM_MASK_PADDING + CA_MASK)) continue;

						float atm = C * exp(-k * d2);
						//                            cout<<"MASK: "<<theDensityMap.ATOM_MASK<<endl;
						float sigmoid_msk = exp(d2 - (ATOM_MASK) * (ATOM_MASK));
						//    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.ATOM_MASK)  );
						//    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.CA_MASK)  );
						float inv_msk = 1 / (1 + sigmoid_msk);
						rhoC(x, y, z) += atm;
						inv_rho_mask0(x, y, z) *= (1 - inv_msk);

						//                            if ( d2 <= (theDensityMap.ATOM_MASK_SQ) ) {
						//                                mask2(x,y,z) = 1.0; // problem?
						//                                if ( d2 <= (theDensityMap.ATOM_DENS_SQ) ) {
						//                                    float atm = C*exp(-k*d2);
						//                                    rhoC(x,y,z) += atm;
						//                                }
						//                            }
					}
				}
			}
			//            cout<<" "<<tmp_i<<" "<<endl;                                       
		}
	}

	//for (int i = 0; i < lig1.num_atoms; i++)
	//{
	//	//        cout<<"YYYY: "<<pointsBx[i].elt_<<endl;;
	//	//        if(pointsBx[i].elt_ == "CA")
	//	//        {
	//	elt_i = lig1.atom_names[i];
	//	elt_i = elt_i[0];
	//	if (lig1.atom_types[i] != "H")
	//	{
	//		atm_idx = vector<float>(3, 0.0);
	//		vector<float> cartX1(3, 0.0);
	//		vector<float> fracX1;
	//		//elt_i = PB[i].elt_;
	//		//elt_i = elt_i[0];
	//		OneGaussianScattering sig_j = get_A(elt_i);
	//		float k = sig_j.k(effectiveB);
	//		float C = sig_j.C(k);
	//		if (C < 1e-6) continue;
	//		//cartX1 = PB[i].x_;
	//		cartX1[0] = lig1.x[i];
	//		cartX1[1] = lig1.y[i];
	//		cartX1[2] = lig1.z[i];
	//		MatrixTimesTransVector(c2f, cartX1, fracX1);
	//		// the location of atom in grid ?
	//		atm_idx[0] = (double(fracX1[0] * grid[0] - origin[0] + 1));
	//		atm_idx[1] = (double(fracX1[1] * grid[1] - origin[1] + 1));
	//		atm_idx[2] = (double(fracX1[2] * grid[2] - origin[2] + 1));
	//		//        atm_idx[tmp_i][0] = fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1 ;
	//		//        atm_idx[tmp_i][1] = fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1 ;
	//		//        atm_idx[tmp_i][2] = fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1 ;
	//		//        cout<<"RRRR: "<<atm_idx[tmp_i][0]<<" "<<atm_idx[tmp_i][1]<<" "<<atm_idx[tmp_i][2]<<endl; 
	//		for (int z = 1; z <= Denz; z++)
	//		{
	//			atm_jx[2] = z;
	//			del_ijx[2] = (atm_idx[2] - atm_jx[2]) / grid[2];
	//			del_ijx[0] = del_ijx[1] = 0.0;
	//			vector<float> frac_tmpz;
	//			MatrixTimesTransVector(f2c, del_ijx, frac_tmpz);
	//			if (square_len(frac_tmpz) > (ATOM_MASK_PADDING + CA_MASK) * (ATOM_MASK_PADDING + CA_MASK)) continue;
	//			for (int y = 1; y <= Deny; y++)
	//			{
	//				atm_jx[1] = y;
	//				del_ijx[1] = (atm_idx[1] - atm_jx[1]) / grid[1];
	//				del_ijx[0] = 0.0;
	//				vector<float> frac_tmpy;
	//				MatrixTimesTransVector(f2c, del_ijx, frac_tmpy);
	//				if (square_len(frac_tmpy) > (ATOM_MASK_PADDING + CA_MASK) * (ATOM_MASK_PADDING + CA_MASK)) continue;
	//				for (int x = 1; x <= Denx; x++)
	//				{
	//					atm_jx[0] = x;
	//					del_ijx[0] = (atm_idx[0] - atm_jx[0]) / grid[0];
	//					vector<float> cart_del_ij2;
	//					MatrixTimesTransVector(f2c, del_ijx, cart_del_ij2);
	//					float d2 = square_len(cart_del_ij2);
	//					if (d2 > (ATOM_MASK_PADDING + CA_MASK) * (ATOM_MASK_PADDING + CA_MASK)) continue;
	//					float atm = C * exp(-k * d2);
	//					//                            cout<<"MASK: "<<theDensityMap.ATOM_MASK<<endl;
	//					float sigmoid_msk = exp(d2 - (ATOM_MASK) * (ATOM_MASK));
	//					//    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.ATOM_MASK)  );
	//					//    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.CA_MASK)  );
	//					float inv_msk = 1 / (1 + sigmoid_msk);
	//					rhoC(x, y, z) += atm;
	//					inv_rho_mask0(x, y, z) *= (1 - inv_msk);
	//					//                            if ( d2 <= (theDensityMap.ATOM_MASK_SQ) ) {
	//					//                                mask2(x,y,z) = 1.0; // problem?
	//					//                                if ( d2 <= (theDensityMap.ATOM_DENS_SQ) ) {
	//					//                                    float atm = C*exp(-k*d2);
	//					//                                    rhoC(x,y,z) += atm;
	//					//                                }
	//					//                            }
	//				}
	//			}
	//		}
	//		//            cout<<" "<<tmp_i<<" "<<endl;                                       
	//	}
	//}

	float sumC_i2 = 0.0, sumO_i2 = 0.0, sumCO_i2 = 0.0, vol_i2 = 0.0, CC_i2 = 0.0;
	float sumO2_i2 = 0.0, sumC2_i2 = 0.0, varC_i2 = 0.0, varO_i2 = 0.0;
	float clc_x2 = 0.0, obs_x2 = 0.0, eps_x2 = 0.0;
	for (int x = 0; x < Denxyz; ++x) {
		//                tmp_den = 0.0; 
		//                tmp_rhc=0.0;
		//              if(theDensityMap.density[x]>(2.0/3.0)*max_den) tmp_den=theDensityMap.density[x];
		//                if(rhoC2[x]>(2.0/3.0)*max_rhc) tmp_rhc=rhoC2[x];
		//                if(theDensityMap.density[x]<(3.0/5.0)*max_den) continue;
				// fetch this point
		clc_x2 = rhoC[x];
		//                clc_x2 = tmp_rhc;
		obs_x2 = density[x];
		//                obs_x2 = tmp_den;
		eps_x2 = 1 - inv_rho_mask0[x]; //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
		//        eps_x2 = 1.0;

				// SMOOTHED
		sumCO_i2 += eps_x2 * clc_x2 * obs_x2;
		sumO_i2 += eps_x2 * obs_x2;
		sumO2_i2 += eps_x2 * obs_x2 * obs_x2;
		sumC_i2 += eps_x2 * clc_x2;
		sumC2_i2 += eps_x2 * clc_x2 * clc_x2;
		vol_i2 += eps_x2;
	}
	varC_i2 = (sumC2_i2 - sumC_i2 * sumC_i2 / vol_i2);
	varO_i2 = (sumO2_i2 - sumO_i2 * sumO_i2 / vol_i2);
	if (varC_i2 == 0 || varO_i2 == 0 || vol_i2 == 0) {
		CC_i2 = 0;
	}
	else {
		CC_i2 = (sumCO_i2 - sumC_i2 * sumO_i2 / vol_i2) / sqrt(varC_i2 * varO_i2);
	}
	return CC_i2;
}

void ElectronDensity::diffmap(DOCKMol lig)
{
	vector<float> del_ijx(3, 0.0);
	vector<float> atm_jx(3, 0.0);
	string elt_i;
	int pnum = 1;
	vector<float> atm_idx;
	ObjexxFCL::FArray3D< float > rhoC;
	ObjexxFCL::FArray3D< float > inv_rho_mask0;
	int Denx = density.u1();
	int Deny = density.u2();
	int Denz = density.u3();
	int Denxyz = Denx * Deny * Denz;
	rhoC.dimension(Denx, Deny, Denz, 0.0);
	inv_rho_mask0.dimension(Denx, Deny, Denz, 1.0);

	//    cout<<"ATOM_MASK_PADDING,CA_MASK: "<<ATOM_MASK_PADDING<<" "<<CA_MASK<<endl;

	for (int i = 0; i < lig.num_atoms; i++)
	{
		//        cout<<"YYYY: "<<pointsBx[i].elt_<<endl;;
		//        if(pointsBx[i].elt_ == "CA")
		//        {
		elt_i = lig.atom_names[i];
		elt_i = elt_i[0];
		{
			atm_idx = vector<float>(3, 0.0);
			vector<float> cartX1(3, 0.0);
			vector<float> fracX1;
			//elt_i = PB[i].elt_;
			//elt_i = elt_i[0];
			OneGaussianScattering sig_j = get_A(elt_i);
			float k = sig_j.k(effectiveB);
			float C = sig_j.C(k);
			if (C < 1e-6) continue;

			//cartX1 = PB[i].x_;
			cartX1[0] = lig.x[i];
			cartX1[1] = lig.y[i];
			cartX1[2] = lig.z[i];
			MatrixTimesTransVector(c2f, cartX1, fracX1);
			// the location of atom in grid ?
			atm_idx[0] = (double(fracX1[0] * grid[0] - origin[0] + 1));
			atm_idx[1] = (double(fracX1[1] * grid[1] - origin[1] + 1));
			atm_idx[2] = (double(fracX1[2] * grid[2] - origin[2] + 1));
			//        atm_idx[tmp_i][0] = fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1 ;
			//        atm_idx[tmp_i][1] = fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1 ;
			//        atm_idx[tmp_i][2] = fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1 ;
			//        cout<<"RRRR: "<<atm_idx[tmp_i][0]<<" "<<atm_idx[tmp_i][1]<<" "<<atm_idx[tmp_i][2]<<endl; 
			for (int z = 1; z <= Denz; z++)
			{
				atm_jx[2] = z;
				del_ijx[2] = (atm_idx[2] - atm_jx[2]) / grid[2];
				del_ijx[0] = del_ijx[1] = 0.0;
				vector<float> frac_tmpz;
				MatrixTimesTransVector(f2c, del_ijx, frac_tmpz);
				if (square_len(frac_tmpz) > (ATOM_MASK_PADDING + CA_MASK) * (ATOM_MASK_PADDING + CA_MASK)) continue;
				for (int y = 1; y <= Deny; y++)
				{
					atm_jx[1] = y;
					del_ijx[1] = (atm_idx[1] - atm_jx[1]) / grid[1];
					del_ijx[0] = 0.0;
					vector<float> frac_tmpy;
					MatrixTimesTransVector(f2c, del_ijx, frac_tmpy);
					if (square_len(frac_tmpy) > (ATOM_MASK_PADDING + CA_MASK) * (ATOM_MASK_PADDING + CA_MASK)) continue;
					for (int x = 1; x <= Denx; x++)
					{
						atm_jx[0] = x;
						del_ijx[0] = (atm_idx[0] - atm_jx[0]) / grid[0];
						vector<float> cart_del_ij2;
						MatrixTimesTransVector(f2c, del_ijx, cart_del_ij2);
						float d2 = square_len(cart_del_ij2);
						if (d2 > (ATOM_MASK_PADDING + CA_MASK) * (ATOM_MASK_PADDING + CA_MASK)) continue;

						float atm = C * exp(-k * d2);
						//                            cout<<"MASK: "<<theDensityMap.ATOM_MASK<<endl;
						float sigmoid_msk = exp(d2 - (ATOM_MASK) * (ATOM_MASK));
						//    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.ATOM_MASK)  );
						//    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.CA_MASK)  );
						float inv_msk = 1 / (1 + sigmoid_msk);
						rhoC(x, y, z) += atm;
						inv_rho_mask0(x, y, z) *= (1 - inv_msk);

						//                            if ( d2 <= (theDensityMap.ATOM_MASK_SQ) ) {
						//                                mask2(x,y,z) = 1.0; // problem?
						//                                if ( d2 <= (theDensityMap.ATOM_DENS_SQ) ) {
						//                                    float atm = C*exp(-k*d2);
						//                                    rhoC(x,y,z) += atm;
						//                                }
						//                            }
					}
				}
			}
			//            cout<<" "<<tmp_i<<" "<<endl;                                       
		}
	}

	float sumC_i2 = 0.0, sumO_i2 = 0.0, sumCO_i2 = 0.0, vol_i2 = 0.0, CC_i2 = 0.0;
	float sumO2_i2 = 0.0, sumC2_i2 = 0.0, varC_i2 = 0.0, varO_i2 = 0.0;
	float clc_x2 = 0.0, obs_x2 = 0.0, eps_x2 = 0.0;
	for (int x = 0; x < Denxyz; ++x) {
		//                tmp_den = 0.0; 
		//                tmp_rhc=0.0;
		//              if(theDensityMap.density[x]>(2.0/3.0)*max_den) tmp_den=theDensityMap.density[x];
		//                if(rhoC2[x]>(2.0/3.0)*max_rhc) tmp_rhc=rhoC2[x];
		//                if(theDensityMap.density[x]<(3.0/5.0)*max_den) continue;
				// fetch this point
		clc_x2 = rhoC[x];
		//                clc_x2 = tmp_rhc;
		obs_x2 = density[x];
		//                obs_x2 = tmp_den;
		cout << "clc:  " << clc_x2 << "  obs:  " << obs_x2 << endl;

		if (rhoC[x] != 0)
		{
			density[x] = 0;
		}
		eps_x2 = 1 - inv_rho_mask0[x]; //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
		//        eps_x2 = 1.0;

				// SMOOTHED
		sumCO_i2 += eps_x2 * clc_x2 * obs_x2;
		sumO_i2 += eps_x2 * obs_x2;
		sumO2_i2 += eps_x2 * obs_x2 * obs_x2;
		sumC_i2 += eps_x2 * clc_x2;
		sumC2_i2 += eps_x2 * clc_x2 * clc_x2;
		vol_i2 += eps_x2;
	}
	varC_i2 = (sumC2_i2 - sumC_i2 * sumC_i2 / vol_i2);
	varO_i2 = (sumO2_i2 - sumO_i2 * sumO_i2 / vol_i2);
	if (varC_i2 == 0 || varO_i2 == 0 || vol_i2 == 0) {
		CC_i2 = 0;
	}
	else {
		CC_i2 = (sumCO_i2 - sumC_i2 * sumO_i2 / vol_i2) / sqrt(varC_i2 * varO_i2);
	}
	cout << "CC:  " << CC_i2 << endl;
	return ;
}


//void ElectronDensity::calcRhoCx(
//	vector<vector<molecule> > const &pose,
//	float highres_limit,
//	ObjexxFCL::FArray3D< float > &rhoC,
//	ObjexxFCL::FArray3D< float > &mask,
//	float fixed_mask_B /* = -1 */,
//	float B_upper_limit /* = 600 */,
//	float force_mask /*=-1*/ )
//	{
//	// get rho_c
//	rhoC.dimension(density.u1() , density.u2() , density.u3());
//	mask.dimension(density.u1() , density.u2() , density.u3());
//	rhoC = 0.0;
//	mask = 0.0;
//	bool use_Bs = pose_has_nonzero_Bs( pose );
//	if ( !use_Bs ) {
//		cout << "Input pose has no nonzero B factors ... setting to baseline" << std::endl;
//	}
//
//	vector<float> per_atom_ks;
//	vector<float> per_atom_Cs;
//	vector<float> ATOM_MASK_SQS;
//	vector<vector<float> > pose_grid;
//	int tmp_p=-1;
//	for(int i=0;i<pose.size();i++)
//	{
//		vector<molecule> chain_x=pose[i];
//		for(int j=0;j<chain_x.size();j++)
//		{
//			molecule residue_x = chain_x[j];
//			for(int t=0; t<residue_x.atoms.size();t++)
//			{
//				tmp_p = tmp_p+1;
//				matom atom_x = residue_x.atoms[t];
//				std::string elt_i = to_string(atom_x.atomname_left[0]); //name2eltx()
////				cout<<"elt_i: "<<elt_i<<endl;
//				OneGaussianScattering sig_j = get_A( elt_i );
//				per_atom_ks.push_back(sig_j.k( atom_x.temfactor, B_upper_limit ));
//				if ( use_Bs ) {
//					per_atom_ks[tmp_p] = std::min ( per_atom_ks[tmp_p], (float) (4*M_PI*M_PI/minimumB ));
//				} else {
//					per_atom_ks[tmp_p] = std::min ( per_atom_ks[tmp_p], (float) (4*M_PI*M_PI/effectiveB ));
//				}
//				per_atom_Cs.push_back(sig_j.C( per_atom_ks[tmp_p] ));
//
//				float C = per_atom_Cs[tmp_p], k = per_atom_ks[tmp_p];
//				if ( C < 1e-6 ) 
//				{
//					ATOM_MASK_SQS.push_back(0.0);
//					vector<float> satom(3,0.0);
//					satom[0]=atom_x.X;
//					satom[1]=atom_x.Y;
//					satom[2]=atom_x.Z;
//					pose_grid.push_back(satom);
//					continue;
//				}
//				ATOM_MASK_SQS.push_back((1.0/k) * ( std::log( C ) - std::log(1e-4) ));  // 1e-4 is the density value where the mask goes to 0
//
//				// for b factor minimization, we don't want the mask to change.
//				// so we fix the mask at a high B value
////				cout<<"fixed_mask_B: "<<fixed_mask_B<<endl;
//				if ( fixed_mask_B>0.0000001 ) {
//					float mK = sig_j.k( fixed_mask_B, B_upper_limit );
//					float mC = sig_j.C( mK );
//					ATOM_MASK_SQS[tmp_p] = (1.0/mK) * (std::log( mC ) - std::log(1e-4));
//				}
//
//				// force mask
////				cout<<"force_mask: "<<force_mask<<endl;
//				if ( force_mask>0.0000001 ) {
//					ATOM_MASK_SQS[tmp_p] = force_mask*force_mask;
//				}
//
//				// TO DO: OPTIONALLY control periodic/nonperiodic boundaries
//				//atm_idx[0] = pos_mod (fracX[0]*grid[0] - origin[0] + 1 , (double)grid[0]);
//				//atm_idx[1] = pos_mod (fracX[1]*grid[1] - origin[1] + 1 , (double)grid[1]);
//				//atm_idx[2] = pos_mod (fracX[2]*grid[2] - origin[2] + 1 , (double)grid[2]);
//				vector<float> cartX(3,0.0);
//				 cartX[0]= atom_x.X; // - getTransform();
//				 cartX[1]= atom_x.Y; // - getTransform();
//				 cartX[2]= atom_x.Z; // - getTransform();
//			//	vector<float> fracX = c2f*cartX;
//				vector<float> fracX;
//				MatrixTimesTransVector(c2f,cartX,fracX);
//				vector<float> atm_idx(3,fracX[0]*grid[0] - origin[0] + 1);
//				atm_idx[1] = fracX[1]*grid[1] - origin[1] + 1;
//				atm_idx[2] = fracX[2]*grid[2] - origin[2] + 1;
//				pose_grid.push_back( atm_idx );				
////				tmp_p = tmp_p+1;
//			}
//		}
//	}
//
//	// PASS 1: rho_calc
//	int tmp_pp=-1;
//	for(int i=0;i<pose.size();i++)
//	{
//		vector<molecule> chain_x=pose[i];
//		for(int j=0;j<chain_x.size();j++)
//		{
//			molecule residue_x = chain_x[j];
//			for(int t=0; t<residue_x.atoms.size();t++)
//			{
//				tmp_pp = tmp_pp +1;
//				matom atom_x = residue_x.atoms[t];
//				std::cout << "\rrho_calc " << tmp_pp << " of " << std::flush;	
//				float C = per_atom_Cs[tmp_pp], k = per_atom_ks[tmp_pp];
//				if ( C < 1e-6 ) continue;	
//				
//				// dist cutoff for mask
//				float ATOM_MASK_SQ = ATOM_MASK_SQS[tmp_pp];
//
//				// dist cutoff for density
//				float ATOM_DENS_SQ = (1.0/k) * (std::log( C ) - std::log(1e-4));
//
//				vector<float> atm_j(3,0.0), del_ij(3,0.0), atm_idx;
//				atm_idx =  pose_grid[tmp_pp];
//
//				if ( atm_idx[0]<1.0 || atm_idx[0]>grid[0] ) continue;
//				if ( atm_idx[1]<1.0 || atm_idx[1]>grid[1] ) continue;
//				if ( atm_idx[2]<1.0 || atm_idx[2]>grid[2] ) continue;
//
//				for(int z=1;z<=density.u3();z++)
//				{
//					atm_j[2] = z;
//					del_ij[2] =(atm_idx[2]-atm_j[2])/grid[2];
//					del_ij[0] = del_ij[1] = 0.0;
//					vector<float> frac_tmpz;
//					MatrixTimesTransVector(f2c,del_ij,frac_tmpz);
//					if(square_len(frac_tmpz)> (ATOM_MASK_SQ) ) continue;					
//					for(int y=1;y<=density.u2();y++)
//					{
//						atm_j[1] = y;
//						del_ij[1] = (atm_idx[1] - atm_j[1])/grid[1];
//						del_ij[0] = 0.0;
//						vector<float> frac_tmpy;
//						MatrixTimesTransVector(f2c,del_ij,frac_tmpy);
//						if(square_len(frac_tmpy)> (ATOM_MASK_SQ) ) continue;						
//						for(int x=1;x<=density.u1();x++)
//						{
//							atm_j[0] = x;
//							del_ij[0] = (atm_idx[0] - atm_j[0])/grid[0];
//							vector<float> cart_del_ij;
//							MatrixTimesTransVector(f2c,del_ij,cart_del_ij);
//							float d2 = square_len(cart_del_ij);
//						
//							if ( d2 <= (ATOM_MASK_SQ) ) {
//								mask(x,y,z) = 1.0; // problem?
//								if ( d2 <= (ATOM_DENS_SQ) ) {
//									float atm = C*exp(-k*d2);
//									rhoC(x,y,z) += atm;
//								}
//							}
//						}
//
//					}
//				}				
////				tmp_pp = tmp_pp +1;
//			}
//		}
//	}
//	cout<<endl;
//
////	if ( highres_limit == 0 ) return;
//	if ( highres_limit <= 0.0000001 ) return;
//
//	// bandlimit mask at 'radius'
//	ObjexxFCL::FArray3D< std::complex<double> > Fmask;
//	fourier::fft3(mask, Fmask);
//	//int H,K,L;
//	for ( int z=1; z<=(int)grid[2]; ++z ) {
//		int H = (z < (int)grid[2]/2) ? z-1 : z-grid[2] - 1;
//		for ( int y=1; y<=(int)grid[1]; ++y ) {
//			int K = (y < (int)grid[1]/2) ? y-1 : y-grid[1]-1;
//			for ( int x=1; x<=(int)grid[0]; ++x ) {
//				int L = (x < (int)grid[0]/2) ? x-1 : x-grid[0]-1;
//				float S2c =  S2(H,K,L);
//
//				// exp fade
//				float scale = exp(-S2c*(highres_limit*highres_limit));
//				Fmask(x,y,z) *= scale;
//			}
//		}
//	}
//	fourier::ifft3(Fmask, mask);
// }

//void ElectronDensity::calcRhoCy(
//	vector<vector<molecule> > const &pose,
//	float highres_limit,
//	ObjexxFCL::FArray3D< float > &rhoC,
//	ObjexxFCL::FArray3D< float > &mask,
//	float fixed_mask_B /* = -1 */,
//	float B_upper_limit /* = 600 */,
//	float force_mask /*=-1*/ )
//	{
//	// get rho_c
//	rhoC.dimension(density.u1() , density.u2() , density.u3());
//	mask.dimension(density.u1() , density.u2() , density.u3());
//	rhoC = 0.0;
//	mask = 0.0;
//	bool use_Bs = pose_has_nonzero_Bs( pose );
//	if ( !use_Bs ) {
//		cout << "Input pose has no nonzero B factors ... setting to baseline" << std::endl;
//	}
//
//	// PASS 1: rho_calc
//	for(int i=0;i<pose.size();i++)
//	{
//		vector<molecule> chain_x=pose[i];
//		for(int j=0;j<chain_x.size();j++)
//		{
//			molecule residue_x = chain_x[j];
//			for(int t=0; t<residue_x.atoms.size();t++)
//			{
//				matom atom_x = residue_x.atoms[t];
//				vector<float> cartX(3,0.0);
//				cartX[0]= atom_x.X;
//				cartX[1]= atom_x.Y;
//				cartX[2]= atom_x.Z;
//				vector<float> fracX ;
//				MatrixTimesTransVector(c2f,cartX,fracX);
//				vector<float> atm_idx(3,fracX[0]*grid[0] - origin[0] + 1);
//				atm_idx[1] = fracX[1]*grid[1] - origin[1] + 1;
//				atm_idx[2] = fracX[2]*grid[2] - origin[2] + 1;
//				vector<float> atm_j(3,0.0), del_ij(3,0.0);
//
//				if ( atm_idx[0]<1.0 || atm_idx[0]>grid[0] ) continue;
//				if ( atm_idx[1]<1.0 || atm_idx[1]>grid[1] ) continue;
//				if ( atm_idx[2]<1.0 || atm_idx[2]>grid[2] ) continue;
//
//				std::string elt_i = to_string(atom_x.atomname_left[0]); // name2eltx
//				OneGaussianScattering sig_j = get_A( elt_i );
//				float per_atom_ks = sig_j.k( atom_x.temfactor, B_upper_limit );
//				if ( use_Bs ) {
//					per_atom_ks = std::min ( per_atom_ks, (float) (4*M_PI*M_PI/minimumB ));
//				} else {
//					per_atom_ks = std::min ( per_atom_ks, (float) (4*M_PI*M_PI/effectiveB ));
//				}				
//				float C= sig_j.C( per_atom_ks);
//				float k = per_atom_ks;
//
//				float ATOM_MASK_SQ=0.0;
//				ATOM_MASK_SQ = (1.0/k) * ( std::log( C ) - std::log(1e-4) );
//				// dist cutoff for density
//				float ATOM_DENS_SQ = (1.0/k) * (std::log( C ) - std::log(1e-4));
//				// for b factor minimization, we don't want the mask to change.
//				// so we fix the mask at a high B value				
//				if ( fixed_mask_B>0 ) {
//					float mK = sig_j.k( fixed_mask_B, B_upper_limit );
//					float mC = sig_j.C( mK );
//					ATOM_MASK_SQ = (1.0/mK) * (std::log( mC ) - std::log(1e-4));
//				}
//				// force mask
//				if ( force_mask>0 ) {
//					ATOM_MASK_SQ = force_mask*force_mask;
//				}
//
//				for(int z=1;z<=density.u3();z++)
//				{
//					atm_j[2] = z;
//					del_ij[2] =(atm_idx[2]-atm_j[2])/grid[2];
//					del_ij[0] = del_ij[1] = 0.0;
//					vector<float> frac_tmpz;
//					MatrixTimesTransVector(f2c,del_ij,frac_tmpz);
//					if(square_len(frac_tmpz)> (ATOM_MASK_SQ) ) continue;					
//					for(int y=1;y<=density.u2();y++)
//					{
//						atm_j[1] = y;
//						del_ij[1] = (atm_idx[1] - atm_j[1])/grid[1];
//						del_ij[0] = 0.0;
//						vector<float> frac_tmpy;
//						MatrixTimesTransVector(f2c,del_ij,frac_tmpy);
//						if(square_len(frac_tmpy)> (ATOM_MASK_SQ) ) continue;						
//						for(int x=1;x<=density.u1();x++)
//						{
//							atm_j[0] = x;
//							del_ij[0] = (atm_idx[0] - atm_j[0])/grid[0];
//							vector<float> cart_del_ij;
//							MatrixTimesTransVector(f2c,del_ij,cart_del_ij);
//							float d2 = square_len(cart_del_ij);
//						
//							if ( d2 <= (ATOM_MASK_SQ) ) {
//								mask(x,y,z) = 1.0; // problem?
//								if ( d2 <= (ATOM_DENS_SQ) ) {
//									float atm = C*exp(-k*d2);
//									rhoC(x,y,z) += atm;
//								}
//							}
//						}
//
//					}
//				}				
//
//			}
//		}
//	}
//	cout<<endl;
//
//	if ( highres_limit == 0 ) return;
//
//	// bandlimit mask at 'radius'
//	ObjexxFCL::FArray3D< std::complex<double> > Fmask;
//	fourier::fft3(mask, Fmask);
//	//int H,K,L;
//	for ( int z=1; z<=(int)grid[2]; ++z ) {
//		int H = (z < (int)grid[2]/2) ? z-1 : z-grid[2] - 1;
//		for ( int y=1; y<=(int)grid[1]; ++y ) {
//			int K = (y < (int)grid[1]/2) ? y-1 : y-grid[1]-1;
//			for ( int x=1; x<=(int)grid[0]; ++x ) {
//				int L = (x < (int)grid[0]/2) ? x-1 : x-grid[0]-1;
//				float S2c =  S2(H,K,L);
//
//				// exp fade
//				float scale = exp(-S2c*(highres_limit*highres_limit));
//				Fmask(x,y,z) *= scale;
//			}
//		}
//	}
//	fourier::ifft3(Fmask, mask);
// }


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

	init();
//	cout<<"ATOM_MASK: "<<ATOM_MASK<<endl;
//	cout<<"111211"<<endl;
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

	float max_den = -10000.0;
	float min_den = 10000.0;
	for(int x=1;x<=density.u1();x++)
	{
		for(int y=1;y<=density.u2();y++)
		{
			for(int z=1;z<=density.u3();z++)
			{
				float den_curr = density(x,y,z);
				if(den_curr > max_den) max_den = den_curr;
				if(den_curr < min_den) min_den = den_curr;
			}
		}
	}

	for(int x=1;x<=density.u1();x++)
	{
		for(int y=1;y<=density.u2();y++)
		{
			for(int z=1;z<=density.u3();z++)
			{
				density(x,y,z) = float(density(x,y,z)-min_den)/float(max_den-min_den);
				//cout << density(x, y, z) << endl;
			}
		}
	}

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
/*
	if( gridSpacing > 0.1)
	{	
		float b_del_grid[3];
		b_del_grid[0] = cellDimensions[0]/float(grid[0]);
		b_del_grid[1] = cellDimensions[1]/float(grid[1]);
		b_del_grid[2] = cellDimensions[2]/float(grid[2]);
		float N_del_grid[3];
		N_del_grid[0] = b_del_grid[0] * gridSpacing;
		N_del_grid[1] = b_del_grid[1] * gridSpacing;
		N_del_grid[2] = b_del_grid[2] * gridSpacing;
		// 'alloc' changes ordering of "extent"
		ObjexxFCL::FArray3D< float > density_5;
		vol_xsize = ceil (float (extent[0])/gridSpacing);
		vol_ysize = ceil (float (extent[1])/gridSpacing);
		vol_zsize = ceil (float (extent[2])/gridSpacing);
		cout<<"YYY: "<<vol_xsize<<" "<<vol_ysize<<" "<<vol_zsize <<endl;
		density_5.dimension( vol_xsize,vol_ysize,vol_zsize,0.0);
		for(int t=0;t<vol_xsize*vol_ysize*vol_zsize;t++)
		{
			density_5[t]=0.0;
		}		
//		density_5.dimension( vol_xsize,vol_ysize,vol_zsize,0.0);
//		cout<<"YYY: "<<vol_xsize<<" "<<vol_ysize<<" "<<vol_zsize <<endl;
//		int i=1,j,k;
//		for( i=1;i<=vol_xsize;i++)
//		{
//			for( j=1;j<=vol_ysize;j++)
//			{
//				for( k=1;k<=vol_zsize;k++)
//				{
//					cout<<"i,j,k: "<<density_5(i,j,k)<<endl;;//i<<" "<<j<<" "<<k<<endl;				
//				}
//			}
//			
//		} 
		
//		cout<<"gridSpacing: "<<gridSpacing<<endl;
//		cout<<"YYY: "<<vol_xsize<<" "<<vol_ysize<<" "<<vol_zsize <<endl;
//		cout<<"densityxyz: "<<density.u1()<<" "<<density.u2()<<" "<<density.u3()<<endl;
//		cout<<endl;
		int sc_g= int(gridSpacing);
		int ixx=1,iyy=1,izz=1;
		int iz = 1,iy=1,ix=1;
		for(ixx=1,ix=1;ixx<=extent[0];)
		{
			for(iyy=1,iy=1;iyy<=extent[1];)
			{
				for(izz=1,iz=1;izz<=extent[2];)
				{
//					cout<<"ix,iy,iz, ixx,iyy,izz: "<<ix<<" "<<iy<<" "<<iz<<" "<<ixx<<" "<<iyy<<" "<<izz<<endl;
					density_5(ix,iy,iz) = density(ixx,iyy,izz);
					izz = izz + sc_g;
					iz= iz + 1;
				}
				iy = iy +1;
				iyy = iyy + sc_g;
			}
			ix = ix + 1;
			ixx = ixx + sc_g;
		}

//		density.clear();
//		density.dimension( vol_xsize,vol_ysize,vol_zsize);
//		density = density_5;
		bool fftshift =false;
		density.dimension( density_5.u1(),density_5.u2(),density_5.u3() );
	//	if ( fftshift ) origin -= grid/2;
		for ( int i=1; i<=(int)density_5.u1(); ++i ) {
			int fi = (int)(fftshift ? pos_mod( i-(density_5.u1()/2)-1 , density_5.u1())+1 : i);
			for ( int j=1; j<=(int)density_5.u2(); ++j ) {
				int fj = (int)(fftshift ? pos_mod( j-(density_5.u2()/2)-1 , density_5.u2())+1 : j);
				for ( int k=1; k<=(int)density_5.u3(); ++k ) {
					int fk = (int)(fftshift ? pos_mod( k-(density_5.u3()/2)-1 , density_5.u3())+1 : k);
					density(fi,fj,fk) = (float)density_5(i,j,k);
				}
			}
		}		

		grid[0] = density_5.u1();
		grid[1] = density_5.u2();
		grid[2] = density_5.u3();
//		extent[0] =  density_5.u1();
//		extent[1] =  density_5.u2();
//		extent[2] =  density_5.u3();		
		cellDimensions[0] = grid[0]*N_del_grid[0];
		cellDimensions[1] = grid[1]*N_del_grid[1];
		cellDimensions[2] = grid[2]*N_del_grid[2]; 
		this->grid[0] = density_5.u1();
		this->grid[1] = density_5.u2();
		this->grid[2] = density_5.u3();
		this->cellDimensions[0] = grid[0]*N_del_grid[0];
		this->cellDimensions[1] = grid[1]*N_del_grid[1];
		this->cellDimensions[2] = grid[2]*N_del_grid[2]; 
//		this->origin[0] =(float) origin[xyz2crs[0]]*gridSpacing;
//		this->origin[1] =(float) origin[xyz2crs[1]]*gridSpacing;
//		this->origin[2] =(float) origin[xyz2crs[2]]*gridSpacing;		

		cout<<"xxxx"<<endl;
	}  */
//	this->grid.push_back(grid[0]); // error
//	this->grid.push_back(grid[1]);
//	this->grid.push_back(grid[2]);

	///////////////////////////////////
	/// POST PROCESSING
	// expand to unit cell
/*	cellDimensions[0] = cellDimensions[0] /this->grid[0];
	cellDimensions[1] = cellDimensions[1] /this->grid[1];
	cellDimensions[2] = cellDimensions[2] /this->grid[2]; */

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
	cout<<"force: "<<force_apix_on_map_load_<<endl;
	if ( force_apix_on_map_load_ > (0.0000001) ) {
		set_voxel_spacing(force_apix_on_map_load_);
	}	
//	if ( gridSpacing > 0.0000001 ) {
//		set_voxel_spacingx( gridSpacing );
//	}	
//	if ( force_apix_on_map_load_ > 0 ) {
//		set_voxel_spacing( force_apix_on_map_load_ );
//	}

	// resample the map
	if ( gridSpacing > 0.0000001 ) {
		//this->resize( gridSpacing );
	} 
	cout<<"gridSpacing: "<<gridSpacing<<endl;	

	cout << density.u1() << endl;
	//max_den = -10000.0;
	//min_den = 10000.0;
	//for (int x = 1; x <= density.u1(); x++)
	//{
	//	for (int y = 1; y <= density.u2(); y++)
	//	{
	//		for (int z = 1; z <= density.u3(); z++)
	//		{
	//			float den_curr = density(x, y, z);
	//			if (den_curr > max_den) max_den = den_curr;
	//			if (den_curr < min_den) min_den = den_curr;
	//		}
	//	}
	//}

	//for (int x = 1; x <= density.u1(); x++)
	//{
	//	for (int y = 1; y <= density.u2(); y++)
	//	{
	//		for (int z = 1; z <= density.u3(); z++)
	//		{
	//			density(x, y, z) = float(density(x, y, z) - min_den) / float(max_den - min_den);
	//			//cout << density(x, y, z) << endl;
	//		}
	//	}
	//}

//	if ( gridSpacing > 0 ) {
//		this->resize( gridSpacing );
//	}

	// Resolution logic is a bit convoluted
	//   We define a "minimum" resolution and an "effective" resolution
	//   The minimum resolution is the absolute minimum resolution component we can represent (grid_sampling/2)
	//   Effective resolution is the "actual" resolution of the data, and comes from -map_reso (if provided)
	//      or is guessed at grid_sampling

	float max_del_grid = std::max( cellDimensions[0]/((double)this->grid[0]) , cellDimensions[1]/((double)this->grid[1]) );
	max_del_grid = std::max( max_del_grid ,(float) (cellDimensions[2]/((double)this->grid[2])) );
	this->max_del_grid = max_del_grid;
//	cout<<"max "<<max_del_grid<<endl;
//	float max_del_grid = std::max( cellDimensions[0] , cellDimensions[1] );
//	max_del_grid = std::max( max_del_grid ,(float) (cellDimensions[2]) );


	minimumB = 16*max_del_grid*max_del_grid;

	cout<<"reso: "<<reso<<endl;
//	if ( reso == 0 ) {
	if ( reso < 0.0000001 ) {
		max_del_grid *= 1.5;
	} else {
		max_del_grid = std::max( max_del_grid, reso/2 );
	}

	cout << "Effective resolution = " << 2*max_del_grid << std::endl;
	effectiveB = 16*max_del_grid*max_del_grid;
//	cout<<"max "<<max_del_grid<<"  "<< effectiveB<<endl;
	// mask width logic is also convoluted
	// make sure mask extends >= 3 carbon stdevs at effectiveB
//	float mask_min = sqrt( effectiveB / (2*M_PI*M_PI) );
	float mask_min = 3.0 * sqrt( effectiveB / (2*M_PI*M_PI) );
//	if ( basic::options::option[ basic::options::OptionKeys::edensity::atom_mask_min ].user() ) {
//		mask_min = basic::options::option[ basic::options::OptionKeys::edensity::atom_mask_min ];
//	}

//	cout<<"ATOM_MASK: "<<ATOM_MASK<<endl;
//	cout<<"mask_min: "<<mask_min<<endl;
	if ( ATOM_MASK < mask_min ) {
		ATOM_MASK = mask_min;
	}
	cout<<"ATOM_MASK: "<<ATOM_MASK<<endl;
//	mask_min = sqrt(2.0) * std::max( (float)(2.4+1.6*reso) , reso ) / M_PI;
	mask_min = 3.0 * sqrt(2.0) * std::max( (float)(2.4+1.6*reso) , reso ) / M_PI;
//	cout<<"CA_MASK: "<<CA_MASK<<endl;
	cout<<"mask_min: "<<mask_min<<endl;
	if ( CA_MASK < mask_min ) {
		CA_MASK = mask_min;
	}
	cout<<"CA_MASK: "<<CA_MASK<<endl;
	// set map resolution
	this->reso = 2*max_del_grid;

	// update/clear derived data
	density_change_trigger();
	test_density=density;


	//ElectronDensity(density,1.0, vector<float>(3,0.0), false ).writeMRC( "diffmap.mrc");
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
//	ObjexxFCL::FArray3D< double > xdensity(density.u1(),density.u2(),density.u3());
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

	resample( density, newDensity, newDims );
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
	//cout << "AAAA" << cellDimensions[0] << "  " << endl;
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

bool ElectronDensity::writeMRC_test(std::string mapfilename) {
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
	ori_int[0] = use_altorigin?  (int)std::floor( origin[0] ):0;
	ori_int[1] = use_altorigin?  (int)std::floor( origin[1] ):0;
	ori_int[2] = use_altorigin?  (int)std::floor( origin[2] ):0;
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
//		idx2cart (vector<float>(3,1), origin_realspace );
		ori_float[0] = (float)( origin[0] );
		ori_float[1] = (float)( origin[1] );
		ori_float[2] = (float)( origin[2] );
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
				buff_f = (float) test_density(coord[0],coord[1],coord[2]);
				outx.write(reinterpret_cast <char*>(&buff_f), sizeof(float));
			}
		}
	}

	return true;
}

void
ElectronDensity::set_voxel_spacingx( float apix ) {
	// 1) scale origin
	//origin[0] *= (apix[0] * this->grid[0]) / cellDimensions[0];
	//origin[1] *= (apix[1] * this->grid[1]) / cellDimensions[1];
	//origin[2] *= (apix[2] * this->grid[2]) / cellDimensions[2];

	cellDimensions[0] = apix * this->grid[0];
	cellDimensions[1] = apix * this->grid[1];
	cellDimensions[2] = apix * this->grid[2];

	// 2) scale cell dimensions
//	cellDimensions[0] = apix * this->grid[0];
//	cellDimensions[1] = apix * this->grid[1];
//	cellDimensions[2] = apix * this->grid[2];

	// 3) update f2c,c2f
	this->computeCrystParams();

//	cout << "Forcing apix to " << apix[0] << "," << apix[1] << "," << apix[2] << std::endl;
	//TR << "new origin:" << this->origin[0] << "," << this->origin[1] << "," << this->origin[2] << std::endl;

	//fpd no need for full density change trigger, just invalidate a few things
	// Fdensity and coeffs_density_ are still valid
	fastdens_score.clear();
	fastgrid = vector<int>(3,0);
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
//				int oldMinMult = MINMULT[j-1];
//				while ( std::fmod(MINMULT[j-1] , denom) > 1e-6 ) MINMULT[j-1] += oldMinMult;
				int oldMinMult = MINMULT[j];
				while ( std::fmod(MINMULT[j] , denom) > 1e-6 ) MINMULT[j] += oldMinMult;				
			}
		}
	}

	return;
}

void ElectronDensity::select_points(DOCKMol &pose, ObjexxFCL::FArray3D< float > densdata) {
//	ObjexxFCL::FArray3D< float > const & densdata = core::scoring::electron_density::getDensityMap().get_data();
	ObjexxFCL::FArray3D< std::complex<double> > Fdens, Frot;
	fourier::fft3(densdata, Fdens);

	// make rotationally averaged pose map
//	utility::vector1< core::Real > pose_1dspec;
	vector<float> pose_1dspec;
	get_spectrum( pose, pose_1dspec);
//	if ( points_defined_ ) return; // points were predefined, don't change
//	vector<vector<float>>  points_to_search_;
	points_to_search_.clear();

	ObjexxFCL::FArray3D< double > rot;
	rot.dimension( densdata.u1(), densdata.u2(), densdata.u3() );
	map_from_spectrum( pose_1dspec, rot);
	fourier::fft3(rot, Frot);

	// FFT convolution
	int npts = densdata.u1()*densdata.u2()*densdata.u3();
	for ( int i=0; i< npts; ++i ) {
		Frot[i] = Fdens[i] * std::conj( Frot[i] );  // reuse
	}
	fourier::ifft3(Frot, rot);

	// mask asu
//	if ( symminfo_.enabled() ) {
//		symminfo_.mask_asu( rot , core::scoring::electron_density::getDensityMap(), 0);
//	}

	// sort points
//	utility::vector1< std::pair< numeric::xyzVector< core::Real >, core::Real > > point_score_pairs;
	int gridStep_=2;
	vector<std::pair< vector< float >, float > > point_score_pairs;
	for ( int z=1; z<=(int)densdata.u3(); z+=gridStep_ ) {
		for ( int y=1; y<=(int)densdata.u2(); y+=gridStep_ ) {
			for ( int x=1; x<=(int)densdata.u1(); x+=gridStep_ ) {
			//	numeric::xyzVector< core::Real > x_idx(x,y,z);
				vector<float> x_idx(3,0.0);
				x_idx[0] = x;
				x_idx[1] = y;
				x_idx[2] = z;
				float dens_value = -rot(x,y,z);
				std::pair< vector<float>, float> point_score_pair = std::make_pair(x_idx, dens_value);
				point_score_pairs.push_back(point_score_pair);
			}
		}
	}
	float minDistNative = 1e4;
	float point_radius_ = 3.0;
	int topNtrans_ = 5000;
	std::sort(point_score_pairs.begin(), point_score_pairs.end(), PointScoreComparator());
	for ( int i=0; i<point_score_pairs.size(); i++ ) {
		bool hasneighbor = false;
		vector<float> x_idx = point_score_pairs[i].first;
		vector<float> x_cart(3,0.0);
		x_cart[0] = x_idx[0];
		x_cart[1] = x_idx[1];
		x_cart[2] = x_idx[2];
		idx2cart( x_idx, x_cart );
		for ( int j=0; j< points_to_search_.size(); j++ ) {
			vector<float> x_idx_stored = points_to_search_[j];
			vector<float> x_cart_stored(3,0.0);
			x_cart_stored[0] = x_idx_stored[0];
			x_cart_stored[1] = x_idx_stored[1];
			x_cart_stored[2] = x_idx_stored[2];
//			(x_idx_stored[0],x_idx_stored[1],x_idx_stored[2]);
			idx2cart( x_idx_stored, x_cart_stored );
			vector<float> tmpB(3,0.0);
			tmpB[0] = x_cart[0] - x_cart_stored[0];
			tmpB[1] = x_cart[1] - x_cart_stored[1];
			tmpB[2] = x_cart[2] - x_cart_stored[2];
			float distance = sqrt(square_len(tmpB));
			if ( distance < point_radius_ ) hasneighbor = true;
		}
		if ( !hasneighbor ) {
			points_to_search_.push_back(point_score_pairs[i].first);
		//	if ( native_ ) {
		//		numeric::xyzVector< core::Real > x_cart;
		//		core::scoring::electron_density::getDensityMap().idx2cart( x_idx, x_cart );
		//		core::Real distNative = symminfo_.min_symm_dist2(x_cart, native_com_);
		//		minDistNative = std::min( minDistNative, distNative );
		//	}
		}  
		if ( points_to_search_.size() >= topNtrans_ ) break;
	}


/*	if ( basic::options::option[ basic::options::OptionKeys::edensity::debug ]() ) {
		core::scoring::electron_density::ElectronDensity mapdmp = core::scoring::electron_density::getDensityMap();
		mapdmp.set_data(rot);
		mapdmp.writeMRC( "filter.mrc" );
		//dump the points
		std::ofstream outpoints;
		outpoints.open("selectedpoints.txt", std::ofstream::app);
		for ( Size i=1; i<=points_to_search_.size(); i++ ) {
			numeric::xyzVector< core::Real > x_cart;
			numeric::xyzVector< core::Real > x_idx = points_to_search_[i];
			core::scoring::electron_density::getDensityMap().idx2cart( x_idx, x_cart );
			outpoints << "ATOM " << utility::to_string(i) << " " << x_cart[0] << " " << x_cart[1] << " " << x_cart[2] << std::endl;
		}
		outpoints.close();
	}  */

//	if ( native_ ) TR << "Closest point to native: " << std::sqrt(minDistNative) << std::endl;
}

void ElectronDensity::get_spectrum(DOCKMol &pose, vector<float> &pose_1dspec) {
	vector<float> com(3,0.0);
	float extent = get_radius( pose, com );

	// grid spacing == delR
	int delR_=2;
	int ngrid = int (std::ceil( float(extent / float (delR_)) + 2.0));
	vector<float>().swap(pose_1dspec);
	pose_1dspec=vector<float>(ngrid, 0.0);
	vector<float> pose_1dspec_for_nrsteps = pose_1dspec;

	float massSum=0.0;
/*	bool convolute_single_residue_ = false;
	// for fine point selection... generate pose_1dspec based on the middle residue of the pose.
	if ( convolute_single_residue_ == true ) {
	//	int midres = (pose.size()+1)/2;
		int midres= int((pose.size()/3 -1)/2) ;
	//	while ( pose.residue(midres).is_virtual_residue() ) {
	//		++midres;
	//	} 
	//	core::conformation::Residue const & rsd( pose.residue(midres) );
	//	core::conformation::Atom const & residue_CA( rsd.atom(2) );
		for ( int j=3*midres; j<= 3*midres+2; ++j ) {
	//		core::conformation::Atom const & atom( rsd.atom(j) );
			float binnum = ( pose[j].x_-residue_CA.xyz() ).length() / delR_ + 1.0;
			core::Real binnum = ( atom.xyz()-residue_CA.xyz() ).length() / delR_ + 1.0;
			core::Real fpart = binnum - std::floor(binnum);
			core::Size binint = (core::Size) std::floor(binnum);
			pose_1dspec[binint] += (1-fpart);
			pose_1dspec[binint+1] += (fpart);
		}
	} else { */
		// for coarse point selection... generate pose_1dspec based on whole pose
		for ( int i=0; i<pose.num_atoms; ++i ) {
	//		core::conformation::Residue const & rsd( pose.residue(i) );
	//		if ( rsd.aa() == core::chemical::aa_vrt ) continue;
	//		for ( core::Size j=1; j<= rsd.nheavyatoms(); ++j ) {
	//			core::conformation::Atom const & atom( rsd.atom(j) );
	//			core::Real binnum = ( atom.xyz()-com ).length() / delR_ + 1.0;
			vector<float> tmpC(3,0.0);
			tmpC[0] = com[0] - pose.x[i];
			tmpC[1] = com[1] - pose.y[i];
			tmpC[2] = com[2] - pose.z[i];
			float binnum = sqrt(square_len( tmpC )) / delR_ + 1.0;
			float fpart = binnum - std::floor(binnum);
			int binint = (int) std::floor(binnum);
			pose_1dspec[binint] += (1-fpart);
			pose_1dspec[binint+1] += (fpart);
	//		}
		}
//	}

	// this is to calculate massSum via full pose for nRsteps_
	for ( int i=0; i<pose.num_atoms; ++i ) {
//		core::conformation::Residue const & rsd( pose.residue(i) );
//		if ( rsd.aa() == core::chemical::aa_vrt ) continue;
//		for ( core::Size j=1; j<= rsd.nheavyatoms(); ++j ) {
//			core::conformation::Atom const & atom( rsd.atom(j) );
			vector<float> tmpC(3,0.0);
			tmpC[0] = com[0] - pose.x[i];
			tmpC[1] = com[1] - pose.y[i];
			tmpC[2] = com[2] - pose.z[i];
			float binnum = sqrt(square_len( tmpC))/ delR_ + 1.0;
			float fpart = binnum - std::floor(binnum);
			int binint = (int) std::floor(binnum);
			pose_1dspec_for_nrsteps[binint] += (1-fpart);
			pose_1dspec_for_nrsteps[binint+1] += (fpart);
			massSum += 1;
//		}
	}

	// now setting nRsteps_
	float  fracDens=0.7; // choose radius covering this fraction of density mass
	float running_total=0.0;
	for ( int i=0; i<(int)ngrid; ++i ) {
		running_total += pose_1dspec_for_nrsteps[i] / massSum;
		if ( running_total > fracDens ) {
			nRsteps_ = i;
			break;
		}
	}

	running_total=0.0;
	for ( int i=0; i<(int)ngrid; ++i ) {
		running_total += pose_1dspec[i] / massSum;
		pose_1dspec[i] /= (i*i); // normalize
//		  cout<< "spectrum " << i << ": " << pose_1dspec[i] << " " << pose_1dspec[i]*i*i/massSum << " " << running_total << std::endl;
	}
}


void ElectronDensity::map_from_spectrum( vector<float> &pose_1dspec, ObjexxFCL::FArray3D< double > &rot) {
	// ugly
	int delR_ =2;
	for ( int z=1; z<=(int)rot.u3(); ++z ) {
		for ( int y=1; y<=(int)rot.u2(); ++y ) {
			for ( int x=1; x<=(int)rot.u1(); ++x ) {
				vector<float> idxX(3,0.0), cartX(3,0.0);
				idxX[0] = x<=rot.u1()/2 ? x-1 : x-1-rot.u1();
				idxX[1] = y<=rot.u2()/2 ? y-1 : y-1-rot.u2();
				idxX[2] = z<=rot.u3()/2 ? z-1 : z-1-rot.u3();

				vector<float> fracX(3,0.0);
				fracX[0] = ( (float) idxX[0] ) / grid[0];
				fracX[1] = ( (float) idxX[1] ) / grid[1];
				fracX[2] = ( (float) idxX[2] ) / grid[2];
			//	cartX = f2c*fracX;
				MatrixTimesTransVector(f2c,fracX,cartX);
			//	idxoffset2cart( idxX, cartX );
				float d = sqrt(square_len(cartX)) / delR_;
				float fpart = d - std::floor(d);
				int dint = (int) std::floor(d) + 1;
				if ( dint<pose_1dspec.size() ) { // last entry is always 0 so this check is valid
					rot(x,y,z) = (1-fpart)*pose_1dspec[dint] + (fpart)*pose_1dspec[dint+1];
				} else {
					rot(x,y,z) = 0.0;
				}
			}
		}
	}

//	if ( basic::options::option[ basic::options::OptionKeys::edensity::debug ]() ) {
//		core::scoring::electron_density::ElectronDensity(rot,1.0, numeric::xyzVector< core::Real >(0,0,0), false).writeMRC( "spectrum.mrc" );
//	}
}

float ElectronDensity::get_radius(DOCKMol & pose, vector<float> &com ){
	vector<float> centerCA(3,0.0);
	com = vector<float> (3,0.0);
	int nAtms=0;
	for ( int i=0; i< pose.num_atoms; ++i ) {
	//	core::conformation::Residue const & rsd( pose.residue(i) );
	//	if ( rsd.aa() == core::chemical::aa_vrt ) continue;
	//	for ( core::Size j=1; j<= rsd.nheavyatoms(); ++j ) {
	//		core::conformation::Atom const & atom( rsd.atom(j) );
			com[0] += pose.x[i];
			com[1] += pose.y[i];
			com[2] += pose.z[i];
		//	if ( i==(pose.size()+1)/2  ) {
			if ( i== 3*int((pose.num_atoms /3 -1)/2) +1) {
				int ii= 3*int((pose.num_atoms /3 -1)/2) +1;
				centerCA[0] = pose.x[ii];
				centerCA[1] = pose.y[ii];
				centerCA[2] = pose.z[ii];
			}
			nAtms++;
	//	}
	}
	com[0] = com[0]/nAtms;
	com[1] = com[1]/nAtms;
	com[2] = com[2]/nAtms;

//	if ( center_on_middle_ca_ ) {
//		com = centerCA;
//	}

	float maxSize = 0;
	for ( int i=0; i< pose.num_atoms; ++i ) {
//		core::conformation::Residue const & rsd( pose.residue(i) );
//		if ( rsd.aa() == core::chemical::aa_vrt ) continue;
//		for ( core::Size j=1; j<= rsd.nheavyatoms(); ++j ) {
//			core::conformation::Atom const & atom( rsd.atom(j) );
		vector<float> tmpCA(3,0.0);
		tmpCA[0] = com[0] - pose.x[i];
		tmpCA[1] = com[1] - pose.y[i];
		tmpCA[2] = com[2] - pose.z[i];
		maxSize = std::max( maxSize, square_len(tmpCA));	
//			maxSize = std::max( maxSize, (com-atom.xyz()).length_squared());
//		}
	}
	maxSize = std::sqrt( maxSize ); // this the radius
	return maxSize;
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
//	this->f2c[2][2] = cellDimensions[2] *(sqrt(1.0-ca*ca-cb*cb-cg*cg+2.0*ca*cb*cg))/sg;

	float D = f2c[0][0] * ( (f2c[1][1]*f2c[2][2]) - (f2c[2][1]*f2c[1][2]) ) -
		      f2c[0][1] * ( (f2c[1][0]*f2c[2][2]) - (f2c[2][0]*f2c[1][2]) ) +
		 	  f2c[0][2] * ( (f2c[1][0]*f2c[2][1]) - (f2c[2][0]*f2c[1][1]) );
//	cout<<" f2c: "<<f2c[0][0]<<" "<<f2c[0][1]<<" "<<f2c[0][2]<<" "<<endl;
	if(D == 0 ) {
		cout<< "Warning, Invalid crystal cell dimensions. " << endl;
		return;
	}
	
	float DD = cellDimensions[0] *cellDimensions[1] *cellDimensions[2] *(sqrt(1.0-ca*ca-cb*cb-cg*cg+2.0*ca*cb*cg));
	if(sg ==0||DD==0)
	{
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
/*	this->c2f[0][0] = 1.0/(cellDimensions[0]);
	this->c2f[0][1] = -cg/(cellDimensions[0]*sg);
	this->c2f[0][2] = cellDimensions[1]*cellDimensions[2]*(ca*cg-cb)/(DD*sg);
	this->c2f[1][0] = 0;
	this->c2f[1][1] = 1.0/(cellDimensions[1]*sg);
	this->c2f[1][2] = cellDimensions[0]*cellDimensions[2]*(cb*cg-ca)/(DD*sg);
	this->c2f[2][0] = 0;
	this->c2f[2][1] = 0;
	this->c2f[2][2] = cellDimensions[0]*cellDimensions[1]*sg/(DD);	 */
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

