void ElectronDensity::calcRhoCy(
	Model const &pose,
	float highres_limit,
	ObjexxFCL::FArray3D< float > &rhoC,
	ObjexxFCL::FArray3D< float > &mask,
	float fixed_mask_B /* = -1 */,
	float B_upper_limit /* = 600 */,
	float force_mask /*=-1*/ ){
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
			for(int t=0; t<pose.chains[i].residues[j].atoms.size();t++)
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
					per_atom_ks = std::min ( per_atom_ks, 4*M_PI*M_PI/minimumB );
				} else {
					per_atom_ks = std::min ( per_atom_ks, 4*M_PI*M_PI/effectiveB );
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