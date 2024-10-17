// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author

#ifndef INCLUDED_xray_scattering_h
#define INCLUDED_xray_scattering_h

//#include <core/types.hh>

#include <cmath>
#include <string>


#ifdef WIN32

#define _USE_MATH_DEFINES

#include <math.h>

#endif

//#define M_PI 3.1415926
// C++ Headers



////////////
// 1-gaussian real-space approximate scattering
class OneGaussianScattering {
public:
	OneGaussianScattering() {
		weight_ = 0;
		sigma_ = 3.0;
	}
	OneGaussianScattering(float w, float s) {
		sigma_ = s;
		weight_ = w;
	}

	inline float B( float k ) const {
		float s = M_PI*M_PI/k;
		float sigma_eff = sigma_;
		float B = ( 4*( s - sigma_eff ) );

		// smooth to flat at B==0
		if ( B < 0 ) B = 0;
		else if ( B<10 ) B = sqrt(10*B);
		return ( B );
	}

	inline float k( float B , float lim=600) const {
		float sigma_eff = sigma_;
		float B_eff = B;
		if ( B<0 ) B_eff = 0;
		else if ( B<1 ) B_eff = B*B;
		else if ( B>lim-10.0 && B<=lim ) B_eff = lim-(1.0/10.0)*(lim-B)*(lim-B);
		else if ( B>lim ) B_eff = lim;
		float s = sigma_eff + B_eff/4;
		float k = M_PI*M_PI/s;
		return k;
	}

	// calculate dK/dB at a given resolution
	inline float dk( float B , float lim=600) const {
		float sigma_eff = sigma_;

		float B_eff = B;
		if ( B<0 || B>lim ) return 0;
		else if ( B<1 ) B_eff = B*B;
		else if ( B>lim-10.0 ) B_eff = lim-(1.0/10.0)*(lim-B)*(lim-B);

		float s = sigma_eff + B_eff/4;
		float dkdb = -M_PI*M_PI/(4*s*s);

		if ( B<1 ) dkdb *= 2*B;
		if ( B>lim-10.0 ) dkdb *= (1.0/5.0)*(lim-B);

		return dkdb;
	}

	// rho = C*exp(-k*X^2)
	// get scale factor, given k
	inline float C( float k ) const {
		float C = pow(k, 1.5);
		return ( C*weight_ );
	}

	inline int a( ) const {
		return ( weight_ );
	}

private:
	float sigma_;
	float weight_;
};


////////////
// 4-gaussian approximation of scattering coefficients
class KromerMann {
public:
	KromerMann() {
		a_[0] = a_[1] = a_[2] = a_[3] = 0;
		b_[0] = b_[1] = b_[2] = b_[3] = 0;
		c_ = 0;
	}
	KromerMann(float c, float a1, float a2, float a3, float a4, float b1, float b2, float b3, float b4) {
		a_[0] = a1; a_[1] = a2; a_[2] = a3; a_[3] = a4;
		b_[0] = b1; b_[1] = b2; b_[2] = b3; b_[3] = b4;
		c_ = c;
	}

	// scattering at reciprocal space distance^2
	inline float f0( float S2 ) {
		return (c_ + a_[0]*exp(b_[0]*S2/4) + a_[1]*exp(b_[1]*S2/4) + a_[2]*exp(b_[2]*S2/4) + a_[3]*exp(b_[3]*S2/4) );
	}

private:
	float a_[4];
	float b_[4];
	float c_;
};

// weight from scattering factors
OneGaussianScattering get_A( std::string elt, bool cryoem_scatterers=true);

// weight from scattering factors
KromerMann get_km( std::string elt );



/*
bool factorsLTE5(int X);
bool factorsLTE19(int X);
int findSampling(double MINSMP, int NMUL);
int findSampling5(double MINSMP, int NMUL);

*/

#endif
