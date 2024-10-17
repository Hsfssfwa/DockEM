#include "CaDensity.h"
using namespace std;

/*
// map atom names to elements
//   loosely based on openbabel logic
std::string
name2elt( std::string line ) {
	std::string atmid = line.substr(12,4);
	while ( !atmid.empty() && atmid[0] == ' ' ) atmid = atmid.substr(1,atmid.size()-1);
	while ( !atmid.empty() && atmid[atmid.size()-1] == ' ' ) atmid = atmid.substr(0,atmid.size()-1);
	std::string resname = line.substr(17,3);
	while ( !resname.empty() && resname[0] == ' ' ) resname = resname.substr(1,resname.size()-1);
	while ( !resname.empty() && resname[resname.size()-1] == ' ' ) resname = resname.substr(0,resname.size()-1);

	std::string type;

	if ( line.substr(0,4) == "ATOM" ) {
		type = atmid.substr(0,2);
		if ( isdigit(type[0]) ) {
			// sometimes non-standard files have, e.g 11HH
			if ( !isdigit(type[1]) ) type = atmid.substr(1,1);
			else type = atmid.substr(2,1);
		} else if ( (line[12] == ' ' && type!="Zn" && type!="Fe" && type!="ZN" && type!="FE")
				|| isdigit(type[1]) ) {
			type = atmid.substr(0,1);     // one-character element
		}

		if ( resname.substr(0,2) == "AS" || resname[0] == 'N' ) {
			if ( atmid == "AD1" ) type = "O";
			if ( atmid == "AD2" ) type = "N";
		}
		if ( resname.substr(0,3) == "HIS" || resname[0] == 'H' ) {
			if ( atmid == "AD1" || atmid == "AE2" ) type = "N";
			if ( atmid == "AE1" || atmid == "AD2" ) type = "C";
		}
		if ( resname.substr(0,2) == "GL" || resname[0] == 'Q' ) {
			if ( atmid == "AE1" ) type = "O";
			if ( atmid == "AE2" ) type = "N";
		}
		if ( atmid.substr(0,2) == "HH" ) { // ARG
			type = "H";
		}
		if ( atmid.substr(0,2) == "HD" || atmid.substr(0,2) == "HE" || atmid.substr(0,2) == "HG" ) {
			type = "H";
		}
	} else {
		if ( isalpha(atmid[0]) ) {
			if ( atmid.size() > 2 && (atmid[2] == '\0' || atmid[2] == ' ') ) {
				type = atmid.substr(0,2);
			} else if ( atmid[0] == 'A' ) { // alpha prefix
				type = atmid.substr(1, atmid.size() - 1);
			} else {
				type = atmid.substr(0,1);
			}
		} else if ( atmid[0] == ' ' ) {
			type = atmid.substr(1,1); // one char element
		} else {
			type = atmid.substr(1,2);
		}

		if ( atmid == resname ) {
			type = atmid;
			if ( type.size() == 2 ) type[1] = toupper(type[1]);
		} else if ( resname == "ADR" || resname == "COA" || resname == "FAD" ||
				resname == "GPG" || resname == "NAD" || resname == "NAL" ||
				resname == "NDP" || resname == "ABA" )  {
			if ( type.size() > 1 ) type = type.substr(0,1);
		} else if ( isdigit(type[0]) ) {
			type = type.substr(1,1);
		} else if ( type.size() > 1 && isdigit(type[1]) ) {
			type = type.substr(0,1);
		} else if ( type.size() > 1 && isalpha(type[1]) ) {
			if ( type[0] == 'O' && type[1] == 'H' ) {
				type = type.substr(0,1); // no "Oh" element (e.g. 1MBN)
			} else if ( islower(type[1]) ) {
				type[1] = toupper(type[1]);
			}
		}
	}
	return type;
} */

/*
// quick and dirty PDB read where we only care about heavyatom locations and atom ids
void
readPDBcoords(std::string filename, poseCoords &atmlist) {
	std::ifstream inpdb(filename.c_str());
//	cout<< filename.c_str()<<endl;
	std::string buf;

	while ( std::getline(inpdb, buf ) ) {
		if ( buf.substr(0,4)!="ATOM" && buf.substr(0,6)!="HETATM" ) continue;
		poseCoord atom_i;

		atom_i.x_.push_back(atof(buf.substr(30,8).c_str()));
		atom_i.x_.push_back(atof(buf.substr(38,8).c_str()));
		atom_i.x_.push_back(atof(buf.substr(46,8).c_str()));
		atom_i.B_ = atof(buf.substr(60,6).c_str());

		atom_i.elt_ = name2elt( buf ); // horrible hacky logic mapping name->elt (could use PDB fields 76-77 if on by default)
		if ( atom_i.elt_ == "H" ) continue;

		atmlist.push_back( atom_i );
	}
} */

int main(int argc, char* argv[])
{
	if(argc<6)
	{
		cout<<"the argument is too smaller"<<endl;
	}else
	{
		Model protein;
		char inputPDB[300];   
		string inputMRC;
		float MRC_reso;
		float mapsampling = 0.0;
//		float grid_spacing = 2.0;
		char outputPDB[300];





		strcpy(inputPDB,argv[1]);   // input PDB structure
//		strcpy(inputMRC,argv[2]);  // input density map
		inputMRC = argv[2];
//		cout<< " inputMRC: "<<inputMRC<<endl;
//		strcpy(MRC_reso,argv[3]); // denstiy map resolution
		MRC_reso = atof(argv[3]);
//		cout<< " MRC_reso: "<<MRC_reso<<endl;
//		strcpy(mapsampling,argv[4]); // desnsity map sampling;
		mapsampling = atof(argv[4]);
//		cout<< " aaaa"<<endl;
		strcpy(outputPDB,argv[5]);  //  output PDB structure
//		cout<<"  outputPDB:"<< outputPDB <<endl;

		// read PDB
		poseCoords pose;
		Model pose1;
		cout<< "Read PDB file : "<<inputPDB<<endl;
//		readpdbstructure(inputPDB,protein);
		readPDBcoords(inputPDB,pose);
		readpdbstructurex(inputPDB,pose1);
		cout<< " bbbb "<<endl;
/*		for(int i=0;i<5;i++)
		{
			cout<<" coord: "<< pose[i].x_[0]<<"  "<<pose[i].x_[1]<<" "<<pose[i].x_[2]<<"  B: "<< pose[i].B_ <<"  elt: "<<pose[i].elt_<<endl;
		}  */

		int nresobins = 200;  //  # resolution bins for statistics
		vector<float> resobins, mapI, maskedmapI, modelI, maskI;
		vector<int> resobin_counts;
		vector<float> perResCC, perResStrain;

		vector<float> mapmapFSC, maskedMapMapFSC;
		vector<float> modelmapFSC, maskedModelMapFSC;

		bool truncate_map = false ;
		bool maskonly = false ;
		bool cutonly = false ;
		bool bin_squared =  false;
		bool mask =false ;
		string alt_mapfile=" ";

		// resolution limits for analysis
		float hires = 0.0;  // high res limit
		float lowres = 1000.0;  // low res limit
		float truncate_hires = 0.0;  // high res truncation
		float truncate_lowres = 1000.0;  // low res truncation
		float mask_resolution = 0.0;   //  radius for masking
		bool perres =false;  //dump extra output

		ObjexxFCL::FArray3D< float > rhoC, rhoMask, rhoOmask, rhoO2mask;
		ObjexxFCL::FArray3D< std::complex<float> > FrhoC, FrhoMask, FrhoCmask, FrhoOmask, FrhoO2mask;
		ObjexxFCL::FArray3D< std::complex<float> > FrhoO, FrhoO2;	

		// read Density map
		ElectronDensity theDensityMap;
		theDensityMap.preso= MRC_reso;
		theDensityMap.pATOM_MASK =4.0;
		theDensityMap.pCA_MASK =6.0;
		theDensityMap.pforce_apix_on_map_load_=0.0;
		theDensityMap.pWINDOW_ =1;
		theDensityMap.pscore_window_context_ = false ;
		theDensityMap.premap_symm_ = false ;
		theDensityMap.pde_edensity = false;
		theDensityMap.pnkbins_ = 0 ;   
/*		theDensityMap.reso= MRC_reso;
		theDensityMap.ATOM_MASK =0.0;
		theDensityMap.CA_MASK =0.0;
		theDensityMap.WINDOW_ =1;
		theDensityMap.score_window_context_ = false ;
		theDensityMap.remap_symm_ = false ;
		theDensityMap.force_apix_on_map_load_ =0.0;
		theDensityMap.nkbins_ = 0 ;   */
		theDensityMap.readMRCandResize(inputMRC,MRC_reso,mapsampling);
		theDensityMap.calcRhoC(pose,4.0,rhoC,rhoMask);
		float CC;
		CC=theDensityMap.getRSCC(rhoC,rhoMask);
		cout<<"CC: "<<CC<<endl;
//		ObjexxFCL::FArray3D< float > Deny;
//		Deny=theDensityMap.get_data();
//		cout<<"Deny: "<<" "<<Deny.u1()<<" "<<Deny.u2()<<" "<<Deny.u3()<<endl;
/*		for(int i=1;i<=Deny.u1();i++)
		{
			for(int j=1;j<=Deny.u2();j++)
			{
				for(int k=1;k<=Deny.u3();k++)
				{
					cout<<Deny(i,j,k)<<" ";
				}
			}
		}
		cout<<endl;   */


	}
}

///***** complie
/// /opt/rh/devtoolset-7/root/usr/bin/g++ UseDensity.cpp CaDensity.cpp SplineInterp.cpp GeometryTools.cpp readpdb.cpp -o UseDensity -I/nfs/amino-home/zbiao/LibC/tt/CaDensity -L/nfs/amino-home/zbiao/LibC/tt/CaDensity/fourier -lfourier

/// ./UseDensity 1a0fA.pdb 1a0fA.mrc 5.0 0.0 11.pdb
