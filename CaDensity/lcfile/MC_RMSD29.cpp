#include "MC.h"
#include "GetVdwRadius.h"
#include <time.h>
#include "Rotbuilder.h"
#include "CaDensity.h"

using namespace std;
std::string
name2eltxx( std::string line ) {
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
}

vector<string> string_splitxx(string const & in, char splitchar)
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

void writePDBcoordss(std::string filename, poseCoords &atmlist)
{
    std::ifstream inpdb(filename.c_str());
//  cout<< filename.c_str()<<endl;
    std::string buf;
    vector<string> tmp_str;
    tmp_str = string_splitxx(filename,'.');
    string filenamex = tmp_str[0] + "_out.pdb";
    std::ofstream outpdb;
    outpdb.open(filenamex.c_str());
    int tmpx=0;

    while ( std::getline(inpdb, buf ) ) {
//        if ( buf.substr(0,4)!="ATOM" && buf.substr(0,6)!="HETATM" ) continue;
    	if(buf.length()<54) continue;
        if ( buf.substr(0,4) =="ATOM" && buf.substr(13,1)!="H") {

	        poseCoord atom_i;
	        string tmp_A,tmp_B;

	        tmp_A = buf.substr(0,30);
	        atom_i.x_.push_back(atof(buf.substr(30,8).c_str()));
	        atom_i.x_.push_back(atof(buf.substr(38,8).c_str()));
	        atom_i.x_.push_back(atof(buf.substr(46,8).c_str()));
	//        tmp_B = buf.substr(54,6);
	//        atom_i.B_ = atof(buf.substr(60,6).c_str());
	    //    cout<<"YYY"<<endl;


	//        atom_i.elt_ = name2elt( buf ); // horrible hacky logic mapping name->elt (could use PDB fields 76-77 if on by default)
	//        if ( atom_i.elt_ == "H" ) continue;

	        outpdb<<setw(30)<<tmp_A;
	        outpdb<<setw(8)<<setiosflags(ios::fixed)<<setprecision(3)<<atmlist[tmpx].x_[0];
	        outpdb<<setw(8)<<setiosflags(ios::fixed)<<setprecision(3)<<atmlist[tmpx].x_[1];
	        outpdb<<setw(8)<<setiosflags(ios::fixed)<<setprecision(3)<<atmlist[tmpx].x_[2];
	//        outpdb<<setw(6)<<tmp_B;
	//        outpdb<<setw(6)<<setiosflags(ios::fixed)<<setprecision(2)<<atmlist[tmpx].B_;
	        outpdb<<endl;
	        tmpx = tmpx + 1;        
	//        atmlist.push_back( atom_i );
    	}
    }
    outpdb.close();
}

// quick and dirty PDB read where we only care about heavyatom locations and atom ids
void
readPDBcoordss(std::string filename, poseCoords &atmlist) {
    std::ifstream inpdb(filename.c_str());
//  cout<< filename.c_str()<<endl;
    std::string buf;

    int nn=0;
    while ( !inpdb.eof() && std::getline(inpdb, buf) ) {
//    	cout<<nn << " ";
    	if(buf.length()<54) continue;
//    	if ( buf.substr(0,4) =="ATOM" && buf.substr(12,4) ==" CA " ) {
    	if ( buf.substr(0,4) =="ATOM" && buf.substr(13,1) !="H" && buf.substr(13,3) !="OXT") {
	        poseCoord atom_i;

	        atom_i.x_.push_back(atof(buf.substr(30,8).c_str()));
	        atom_i.x_.push_back(atof(buf.substr(38,8).c_str()));
	        atom_i.x_.push_back(atof(buf.substr(46,8).c_str()));
	        atom_i.B_ = 0.0;// atof(buf.substr(60,6).c_str())
	        atom_i.elt_ = buf.substr(12,4);

//	        nn=nn+1;
//	        cout<< nn<<" ";
	        atmlist.push_back( atom_i );
    	}
    }
}

vector<float> Generate_random_point(float A,float B,float r)
{
	vector<float> coordx;
	float u=A+(B-A)*rand()/(RAND_MAX+1.0) ; // [a,b] random number,(rand()%(b-a+1))+a
	float v=A+(B-A)*rand()/(RAND_MAX+1.0);
	float theta = 2.0*MX_PI*u;
	float phi=acos(2*v-1.0);
	float x=r*sin(phi)*sin(theta);
	float y=r*cos(phi)*sin(theta);
	float z=r*cos(theta);
	coordx.push_back(x);
	coordx.push_back(y);
	coordx.push_back(z);

	return coordx;
}

float randomAB(float A,float B)
{
	return A+(B-A)*rand()/(RAND_MAX);
}

float Distance_point(vector<float> A,vector<float> B)
{
	float dis=0.0;
	for(int i=0;i<A.size();i++)
	{
		dis = dis + (A[i]-B[i])*(A[i]-B[i]);
	}
	dis=sqrt(dis);
	return dis;
}
float Distance_pointy(vector<vector<float> > A,vector<vector<float> > B)
{

}
float Distance_pointx(vector<vector<float> > A,vector<vector<float> > B)
{
	float dis=0.0;
	float a1 = sqrt((A[0][0]-A[1][0])*(A[0][0]-A[1][0])+(A[0][1]-A[1][1])*(A[0][1]-A[1][1])+(A[0][2]-A[1][2])*(A[0][2]-A[1][2]));
	float a2 = sqrt((A[0][0]-A[2][0])*(A[0][0]-A[2][0])+(A[0][1]-A[2][1])*(A[0][1]-A[2][1])+(A[0][2]-A[2][2])*(A[0][2]-A[2][2]));
	float a3 = sqrt((A[0][0]-A[3][0])*(A[0][0]-A[3][0])+(A[0][1]-A[3][1])*(A[0][1]-A[3][1])+(A[0][2]-A[3][2])*(A[0][2]-A[3][2]));
	float b1 = sqrt((B[0][0]-B[1][0])*(B[0][0]-B[1][0])+(B[0][1]-B[1][1])*(B[0][1]-B[1][1])+(B[0][2]-B[1][2])*(B[0][2]-B[1][2]));
	float b2 = sqrt((B[0][0]-B[2][0])*(B[0][0]-B[2][0])+(B[0][1]-B[2][1])*(B[0][1]-B[2][1])+(B[0][2]-B[2][2])*(B[0][2]-B[2][2]));
	float b3 = sqrt((B[0][0]-B[3][0])*(B[0][0]-B[3][0])+(B[0][1]-B[3][1])*(B[0][1]-B[3][1])+(B[0][2]-B[3][2])*(B[0][2]-B[3][2]));

/*	for(int i=0;i<A.size();i++)
	{
		for(int j=0;j<A[i].size();j++)
		{
//			dis = dis + (A[i][j]-B[i][j])*(A[i][j]-B[i][j]);	
			dis = dis + (a1-b1)*(a1-b1)+		
		}
		
	} */
	dis = (a1-b1)*(a1-b1) + (a2-b2)*(a2-b2) + (a3-b3)*(a3-b3);
	dis=sqrt(dis);
	return dis;
} 


vector<float> calculateCoordinatesx(vector<float> refA,vector<float> refB,vector<float> refC,float L,float ang,float di)
{
	vector<float> AV(3,0.0),BV(3,0.0),CV(3,0.0);
	AV=refA;
	BV=refB;
	CV=refC;

	vector<float> CA(3,0.0),CB(3,0.0); 
	CA[0]=AV[0]-CV[0];
	CA[1]=AV[1]-CV[1];
	CA[2]=AV[2]-CV[2];
	CB[0]=BV[0]-CV[0];
	CB[1]=BV[1]-CV[1];
	CB[2]=BV[2]-CV[2];

	float AX = CA[0];
	float AY = CA[1];
	float AZ = CA[2];

	float BX = CB[0];
	float BY = CB[1];
	float BZ = CB[2];

	float A = (AY*BZ)-(AZ*BY);
	float B = (AZ*BX)-(AX*BZ);
	float G = (AX*BY)-(AY*BX);

	float F = sqrt(BX*BX + BY*BY + BZ*BZ) * L * cos(ang * (3.1415926/180.0));
//	float F = sqrt(BX*BX + BY*BY + BZ*BZ) * L * cos(ang);

	float constx = sqrt((pow(B*BZ-BY*G,2))*(-(F*F)*(A*A+B*B+G*G)+(B*B*(BX*BX+BZ*BZ) + A*A*(BY*BY+BZ*BZ)- (2*A*BX*BZ*G) + (BX*BX+ BY*BY)*G*G - (2*B*BY)*(A*BX+BZ*G))*L*L));
	float denom = (B*B)*(BX*BX+BZ*BZ)+ (A*A)*(BY*BY+BZ*BZ) - (2*A*BX*BZ*G) + (BX*BX+BY*BY)*(G*G) - (2*B*BY)*(A*BX+BZ*G);

	float X= ((B*B*BX*F)-(A*B*BY*F)+(F*G)*(-A*BZ+BX*G)+constx)/denom;

	float Y=0.0,Z=0.0;
    if((B==0 || BZ==0) && (BY==0 || G==0))
    {
        float const1=sqrt( G*G*(-A*A*X*X+(B*B+G*G)*(L-X)*(L+X)));
        Y= ((-A*B*X)+const1)/(B*B+G*G);
        Z= -(A*G*G*X+B*const1)/(G*(B*B+G*G));
    }
    else
    {
        Y= ((A*A*BY*F)*(B*BZ-BY*G)+ G*( -F*pow(B*BZ-BY*G,2) + BX*constx) - A*( B*B*BX*BZ*F- B*BX*BY*F*G + BZ*constx)) / ((B*BZ-BY*G)*denom);
        Z= ((A*A*BZ*F)*(B*BZ-BY*G) + (B*F)*pow(B*BZ-BY*G,2) + (A*BX*F*G)*(-B*BZ+BY*G) - B*BX*constx + A*BY*constx) / ((B*BZ-BY*G)*denom);
	}
	vector<float> D(3,0.0);
	D[0] = X + CV[0];
	D[1] = Y + CV[1];
	D[2] = Z + CV[2];

	float temp=0.0;

	temp=Points2Dihedral(AV,BV,CV,D) * (180.0/3.1415926);
//	cout<<"temp: "<<temp<<endl;

	di = di - temp;
	vector<float> Dx(3,0.0);
	Dx[0] = D[0] - BV[0];
	Dx[1] = D[1] - BV[1];
	Dx[2] = D[2] - BV[2];

	di = di * 3.1415926/180.0;
	GroupRotationt(CV,BV,di,Dx);
	D[0] = Dx[0] + BV[0];
	D[1] = Dx[1] + BV[1];
	D[2] = Dx[2] + BV[2];

	return D;
}

//Model REMC_sampley(Model posex, Model posey)
Model REMC_sampley(Model posex, ElectronDensity theDensityMap)
{
//	cout<<" 222 "<<endl;
//	int num=0;
	Model fin_pose;
	fin_pose = posex;
//	vector<Rot> fin_vec;
//	vector<Rot> beg_vec;
//	float min_dis=100.0;
//	vector<vector<float> > point;
//	vector<float> old_v(3,0.0),old_v1(3,0.0);
//	vector<vector<float> > old_v,old_v1,old_tmp,new_v;
//	old_v=vector<vector<float> >(4,vector<float>(3,0.0));
//	old_v1=vector<vector<float> >(4,vector<float>(3,0.0));
//	old_tmp = vector<vector<float> >(4,vector<float>(3,0.0));
//	new_v = vector<vector<float> >(4,vector<float>(3,0.0));
//	vector<float> new_v(3,0.0);
	float KT=0.005;
//	vector<vector<Rot> > frag_pre(7),frag_nat(7);


	
	vector<Rot> pointsx;
	vector<vector<float> > pointsABx; // all backbone atom
	vector<poseCoord> pointsBx;
	vector<vector<float> > pointsABx0; 
	vector<string> pointsrenm;
//	vector<vector<float> > pointsx,pointsy;
//	int pnum = posex.size();

//	cout<<"33 "<<endl;
//	cout<< posex.size()<<endl;
//	cout<< posey.size()<<endl;
//	for(int j=0;j<pnum;j++)
//	{
//		pointsx.push_back(posex[j].x_);
//		pointsy.push_back(posey[j].x_);
//	}
	int pnum =0,pnumx=0;
//	cout<< " 1111"<<endl;
	for(int i=0;i<posex.chains.size();i++)
	{
//		cout<< "chainx size: "<<posex.chains.size()<<endl;
//		cout<< "chainy size: "<<posey.chains.size()<<endl;
		Chain Chanx = posex.chains[i];
//		Chain Chany = posey.chains[i];
		for(int j=0;j<Chanx.residues.size();j++)
		{
//			cout<< "residuex size: "<< Chanx.residues.size()<<endl;
//			cout<< "residuey size: "<< Chany.residues.size()<<endl;
			Residue Resdx = Chanx.residues[j];
	//		Residue Resdy = Chany.residues[j];
			Rot rotx;
//			cout<<"Res: "<<Resdx.atoms.size()<<endl;
			pointsrenm.push_back(Resdx.resname);
//			pointsrenmy.push_back(Resdy.resname);
			for(int t=0;t<Resdx.atoms.size();t++)
			{
//				cout<< "atomx size: "<<Resdx.atoms.size()<<endl;
				Atom Atmx=Resdx.atoms[t];
				if(Atmx.atona != " CA " && Atmx.atona != " C  " && Atmx.atona != " N  " && Atmx.atona != " O  ")
				{
					rotx.res_atm.push_back(Atmx.xyzVect);
					vector<string> atmnmx = string_splitxx(Atmx.atona,' ');
					string atmnmy = atmnmx[0];
					rotx.resatmnm.push_back(atmnmy);
				}				
				if(Atmx.atona == " CA ") {rotx.ca = Atmx.xyzVect;pointsABx.push_back(Atmx.xyzVect);pnumx = pnumx +1;}
				if(Atmx.atona == " C  ") {rotx.cc = Atmx.xyzVect;pointsABx.push_back(Atmx.xyzVect);pnumx = pnumx +1;}
				if(Atmx.atona == " N  ") {rotx.cn = Atmx.xyzVect;pointsABx.push_back(Atmx.xyzVect);pnumx = pnumx +1;}
				if(Atmx.atona == " O  ") {rotx.co = Atmx.xyzVect;}
				if(Atmx.atona == " CB ") {rotx.cb = Atmx.xyzVect;}
				if(Atmx.atona == " CG ") {rotx.cg = Atmx.xyzVect;}
				if(Atmx.atona == " CG1") {rotx.cg1 = Atmx.xyzVect;}
				if(Atmx.atona == " CG2") {rotx.cg2 = Atmx.xyzVect;}
				if(Atmx.atona == " OG ") {rotx.og  = Atmx.xyzVect;}
				if(Atmx.atona == " SG ") {rotx.sg = Atmx.xyzVect;}
				if(Atmx.atona == " CD ") {rotx.cd = Atmx.xyzVect;}
				if(Atmx.atona == " SD ") {rotx.sd = Atmx.xyzVect;}
				if(Atmx.atona == " CD1") {rotx.cd1 = Atmx.xyzVect;}
				if(Atmx.atona == " CD2") {rotx.cd2 = Atmx.xyzVect;}
				if(Atmx.atona == " ND2") {rotx.nd2 = Atmx.xyzVect;}
				if(Atmx.atona == " ND1") {rotx.nd1 = Atmx.xyzVect;}
				if(Atmx.atona == " OD2") {rotx.od2 = Atmx.xyzVect;}
				if(Atmx.atona == " OE2") {rotx.oe2 = Atmx.xyzVect;}
				if(Atmx.atona == " CE ") {rotx.ce = Atmx.xyzVect;}
				if(Atmx.atona == " NE2") {rotx.ne2 = Atmx.xyzVect;}
				if(Atmx.atona == " NE ") {rotx.ne = Atmx.xyzVect;}
				if(Atmx.atona == " NZ ") {rotx.nz = Atmx.xyzVect;}
				if(Atmx.atona == " CZ ") {rotx.cz = Atmx.xyzVect;}					
	/*			if(Atmy.atona == " CB ") roty.cb = Atmy.xyzVect;					
				if(Atmx.atona == " CG ") rotx.cg = Atmy.xyzVect;
				if(Atmy.atona == " CG ") roty.cg = Atmy.xyzVect;
				if(Atmx.atona == " CD ") rotx.cd = Atmy.xyzVect;
				if(Atmy.atona == " CD ") roty.cd = Atmy.xyzVect;
				if(Atmx.atona == " CE ") rotx.ce = Atmy.xyzVect;
				if(Atmy.atona == " CE ") roty.ce = Atmy.xyzVect;
				if(Atmx.atona == " CZ ") rotx.cz = Atmy.xyzVect;
				if(Atmy.atona == " CZ ") roty.cz = Atmy.xyzVect;		*/					

			}
	/*		for(int t=0;t<Resdy.atoms.size();t++)
			{
//				cout<< "atomy size: "<<Resdy.atoms.size()<<endl;
				Atom Atmy=Resdy.atoms[t];
				if(Atmy.atona != " CA " && Atmy.atona != " C  " && Atmy.atona != " N  " && Atmy.atona != " O  ")
				{
					roty.res_atm.push_back(Atmy.xyzVect);
				}				
				if(Atmy.atona == " CA ") {roty.ca = Atmy.xyzVect;pointsABy.push_back(Atmy.xyzVect);}
				if(Atmy.atona == " C  ") {roty.cc = Atmy.xyzVect;pointsABy.push_back(Atmy.xyzVect);}
				if(Atmy.atona == " N  ") {roty.cn = Atmy.xyzVect;pointsABy.push_back(Atmy.xyzVect);}
				if(Atmy.atona == " O  ") {roty.co = Atmy.xyzVect;}
				if(Atmy.atona == " CB ") {roty.cb = Atmy.xyzVect;}
				if(Atmy.atona == " CG ") {roty.cg = Atmy.xyzVect;}
				if(Atmy.atona == " CG1") {roty.cg1 = Atmy.xyzVect;}
				if(Atmy.atona == " CG2") {roty.cg2 = Atmy.xyzVect;}
				if(Atmy.atona == " OG ") {roty.og  = Atmy.xyzVect;}
				if(Atmy.atona == " SG ") {roty.sg = Atmy.xyzVect;}
				if(Atmy.atona == " CD ") {roty.cd = Atmy.xyzVect;}
				if(Atmy.atona == " SD ") {roty.sd = Atmy.xyzVect;}
				if(Atmy.atona == " CD1") {roty.cd1 = Atmy.xyzVect;}
				if(Atmy.atona == " CD2") {roty.cd2 = Atmy.xyzVect;}
				if(Atmy.atona == " ND2") {roty.nd2 = Atmy.xyzVect;}
				if(Atmy.atona == " ND1") {roty.nd1 = Atmy.xyzVect;}
				if(Atmy.atona == " OD2") {roty.od2 = Atmy.xyzVect;}
				if(Atmy.atona == " OE2") {roty.oe2 = Atmy.xyzVect;}
				if(Atmy.atona == " CE ") {roty.ce = Atmy.xyzVect;}
				if(Atmy.atona == " NE2") {roty.ne2 = Atmy.xyzVect;}
				if(Atmy.atona == " NE ") {roty.ne = Atmy.xyzVect;}
				if(Atmy.atona == " NZ ") {roty.nz = Atmy.xyzVect;}
				if(Atmy.atona == " CZ ") {roty.cz = Atmy.xyzVect;}				
			}			*/
			pointsx.push_back(rotx);
//			pointsy.push_back(roty);
			pnum = pnum +1;	
//			cout<< pnum<<endl;		
		}
	}
	cout<<"pnum,pnumx: "<<pnum<<" "<<pnumx<<endl;
	pointsABx0 = pointsABx;

	for(int ttt=0;ttt<1;ttt++)
	{
	vector<float> axyz;
  	axyz=vector<float>(3,0.0);
  	float ang; 
  	float angle_rotate;	
  	float old_dE1=0.0;
  	float new_dE1=0.0;
  	float old_dE2=0.0;
  	float new_dE2=0.0; 
  	float old_vwd1=0.0;
  	float old_vwd2=0.0;
  	float new_vwd1=0.0;
  	float new_vwd2=0.0; 
 	float new_dE=0.0;
  	float old_dE=0.0;  		
  	float dE=0.0;
  	vector<vector<float> > fin_mat;
  	vector<float> coord_change(3,0.0);

//	int rand_point = rand()%(pnumx);
//	axyz[0]=pointsx[rand_point].ca[0];
//	axyz[1]=pointsx[rand_point].ca[1];
//	axyz[2]=pointsx[rand_point].ca[2];  	
	for(int i=0;i<4000;i++)
	{
	  	// Rotation axis 
		float asin_theta=2.0*RandomDoubleX(0.0,1.0)-1.0;
		float acos_theta=sqrt(1.0-asin_theta*asin_theta);
		float apha=2.0*PI*RandomDoubleX(0.0,1.0);
		float awx=acos_theta*cos(apha);
		float awy=acos_theta*sin(apha);
		float awz=asin_theta;
		// Translation Vector
		float t0=0.3;
		float t1=(RandomDoubleX(0.0,1.0)*2.0-1.0)*t0+0.0;
		float t2=(RandomDoubleX(0.0,1.0)*2.0-1.0)*t0+0.0;
		float t3=(RandomDoubleX(0.0,1.0)*2.0-1.0)*t0+0.0;

		// Rotation matrix
		ang=1.0;
		angle_rotate=(2.0*RandomDoubleX(0.0,1.0)-1.0)*ang; // rotate angle
		float asin=sin(angle_rotate);
		float acos=cos(angle_rotate);
		float u[3][3];
		u[0][0]=acos+awx*awx*(1.0-acos);
		u[0][1]=awx*awy*(1.0-acos)-awz*asin;
		u[0][2]=awx*awz*(1.0-acos)+awy*asin;
		u[1][0]=awx*awy*(1.0-acos)+awz*asin;
		u[1][1]=acos+awy*awy*(1.0-acos);
		u[1][2]=awy*awz*(1.0-acos)-awx*asin;
		u[2][0]=awx*awz*(1.0-acos)-awy*asin;
		u[2][1]=awy*awz*(1.0-acos)+awx*asin;
		u[2][2]=acos+awz*awz*(1.0-acos);			
		// Rotation points
//		int rand_point = rand()%(pnumx);
//		int rand_point = rand()%(pnum);
//		axyz[0]=pointsx[rand_point].ca[0];
//		axyz[1]=pointsx[rand_point].ca[1];
//		axyz[2]=pointsx[rand_point].ca[2];
		int rand_point = rand()%(pnumx);
		axyz[0]=pointsABx[rand_point][0];
		axyz[1]=pointsABx[rand_point][1];
		axyz[2]=pointsABx[rand_point][2];

		theDensityMap.calcRhoCx(pose,4.0,rhoC,rhoMask);
		old_dE=theDensityMap.getRSCCX(rhoC,rhoMask);
		cout<<"old_dE: "<<old_dE<<endl;			
/*	  	for(int j=0;j<rand_point;j++)
	  	{
	  		old_dE1 = old_dE1 + Distance_point(pointsABx[j],pointsABy[j]);
	  	}

	  	old_dE2 =0.0;
	  	for(int j=rand_point;j<pnumx;j++)
	  	{
	  		old_dE2 = old_dE2 + Distance_point(pointsABx[j],pointsABy[j]);
	  	}		  			
	  	old_dE = old_dE1 + old_dE2; */
//	  	cout<<"old_dE,vwd: "<<old_dE<<" "<<GetVdwEg()<<endl;
		if(rand_point>=int(float(pnumx)/2.0))  // fix N terminal and rotate C terminal
		{
//			vector<vector<float> >(0,vector<float>()).swap(fin_mat);
			fin_mat.clear();
			for(int j=rand_point;j<pnumx;j++)
			{
				fin_mat.push_back(pointsABx[j]);
			}
			int tmp3 =0;
	  		for(int j=rand_point;j<pnumx;j++)
	  		{
				fin_mat[tmp3][0]=t1+axyz[0]+(pointsABx[j][0]-axyz[0])*u[0][0]+(pointsABx[j][1]-axyz[1])*u[0][1]+(pointsABx[j][2]-axyz[2])*u[0][2];
				fin_mat[tmp3][1]=t2+axyz[1]+(pointsABx[j][0]-axyz[0])*u[1][0]+(pointsABx[j][1]-axyz[1])*u[1][1]+(pointsABx[j][2]-axyz[2])*u[1][2];
				fin_mat[tmp3][2]=t3+axyz[2]+(pointsABx[j][0]-axyz[0])*u[2][0]+(pointsABx[j][1]-axyz[1])*u[2][1]+(pointsABx[j][2]-axyz[2])*u[2][2];	  			
	  			tmp3= tmp3 +1 ;
	  		}
			new_dE2 = 0.0;
			int tmp1=0;
	  		for(int j=rand_point;j<pnumx;j++)
	  		{
	  			new_dE2 = new_dE2 + Distance_point(fin_mat[tmp1],pointsABy[j]);
	  			tmp1 = tmp1 + 1;
	  		}
	  		new_dE = old_dE1 + new_dE2;
//	  		cout<<"new_dE,vwd: "<<old_dE<<" "<<vwd<<endl;
	  		dE = new_dE - old_dE;
	  		if(new_dE<old_dE)
	  		{
	  			cout<<"EEEE"<<endl;	
	  			for(int j=0;j<rand_point;j++)
	  			{
		  			// translate CN CA CC CO
		  			pointsABx[j][0] = t1 + pointsABx[j][0];
		  			pointsABx[j][1] = t2 + pointsABx[j][1];
		  			pointsABx[j][2] = t3 + pointsABx[j][2];		  				
	  			}
	  			int tmp2 =0;
	  			for(int j=rand_point;j<pnumx;j++)
	  			{
		    		// rotate N
			        fin_mat[tmp2][0]=t1+axyz[0]+(pointsABx[j][0]-axyz[0])*u[0][0]+(pointsABx[j][1]-axyz[1])*u[0][1]+(pointsABx[j][2]-axyz[2])*u[0][2];
			        fin_mat[tmp2][1]=t2+axyz[1]+(pointsABx[j][0]-axyz[0])*u[1][0]+(pointsABx[j][1]-axyz[1])*u[1][1]+(pointsABx[j][2]-axyz[2])*u[1][2];
			        fin_mat[tmp2][2]=t3+axyz[2]+(pointsABx[j][0]-axyz[0])*u[2][0]+(pointsABx[j][1]-axyz[1])*u[2][1]+(pointsABx[j][2]-axyz[2])*u[2][2];
			        tmp2 = tmp2 +1;
	  			}
	  			tmp2 = 0;
				for(int j=rand_point;j<pnumx;j++)
				{
				    pointsABx[j] = fin_mat[tmp2];
				    tmp2 = tmp2 +1;
				}		  			
	  		}
	  		else
	  		{
//	  			cout<<"EEER"<<endl;	
	  			float tmpx=rand()/double(RAND_MAX);
	  			float mc_v = exp(-dE/(KT));
	  			if(tmpx< mc_v)
	  			{
		  			for(int j=0;j<rand_point;j++)
		  			{
			  			// translate CN CA CC CO
			  			pointsABx[j][0] = t1 + pointsABx[j][0];
			  			pointsABx[j][1] = t2 + pointsABx[j][1];
			  			pointsABx[j][2] = t3 + pointsABx[j][2];
	  				
		  			}
		  			int tmp2 =0;
		  			for(int j=rand_point;j<pnumx;j++)
		  			{
			    		// rotate N
				        fin_mat[tmp2][0]=t1+axyz[0]+(pointsABx[j][0]-axyz[0])*u[0][0]+(pointsABx[j][1]-axyz[1])*u[0][1]+(pointsABx[j][2]-axyz[2])*u[0][2];
				        fin_mat[tmp2][1]=t2+axyz[1]+(pointsABx[j][0]-axyz[0])*u[1][0]+(pointsABx[j][1]-axyz[1])*u[1][1]+(pointsABx[j][2]-axyz[2])*u[1][2];
				        fin_mat[tmp2][2]=t3+axyz[2]+(pointsABx[j][0]-axyz[0])*u[2][0]+(pointsABx[j][1]-axyz[1])*u[2][1]+(pointsABx[j][2]-axyz[2])*u[2][2];
				        tmp2 = tmp2 +1;	  				
		  			}
		  			tmp2 = 0;
					for(int j=rand_point;j<pnumx;j++)
					{
					    pointsABx[j] = fin_mat[tmp2];
					    tmp2 = tmp2 +1;
					}	
	  			}
	  		}					
		}
		else
		{
//			vector<Rot>().swap(fin_mat);
			fin_mat.clear();
			for(int j=0;j<rand_point;j++)
			{
				fin_mat.push_back(pointsABx[j]);
			}
			int tmp3 =0;
	  		for(int j=0;j<rand_point;j++)
	  		{
				fin_mat[tmp3][0]=t1+axyz[0]+(pointsABx[j][0]-axyz[0])*u[0][0]+(pointsABx[j][1]-axyz[1])*u[0][1]+(pointsABx[j][2]-axyz[2])*u[0][2];
				fin_mat[tmp3][1]=t2+axyz[1]+(pointsABx[j][0]-axyz[0])*u[1][0]+(pointsABx[j][1]-axyz[1])*u[1][1]+(pointsABx[j][2]-axyz[2])*u[1][2];
				fin_mat[tmp3][2]=t3+axyz[2]+(pointsABx[j][0]-axyz[0])*u[2][0]+(pointsABx[j][1]-axyz[1])*u[2][1]+(pointsABx[j][2]-axyz[2])*u[2][2];	  			
	  			tmp3 =tmp3 +1;
	  		}
			new_dE1 = 0.0;
			int tmp1=0;
	  		for(int j=0;j<rand_point;j++)
	  		{
	  			new_dE1 = new_dE1 + Distance_point(fin_mat[tmp1],pointsABy[j]);
	  			tmp1 = tmp1 + 1;
	  		}
	  		new_dE = new_dE1 + old_dE2;
	  		dE = new_dE - old_dE;
	  		if(new_dE<old_dE)
	  		{
	  			cout<<"EEEE"<<endl;	
	  			for(int j=rand_point;j<pnumx;j++)
	  			{
		  			// translate CN CA CC CO
		  			pointsABx[j][0] = t1 + pointsABx[j][0];
		  			pointsABx[j][1] = t2 + pointsABx[j][1];
		  			pointsABx[j][2] = t3 + pointsABx[j][2];		  				
	  			}
	  			int tmp2 =0;
	  			for(int j=0;j<rand_point;j++)
	  			{
		    		// rotate N
			        fin_mat[tmp2][0]=t1+axyz[0]+(pointsABx[j][0]-axyz[0])*u[0][0]+(pointsABx[j][1]-axyz[1])*u[0][1]+(pointsABx[j][2]-axyz[2])*u[0][2];
			        fin_mat[tmp2][1]=t2+axyz[1]+(pointsABx[j][0]-axyz[0])*u[1][0]+(pointsABx[j][1]-axyz[1])*u[1][1]+(pointsABx[j][2]-axyz[2])*u[1][2];
			        fin_mat[tmp2][2]=t3+axyz[2]+(pointsABx[j][0]-axyz[0])*u[2][0]+(pointsABx[j][1]-axyz[1])*u[2][1]+(pointsABx[j][2]-axyz[2])*u[2][2];
			        tmp2 = tmp2 +1;	  				
	  			}
	  			tmp2 = 0;
				for(int j=0;j<rand_point;j++)
				{
				    pointsABx[j] = fin_mat[tmp2];
				    tmp2 = tmp2 +1;
				}		  			
	  		}
	  		else
	  		{
	  			float tmpx=rand()/double(RAND_MAX);
	  			float mc_v = exp(-dE/(KT));
	  			if(tmpx< mc_v)
	  			{
		  			for(int j=rand_point;j<pnumx;j++)
		  			{
			  			// translate CN CA CC CO
			  			pointsABx[j][0] = t1 + pointsABx[j][0];
			  			pointsABx[j][1] = t2 + pointsABx[j][1];
			  			pointsABx[j][2] = t3 + pointsABx[j][2];		  				
		  			}
		  			int tmp2 =0;
		  			for(int j=0;j<rand_point;j++)
		  			{
					    		// rotate N
						fin_mat[tmp2][0]=t1+axyz[0]+(pointsABx[j][0]-axyz[0])*u[0][0]+(pointsABx[j][1]-axyz[1])*u[0][1]+(pointsABx[j][2]-axyz[2])*u[0][2];
						fin_mat[tmp2][1]=t2+axyz[1]+(pointsABx[j][0]-axyz[0])*u[1][0]+(pointsABx[j][1]-axyz[1])*u[1][1]+(pointsABx[j][2]-axyz[2])*u[1][2];
						fin_mat[tmp2][2]=t3+axyz[2]+(pointsABx[j][0]-axyz[0])*u[2][0]+(pointsABx[j][1]-axyz[1])*u[2][1]+(pointsABx[j][2]-axyz[2])*u[2][2];
						tmp2 = tmp2 +1;	  				
		  			}
		  			tmp2 = 0;
					for(int j=0;j<rand_point;j++)
					{
					    pointsABx[j] = fin_mat[tmp2];
					    tmp2 = tmp2 +1;
					}
	  			}
	  		}						  					
		}

	}
	cout<< "EEEEEEE"<<endl;
//		vector<Rot> tmp_mat;
	
		vector<vector<float> > tmp_mat;
		int randp= 0;
		for(int tttt=0;tttt<4000;tttt++)
		{ 
			vector<float> min_dE(pnumx,0.0);
			float thrshd = 0.001;
			float max_ca = thrshd;
			float min_ca = max_ca;
			int max_ca_index = -1,min_ca_index=-1;
			for(int i=0;i<pnumx;i++)
			{
				min_dE[i] = Distance_point(pointsABx[i],pointsABy[i]);//pointsx[i].ca[0] - pointsy[i].ca[0] 
				if(min_dE[i]>max_ca)
				{
					max_ca = min_dE[i];
					max_ca_index = i;
				}
				if(min_dE[i]<min_ca)
				{
					min_ca = min_dE[i];
					min_ca_index = i;
				}
		//		cout<< min_dE[i]<<" ";		
			}
			cout<< "pnumx: "<<pnumx<<endl; 	
			cout<<"max_ca,min_ca: "<<max_ca<<" "<<min_ca<<endl;
			cout<<"max_ca_index,min_ca_index: "<<max_ca_index<<" "<<min_ca_index<<endl; 
			if(max_ca_index<0) randp =rand()%(pnumx);
			else randp = max_ca_index;
		  	for(int i=0;i<1000;i++)
		  	{
			//	int rand_a = int(RandomDoubleX(0.0,float(coord_num-1)));
		  	//	vector<Rot>().swap(tmp_mat);
		  		tmp_mat.clear();
//	  			int randp= rand()%(pnumx);
		  		int rand_a=randp;
		  		int rand_b=randp;
		  		int randt = rand()%(9) + 1;
		  	//	int randt = 2;
		  		if(randt==0) continue;
		  		for(int j=0;j<randt;j++)
		  		{
		  			if((randp+j)>=(pnumx-1))
		  			{
		  				rand_b = randp + j;
		  				break;
		  			}
		  			rand_b = randp + j;
		  		}
		  		for(int j=0;j<randt;j++)
		  		{
		  			if((randp-j)<=0)
		  			{
		  				rand_a = randp - j;
		  				break;
		  			}
		  			rand_a = randp-j;
		  		}
		 // 		cout<<"TTTTTTT"<<endl;
		//		cout<< " rand: "<< rand_a<<" "<< rand_b<<endl;;
				if( rand_a == rand_b) continue;
				if(rand_a>rand_b)
				{
					int tmp_rand = rand_b;
					rand_b = rand_a;
					rand_a = tmp_rand;
				}
		//		cout<< " rand a b "<< rand_a<<" "<<rand_b<<endl;
				old_dE = 0.0;
				for(int j=rand_a;j<=rand_b;j++)
				{
					tmp_mat.push_back(pointsABx[j]);
					old_dE = old_dE + Distance_point(pointsABx[j],pointsABy[j]);
				}

				// Rotate angle
				ang=1.0;
				angle_rotate = (2.0*RandomDoubleX(0.0,1.0)-1.0)*ang; // rotate angle
			//	GroupRotationx(tmp_mat[0].ca,tmp_mat[tmp_mat.size()-1].ca,angle_rotate,tmp_mat);
			//	GroupRotationy(tmp_mat[0].ca,tmp_mat[tmp_mat.size()-1].ca,angle_rotate,tmp_mat);
				GroupRotation(tmp_mat[0],tmp_mat[tmp_mat.size()-1],angle_rotate,tmp_mat);
				int randxx = rand()%(rand_b-rand_a+1)+rand_a;
				int inta = randxx%3;
				if(inta == 0)
				{
					vector<float> tmpcc(3,0.0);
					int randyy=randxx-rand_a;  // randxx and randyy is the same point
					float deta0 = (RandomDoubleX(0.0,1.0)*2.0-1.0)*0.01+0.0;
					float deta1 = (RandomDoubleX(0.0,1.0)*2.0-1.0)*0.01+0.0;
					float deta2 = (RandomDoubleX(0.0,1.0)*2.0-1.0)*0.01+0.0;
					tmpcc[0] = tmp_mat[randyy][0] + deta0;
					tmpcc[1] = tmp_mat[randyy][1] + deta1;
					tmpcc[2] = tmp_mat[randyy][2] + deta2;
					float dis_s0=0.0,dis_s1=0.0;
					if(randyy+1<tmp_mat.size()) dis_s0 = Distance_point(tmpcc,tmp_mat[randyy+1]);
					if(randyy==tmp_mat.size()) dis_s0 = Distance_point(tmpcc,pointsABx[randxx+1]); // because inta =0, so randxx can't be the last atom
					if(randyy-1>=0) dis_s1 = Distance_point(tmpcc,tmp_mat[randyy-1]);
					if(randxx==0) 
					{
						dis_s1 = 1.33000;
					} 
					else
					{
						if(randyy==0) dis_s1 = Distance_point(tmpcc,pointsABx[randxx-1]);  // the situation of boundary for N terminal
					}
			//		if(dis_s0>1.44906 && dis_s0<1.47186 && dis_s1<(1.33883) && dis_s1>(1.32169)) tmp_mat[randyy] = tmpcc;
					if(dis_s0>1.43001 && dis_s0<1.49999 && dis_s1<(1.33999) && dis_s1>(1.31001)) tmp_mat[randyy] = tmpcc;
				}
				if(inta == 1)
				{
					vector<float> tmpcc(3,0.0);
					int randyy=randxx-rand_a;  // randxx and randyy is the same point
					float deta0 = (RandomDoubleX(0.0,1.0)*2.0-1.0)*0.01+0.0;
					float deta1 = (RandomDoubleX(0.0,1.0)*2.0-1.0)*0.01+0.0;
					float deta2 = (RandomDoubleX(0.0,1.0)*2.0-1.0)*0.01+0.0;
					tmpcc[0] = tmp_mat[randyy][0] + deta0;
					tmpcc[1] = tmp_mat[randyy][1] + deta1;
					tmpcc[2] = tmp_mat[randyy][2] + deta2;
					float dis_s0=0.0,dis_s1=0.0;
					if(randyy+1<tmp_mat.size()) dis_s0 = Distance_point(tmpcc,tmp_mat[randyy+1]);
					if(randyy==tmp_mat.size()) dis_s0 = Distance_point(tmpcc,pointsABx[randxx+1]); // because inta =1, so randxx can't be the last atom
					if(randyy-1>=0) dis_s1 = Distance_point(tmpcc,tmp_mat[randyy-1]);
					if(randyy==0) dis_s1 = Distance_point(tmpcc,pointsABx[randxx-1]);  // because inta =1, so randxx can't be the first atom
			//		if(dis_s0>1.51328 && dis_s0<1.53698 && dis_s1<(1.47186) && dis_s1>(1.44906)) tmp_mat[randyy] = tmpcc;
					if(dis_s0>1.51001 && dis_s0<1.55999 && dis_s1<(1.49999) && dis_s1>(1.43001)) tmp_mat[randyy] = tmpcc;					
				}
				if(inta == 2)
				{
					vector<float> tmpcc(3,0.0);
					int randyy=randxx-rand_a; // randxx and randyy is the same point
					float deta0 = (RandomDoubleX(0.0,1.0)*2.0-1.0)*0.008+0.0;
					float deta1 = (RandomDoubleX(0.0,1.0)*2.0-1.0)*0.008+0.0;
					float deta2 = (RandomDoubleX(0.0,1.0)*2.0-1.0)*0.008+0.0;
					tmpcc[0] = tmp_mat[randyy][0] + deta0;
					tmpcc[1] = tmp_mat[randyy][1] + deta1;
					tmpcc[2] = tmp_mat[randyy][2] + deta2;
					float dis_s0=0.0,dis_s1=0.0;
					if(randyy+1<tmp_mat.size()) dis_s0 = Distance_point(tmpcc,tmp_mat[randyy+1]);
					if(randxx==(pnumx-1)) 
					{
						dis_s0 = 1.32469;
					} 
					else
					{
						if(randyy==tmp_mat.size()) dis_s0 = Distance_point(tmpcc,pointsABx[randxx+1]);  // the situation of boundary for C terminal
					}					
					if(randyy-1>=0) dis_s1 = Distance_point(tmpcc,tmp_mat[randyy-1]);
					if(randyy==0) dis_s1 = Distance_point(tmpcc,pointsABx[randxx-1]);  // because inta =2, so randxx can't be the first atom
			//		if(dis_s0>1.32169 && dis_s0<1.33883 && dis_s1<(1.53698) && dis_s1>(1.51328)) tmp_mat[randyy] = tmpcc;
					if(dis_s0>1.31001 && dis_s0<1.33999 && dis_s1<(1.55999) && dis_s1>(1.51001)) tmp_mat[randyy] = tmpcc;					
				}												
				new_dE = 0.0;
				int tmp_tm=0;
				for(int j=rand_a;j<=rand_b;j++)
				{
					new_dE = new_dE + Distance_point(tmp_mat[tmp_tm],pointsABy[j]);
					tmp_tm = tmp_tm +1;
				}
				dE = new_dE - old_dE;
				if(new_dE < old_dE)
				{
		//			cout<< tmp_count<<endl;
					int tmp_tmx=0;
			//		point_cn = points_pre3[rand_a].cn;
			//		point_cc = points_pre3[rand_b].cc;
			//		point_co = points_pre3[rand_b].co;
					for(int j=rand_a;j<=rand_b;j++)
					{
						pointsABx[j] = tmp_mat[tmp_tmx];
						tmp_tmx = tmp_tmx +1;
					}
			//		points_pre3[rand_a].cn = point_cn;
			//		points_pre3[rand_b].cc = point_cc;
			//		points_pre3[rand_b].co = point_co;
		//			tmp_count = tmp_count +1;	
			 //   	if(new_dE<fin_dE)
			 //   	{
			  //  		vector<Rot>().swap(fin_frag);
			 //   		fin_frag = points_pre3;
			 //   		fin_dE = new_dE;
			 //   	}
				}
				else
				{
				//	float tmpx=randomAB(0.0,1.0);
					float tmpx=rand()/double(RAND_MAX);
					float mc_v = exp(-dE/(KT));
					if(tmpx < mc_v)
					{
						int tmp_tmx=0;
				//		point_cn = points_pre3[rand_a].cn;
				//		point_cc = points_pre3[rand_b].cc;
				//		point_co = points_pre3[rand_b].co;
						for(int j=rand_a;j<=rand_b;j++)
						{
							pointsABx[j] = tmp_mat[tmp_tmx];
							tmp_tmx = tmp_tmx +1;
						}
					}		
				}
			//	tmp_mat.clear();
			//	vector<Rot>().swap(tmp_mat);
				tmp_mat.clear();		
		  	}
	  	}	
	}
//	vector<Rotm> CBp(pnum,vector<float>(3,0.0)),CBpx(pnum,vector<float>(3,0.0));
	vector<Rotm> CBp(pnum),CBpy(pnum);
	for(int i=0;i<pnum;i++)
	{
		cout<<pointsrenm[i]<<endl;
		CBp[i] =CalculateRotm(pointsrenm[i],pointsABx[3*i],pointsABx[3*i+1],pointsABx[3*i+2]);
//		CBpy[i]=CalculateRotm(pointsrenmy[i],pointsABy[3*i],pointsABy[3*i+1],pointsABy[3*i+2]);
		//for CO atom
		cout<<"0X"<<endl;
		for(int j=0;j<2000;j++)
		{
			float dis_old = Distance_point(CBp[i].Ocd,pointsy[i].co);
			float ang=0.5;
			float angle_rotate = (2.0*RandomDoubleX(0.0,1.0)-1.0)*ang; // rotate angle
			vector<float> tmp_matx=CBp[i].Ocd;
			GroupRotationt(pointsABx[3*i+1],pointsABx[3*i+2],angle_rotate,tmp_matx);
			float dis_new = Distance_point(tmp_matx,pointsy[i].co);
			float dE = dis_new- dis_old;
			if(dis_new<dis_old)
			{
				CBp[i].Ocd = tmp_matx;
			}
			else
			{
				float tmpx=rand()/double(RAND_MAX);
				float mc_v = exp(-dE/(KT));
				if(tmpx < mc_v)
				{
					CBp[i].Ocd = tmp_matx;
//					cout<<"111"<<endl;
				}
			}			
		}		
		cout<<"X"<<endl;
		// for side chain atom
		// for CB
		if(strcmp(pointsrenm[i].c_str(),"G")==0 || strcmp(pointsrenm[i].c_str(),"P")==0) continue;
		// for CA-CB
		if(strcmp(pointsrenm[i].c_str(),"A")==0) continue;
		if(strcmp(pointsrenm[i].c_str(),"C")==0)
		{
			if(pointsy[i].res_atm.size()< 2) continue;
			if(pointsy[i].sg.size()<2) continue;
			for(int j=0;j<2000;j++)
			{
				float dis_old = Distance_point(CBp[i].Rotcd[1],pointsy[i].sg);
				float ang=0.5;
				float angle_rotate = (2.0*RandomDoubleX(0.0,1.0)-1.0)*ang; // rotate angle
				vector<vector<float> > tmp_matx = CBp[i].Rotcd;
				GroupRotation(pointsABx[3*i+1],CBp[i].Rotcd[0],angle_rotate,tmp_matx);
				float dis_new = Distance_point(tmp_matx[1],pointsy[i].sg);
				float dE = dis_new- dis_old;
				if(dis_new<dis_old)
				{
					for(int t=0;t<CBp[i].Rotcd.size();t++)
					{
						CBp[i].Rotcd[t]=tmp_matx[t];
					}
				}
				else
				{
					float tmpx=rand()/double(RAND_MAX);
					float mc_v = exp(-dE/(KT));
					if(tmpx < mc_v)
					{
						for(int t=0;t<CBp[i].Rotcd.size();t++)
						{
							CBp[i].Rotcd[t]=tmp_matx[t];
						}					
					}				
				}
			}			
		}
		if(strcmp(pointsrenm[i].c_str(),"C")==0) continue;
		if(strcmp(pointsrenm[i].c_str(),"S")==0)
		{
			if(pointsy[i].res_atm.size()< 2) continue;
			if(pointsy[i].og.size()<2) continue;
			for(int j=0;j<2000;j++)
			{
				float dis_old = Distance_point(CBp[i].Rotcd[1],pointsy[i].og);
				float ang=0.5;
				float angle_rotate = (2.0*RandomDoubleX(0.0,1.0)-1.0)*ang; // rotate angle
				vector<vector<float> > tmp_matx = CBp[i].Rotcd;
				GroupRotation(pointsABx[3*i+1],CBp[i].Rotcd[0],angle_rotate,tmp_matx);
				float dis_new = Distance_point(tmp_matx[1],pointsy[i].og);
				float dE = dis_new- dis_old;
				if(dis_new<dis_old)
				{
					for(int t=0;t<CBp[i].Rotcd.size();t++)
					{
						CBp[i].Rotcd[t]=tmp_matx[t];
					}
				}
				else
				{
					float tmpx=rand()/double(RAND_MAX);
					float mc_v = exp(-dE/(KT));
					if(tmpx < mc_v)
					{
						for(int t=0;t<CBp[i].Rotcd.size();t++)
						{
							CBp[i].Rotcd[t]=tmp_matx[t];
						}					
					}				
				}
			}				
		}
		cout<<"Xx"<<endl;
		if(strcmp(pointsrenm[i].c_str(),"T")==0)
		{
			if(pointsy[i].res_atm.size()< 3) continue;
			if(pointsy[i].cg2.size()<2) continue;
			for(int j=0;j<2000;j++)
			{
				float dis_old = Distance_point(CBp[i].Rotcd[2],pointsy[i].cg2);
				float ang=0.5;
				float angle_rotate = (2.0*RandomDoubleX(0.0,1.0)-1.0)*ang; // rotate angle
				vector<vector<float> > tmp_matx = CBp[i].Rotcd;
				GroupRotation(pointsABx[3*i+1],CBp[i].Rotcd[0],angle_rotate,tmp_matx);
				float dis_new = Distance_point(tmp_matx[2],pointsy[i].cg2);
				float dE = dis_new- dis_old;
				if(dis_new<dis_old)
				{
					for(int t=0;t<CBp[i].Rotcd.size();t++)
					{
						CBp[i].Rotcd[t]=tmp_matx[t];
					}
				}
				else
				{
					float tmpx=rand()/double(RAND_MAX);
					float mc_v = exp(-dE/(KT));
					if(tmpx < mc_v)
					{
						for(int t=0;t<CBp[i].Rotcd.size();t++)
						{
							CBp[i].Rotcd[t]=tmp_matx[t];
						}					
					}				
				}
			}				
		}	
		cout<< "10XX"<<endl;
		if(strcmp(pointsrenm[i].c_str(),"I")==0)
		{
			if(pointsy[i].res_atm.size()< 2) continue;
			if(pointsy[i].cg1.size()<2) continue;
			for(int j=0;j<2000;j++)
			{
				float dis_old = Distance_point(CBp[i].Rotcd[1],pointsy[i].cg1);
				float ang=0.5;
				float angle_rotate = (2.0*RandomDoubleX(0.0,1.0)-1.0)*ang; // rotate angle
				vector<vector<float> > tmp_matx = CBp[i].Rotcd;
				GroupRotation(pointsABx[3*i+1],CBp[i].Rotcd[0],angle_rotate,tmp_matx);
				float dis_new = Distance_point(tmp_matx[1],pointsy[i].cg1);
				float dE = dis_new- dis_old;
				if(dis_new<dis_old)
				{
					for(int t=0;t<CBp[i].Rotcd.size();t++)
					{
						CBp[i].Rotcd[t]=tmp_matx[t];
					}
				}
				else
				{
					float tmpx=rand()/double(RAND_MAX);
					float mc_v = exp(-dE/(KT));
					if(tmpx < mc_v)
					{
						for(int t=0;t<CBp[i].Rotcd.size();t++)
						{
							CBp[i].Rotcd[t]=tmp_matx[t];
						}					
					}				
				}
			}
			if(pointsy[i].res_atm.size()< 4) continue;
			if(pointsy[i].cd1.size()<2) continue;
			for(int j=0;j<2000;j++)
			{
				float dis_old = Distance_point(CBp[i].Rotcd[3],pointsy[i].cd1);
				float ang=0.5;
				float angle_rotate = (2.0*RandomDoubleX(0.0,1.0)-1.0)*ang; // rotate angle
				vector<vector<float> > tmp_matx = CBp[i].Rotcd;
				GroupRotation(CBp[i].Rotcd[0],CBp[i].Rotcd[1],angle_rotate,tmp_matx);
				float dis_new = Distance_point(tmp_matx[3],pointsy[i].cd1);
				float dE = dis_new- dis_old;
				if(dis_new<dis_old)
				{
					for(int t=0;t<CBp[i].Rotcd.size();t++)
					{
						CBp[i].Rotcd[t]=tmp_matx[t];
					}
				}
				else
				{
					float tmpx=rand()/double(RAND_MAX);
					float mc_v = exp(-dE/(KT));
					if(tmpx < mc_v)
					{
						for(int t=0;t<CBp[i].Rotcd.size();t++)
						{
							CBp[i].Rotcd[t]=tmp_matx[t];
						}					
					}				
				}
			}							
		}
		cout<< "9XX"<<endl;
		if(strcmp(pointsrenm[i].c_str(),"V")==0)
		{
			if(pointsy[i].res_atm.size()< 2) continue;
			if(pointsy[i].cg1.size()<2) continue;
			for(int j=0;j<2000;j++)
			{
				float dis_old;
				if(pointsy[i].cg1.size()>0) 
				{
					dis_old = Distance_point(CBp[i].Rotcd[1],pointsy[i].cg1);
				}
				else 
				{
					break;
				}
				float ang=0.5;
				float angle_rotate = (2.0*RandomDoubleX(0.0,1.0)-1.0)*ang; // rotate angle
				vector<vector<float> > tmp_matx = CBp[i].Rotcd;
				GroupRotation(pointsABx[3*i+1],CBp[i].Rotcd[0],angle_rotate,tmp_matx);
				float dis_new = Distance_point(tmp_matx[1],pointsy[i].cg1);
				float dE = dis_new- dis_old;
				if(dis_new<dis_old)
				{
					for(int t=0;t<CBp[i].Rotcd.size();t++)
					{
						CBp[i].Rotcd[t]=tmp_matx[t];
					}
				}
				else
				{
					float tmpx=rand()/double(RAND_MAX);
					float mc_v = exp(-dE/(KT));
					if(tmpx < mc_v)
					{
						for(int t=0;t<CBp[i].Rotcd.size();t++)
						{
							CBp[i].Rotcd[t]=tmp_matx[t];
						}					
					}				
				}
			}				
		}		
		cout<< "8XX"<<endl;			
		if(strcmp(pointsrenm[i].c_str(),"S")==0 || strcmp(pointsrenm[i].c_str(),"T")==0 || strcmp(pointsrenm[i].c_str(),"I")==0 || strcmp(pointsrenm[i].c_str(),"V")==0) continue;
		if(pointsy[i].res_atm.size()< 2) continue;
		if(pointsy[i].cg.size()<2) continue;
		for(int j=0;j<2000;j++)
		{
			float dis_old = Distance_point(CBp[i].Rotcd[1],pointsy[i].cg);
			float ang=0.5;
			float angle_rotate = (2.0*RandomDoubleX(0.0,1.0)-1.0)*ang; // rotate angle
			vector<vector<float> > tmp_matx = CBp[i].Rotcd;
			GroupRotation(pointsABx[3*i+1],CBp[i].Rotcd[0],angle_rotate,tmp_matx);
			float dis_new = Distance_point(tmp_matx[1],pointsy[i].cg);
			float dE = dis_new- dis_old;
			if(dis_new<dis_old)
			{
				for(int t=0;t<CBp[i].Rotcd.size();t++)
				{
					CBp[i].Rotcd[t]=tmp_matx[t];
				}
			}
			else
			{
				float tmpx=rand()/double(RAND_MAX);
				float mc_v = exp(-dE/(KT));
				if(tmpx < mc_v)
				{
					for(int t=0;t<CBp[i].Rotcd.size();t++)
					{
						CBp[i].Rotcd[t]=tmp_matx[t];
					}					
				}				
			}
		}
		cout<< "7XX"<<endl;
//		if(CBp[i].Rotnm.size()< 2) continue;
		if(strcmp(pointsrenm[i].c_str(),"R")==0 || strcmp(pointsrenm[i].c_str(),"Q")==0  || strcmp(pointsrenm[i].c_str(),"K")==0 || strcmp(pointsrenm[i].c_str(),"E")==0 )
		{
			if(pointsy[i].res_atm.size()< 3) continue;
			if(pointsy[i].cd.size()<2) continue;
			for(int j=0;j<2000;j++)
			{
				float dis_old = Distance_point(CBp[i].Rotcd[2],pointsy[i].cd);
				float ang=0.5;
				float angle_rotate = (2.0*RandomDoubleX(0.0,1.0)-1.0)*ang; // rotate angle
				vector<vector<float> > tmp_matx = CBp[i].Rotcd;
				GroupRotation(CBp[i].Rotcd[0],CBp[i].Rotcd[1],angle_rotate,tmp_matx);
				float dis_new = Distance_point(tmp_matx[2],pointsy[i].cd);
				float dE = dis_new- dis_old;
				if(dis_new<dis_old)
				{
					for(int t=0;t<CBp[i].Rotcd.size();t++)
					{
						CBp[i].Rotcd[t]=tmp_matx[t];
					}
				}
				else
				{
					float tmpx=rand()/double(RAND_MAX);
					float mc_v = exp(-dE/(KT));
					if(tmpx < mc_v)
					{
						for(int t=0;t<CBp[i].Rotcd.size();t++)
						{
							CBp[i].Rotcd[t]=tmp_matx[t];
						}					
					}				
				}
			}			
		}
		cout<< "2XX"<<endl;
		if(strcmp(pointsrenm[i].c_str(),"Y")==0 || strcmp(pointsrenm[i].c_str(),"W")==0 || strcmp(pointsrenm[i].c_str(),"L")==0 || strcmp(pointsrenm[i].c_str(),"F")==0 )
		{
			if(pointsy[i].res_atm.size()< 3) continue;
			if(pointsy[i].cd1.size()<2) continue;
			for(int j=0;j<2500;j++)
			{
				float dis_old = Distance_point(CBp[i].Rotcd[2],pointsy[i].cd1);
				float ang=0.5;
				float angle_rotate = (2.0*RandomDoubleX(0.0,1.0)-1.0)*ang; // rotate angle
				vector<vector<float> > tmp_matx = CBp[i].Rotcd;
				GroupRotation(CBp[i].Rotcd[0],CBp[i].Rotcd[1],angle_rotate,tmp_matx);
				float dis_new = Distance_point(tmp_matx[2],pointsy[i].cd1);
				float dE = dis_new- dis_old;
				if(dis_new<dis_old)
				{
					for(int t=0;t<CBp[i].Rotcd.size();t++)
					{
						CBp[i].Rotcd[t]=tmp_matx[t];
					}
				}
				else
				{
					float tmpx=rand()/double(RAND_MAX);
					float mc_v = exp(-dE/(KT));
					if(tmpx < mc_v)
					{
						for(int t=0;t<CBp[i].Rotcd.size();t++)
						{
							CBp[i].Rotcd[t]=tmp_matx[t];
						}					
					}				
				}
			}
		}
		cout<< "3XX"<<endl;
		if(strcmp(pointsrenm[i].c_str(),"N")==0)
		{
			if(pointsy[i].res_atm.size()< 4) continue;
			if(pointsy[i].nd2.size()<2) continue;
			for(int j=0;j<2000;j++)
			{
				float dis_old = Distance_point(CBp[i].Rotcd[3],pointsy[i].nd2);
				float ang=0.5;
				float angle_rotate = (2.0*RandomDoubleX(0.0,1.0)-1.0)*ang; // rotate angle
				vector<vector<float> > tmp_matx = CBp[i].Rotcd;
				GroupRotation(CBp[i].Rotcd[0],CBp[i].Rotcd[1],angle_rotate,tmp_matx);
				float dis_new = Distance_point(tmp_matx[3],pointsy[i].nd2);
				float dE = dis_new- dis_old;
				if(dis_new<dis_old)
				{
					for(int t=0;t<CBp[i].Rotcd.size();t++)
					{
						CBp[i].Rotcd[t]=tmp_matx[t];
					}
				}
				else
				{
					float tmpx=rand()/double(RAND_MAX);
					float mc_v = exp(-dE/(KT));
					if(tmpx < mc_v)
					{
						for(int t=0;t<CBp[i].Rotcd.size();t++)
						{
							CBp[i].Rotcd[t]=tmp_matx[t];
						}					
					}				
				}
			}
		}			
		cout<< "5XX"<<endl;			
		if(strcmp(pointsrenm[i].c_str(),"H")==0)
		{
			if(pointsy[i].res_atm.size()< 4) continue;
			if(pointsy[i].cd2.size()<2) continue;
			for(int j=0;j<2000;j++)
			{
				float dis_old = Distance_point(CBp[i].Rotcd[3],pointsy[i].cd2);
				float ang=0.5;
				float angle_rotate = (2.0*RandomDoubleX(0.0,1.0)-1.0)*ang; // rotate angle
				vector<vector<float> > tmp_matx = CBp[i].Rotcd;
				GroupRotation(CBp[i].Rotcd[0],CBp[i].Rotcd[1],angle_rotate,tmp_matx);
				float dis_new = Distance_point(tmp_matx[3],pointsy[i].cd2);
				float dE = dis_new- dis_old;
				if(dis_new<dis_old)
				{
					for(int t=0;t<CBp[i].Rotcd.size();t++)
					{
						CBp[i].Rotcd[t]=tmp_matx[t];
					}
				}
				else
				{
					float tmpx=rand()/double(RAND_MAX);
					float mc_v = exp(-dE/(KT));
					if(tmpx < mc_v)
					{
						for(int t=0;t<CBp[i].Rotcd.size();t++)
						{
							CBp[i].Rotcd[t]=tmp_matx[t];
						}					
					}				
				}
			}
		}
		cout<< "6XX"<<endl;
		if(strcmp(pointsrenm[i].c_str(),"D")==0)
		{
			if(pointsy[i].res_atm.size()< 4) continue;
			if(pointsy[i].od2.size()<2) continue;
			for(int j=0;j<2000;j++)
			{
				float dis_old = Distance_point(CBp[i].Rotcd[3],pointsy[i].od2);
				float ang=0.5;
				float angle_rotate = (2.0*RandomDoubleX(0.0,1.0)-1.0)*ang; // rotate angle
				vector<vector<float> > tmp_matx = CBp[i].Rotcd;
				GroupRotation(CBp[i].Rotcd[0],CBp[i].Rotcd[1],angle_rotate,tmp_matx);
				float dis_new = Distance_point(tmp_matx[3],pointsy[i].od2);
				float dE = dis_new- dis_old;
				if(dis_new<dis_old)
				{
					for(int t=0;t<CBp[i].Rotcd.size();t++)
					{
						CBp[i].Rotcd[t]=tmp_matx[t];
					}
				}
				else
				{
					float tmpx=rand()/double(RAND_MAX);
					float mc_v = exp(-dE/(KT));
					if(tmpx < mc_v)
					{
						for(int t=0;t<CBp[i].Rotcd.size();t++)
						{
							CBp[i].Rotcd[t]=tmp_matx[t];
						}					
					}				
				}
			}
		}
		cout<< "4XX"<<endl;
		if(strcmp(pointsrenm[i].c_str(),"M")==0)
		{
			if(pointsy[i].res_atm.size()< 3) continue;
			if(pointsy[i].sd.size()<2) continue;
			for(int j=0;j<2000;j++)
			{
				float dis_old = Distance_point(CBp[i].Rotcd[2],pointsy[i].sd);
				float ang=0.5;
				float angle_rotate = (2.0*RandomDoubleX(0.0,1.0)-1.0)*ang; // rotate angle
				vector<vector<float> > tmp_matx = CBp[i].Rotcd;
				GroupRotation(CBp[i].Rotcd[0],CBp[i].Rotcd[1],angle_rotate,tmp_matx);
				float dis_new = Distance_point(tmp_matx[2],pointsy[i].sd);
				float dE = dis_new- dis_old;
				if(dis_new<dis_old)
				{
					for(int t=0;t<CBp[i].Rotcd.size();t++)
					{
						CBp[i].Rotcd[t]=tmp_matx[t];
					}
				}
				else
				{
					float tmpx=rand()/double(RAND_MAX);
					float mc_v = exp(-dE/(KT));
					if(tmpx < mc_v)
					{
						for(int t=0;t<CBp[i].Rotcd.size();t++)
						{
							CBp[i].Rotcd[t]=tmp_matx[t];
						}					
					}				
				}
			}
			if(pointsy[i].res_atm.size()< 4) continue;
			if(pointsy[i].ce.size()<2) continue;
			for(int j=0;j<2000;j++)
			{
				float dis_old = Distance_point(CBp[i].Rotcd[3],pointsy[i].ce);
				float ang=0.5;
				float angle_rotate = (2.0*RandomDoubleX(0.0,1.0)-1.0)*ang; // rotate angle
				vector<vector<float> > tmp_matx = CBp[i].Rotcd;
				GroupRotationid(CBp[i].Rotcd[1],CBp[i].Rotcd[2],angle_rotate,tmp_matx,2);
				float dis_new = Distance_point(tmp_matx[3],pointsy[i].ce);
				float dE = dis_new- dis_old;
				if(dis_new<dis_old)
				{
					for(int t=0;t<CBp[i].Rotcd.size();t++)
					{
						CBp[i].Rotcd[t]=tmp_matx[t];
					}
				}
				else
				{
					float tmpx=rand()/double(RAND_MAX);
					float mc_v = exp(-dE/(KT));
					if(tmpx < mc_v)
					{
						for(int t=0;t<CBp[i].Rotcd.size();t++)
						{
							CBp[i].Rotcd[t]=tmp_matx[t];
						}					
					}				
				}
			}
		}
		cout<< "1XX"<<endl;
		if(strcmp(pointsrenm[i].c_str(),"R")==0)
		{
			if(pointsy[i].res_atm.size()< 4) continue;
			if(pointsy[i].ne.size()<2) continue;
			for(int j=0;j<2000;j++)
			{
				float dis_old = Distance_point(CBp[i].Rotcd[3],pointsy[i].ne);
				float ang=0.5;
				float angle_rotate = (2.0*RandomDoubleX(0.0,1.0)-1.0)*ang; // rotate angle
				vector<vector<float> > tmp_matx = CBp[i].Rotcd;
				GroupRotationid(CBp[i].Rotcd[1],CBp[i].Rotcd[2],angle_rotate,tmp_matx,2);
				float dis_new = Distance_point(tmp_matx[3],pointsy[i].ne);
				float dE = dis_new- dis_old;
				if(dis_new<dis_old)
				{
					for(int t=0;t<CBp[i].Rotcd.size();t++)
					{
						CBp[i].Rotcd[t]=tmp_matx[t];
					}
				}
				else
				{
					float tmpx=rand()/double(RAND_MAX);
					float mc_v = exp(-dE/(KT));
					if(tmpx < mc_v)
					{
						for(int t=0;t<CBp[i].Rotcd.size();t++)
						{
							CBp[i].Rotcd[t]=tmp_matx[t];
						}					
					}				
				}
			}
			if(pointsy[i].res_atm.size()< 5) continue;
			if(pointsy[i].cz.size()<2) continue;
			for(int j=0;j<2000;j++)
			{
				float dis_old = Distance_point(CBp[i].Rotcd[4],pointsy[i].cz);					
				float ang=0.5;
				float angle_rotate = (2.0*RandomDoubleX(0.0,1.0)-1.0)*ang; // rotate angle
				vector<vector<float> > tmp_matx = CBp[i].Rotcd;
				GroupRotationid(CBp[i].Rotcd[2],CBp[i].Rotcd[3],angle_rotate,tmp_matx,3);
				float dis_new = Distance_point(tmp_matx[4],pointsy[i].cz);
				float dE = dis_new- dis_old;
				if(dis_new<dis_old)
				{
					for(int t=0;t<CBp[i].Rotcd.size();t++)
					{
						CBp[i].Rotcd[t]=tmp_matx[t];
					}
				}
				else
				{
					float tmpx=rand()/double(RAND_MAX);
					float mc_v = exp(-dE/(KT));
					if(tmpx < mc_v)
					{
						for(int t=0;t<CBp[i].Rotcd.size();t++)
						{
							CBp[i].Rotcd[t]=tmp_matx[t];
						}					
					}				
				}
			}			
		}	
		if(strcmp(pointsrenm[i].c_str(),"Q")==0)
		{
			if(pointsy[i].res_atm.size()< 5) continue;
			if(pointsy[i].ne2.size()<2) continue;
			for(int j=0;j<2000;j++)
			{
				float dis_old = Distance_point(CBp[i].Rotcd[4],pointsy[i].ne2);				
				float ang=0.5;
				float angle_rotate = (2.0*RandomDoubleX(0.0,1.0)-1.0)*ang; // rotate angle
				vector<vector<float> > tmp_matx = CBp[i].Rotcd;
				GroupRotationid(CBp[i].Rotcd[1],CBp[i].Rotcd[2],angle_rotate,tmp_matx,2);
				float dis_new = Distance_point(tmp_matx[4],pointsy[i].ne2);
				float dE = dis_new- dis_old;
				if(dis_new<dis_old)
				{
					for(int t=0;t<CBp[i].Rotcd.size();t++)
					{
						CBp[i].Rotcd[t]=tmp_matx[t];
					}
				}
				else
				{
					float tmpx=rand()/double(RAND_MAX);
					float mc_v = exp(-dE/(KT));
					if(tmpx < mc_v)
					{
						for(int t=0;t<CBp[i].Rotcd.size();t++)
						{
							CBp[i].Rotcd[t]=tmp_matx[t];
						}					
					}				
				}
			}			
		}
		if(strcmp(pointsrenm[i].c_str(),"K")==0)
		{
			if(pointsy[i].res_atm.size()< 4) continue;
			if(pointsy[i].ce.size()<2) continue;
			for(int j=0;j<2000;j++)
			{
				float dis_old = Distance_point(CBp[i].Rotcd[3],pointsy[i].ce);
				float ang=0.5;
				float angle_rotate = (2.0*RandomDoubleX(0.0,1.0)-1.0)*ang; // rotate angle
				vector<vector<float> > tmp_matx = CBp[i].Rotcd;
				GroupRotationid(CBp[i].Rotcd[1],CBp[i].Rotcd[2],angle_rotate,tmp_matx,2);
				float dis_new = Distance_point(tmp_matx[3],pointsy[i].ce);
				float dE = dis_new- dis_old;
				if(dis_new<dis_old)
				{
					for(int t=0;t<CBp[i].Rotcd.size();t++)
					{
						CBp[i].Rotcd[t]=tmp_matx[t];
					}
				}
				else
				{
					float tmpx=rand()/double(RAND_MAX);
					float mc_v = exp(-dE/(KT));
					if(tmpx < mc_v)
					{
						for(int t=0;t<CBp[i].Rotcd.size();t++)
						{
							CBp[i].Rotcd[t]=tmp_matx[t];
						}					
					}				
				}
			}
			if(pointsy[i].res_atm.size()< 5) continue;
			if(pointsy[i].nz.size()<2) continue;
			for(int j=0;j<2000;j++)
			{
				float dis_old = Distance_point(CBp[i].Rotcd[4],pointsy[i].nz);
				float ang=0.5;
				float angle_rotate = (2.0*RandomDoubleX(0.0,1.0)-1.0)*ang; // rotate angle
				vector<vector<float> > tmp_matx = CBp[i].Rotcd;
				GroupRotationid(CBp[i].Rotcd[2],CBp[i].Rotcd[3],angle_rotate,tmp_matx,3);
				float dis_new = Distance_point(tmp_matx[4],pointsy[i].nz);
				float dE = dis_new- dis_old;
				if(dis_new<dis_old)
				{
					for(int t=0;t<CBp[i].Rotcd.size();t++)
					{
						CBp[i].Rotcd[t]=tmp_matx[t];
					}
				}
				else
				{
					float tmpx=rand()/double(RAND_MAX);
					float mc_v = exp(-dE/(KT));
					if(tmpx < mc_v)
					{
						for(int t=0;t<CBp[i].Rotcd.size();t++)
						{
							CBp[i].Rotcd[t]=tmp_matx[t];
						}					
					}				
				}
			}			
		}
		cout<< "XXXX"<<endl;
		if(strcmp(pointsrenm[i].c_str(),"E")==0)
		{
			if(pointsy[i].res_atm.size()< 5) continue;
			if(pointsy[i].oe2.size()<2) continue;
			for(int j=0;j<2000;j++)
			{
				float dis_old = Distance_point(CBp[i].Rotcd[4],pointsy[i].oe2);
				float ang=0.5;
				float angle_rotate = (2.0*RandomDoubleX(0.0,1.0)-1.0)*ang; // rotate angle
				vector<vector<float> > tmp_matx = CBp[i].Rotcd;
				GroupRotationid(CBp[i].Rotcd[1],CBp[i].Rotcd[2],angle_rotate,tmp_matx,2);
				float dis_new = Distance_point(tmp_matx[4],pointsy[i].oe2);
				float dE = dis_new- dis_old;
				if(dis_new<dis_old)
				{
					for(int t=0;t<CBp[i].Rotcd.size();t++)
					{
						CBp[i].Rotcd[t]=tmp_matx[t];
					}
				}
				else
				{
					float tmpx=rand()/double(RAND_MAX);
					float mc_v = exp(-dE/(KT));
					if(tmpx < mc_v)
					{
						for(int t=0;t<CBp[i].Rotcd.size();t++)
						{
							CBp[i].Rotcd[t]=tmp_matx[t];
						}					
					}				
				}
			}			
		}											
/*		//For CB-CG
		if(strcmp(pointsrenm[i].c_str(),"V")==0 || strcmp(pointsrenm[i].c_str(),"T")==0|| strcmp(pointsrenm[i].c_str(),"S")==0) continue;
		if(CBp[i].Rotnm.size()< 2) continue;
		for(int j=0;j<1000;j++)
		{
			float dis_old = Distance_point(CBp[i].Rotcd[2],pointsy[i].cd1);
			float ang=1.0;
			float angle_rotate = (2.0*RandomDoubleX(0.0,1.0)-1.0)*ang; // rotate angle
			vector<vector<float> > tmp_matx = CBp[i].Rotcd;
			GroupRotation(CBp[i].Rotcd[0],CBp[i].Rotcd[1],angle_rotate,tmp_matx);
			float dis_new = Distance_point(tmp_matx[2],CBpy[i].Rotcd[2]);
			float dE = dis_new- dis_old;
			if(dis_new<dis_old)
			{
				for(int t=0;t<CBp[i].Rotcd.size();t++)
				{
					CBp[i].Rotcd[t]=tmp_matx[t];
				}
			}
			else
			{
				float tmpx=rand()/double(RAND_MAX);
				float mc_v = exp(-dE/(KT));
				if(tmpx < mc_v)
				{
					for(int t=0;t<CBp[i].Rotcd.size();t++)
					{
						CBp[i].Rotcd[t]=tmp_matx[t];
					}					
				}				
			}
		}		*/
//		CBp[i]=calculateCoordinatesx(pointsABx[3*i],pointsABx[3*i+2],pointsABx[3*i+1],1.52,109.5,122.5);
/*		float NCA = Distance_point(pointsABx[3*i+0],pointsABx[3*i+1]);
		float NCB = sqrt(abs(1.54*1.54 + NCA*NCA -2.0*1.54*NCA*cos(3.1415926*109.5/180)));
		float angNCBCA=asin((180.0/3.1415926)*(sin(3.1415926*109.5/180.0))*NCA/NCB);
		float aj= 3.1415926*angNCBCA/180.0;
		float bj= 3.1415926*122.5/180.0 ;
		float Rj= NCB;
		CBp[i][0] = -(cos(aj))*pointsABx[3*i+0][0]-(sin(aj))*pointsABx[3*i+0][1]-Rj*cos(aj);
		CBp[i][1] = (sin(aj))*(cos(bj))*pointsABx[3*i+0][0]-(cos(aj))*(cos(bj))*pointsABx[3*i+0][1]-(sin(bj))*pointsABx[3*i+0][2]+Rj*(sin(aj))*cos(bj);
		CBp[i][2] = (sin(aj))*(sin(bj))*pointsABx[3*i+0][0]-(cos(aj))*(sin(bj))*pointsABx[3*i+0][1]+(cos(bj))*pointsABx[3*i+0][2]+Rj*(sin(aj))*sin(bj);		*/
//		CBp[0] = -(cos(3.1415926*angNCBCA/180.0))*pointsABx[3*i+0][0]-(sin(3.1415926*angNCBCA/180.0))*pointsABx[3*i+0][1]-NCB*cos(3.1415926*angNCBCA/180.0);
//		CBp[1] = (sin(3.1415926*angNCBCA/180.0))*(cos(3.1415926*122.5/180.0))*pointsABx[3*i+0][0]-(cos(3.1415926*angNCBCA/180.0))*(cos(3.1415926*122.5/180.0))*pointsABx[3*i+0][1]-(sin(3.1415926*122.5/180.0))*pointsABx[3*i+0][2]+NCB*(sin(3.1415926*angNCBCA/180.0))*cos(3.1415926*122.5/180.0);
//		CBp[2] = (sin(3.1415926*angNCBCA/180.0))*(sin(3.1415926*122.5/180.0))*pointsABx[3*i+0][0]-(cos(3.1415926*angNCBCA/180.0))*(sin(3.1415926*122.5/180.0))*pointsABx[3*i+0][1]+(cos(3.1415926*122.5/180.0))*pointsABx[3*i+0][2]+NCB*(sin(3.1415926*angNCBCA/180.0))*sin(3.1415926*122.5/180.0);
	/*	float xy1=0.0,xy2=0.0,xy3=0.0;
		float CAC= Distance_point(pointsABx[3*i+1],pointsx[3*i+2])
		vector<float> a1(3,0.0),a2(3,0.0);
		// C-CA
		a1[0] = pointsABx[3*i+2][0] - pointsABx[3*i+1][0];
		a1[1] = pointsABx[3*i+2][1] - pointsABx[3*i+1][1];
		a1[2] = pointsABx[3*i+2][2] - pointsABx[3*i+1][2];
		a1[0] = a1[0]/(CAC);	
		a1[1] = a1[1]/(CAC);
		a1[2] = a1[2]/(CAC);
		(xy1-pointsABx[3*i+1][0])*a1[0] + (xy2-pointsABx[3*i+1][1])*a1[1]+(xy3-pointsABx[3*i+1][2])*a1[2] =(1.54)*cos(3.1415926*109.5/180);  

		CBp[i]= pointsABx[3*i].cn;
		vector<float> Nvt(3,0.0);
		Nvt=CBp[i];
		GroupRotationt(pointsABx[3*i+1],pointsABx[3*i+2],120,Nvt);
		float a1= Distance_point(pointsABx[i].cn,pointsABx[i].ca);
		float d1= Distance_point(pointsABx[i].ca,pointsABx[i].cc);
		float e1= Distance_point(pointsABx[i].cn,pointsABx[i].cc);
		float NCAC = (180/3.1415926)*acos((a1*a1+d1*d1-e1*e1)/(2.0*a1*d1));  
		float v0,v1,v2;
		float x1 = pointsABx[3*i][0];
		float y1 = pointsABx[3*i][1];
		float z1 = pointsABx[3*i][2];
		float x2 = pointsABx[3*i+1][0];
		float y2 = pointsABx[3*i+1][1];
		float z2 = pointsABx[3*i+1][2];
		float x3 = pointsABx[3*i+2][0];
		float y3 = pointsABx[3*i+2][1];
		float z3 = pointsABx[3*i+2][2];
		v0 = (y3-y1)*(z3-z1)-(z2-z1)*(y3-y1);
		v1 = (x3-x1)*(z2-z1)-(x2-x1)*(z3-z1);
		v2 = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);
		float D = -(v0*x1+v1*y1+v2*z1);
		v0*xy1 + v1*xy2+v2*xy3 + D = 0;  */


	}
	for(int i=0;i<posex.chains.size();i++)
	{
//		cout<< "chainx size: "<<posex.chains.size()<<endl;
//		cout<< "chainy size: "<<posey.chains.size()<<endl;
		Chain Chanx = posex.chains[i];
//		Chain Chany = posey.chains[i];
		int tmp_nx=0;
		for(int j=0;j<pnum;j++)	
		{
//			cout<<"j: "<<j<<endl;
//			cout<< "residuex size: "<< Chanx.residues.size()<<endl;
//			cout<< "residuey size: "<< Chany.residues.size()<<endl;
			Residue Resdx = Chanx.residues[j];
//			Residue Resdy = Chany.residues[j];
//			Rot rotx,roty;
//			cout<<"Res: "<<Resdx.atoms.size()<<endl;
//			cout<<"res: "<<points_pre3[tmp_nx].res_atm.size()<<endl;
			int tt=0;
			for(int t=0;t<Resdx.atoms.size();t++)
			{
//				cout<<"t: "<<t<<endl;
//				cout<< "atomx size: "<<Resdx.atoms.size()<<endl;
//				cout<< "atomy size: "<<Resdy.atoms.size()<<endl;
				Atom Atmx=Resdx.atoms[t];
//				Atom Atmy=Resdy.heavyatoms[t];
				if(Atmx.atona != " CA " && Atmx.atona != " C  " && Atmx.atona != " N  " && Atmx.atona != " O  ")
				{	
					for(int tj=0;tj<CBp[j].Rotnm.size();tj++)
					{
						vector<string> atmnmx = string_splitxx(Atmx.atona,' ');
						string atmnmy = atmnmx[1];
//						cout<<" atmnmy,CBp: "<<atmnmy<<" "<<CBp[j].Rotnm[tj]<<endl;					
				//		if(atmnmy==CBp[j].Rotnm[tj]) {fin_pose.chains[i].residues[j].atoms[t].xyzVect = CBp[j].Rotcd[tj]; break;}
						if(strcmp(atmnmy.c_str(),CBp[j].Rotnm[tj].c_str())==0) {fin_pose.chains[i].residues[j].atoms[t].xyzVect = CBp[j].Rotcd[tj]; break;}
					}
//					fin_pose.chains[i].residues[j].atoms[t].xyzVect[0] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[0] - pointsABx0[3*j+1][0]+pointsABx[3*j+1][0];
//					fin_pose.chains[i].residues[j].atoms[t].xyzVect[1] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[1] - pointsABx0[3*j+1][1]+pointsABx[3*j+1][1];
//					fin_pose.chains[i].residues[j].atoms[t].xyzVect[2] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[2] - pointsABx0[3*j+1][2]+pointsABx[3*j+1][2];
//					tt=tt+1;
				}									
				if(Atmx.atona == " CA ") 
				{
					fin_pose.chains[i].residues[j].atoms[t].xyzVect = pointsABx[3*j+1];

				}
				if(Atmx.atona == " N  ") 
				{
					fin_pose.chains[i].residues[j].atoms[t].xyzVect = pointsABx[3*j];
				}
				if(Atmx.atona == " C  ") 
				{
					fin_pose.chains[i].residues[j].atoms[t].xyzVect = pointsABx[3*j+2];
				}
				if(Atmx.atona == " O  ") 
				{
					fin_pose.chains[i].residues[j].atoms[t].xyzVect = CBp[j].Ocd;
//					fin_pose.chains[i].residues[j].atoms[t].xyzVect[0] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[0] + pointsABx[3*j+2][0] - pointsABx0[3*j+2][0];
//					fin_pose.chains[i].residues[j].atoms[t].xyzVect[1] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[1] + pointsABx[3*j+2][1] - pointsABx0[3*j+2][1];
//					fin_pose.chains[i].residues[j].atoms[t].xyzVect[2] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[2] + pointsABx[3*j+2][2] - pointsABx0[3*j+2][2];
				}													
			}	
			tmp_nx = tmp_nx + 1; 
//			cout<< "tmP: "<<tmp_nx<<endl;
		}
	}	
/*	for(int i=0;i<posex.chains.size();i++)
	{
//		cout<< "chainx size: "<<posex.chains.size()<<endl;
//		cout<< "chainy size: "<<posey.chains.size()<<endl;
		Chain Chanx = posex.chains[i];
//		Chain Chany = posey.chains[i];
		int tmp_nx=0;
		for(int j=0;j<pointsx.size();j++)	
		{
//			cout<<"j: "<<j<<endl;
//			cout<< "residuex size: "<< Chanx.residues.size()<<endl;
//			cout<< "residuey size: "<< Chany.residues.size()<<endl;
			Residue Resdx = Chanx.residues[j];
//			Residue Resdy = Chany.residues[j];
//			Rot rotx,roty;
//			cout<<"Res: "<<Resdx.atoms.size()<<endl;
//			cout<<"res: "<<points_pre3[tmp_nx].res_atm.size()<<endl;
			int tt=0;
			for(int t=0;t<Resdx.atoms.size();t++)
			{
//				cout<<"t: "<<t<<endl;
//				cout<< "atomx size: "<<Resdx.atoms.size()<<endl;
//				cout<< "atomy size: "<<Resdy.atoms.size()<<endl;
				Atom Atmx=Resdx.atoms[t];
//				Atom Atmy=Resdy.heavyatoms[t];
				if(Atmx.atona != " CA " && Atmx.atona != " C  " && Atmx.atona != " N  " && Atmx.atona != " O  ")
				{
					fin_pose.chains[i].residues[j].atoms[t].xyzVect = pointsx[j].res_atm[tt];
					tt=tt+1;
				}					
				if(Atmx.atona == " CA ") 
				{
					fin_pose.chains[i].residues[j].atoms[t].xyzVect = pointsx[j].ca;
				}
				if(Atmx.atona == " N  ") 
				{
					fin_pose.chains[i].residues[j].atoms[t].xyzVect = pointsx[j].cn;
				}
				if(Atmx.atona == " C  ") 
				{
					fin_pose.chains[i].residues[j].atoms[t].xyzVect = pointsx[j].cc;
				}
				if(Atmx.atona == " O  ") 
				{
					fin_pose.chains[i].residues[j].atoms[t].xyzVect = pointsx[j].co;
				}													
			}	
			tmp_nx = tmp_nx + 1; 
//			cout<< "tmP: "<<tmp_nx<<endl;
		}
	}	*/

	cout<< "XXXXXXXXXX"<<endl;
	return fin_pose;	
//	writePDBcoordss(inputPDB, fin_pose);
}

int main(int argc, char* argv[])
{
	char inputPDB1[600]; 
	string inputMRC;
//	poseCoords pose1,pose2;
	Model pose;

	clock_t start,finish;
	double totaltime;
	start=clock();	

// 	readPDB
	strcpy(inputPDB1,argv[1]);
	readpdbstructurex(inputPDB1,pose);

	poseCoords posey;
	readPDBcoords(inputPDB1,posey);
	cout<< "Read PDB file : "<<inputPDB1<<endl;


	float MRC_reso;
	float mapsampling = 0.0;
	inputMRC = argv[2];
	MRC_reso = atof(argv[3]);
	mapsampling = atof(argv[4]);

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
	theDensityMap.pforce_apix_on_map_load_= 0.0;
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
//	theDensityMap.calcRhoC(posey,4.0,rhoC,rhoMask);
	theDensityMap.calcRhoCx(pose,4.0,rhoC,rhoMask);
//	theDensityMap.calcRhoCy(pose,4.0,rhoC,rhoMask);
	float CC;
//	CC=theDensityMap.getRSCC(rhoC,rhoMask);
	CC=theDensityMap.getRSCCX(rhoC,rhoMask);
	cout<<"CC: "<<CC<<endl;	


//	strcpy(inputPDB2,argv[2]);
//	readPDBcoordss(inputPDB1,pose1);
//	cout<<"000  "<<endl;
//	readPDBcoordss(inputPDB2,pose2);
//	readpdbstructurex(inputPDB2,pose2);
//	cout<< "11 "<<endl;
//	cout<<" 11 "<<endl;
	Model fin_p,fin_px,fin_py;	
	cout<<" 11 "<<endl;
	fin_p=REMC_sampley(pose,theDensityMap);
	cout<< "444"<<endl;
//	fin_px=REMC_samplex(fin_p,pose2);
//	fin_py=REMC_samplexy(fin_px,pose2);
///	writePDBcoordss(inputPDB1, fin_p);
	string outname;
	vector<string> tmp_str,tmp_strx;
	tmp_str = string_splitx(argv[1],'.');
	tmp_strx = string_splitx(tmp_str[0],'/');

	if(tmp_strx.size()>0)
	{
		outname = tmp_strx[tmp_strx.size()-1] + "_out.pdb";
	}
	else
	{
		outname = tmp_strx[0] + "_out.pdb";
	}
	
	writePDBStructure(inputPDB1,fin_p,outname);
   	finish=clock();
   	totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
   	cout<<"\nrunning time: "<<totaltime<<" s！"<<endl;  
   	return 0;
}

// g++ GeometryTools.cpp readpdb.cpp MC_RMSD.cpp -o MC_RMSD
// ./MC_RMSD 1a0fA_1.pdb 1a0fA.pdb
// g++ readpdb.cpp GeometryTools.cpp Rotbuilder.cpp GetVdwRadius.cpp MC_RMSD27.cpp -o MC_RMSD27