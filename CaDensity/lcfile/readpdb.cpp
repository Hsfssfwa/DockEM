/*
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <cstring>
#include <iomanip>
#include <vector>
#include <math.h>  */
#include "readpdb.h"

// using namespace std;



float RandomDoubleX(double numx,double numy)
{
    float ran_num=0.0;
 //   srand((unsigned)time(0));
 //   srand((unsigned int)(time(NULL)));
 //   ran_num=(rand()/((float) numy-numx))+numx;
    ran_num = numx + (numy-numx)*rand()/double(RAND_MAX);
    return ran_num;
}  
/* float RandomDoubleX(double numx,double numy)
{
    float ran_num=0.0;
 //   srand((unsigned)time(0));
 //   srand((unsigned int)(time(NULL)));
 //   ran_num=(rand()/((float) numy-numx))+numx;
 //   ran_num = numx + (numy-numx)*rand() / double(RAND_MAX+1.0);
    default_random_engine gen;
    uniform_real_distribution<float> dis(numx,numy);
    ran_num = dis(gen);
    return ran_num;
}  */
 

string AAtoS(string animoacide)
{
    if(animoacide==" ALA"||animoacide=="AALA"){return ("A") ;}
    if(animoacide==" CYS"||animoacide=="ACYS"){return ("C");}
    if(animoacide==" ASP"||animoacide=="AASP"){return ("D");}
    if(animoacide==" GLU"||animoacide=="AGLU"){return ("E");}
    if(animoacide==" PHE"||animoacide=="APHE"){return ("F");}
    if(animoacide==" GLY"||animoacide=="AGLY"){return ("G");}
    if(animoacide==" HIS"||animoacide=="AHIS"){return ("H");}
    if(animoacide==" ILE"||animoacide=="AILE"){return ("I");}
    if(animoacide==" LYS"||animoacide=="ALYS"){return ("K");}
    if(animoacide==" LEU"||animoacide=="ALEU"){return ("L");}
    if(animoacide==" MET"||animoacide=="AMET"){return ("M");}
    if(animoacide==" ASN"||animoacide=="AASN"){return ("N");}
    if(animoacide==" PRO"||animoacide=="APRO"){return ("P");}
    if(animoacide==" GLN"||animoacide=="AGLN"){return ("Q");}
    if(animoacide==" ARG"||animoacide=="AARG"){return ("R");}
    if(animoacide==" SER"||animoacide=="ASER"){return ("S");}
    if(animoacide==" THR"||animoacide=="ATHR"){return ("T");}
    if(animoacide==" VAL"||animoacide=="AVAL"){return ("V");}
    if(animoacide==" TRP"||animoacide=="ATRP"){return ("W");}
    if(animoacide==" TYR"||animoacide=="ATYR"){return ("Y");}
}

string StoAA(string animoacide)
{
    if(animoacide=="A"){return ("ALA") ;}
    if(animoacide=="C"){return ("CYS");}
    if(animoacide=="D"){return ("ASP");}
    if(animoacide=="E"){return ("GLU");}
    if(animoacide=="F"){return ("PHE");}
    if(animoacide=="G"){return ("GLY");}
    if(animoacide=="H"){return ("HIS");}
    if(animoacide=="I"){return ("ILE");}
    if(animoacide=="K"){return ("LYS");}
    if(animoacide=="L"){return ("LEU");}
    if(animoacide=="M"){return ("MET");}
    if(animoacide=="N"){return ("ASN");}
    if(animoacide=="P"){return ("PRO");}
    if(animoacide=="Q"){return ("GLN");}
    if(animoacide=="R"){return ("ARG");}
    if(animoacide=="S"){return ("SER");}
    if(animoacide=="T"){return ("THR");}
    if(animoacide=="V"){return ("VAL");}
    if(animoacide=="W"){return ("TRP");}
    if(animoacide=="Y"){return ("TYR");}
}

/* atomic_detail:0 - CA only, 1-backbone heavy atoms (CA C N O), 2 - all atom
 * allowX: 0 - ATOM, 1 - ATOM and MSE, converting MSE to MET,
 *         2 - all, converting MSE to MET, 3 - all, no convertsion
 * filename: full filename path, stdin if filename=="-"
*/
void readpdbstructure(const char *filename,Model &proteins)
{
    string atomline;
    Atom newatom;
    Residue newresidue;
    Chain newchain;
//    Model proteins;
    string chain_old="as";
    string chain_new="xd";
    string residue_old="ab";
    string residue_new="cd";
    string old_icode="sa";
    int Nresidue_old=11111;
    int Nresidue_new=11112;
    int tgx=0;
    int tgy=0;
    int countx=0;

    ifstream ReadFile;
    ReadFile.open(filename);
    if(!ReadFile.is_open())
    {
        cout<<"error to open the file"<<endl;
    }
    while(!ReadFile.eof())
    {
        getline(ReadFile,atomline);
        if(atomline.substr(0,6)=="ATOM  "&& atomline.substr(13,3) !="OXT")
        {
            chain_new=atomline.substr(21,1);
            residue_new=atomline.substr(16,4);
        Nresidue_new=atoi(atomline.substr(22,4).c_str());

        newatom.recna=atomline.substr(0,6);
        newatom.serial=atoi(atomline.substr(6,5).c_str());
            newatom.atona=atomline.substr(12,4);
        newatom.chara=atomline.substr(16,1);
            newatom.xyzVect.push_back(strtof(atomline.substr(30,8).c_str(),0));
            newatom.xyzVect.push_back(strtof(atomline.substr(38,8).c_str(),0));
            newatom.xyzVect.push_back(strtof(atomline.substr(46,8).c_str(),0));
            if(atomline.size()>60){
                newatom.tempfac=(strtof(atomline.substr(60,6).c_str(),0));
            }
         
            

            // heavy atoms
  //          if(atomline.substr(13,14)=="H")
  //          {
//
//            }
            
            if(Nresidue_new!=Nresidue_old && tgx>0)
            {
                newresidue.resseq=Nresidue_old;
                newresidue.icode=atomline.substr(26,1);
                newresidue.resname=AAtoS(residue_old);
                newchain.residues.push_back(newresidue);
                newresidue.atoms.clear();
//                vector<Atom>(newresidue.atoms).swap(newresidue.atoms);
//                newresidue.atoms.push_back(newatom);
//                newatom.xyzVect.clear();
		countx=countx+1;
            }
            tgx=1;

            // heavy atoms
            if(atomline.substr(13,1) !="H" && atomline.substr(13,3) !="OXT")
            {
                newresidue.heavyatoms.push_back(newatom);
            }   

            newresidue.atoms.push_back(newatom);
            newatom.xyzVect.clear();
            residue_old=residue_new;
        Nresidue_old=Nresidue_new;
        old_icode=atomline.substr(26,1);
/*
            if(residue_new==residue_old)
            {
                newresidue.atoms.push_back(newatom);
                newatom.xyzVect.clear();
                residue_old==residue_new;
            }
            else
            {
                newresidue.resseq=atoi(atomline.substr(22,4).c_str());
                newresidue.icode=atomline.substr(26,1);
                newresidue.resname=residue_new;
                newchain.residues.push_back(newresidue);
                newresidue.atoms.clear();
//                vector<Atom>(newresidue.atoms).swap(newresidue.atoms);
                newresidue.atoms.push_back(newatom);
                newatom.xyzVect.clear();
                residue_old==residue_new;
            }  */

            if(chain_new!=chain_old && tgy>0)
            {
                newchain.chaid=chain_old;
                proteins.chains.push_back(newchain);
                newchain.residues.clear();
//                vector<Residue>(newchain.residues).swap(newchain.residues);                
            }
            chain_old=chain_new;
            tgy=1;            
            
        }
        if(atomline.substr(0,6)=="ENDMDL" ){
        break;
    }
    }
    newresidue.resseq=Nresidue_old;
    newresidue.icode=old_icode;
    newresidue.resname=AAtoS(residue_old);
    newchain.residues.push_back(newresidue);
    newresidue.atoms.clear();
    newchain.chaid=chain_old;
    proteins.chains.push_back(newchain);
    newchain.residues.clear();
    ReadFile.close();
}

void readpdbstructurex(const char *filename,Model &proteins)
{
    string atomline;
    Atom newatom;
    Residue newresidue;
    Chain newchain;
//    Model proteins;
    string chain_old="as";
    string chain_new="xd";
    string residue_old="ab";
    string residue_new="cd";
    string old_icode="sa";
    int Nresidue_old=11111;
    int Nresidue_new=11112;
    int tgx=0;
    int tgy=0;
    int countx=0;

    ifstream ReadFile;
    ReadFile.open(filename);
    if(!ReadFile.is_open())
    {
        cout<<"error to open the file"<<endl;
    }
    while(!ReadFile.eof())
    {
        getline(ReadFile,atomline);
        if(atomline.substr(0,6)=="ATOM  "&& atomline.substr(13,3) !="OXT" && atomline.substr(13,1) !="H")
        {
            chain_new=atomline.substr(21,1);
            residue_new=atomline.substr(16,4);
        Nresidue_new=atoi(atomline.substr(22,4).c_str());

        newatom.recna=atomline.substr(0,6);
        newatom.serial=atoi(atomline.substr(6,5).c_str());
            newatom.atona=atomline.substr(12,4);
        newatom.chara=atomline.substr(16,1);
            newatom.xyzVect.push_back(strtof(atomline.substr(30,8).c_str(),0));
            newatom.xyzVect.push_back(strtof(atomline.substr(38,8).c_str(),0));
            newatom.xyzVect.push_back(strtof(atomline.substr(46,8).c_str(),0));
            if(atomline.size()>60){
                newatom.tempfac=(strtof(atomline.substr(60,6).c_str(),0));
            }
         
            

            // heavy atoms
  //          if(atomline.substr(13,14)=="H")
  //          {
//
//            }
            
            if(Nresidue_new!=Nresidue_old && tgx>0)
            {
                newresidue.resseq=Nresidue_old;
                newresidue.icode=atomline.substr(26,1);
                newresidue.resname=AAtoS(residue_old);
                newchain.residues.push_back(newresidue);
                newresidue.atoms.clear();
//                vector<Atom>(newresidue.atoms).swap(newresidue.atoms);
//                newresidue.atoms.push_back(newatom);
//                newatom.xyzVect.clear();
        countx=countx+1;
            }
            tgx=1;

            // heavy atoms
            if(atomline.substr(13,1) !="H" && atomline.substr(13,3) !="OXT")
            {
                newresidue.heavyatoms.push_back(newatom);
            }   

            newresidue.atoms.push_back(newatom);
            newatom.xyzVect.clear();
            residue_old=residue_new;
        Nresidue_old=Nresidue_new;
        old_icode=atomline.substr(26,1);
/*
            if(residue_new==residue_old)
            {
                newresidue.atoms.push_back(newatom);
                newatom.xyzVect.clear();
                residue_old==residue_new;
            }
            else
            {
                newresidue.resseq=atoi(atomline.substr(22,4).c_str());
                newresidue.icode=atomline.substr(26,1);
                newresidue.resname=residue_new;
                newchain.residues.push_back(newresidue);
                newresidue.atoms.clear();
//                vector<Atom>(newresidue.atoms).swap(newresidue.atoms);
                newresidue.atoms.push_back(newatom);
                newatom.xyzVect.clear();
                residue_old==residue_new;
            }  */

            if(chain_new!=chain_old && tgy>0)
            {
                newchain.chaid=chain_old;
                proteins.chains.push_back(newchain);
                newchain.residues.clear();
//                vector<Residue>(newchain.residues).swap(newchain.residues);                
            }
            chain_old=chain_new;
            tgy=1;            
            
        }
        if(atomline.substr(0,6)=="ENDMDL" ){
            break;
        }
    }
    newresidue.resseq=Nresidue_old;
    newresidue.icode=old_icode;
    newresidue.resname=AAtoS(residue_old);
    newchain.residues.push_back(newresidue);
    newresidue.atoms.clear();
    newchain.chaid=chain_old;
    proteins.chains.push_back(newchain);
    newchain.residues.clear();
    ReadFile.close();

}


void readcoordinatefrompdb(const char *filename,coordtype &coord)
{
        string atomline;
//        vector<vector<float> > coord(3);

        ifstream ReadFile;
 //       ofstream outfile;
  //      outfile.open(outputfilex);
        ReadFile.open(filename);
        if(!ReadFile.is_open())
        {
                cout<<"error to open the file"<<endl;
        }
        while(!ReadFile.eof())
        {
                getline(ReadFile,atomline);
                if(atomline.substr(0,6)=="ATOM  "){
                        coord.recna.push_back(atomline.substr(0,6));
  //                      outfile<<atomline.substr(0,5)<<"   ";
                        coord.serial.push_back(atoi(atomline.substr(6,5).c_str()));
                        coord.atona.push_back(atomline.substr(12,4));
                        coord.chara.push_back(atomline.substr(16,1));
                        coord.resname.push_back(atomline.substr(17,3));
                        coord.chaid.push_back(atomline.substr(21,1));
                        coord.resseq.push_back(atoi(atomline.substr(22,4).c_str()));
                        coord.icode.push_back(atomline.substr(26,1));
                        coord.cx.push_back(strtof(atomline.substr(30,8).c_str(),0));
                        coord.cy.push_back(strtof(atomline.substr(38,8).c_str(),0));
                        coord.cz.push_back(strtof(atomline.substr(46,8).c_str(),0));
                        coord.occ.push_back(strtof(atomline.substr(54,6).c_str(),0));
                        coord.tempfac.push_back(strtof(atomline.substr(60,6).c_str(),0));
                        coord.ele.push_back(atomline.substr(76,2));
                        coord.chag.push_back(atomline.substr(78,2));
              /*          coord[0].push_back(strtof(atomline.substr(31,7).c_str(),0));
                        outfile<<strtof(atomline.substr(31,7).c_str(),0)<<"   ";
                        coord[1].push_back(strtof(atomline.substr(39,7).c_str(),0));
                        outfile<<strtof(atomline.substr(39,7).c_str(),0)<<"   ";
                        coord[2].push_back(strtof(atomline.substr(48,7).c_str(),0));
                        outfile<<strtof(atomline.substr(48,7).c_str(),0)<<"   "<<endl; */
                }
                else if(atomline.substr(0,6)=="ENDMDL"){break;}
                else{
                        coord.costr.push_back(atomline);
                }  
        }
//        outfile.close();
        ReadFile.close();
}
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
}

std::string
name2eltx( std::string atmid ) {
    std::string atmidx = atmid ;
    while ( !atmidx.empty() && atmidx[0] == ' ' ) atmidx = atmidx.substr(1,atmidx.size()-1);
    while ( !atmidx.empty() && atmidx[atmidx.size()-1] == ' ' ) atmidx = atmidx.substr(0,atmidx.size()-1);

    std::string type;
    if(atmidx=="C") type = "C";
    else if(atmidx=="N") type = "N";
    else if(atmidx=="CA") type = "CA";
    else if(atmidx=="O") type = "O";
    else type = atmid[1];
//    cout<< type<<endl;
    return type;
}


// quick and dirty PDB read where we only care about heavyatom locations and atom ids
void
readPDBcoords(std::string filename, poseCoords &atmlist) {
    std::ifstream inpdb(filename.c_str());
//  cout<< filename.c_str()<<endl;
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

void writePDBcoords(std::string filename, poseCoords &atmlist)
{
    std::ifstream inpdb(filename.c_str());
//  cout<< filename.c_str()<<endl;
    std::string buf;
    vector<string> tmp_str;
    tmp_str = string_splitx(filename,'.');
    string filenamex = tmp_str[0] + "_out.pdb";
    std::ofstream outpdb;
    outpdb.open(filenamex.c_str());
    int tmpx=0;

    while ( std::getline(inpdb, buf ) ) {
        if ( buf.substr(0,4)!="ATOM" && buf.substr(0,6)!="HETATM" ) continue;

        poseCoord atom_i;
        string tmp_A,tmp_B;

        tmp_A = buf.substr(0,30);
        atom_i.x_.push_back(atof(buf.substr(30,8).c_str()));
        atom_i.x_.push_back(atof(buf.substr(38,8).c_str()));
        atom_i.x_.push_back(atof(buf.substr(46,8).c_str()));
        tmp_B = buf.substr(54,6);
        atom_i.B_ = atof(buf.substr(60,6).c_str());
    //    cout<<"YYY"<<endl;


        atom_i.elt_ = name2elt( buf ); // horrible hacky logic mapping name->elt (could use PDB fields 76-77 if on by default)
        if ( atom_i.elt_ == "H" ) continue;

        outpdb<<setw(30)<<tmp_A;
        outpdb<<setw(8)<<setiosflags(ios::fixed)<<setprecision(3)<<atmlist[tmpx].x_[0];
        outpdb<<setw(8)<<setiosflags(ios::fixed)<<setprecision(3)<<atmlist[tmpx].x_[1];
        outpdb<<setw(8)<<setiosflags(ios::fixed)<<setprecision(3)<<atmlist[tmpx].x_[2];
        outpdb<<setw(6)<<tmp_B;
        outpdb<<setw(6)<<setiosflags(ios::fixed)<<setprecision(2)<<atmlist[tmpx].B_;
        outpdb<<endl;
        tmpx = tmpx + 1;        
//        atmlist.push_back( atom_i );
    }
    outpdb.close();
}

poseCoords rotatePDB(poseCoords &pose)
{
    poseCoords new_pose;
    int pnum=pose.size();
    new_pose = pose;
    float R=RandomDoubleX(0.0,1.0);

// Rotation axis    
    float asin_theta=1.0-2.0*R;
    float acos_theta=sqrt(1.0-asin_theta*asin_theta);
    float apha=2.0*PI*R;
    float awx=acos_theta*cos(apha);
    float awy=acos_theta*sin(apha);
    float awz=asin_theta;

// Translation Vector
    float t0=0.0;
    float t1=(R * 2.0-1.0)*t0;
    float t2=(R * 2.0-1.0)*t0;
    float t3=(R * 2.0-1.0)*t0;

// Rotation points
//    float *axyz =new float[3];
    float axyz[3];
 /*   for(int i=0;i<3;i++)
        axyz[i]=m[i][0]; */
 //   axyz[0]=m.cx.at(0);
 //   axyz[1]=m.cy.at(0);
 //   axyz[2]=m.cz.at(0);
    axyz[0]=0.0;
    axyz[1]=0.0;
    axyz[2]=0.0;
 /*  for(int i=0;i<pnum;i++)
    {
        axyz[0] = axyz[0] + pose[i].x_[0];
        axyz[1] = axyz[1] + pose[i].x_[1];
        axyz[2] = axyz[2] + pose[i].x_[2];
    }
    axyz[0] = axyz[0]/pnum;
    axyz[1] = axyz[1]/pnum;
    axyz[2] = axyz[2]/pnum;  */

// Rotation matrix
    float angle_rotate=(PI/2)*(2.0 * R-1.0); // rotate angle (-PI,PI)
    float asin=sin(angle_rotate);
    float acos=cos(angle_rotate);
    float u[3][3];
    u[0][0]=acos+awx*awx*(1.0-acos);
    u[0][1]=awx*awy*(1.0-acos)-awz*asin;
    u[0][2]=awx*awz*(1.0-acos)+awy*asin;
    u[1][0]=awy*awz*(1.0-acos)+awz*asin;
    u[1][1]=acos+awy*awy*(1.0-acos);
    u[1][2]=awy*awz*(1.0-acos)-awx*asin;
    u[2][0]=awx*awz*(1.0-acos)-awy*asin;
    u[2][1]=awy*awz*(1.0-acos)+awx*asin;
    u[2][2]=acos+awz*awz*(1.0-acos);

// Translation and Rotation(including translate to the rotation point)
    for(int i=0;i<pnum;i++)
    {
        new_pose[i].x_[0]=t1+axyz[0]+(pose[i].x_[0]-axyz[0])*u[0][0]+(pose[i].x_[1]-axyz[1])*u[0][1]+(pose[i].x_[2]-axyz[2])*u[0][2];
        new_pose[i].x_[1]=t2+axyz[1]+(pose[i].x_[0]-axyz[0])*u[1][0]+(pose[i].x_[1]-axyz[1])*u[1][1]+(pose[i].x_[2]-axyz[2])*u[1][2];
        new_pose[i].x_[2]=t3+axyz[2]+(pose[i].x_[0]-axyz[0])*u[2][0]+(pose[i].x_[1]-axyz[1])*u[2][1]+(pose[i].x_[2]-axyz[2])*u[2][2];
    }
/*    for(int i=0;i<pnum;i++)
    {
        new_pose[i].x_[0]=t1+axyz[0]+(pose[i].x_[0]-pose[0].x_[0])*u[0][0]+(pose[i].x_[1]-pose[0].x_[1])*u[0][1]+(pose[i].x_[2]-pose[0].x_[2])*u[0][2];
        new_pose[i].x_[1]=t2+axyz[1]+(pose[i].x_[0]-pose[0].x_[0])*u[1][0]+(pose[i].x_[1]-pose[0].x_[1])*u[1][1]+(pose[i].x_[2]-pose[0].x_[2])*u[1][2];
        new_pose[i].x_[2]=t3+axyz[2]+(pose[i].x_[0]-pose[0].x_[0])*u[2][0]+(pose[i].x_[1]-pose[0].x_[1])*u[2][1]+(pose[i].x_[2]-pose[0].x_[2])*u[2][2];
    }   */
    return new_pose;
}

void rotatepdb(coordtype &m)
{
    int coord_num=m.cx.size();
    vector<vector<float> > fin_mat(3,vector<float>(coord_num));

// Rotation axis    
    float asin_theta=1.0-2.0*RandomDoubleX(0.0,1.0);
    float acos_theta=sqrt(1.0-asin_theta*asin_theta);
    float apha=2.0*PI*RandomDoubleX(0.0,1.0);
    float awx=acos_theta*cos(apha);
    float awy=acos_theta*sin(apha);
    float awz=asin_theta;

// Translation Vector
    float t0=0.3;
    float t1=(RandomDoubleX(0.0,1.0)*2.0-1.0)*t0;
    float t2=(RandomDoubleX(0.0,1.0)*2.0-1.0)*t0;
    float t3=(RandomDoubleX(0.0,1.0)*2.0-1.0)*t0;

// Rotation points
    float *axyz =new float[3];
 /*   for(int i=0;i<3;i++)
        axyz[i]=m[i][0]; */
    axyz[0]=m.cx.at(0);
    axyz[1]=m.cy.at(0);
    axyz[2]=m.cz.at(0);

// Rotation matrix
    float angle_rotate=(2.0*RandomDoubleX(0.0,1.0)-1.0); // rotate angle
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

// Translation and Rotation(including translate to the rotation point)
    for(int i=0;i<coord_num;i++)
    {
        fin_mat[0][i]=t1+axyz[0]+(m.cx.at(i)-axyz[0])*u[0][0]+(m.cy.at(i)-axyz[1])*u[0][1]+(m.cz.at(i)-axyz[3])*u[0][2];
        fin_mat[1][i]=t2+axyz[1]+(m.cx.at(i)-axyz[0])*u[1][0]+(m.cy.at(i)-axyz[1])*u[1][1]+(m.cz.at(i)-axyz[3])*u[1][2];
        fin_mat[2][i]=t3+axyz[2]+(m.cx.at(i)-axyz[0])*u[2][0]+(m.cy.at(i)-axyz[1])*u[2][1]+(m.cz.at(i)-axyz[3])*u[2][2];
    }
    for(int i=0;i<coord_num;i++)
    {
        m.cx.at(i)=fin_mat[0][i];
        m.cy.at(i)=fin_mat[1][i];
        m.cz.at(i)=fin_mat[2][i];
    }
 //   outputm=fin_mat;
}

void writePDBStructure(const char *filename,Model &protein)
{
        ofstream outfile;
        string fin_name = "fin.pdb";
//        outfile.open(filename);
        outfile.open(fin_name.c_str());
        if(!outfile.is_open())
        {
                cout<<"error to open the file"<<endl;
        }

    int n=1;
	for(int i=0;i<protein.chains.size();i++)
	{
		for(int j=0;j<protein.chains[i].residues.size();j++)
		{
			for(int k=0;k<protein.chains[i].residues[j].atoms.size();k++)
			{
				outfile<<protein.chains[i].residues[j].atoms[k].recna;
                //		outfile<<setw(5)<<protein.chains[i].residues[j].atoms[k].serial;
                        outfile<<setw(5)<<n;
               			outfile<<setw(5)<<protein.chains[i].residues[j].atoms[k].atona;
              			outfile<<setw(1)<<protein.chains[i].residues[j].atoms[k].chara;
                		outfile<<setw(3)<<StoAA(protein.chains[i].residues[j].resname);
                		outfile<<setw(2)<<protein.chains[i].chaid;
                		outfile<<setw(4)<<protein.chains[i].residues[j].resseq;
                		outfile<<setw(1)<<protein.chains[i].residues[j].icode;
                		outfile<<setw(11)<<setiosflags(ios::fixed)<<setprecision(3)<<protein.chains[i].residues[j].atoms[k].xyzVect[0];
                		outfile<<setw(8)<<setiosflags(ios::fixed)<<setprecision(3)<<protein.chains[i].residues[j].atoms[k].xyzVect[1];
                		outfile<<setw(8)<<setiosflags(ios::fixed)<<setprecision(3)<<protein.chains[i].residues[j].atoms[k].xyzVect[2];
//                		outfile<<setw(6)<<setiosflags(ios::fixed)<<setprecision(2)<<coordx.occ[i];
//                		outfile<<setw(6)<<setiosflags(ios::fixed)<<setprecision(2)<<coordx.tempfac[i];
//                		outfile<<setw(12)<<coordx.ele[i];
//                		outfile<<setw(2)<<coordx.chag[i];
                		outfile<<endl;
                        n=n+1;	
			}
		}
	}
	outfile.close();	
}

void writePDBStructure(const char *filename,Model &protein,string outname)
{
        ofstream outfile;
 //       string fin_name = "fin.pdb";
        string fin_name = outname;
//        outfile.open(filename);
        outfile.open(fin_name.c_str());
        if(!outfile.is_open())
        {
                cout<<"error to open the file"<<endl;
        }

    int n=1;
    for(int i=0;i<protein.chains.size();i++)
    {
        for(int j=0;j<protein.chains[i].residues.size();j++)
        {
            for(int k=0;k<protein.chains[i].residues[j].atoms.size();k++)
            {
                outfile<<protein.chains[i].residues[j].atoms[k].recna;
                //      outfile<<setw(5)<<protein.chains[i].residues[j].atoms[k].serial;
                        outfile<<setw(5)<<n;
                        outfile<<setw(5)<<protein.chains[i].residues[j].atoms[k].atona;
                        outfile<<setw(1)<<protein.chains[i].residues[j].atoms[k].chara;
                        outfile<<setw(3)<<StoAA(protein.chains[i].residues[j].resname);
                        outfile<<setw(2)<<protein.chains[i].chaid;
                        outfile<<setw(4)<<protein.chains[i].residues[j].resseq;
                        outfile<<setw(1)<<protein.chains[i].residues[j].icode;
                        outfile<<setw(11)<<setiosflags(ios::fixed)<<setprecision(3)<<protein.chains[i].residues[j].atoms[k].xyzVect[0];
                        outfile<<setw(8)<<setiosflags(ios::fixed)<<setprecision(3)<<protein.chains[i].residues[j].atoms[k].xyzVect[1];
                        outfile<<setw(8)<<setiosflags(ios::fixed)<<setprecision(3)<<protein.chains[i].residues[j].atoms[k].xyzVect[2];
//                      outfile<<setw(6)<<setiosflags(ios::fixed)<<setprecision(2)<<coordx.occ[i];
//                      outfile<<setw(6)<<setiosflags(ios::fixed)<<setprecision(2)<<coordx.tempfac[i];
//                      outfile<<setw(12)<<coordx.ele[i];
//                      outfile<<setw(2)<<coordx.chag[i];
                        outfile<<endl;
                        n=n+1;  
            }
        }
    }
    outfile.close();    
}

void writePDB(const char *filename,coordtype &coordx)
{
    //    string atomline;
        
        ofstream outfile;
        outfile.open(filename);
        if(!outfile.is_open())
        {
                cout<<"error to open the file"<<endl;
        }
//        outfile<<coordx.recna.size();
        for(int i=0;i<coordx.recna.size();i++)
        {
 //               outfile<<111;
 //               outfile.seekp(1150,ios::cur);
                outfile<<coordx.recna[i];
                outfile<<setw(5)<<coordx.serial[i];
                outfile<<setw(5)<<coordx.atona[i];
                outfile<<setw(1)<<coordx.chara[i];
                outfile<<setw(3)<<coordx.resname[i];
                outfile<<setw(2)<<coordx.chaid[i];
                outfile<<setw(4)<<coordx.resseq[i];
                outfile<<setw(1)<<coordx.icode[i];
                outfile<<setw(11)<<setiosflags(ios::fixed)<<setprecision(3)<<coordx.cx[i];
                outfile<<setw(8)<<setiosflags(ios::fixed)<<setprecision(3)<<coordx.cy[i];
                outfile<<setw(8)<<setiosflags(ios::fixed)<<setprecision(3)<<coordx.cz[i];
                outfile<<setw(6)<<setiosflags(ios::fixed)<<setprecision(2)<<coordx.occ[i];
                outfile<<setw(6)<<setiosflags(ios::fixed)<<setprecision(2)<<coordx.tempfac[i];
                outfile<<setw(12)<<coordx.ele[i];
                outfile<<setw(2)<<coordx.chag[i];
                outfile<<endl;
        }
        outfile.close();
}

bool GroupRotationy(const vector<float> &axisA, const vector<float> &axisB, float angle, vector<Rot> &pointB)
{
   if(axisA.size()!=3 || axisB.size()!=3)
   {
      cout<<"Error in GroupRotation()"<<endl;
      return false;
   }
   
   vector<float> axis(3, 0);
   axis[0]=axisB[0]-axisA[0]; 
   axis[1]=axisB[1]-axisA[1];
   axis[2]=axisB[2]-axisA[2]; 
   
   matrix rotmtx;
   if(!RotationMatrixB(axis, deg2rad(angle), rotmtx)) return false;
   
   //vector<float> point_A(4, 0);
   float point_A[3];
   for(int i=0;i<(pointB.size());i++)
   {
      //rotate side chain
      for(int j=0;j<pointB[i].res_atm.size();j++)
      {
         point_A[0] = pointB[i].res_atm[j][0]- axisA[0];
         point_A[1] = pointB[i].res_atm[j][1]- axisA[1];
         point_A[2] = pointB[i].res_atm[j][2]- axisA[2];

         pointB[i].res_atm[j] = axisA;
         pointB[i].res_atm[j][0] += rotmtx[0][0]*point_A[0]+rotmtx[0][1]*point_A[1]+rotmtx[0][2]*point_A[2]+rotmtx[0][3];
         pointB[i].res_atm[j][1] += rotmtx[1][0]*point_A[0]+rotmtx[1][1]*point_A[1]+rotmtx[1][2]*point_A[2]+rotmtx[1][3];
         pointB[i].res_atm[j][2] += rotmtx[2][0]*point_A[0]+rotmtx[2][1]*point_A[1]+rotmtx[2][2]*point_A[2]+rotmtx[2][3];
      }
      // rotate N, CA,C,O
       point_A[0] = pointB[i].cn[0]- axisA[0];
       point_A[1] = pointB[i].cn[1]- axisA[1];
       point_A[2] = pointB[i].cn[2]- axisA[2];
       pointB[i].cn = axisA;
       pointB[i].cn[0] += rotmtx[0][0]*point_A[0]+rotmtx[0][1]*point_A[1]+rotmtx[0][2]*point_A[2]+rotmtx[0][3];
       pointB[i].cn[1] += rotmtx[1][0]*point_A[0]+rotmtx[1][1]*point_A[1]+rotmtx[1][2]*point_A[2]+rotmtx[1][3];
       pointB[i].cn[2] += rotmtx[2][0]*point_A[0]+rotmtx[2][1]*point_A[1]+rotmtx[2][2]*point_A[2]+rotmtx[2][3];    

       point_A[0] = pointB[i].ca[0]- axisA[0];
       point_A[1] = pointB[i].ca[1]- axisA[1];
       point_A[2] = pointB[i].ca[2]- axisA[2];
       pointB[i].ca = axisA;
       pointB[i].ca[0] += rotmtx[0][0]*point_A[0]+rotmtx[0][1]*point_A[1]+rotmtx[0][2]*point_A[2]+rotmtx[0][3];
       pointB[i].ca[1] += rotmtx[1][0]*point_A[0]+rotmtx[1][1]*point_A[1]+rotmtx[1][2]*point_A[2]+rotmtx[1][3];
       pointB[i].ca[2] += rotmtx[2][0]*point_A[0]+rotmtx[2][1]*point_A[1]+rotmtx[2][2]*point_A[2]+rotmtx[2][3];       

       point_A[0] = pointB[i].cc[0]- axisA[0];
       point_A[1] = pointB[i].cc[1]- axisA[1];
       point_A[2] = pointB[i].cc[2]- axisA[2];
       pointB[i].cc = axisA;
       pointB[i].cc[0] += rotmtx[0][0]*point_A[0]+rotmtx[0][1]*point_A[1]+rotmtx[0][2]*point_A[2]+rotmtx[0][3];
       pointB[i].cc[1] += rotmtx[1][0]*point_A[0]+rotmtx[1][1]*point_A[1]+rotmtx[1][2]*point_A[2]+rotmtx[1][3];
       pointB[i].cc[2] += rotmtx[2][0]*point_A[0]+rotmtx[2][1]*point_A[1]+rotmtx[2][2]*point_A[2]+rotmtx[2][3];       

       point_A[0] = pointB[i].co[0]- axisA[0];
       point_A[1] = pointB[i].co[1]- axisA[1];
       point_A[2] = pointB[i].co[2]- axisA[2];
       pointB[i].co = axisA;
       pointB[i].co[0] += rotmtx[0][0]*point_A[0]+rotmtx[0][1]*point_A[1]+rotmtx[0][2]*point_A[2]+rotmtx[0][3];
       pointB[i].co[1] += rotmtx[1][0]*point_A[0]+rotmtx[1][1]*point_A[1]+rotmtx[1][2]*point_A[2]+rotmtx[1][3];
       pointB[i].co[2] += rotmtx[2][0]*point_A[0]+rotmtx[2][1]*point_A[1]+rotmtx[2][2]*point_A[2]+rotmtx[2][3];       
   }    
   return true;
}

bool GroupRotationx(const vector<float> &axisA, const vector<float> &axisB, float angle, vector<Rot> &pointB)
{
   if(axisA.size()!=3 || axisB.size()!=3)
   {
      cout<<"Error in GroupRotation()"<<endl;
      return false;
   }
   
   vector<float> axis(3, 0);
   axis[0]=axisB[0]-axisA[0]; 
   axis[1]=axisB[1]-axisA[1];
   axis[2]=axisB[2]-axisA[2]; 
   
   matrix rotmtx;
   if(!RotationMatrixB(axis, deg2rad(angle), rotmtx)) return false;
   
   //vector<float> point_A(4, 0);
   float point_A[3];
   for(int i=1;i<(pointB.size()-1);i++)
   {
      //rotate side chain
      for(int j=0;j<pointB[i].res_atm.size();j++)
      {
         point_A[0] = pointB[i].res_atm[j][0]- axisA[0];
         point_A[1] = pointB[i].res_atm[j][1]- axisA[1];
         point_A[2] = pointB[i].res_atm[j][2]- axisA[2];

         pointB[i].res_atm[j] = axisA;
         pointB[i].res_atm[j][0] += rotmtx[0][0]*point_A[0]+rotmtx[0][1]*point_A[1]+rotmtx[0][2]*point_A[2];//+rotmtx[0][3];
         pointB[i].res_atm[j][1] += rotmtx[1][0]*point_A[0]+rotmtx[1][1]*point_A[1]+rotmtx[1][2]*point_A[2];//+rotmtx[1][3];
         pointB[i].res_atm[j][2] += rotmtx[2][0]*point_A[0]+rotmtx[2][1]*point_A[1]+rotmtx[2][2]*point_A[2];//+rotmtx[2][3];
      }
      // rotate N, CA,C,O
       point_A[0] = pointB[i].cn[0]- axisA[0];
       point_A[1] = pointB[i].cn[1]- axisA[1];
       point_A[2] = pointB[i].cn[2]- axisA[2];
       pointB[i].cn = axisA;
       pointB[i].cn[0] += rotmtx[0][0]*point_A[0]+rotmtx[0][1]*point_A[1]+rotmtx[0][2]*point_A[2];//+rotmtx[0][3];
       pointB[i].cn[1] += rotmtx[1][0]*point_A[0]+rotmtx[1][1]*point_A[1]+rotmtx[1][2]*point_A[2];//+rotmtx[1][3];
       pointB[i].cn[2] += rotmtx[2][0]*point_A[0]+rotmtx[2][1]*point_A[1]+rotmtx[2][2]*point_A[2];//+rotmtx[2][3];    

       point_A[0] = pointB[i].ca[0]- axisA[0];
       point_A[1] = pointB[i].ca[1]- axisA[1];
       point_A[2] = pointB[i].ca[2]- axisA[2];
       pointB[i].ca = axisA;
       pointB[i].ca[0] += rotmtx[0][0]*point_A[0]+rotmtx[0][1]*point_A[1]+rotmtx[0][2]*point_A[2];//+rotmtx[0][3];
       pointB[i].ca[1] += rotmtx[1][0]*point_A[0]+rotmtx[1][1]*point_A[1]+rotmtx[1][2]*point_A[2];//+rotmtx[1][3];
       pointB[i].ca[2] += rotmtx[2][0]*point_A[0]+rotmtx[2][1]*point_A[1]+rotmtx[2][2]*point_A[2];//+rotmtx[2][3];       

       point_A[0] = pointB[i].cc[0]- axisA[0];
       point_A[1] = pointB[i].cc[1]- axisA[1];
       point_A[2] = pointB[i].cc[2]- axisA[2];
       pointB[i].cc = axisA;
       pointB[i].cc[0] += rotmtx[0][0]*point_A[0]+rotmtx[0][1]*point_A[1]+rotmtx[0][2]*point_A[2];//+rotmtx[0][3];
       pointB[i].cc[1] += rotmtx[1][0]*point_A[0]+rotmtx[1][1]*point_A[1]+rotmtx[1][2]*point_A[2];//+rotmtx[1][3];
       pointB[i].cc[2] += rotmtx[2][0]*point_A[0]+rotmtx[2][1]*point_A[1]+rotmtx[2][2]*point_A[2];//+rotmtx[2][3];       

       point_A[0] = pointB[i].co[0]- axisA[0];
       point_A[1] = pointB[i].co[1]- axisA[1];
       point_A[2] = pointB[i].co[2]- axisA[2];
       pointB[i].co = axisA;
       pointB[i].co[0] += rotmtx[0][0]*point_A[0]+rotmtx[0][1]*point_A[1]+rotmtx[0][2]*point_A[2];//+rotmtx[0][3];
       pointB[i].co[1] += rotmtx[1][0]*point_A[0]+rotmtx[1][1]*point_A[1]+rotmtx[1][2]*point_A[2];//+rotmtx[1][3];
       pointB[i].co[2] += rotmtx[2][0]*point_A[0]+rotmtx[2][1]*point_A[1]+rotmtx[2][2]*point_A[2];//+rotmtx[2][3];       
   }
   // rorate 0 point in C O
    point_A[0] = pointB[0].cc[0]- axisA[0];
    point_A[1] = pointB[0].cc[1]- axisA[1];
    point_A[2] = pointB[0].cc[2]- axisA[2];
    pointB[0].cc = axisA;
    pointB[0].cc[0] += rotmtx[0][0]*point_A[0]+rotmtx[0][1]*point_A[1]+rotmtx[0][2]*point_A[2];//+rotmtx[0][3];
    pointB[0].cc[1] += rotmtx[1][0]*point_A[0]+rotmtx[1][1]*point_A[1]+rotmtx[1][2]*point_A[2];//+rotmtx[1][3];
    pointB[0].cc[2] += rotmtx[2][0]*point_A[0]+rotmtx[2][1]*point_A[1]+rotmtx[2][2]*point_A[2];//+rotmtx[2][3];

    point_A[0] = pointB[0].co[0]- axisA[0];
    point_A[1] = pointB[0].co[1]- axisA[1];
    point_A[2] = pointB[0].co[2]- axisA[2];
    pointB[0].co = axisA;
    pointB[0].co[0] += rotmtx[0][0]*point_A[0]+rotmtx[0][1]*point_A[1]+rotmtx[0][2]*point_A[2];//+rotmtx[0][3];
    pointB[0].co[1] += rotmtx[1][0]*point_A[0]+rotmtx[1][1]*point_A[1]+rotmtx[1][2]*point_A[2];//+rotmtx[1][3];
    pointB[0].co[2] += rotmtx[2][0]*point_A[0]+rotmtx[2][1]*point_A[1]+rotmtx[2][2]*point_A[2];//+rotmtx[2][3];    

    //  rorate finally(pointB.size()-1) point in N
    point_A[0] = pointB[pointB.size()-1].cn[0]- axisA[0];
    point_A[1] = pointB[pointB.size()-1].cn[1]- axisA[1];
    point_A[2] = pointB[pointB.size()-1].cn[2]- axisA[2];
    pointB[pointB.size()-1].cn = axisA;
    pointB[pointB.size()-1].cn[0] += rotmtx[0][0]*point_A[0]+rotmtx[0][1]*point_A[1]+rotmtx[0][2]*point_A[2];//+rotmtx[0][3];
    pointB[pointB.size()-1].cn[1] += rotmtx[1][0]*point_A[0]+rotmtx[1][1]*point_A[1]+rotmtx[1][2]*point_A[2];//+rotmtx[1][3];
    pointB[pointB.size()-1].cn[2] += rotmtx[2][0]*point_A[0]+rotmtx[2][1]*point_A[1]+rotmtx[2][2]*point_A[2];//+rotmtx[2][3];      
   return true;
}

// rotate the amount of point beginning from index 
bool GroupRotationx(const vector<float> &axisA, const vector<float> &axisB, float angle, vector<Rot> &pointB, int index)
{
   if(axisA.size()!=3 || axisB.size()!=3)
   {
      cout<<"Error in GroupRotation()"<<endl;
      return false;
   }
   
   vector<float> axis(3, 0);
   axis[0]=axisB[0]-axisA[0]; 
   axis[1]=axisB[1]-axisA[1];
   axis[2]=axisB[2]-axisA[2]; 
   
   matrix rotmtx;
   if(!RotationMatrixB(axis, deg2rad(angle), rotmtx)) return false;
   
   //vector<float> point_A(4, 0);
   float point_A[3];
   for(int i=index;i<pointB.size();i++)
   {
      //rotate side chain
      for(int j=0;j<pointB[i].res_atm.size();j++)
      {
         point_A[0] = pointB[i].res_atm[j][0]- axisA[0];
         point_A[1] = pointB[i].res_atm[j][1]- axisA[1];
         point_A[2] = pointB[i].res_atm[j][2]- axisA[2];

         pointB[i].res_atm[j] = axisA;
         pointB[i].res_atm[j][0] += rotmtx[0][0]*point_A[0]+rotmtx[0][1]*point_A[1]+rotmtx[0][2]*point_A[2]+rotmtx[0][3];
         pointB[i].res_atm[j][1] += rotmtx[1][0]*point_A[0]+rotmtx[1][1]*point_A[1]+rotmtx[1][2]*point_A[2]+rotmtx[1][3];
         pointB[i].res_atm[j][2] += rotmtx[2][0]*point_A[0]+rotmtx[2][1]*point_A[1]+rotmtx[2][2]*point_A[2]+rotmtx[2][3];
      }
      // rotate N, CA,C,O
       point_A[0] = pointB[i].cn[0]- axisA[0];
       point_A[1] = pointB[i].cn[1]- axisA[1];
       point_A[2] = pointB[i].cn[2]- axisA[2];
       pointB[i].cn = axisA;
       pointB[i].cn[0] += rotmtx[0][0]*point_A[0]+rotmtx[0][1]*point_A[1]+rotmtx[0][2]*point_A[2]+rotmtx[0][3];
       pointB[i].cn[1] += rotmtx[1][0]*point_A[0]+rotmtx[1][1]*point_A[1]+rotmtx[1][2]*point_A[2]+rotmtx[1][3];
       pointB[i].cn[2] += rotmtx[2][0]*point_A[0]+rotmtx[2][1]*point_A[1]+rotmtx[2][2]*point_A[2]+rotmtx[2][3];    

       point_A[0] = pointB[i].ca[0]- axisA[0];
       point_A[1] = pointB[i].ca[1]- axisA[1];
       point_A[2] = pointB[i].ca[2]- axisA[2];
       pointB[i].ca = axisA;
       pointB[i].ca[0] += rotmtx[0][0]*point_A[0]+rotmtx[0][1]*point_A[1]+rotmtx[0][2]*point_A[2]+rotmtx[0][3];
       pointB[i].ca[1] += rotmtx[1][0]*point_A[0]+rotmtx[1][1]*point_A[1]+rotmtx[1][2]*point_A[2]+rotmtx[1][3];
       pointB[i].ca[2] += rotmtx[2][0]*point_A[0]+rotmtx[2][1]*point_A[1]+rotmtx[2][2]*point_A[2]+rotmtx[2][3];       

       point_A[0] = pointB[i].cc[0]- axisA[0];
       point_A[1] = pointB[i].cc[1]- axisA[1];
       point_A[2] = pointB[i].cc[2]- axisA[2];
       pointB[i].cc = axisA;
       pointB[i].cc[0] += rotmtx[0][0]*point_A[0]+rotmtx[0][1]*point_A[1]+rotmtx[0][2]*point_A[2]+rotmtx[0][3];
       pointB[i].cc[1] += rotmtx[1][0]*point_A[0]+rotmtx[1][1]*point_A[1]+rotmtx[1][2]*point_A[2]+rotmtx[1][3];
       pointB[i].cc[2] += rotmtx[2][0]*point_A[0]+rotmtx[2][1]*point_A[1]+rotmtx[2][2]*point_A[2]+rotmtx[2][3];       

       point_A[0] = pointB[i].co[0]- axisA[0];
       point_A[1] = pointB[i].co[1]- axisA[1];
       point_A[2] = pointB[i].co[2]- axisA[2];
       pointB[i].co = axisA;
       pointB[i].co[0] += rotmtx[0][0]*point_A[0]+rotmtx[0][1]*point_A[1]+rotmtx[0][2]*point_A[2]+rotmtx[0][3];
       pointB[i].co[1] += rotmtx[1][0]*point_A[0]+rotmtx[1][1]*point_A[1]+rotmtx[1][2]*point_A[2]+rotmtx[1][3];
       pointB[i].co[2] += rotmtx[2][0]*point_A[0]+rotmtx[2][1]*point_A[1]+rotmtx[2][2]*point_A[2]+rotmtx[2][3];       
   }
   return true;
}
bool GroupTranslationx(const vector<float> &trans, vector<Rot> &pointB)
{
   if(trans.size()!=3)
   {
      cout<<"Error in GroupTranslation\n";
      return false;
   }
   for(int i=0;i<pointB.size();i++)
   {
        // translate sidechain
        for(int j=0;j<pointB[i].res_atm.size();j++)
        {
            pointB[i].res_atm[j][0] += trans[0];
            pointB[i].res_atm[j][1] += trans[1];
            pointB[i].res_atm[j][2] += trans[2];
        }
        // translate N, CA,C,O
       pointB[i].cn[0] += trans[0];
       pointB[i].cn[1] += trans[1];
       pointB[i].cn[2] += trans[2];    

       pointB[i].ca[0] += trans[0];
       pointB[i].ca[1] += trans[1];
       pointB[i].ca[2] += trans[2];       

       pointB[i].cc[0] += trans[0];
       pointB[i].cc[1] += trans[1];
       pointB[i].cc[2] += trans[2];       

       pointB[i].co[0] += trans[0];
       pointB[i].co[1] += trans[1];
       pointB[i].co[2] += trans[2]; 
   }
   return true;
}
bool GroupTranslationx(const vector<float> &trans, vector<Rot> &pointB, int index)
{
   if(trans.size()!=3)
   {
      cout<<"Error in GroupTranslation\n";
      return false;
   }
   for(int i=index;i<pointB.size();i++)
   {
       // translate sidechain
       for(int j=0;j<pointB[i].res_atm.size();j++)
       {
           pointB[i].res_atm[j][0] += trans[0];
           pointB[i].res_atm[j][1] += trans[1];
           pointB[i].res_atm[j][2] += trans[2];
       }
        // translate N, CA,C,O
       pointB[i].cn[0] += trans[0];
       pointB[i].cn[1] += trans[1];
       pointB[i].cn[2] += trans[2];    

       pointB[i].ca[0] += trans[0];
       pointB[i].ca[1] += trans[1];
       pointB[i].ca[2] += trans[2];       

       pointB[i].cc[0] += trans[0];
       pointB[i].cc[1] += trans[1];
       pointB[i].cc[2] += trans[2];       

       pointB[i].co[0] += trans[0];
       pointB[i].co[1] += trans[1];
       pointB[i].co[2] += trans[2]; 
   }
   return true;
}

bool GroupRotation1(const vector<float> &axisA, const vector<float> &axisB, float angle, vector<vector<float> > &pointB)
{
   if(axisA.size()!=3 || axisB.size()!=3)
   {
      cout<<"Error in GroupRotation()"<<endl;
      return false;
   }
   
   vector<float> axis(3, 0);
   axis[0]=axisB[0]-axisA[0]; 
   axis[1]=axisB[1]-axisA[1];
   axis[2]=axisB[2]-axisA[2]; 
   
   matrix rotmtx;
   if(!RotationMatrixB(axis, deg2rad(angle), rotmtx)) return false;
   
   //vector<float> point_A(4, 0);
   float point_A[3];
   for(int i=0; i<pointB.size(); ++i)
   {
      point_A[0]=pointB[i][0]-axisA[0];
      point_A[1]=pointB[i][1]-axisA[1];
      point_A[2]=pointB[i][2]-axisA[2];

      pointB[i]=axisA;
      pointB[i][0]+=rotmtx[0][0]*point_A[0]+rotmtx[0][1]*point_A[1]+rotmtx[0][2]*point_A[2]+rotmtx[0][3];
      pointB[i][1]+=rotmtx[1][0]*point_A[0]+rotmtx[1][1]*point_A[1]+rotmtx[1][2]*point_A[2]+rotmtx[1][3];
      pointB[i][2]+=rotmtx[2][0]*point_A[0]+rotmtx[2][1]*point_A[1]+rotmtx[2][2]*point_A[2]+rotmtx[2][3];
      /*
      point_A[0]=pointB[i][0]-axisA[0];
      point_A[1]=pointB[i][1]-axisA[1];
      point_A[2]=pointB[i][2]-axisA[2];
      point_A[3]=1;
      
      if(!MatrixTimesTransVector(rotmtx, point_A, pointB[i]))
      { 
         return false;
      }
      pointB[i].pop_back();
   
      pointB[i][0]+=axisA[0];
      pointB[i][1]+=axisA[1];
      pointB[i][2]+=axisA[2];
      */
   }
   return true;
}

bool GroupRotationt(const vector<float> &axisA, const vector<float> &axisB, float angle, vector<float> &pointB)
{
   if(axisA.size()!=3 || axisB.size()!=3)
   {
      cout<<"Error in GroupRotation()"<<endl;
      return false;
   }
   
   vector<float> axis(3, 0);
   axis[0]=axisB[0]-axisA[0]; 
   axis[1]=axisB[1]-axisA[1];
   axis[2]=axisB[2]-axisA[2]; 
   
   matrix rotmtx;
   if(!RotationMatrixB(axis, deg2rad(angle), rotmtx)) return false;
   
   //vector<float> point_A(4, 0);
    float point_A[3];
    point_A[0]=pointB[0]-axisA[0];
    point_A[1]=pointB[1]-axisA[1];
    point_A[2]=pointB[2]-axisA[2];

    pointB=axisA;
    pointB[0]+=rotmtx[0][0]*point_A[0]+rotmtx[0][1]*point_A[1]+rotmtx[0][2]*point_A[2]+rotmtx[0][3];
    pointB[1]+=rotmtx[1][0]*point_A[0]+rotmtx[1][1]*point_A[1]+rotmtx[1][2]*point_A[2]+rotmtx[1][3];
    pointB[2]+=rotmtx[2][0]*point_A[0]+rotmtx[2][1]*point_A[1]+rotmtx[2][2]*point_A[2]+rotmtx[2][3];
      /*
      point_A[0]=pointB[i][0]-axisA[0];
      point_A[1]=pointB[i][1]-axisA[1];
      point_A[2]=pointB[i][2]-axisA[2];
      point_A[3]=1;
      
      if(!MatrixTimesTransVector(rotmtx, point_A, pointB[i]))
      { 
         return false;
      }
      pointB[i].pop_back();
   
      pointB[i][0]+=axisA[0];
      pointB[i][1]+=axisA[1];
      pointB[i][2]+=axisA[2];
      */
   return true;
}

bool GroupRotationid(const vector<float> &axisA, const vector<float> &axisB, float angle, vector<vector<float> > &pointB, int index)
{
   if(axisA.size()!=3 || axisB.size()!=3)
   {
      cout<<"Error in GroupRotation()"<<endl;
      return false;
   }
   
   vector<float> axis(3, 0);
   axis[0]=axisB[0]-axisA[0]; 
   axis[1]=axisB[1]-axisA[1];
   axis[2]=axisB[2]-axisA[2]; 
   
   matrix rotmtx;
   if(!RotationMatrixB(axis, deg2rad(angle), rotmtx)) return false;
   
   float point_A[3];
   int i, j, m, n;
   for(j=index; j<pointB.size(); j++)
   {
      i=j;
      point_A[0]=pointB[i][0]-axisA[0];
      point_A[1]=pointB[i][1]-axisA[1];
      point_A[2]=pointB[i][2]-axisA[2];

      pointB[i]=axisA;
      pointB[i][0]+=rotmtx[0][0]*point_A[0]+rotmtx[0][1]*point_A[1]+rotmtx[0][2]*point_A[2]+rotmtx[0][3];
      pointB[i][1]+=rotmtx[1][0]*point_A[0]+rotmtx[1][1]*point_A[1]+rotmtx[1][2]*point_A[2]+rotmtx[1][3];
      pointB[i][2]+=rotmtx[2][0]*point_A[0]+rotmtx[2][1]*point_A[1]+rotmtx[2][2]*point_A[2]+rotmtx[2][3];
   
   }
   return true;
}

bool GroupRotationp(const vector<float> &axisA, const vector<float> &axisB, float angle, vector<poseCoord> &pointB)
{
   if(axisA.size()!=3 || axisB.size()!=3)
   {
      cout<<"Error in GroupRotation()"<<endl;
      return false;
   }
   
   vector<float> axis(3, 0);
   axis[0]=axisB[0]-axisA[0]; 
   axis[1]=axisB[1]-axisA[1];
   axis[2]=axisB[2]-axisA[2]; 
   
   matrix rotmtx;
   if(!RotationMatrixB(axis, deg2rad(angle), rotmtx)) return false;
   
   //vector<float> point_A(4, 0);
   float point_A[3];
   for(int i=0; i<pointB.size(); ++i)
   {
      point_A[0]=pointB[i].x_[0]-axisA[0];
      point_A[1]=pointB[i].x_[1]-axisA[1];
      point_A[2]=pointB[i].x_[2]-axisA[2];

      pointB[i].x_=axisA;
      pointB[i].x_[0]+=rotmtx[0][0]*point_A[0]+rotmtx[0][1]*point_A[1]+rotmtx[0][2]*point_A[2]+rotmtx[0][3];
      pointB[i].x_[1]+=rotmtx[1][0]*point_A[0]+rotmtx[1][1]*point_A[1]+rotmtx[1][2]*point_A[2]+rotmtx[1][3];
      pointB[i].x_[2]+=rotmtx[2][0]*point_A[0]+rotmtx[2][1]*point_A[1]+rotmtx[2][2]*point_A[2]+rotmtx[2][3];
      /*
      point_A[0]=pointB[i][0]-axisA[0];
      point_A[1]=pointB[i][1]-axisA[1];
      point_A[2]=pointB[i][2]-axisA[2];
      point_A[3]=1;
      
      if(!MatrixTimesTransVector(rotmtx, point_A, pointB[i]))
      { 
         return false;
      }
      pointB[i].pop_back();
   
      pointB[i][0]+=axisA[0];
      pointB[i][1]+=axisA[1];
      pointB[i][2]+=axisA[2];
      */
   }
   return true;
}

/*
int main(int argc, char** argv)
{
        if(argc<3)
        {
                cout<<"the argument is too smaller"<<endl;
        }
        else
        {
		Model protein;
                char inputfile[30];
                char outputfile[30];
                strcpy(inputfile,argv[1]);
                strcpy(outputfile,argv[2]);
		readpdbstructure(inputfile,protein);
//		cout<<protein.chains.size()<<"  "<<protein.chains[0].residues.size();//<<"  "<<protein.chains[0].residues[0].atoms.size();
		writePDBStructure(outputfile,protein);
		
   //             coordtype coor;//cood[0] the string before the ATOM,cood[1]--cood[3] ATOM coordinate,coor[4]
   //             char inputfile[30];
   //             char outputfile[30];
   //             strcpy(inputfile,argv[1]);
   //             strcpy(outputfile,argv[2]);
   //             readcoordinatefrompdb(inputfile,coor);
   //             rotatepdb(coor);
   //             writePDB(outputfile,coor); 
        }
        return 0;
} */

