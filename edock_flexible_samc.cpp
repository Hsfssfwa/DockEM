#include <sstream> 
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <random>

#include "readligand.cpp"
#include "output.h"
#include "match.h"
#include "energy.h"
#include "translate_rotate.h"
#include "MC_simulation.h"
#include "REMC.h"
#include "precon.h"
#include "initialmol.h"
#include "flex_dock.h"
#include "generaterandnum.h"
#include "CaDensityx.h"
using namespace std;

/*example:
docking outfile tmin tmax grid sphere 
edock_MC re1 2 100
*/


int main(int argc, char ** argv){

    random_device rd;
    mt19937 gen(rd());
    //uniform_real_distribution<> randgeneration(-1.0, 1.0);
    //uniform_real_distribution<> rand0to1(0, 1.0);
    DOCKMol mol;
    DOCKMol receptor;
    Orient c_orient;
    AMBER_TYPER amber;
    Energy_Score c_nrg;
    
    //const char * ligand_in = "lig_charge.mol2";
    //const char * receptor_mol = "af2_protein.mol2";
    //string bindingsitesfile="./bindingsite_pre.dat";
    
    char * receptor_mol = argv[1];
    char * ligand_in = argv[2];
    string bindingsitesfile = argv[3];
  
    string outfile=string(argv[4]); //docking results 
    // dock6 grid app, save the protein into the grid file for calcuating the energy
    string grid_file = argv[5]; // dock6 grid app, save the protein into the grid file for calcuating the energy
    char * receptor_in = argv[6]; // is a .sph file, for binding sites
    //
    /***************************************EM*************************************************/
    string inputMRC = argv[7];
    float resolution = atof(argv[8]);
    float sampleing = atof(argv[9]);

    ElectronDensity theDensityMap;
    theDensityMap.preso = resolution;
    theDensityMap.pATOM_MASK = 7.0;
    //    theDensityMap.pATOM_MASK =7.0;
    theDensityMap.pCA_MASK = 7.0;
    theDensityMap.pforce_apix_on_map_load_ = 0.0;
    theDensityMap.pWINDOW_ = 1;
    theDensityMap.pscore_window_context_ = false;
    theDensityMap.premap_symm_ = false;
    theDensityMap.pde_edensity = false;
    theDensityMap.pnkbins_ = 0;

    theDensityMap.readMRCandResize(inputMRC, resolution, sampleing);
    cout << "finished load map" << endl;
    /*******************************************************************************************/

    //char * smin=argv[4]; // min tempreture for remc
    //char * smax=argv[5]; // max tempreture for remc
    //char * vdw_weight_c = argv[6]; //user input vdw weight
    //char * bindsite_weight_c = argv[7];//user input binding sites weight
    //float Tmin=atof(smin);
    //float Tmax=atof(smax);
    float Tmin=1.0;
    float Tmax=60.0;
    float cutoff = 100.0; // total energy cutoff
    //float vdw_weight = atof(vdw_weight_c);
    float vdw_weight = 0.001;
    //float  bindsite_weight = atof(bindsite_weight_c);
    float  bindsite_weight = 1;
    c_nrg.vdw_scale = vdw_weight;
    int initialnum=40; //REMC number
    int swapnumber = 100;
    int MC_steps = 201;
    bool flexible_flag =1; //whethe do flexible docking
    //bool flexible_flag = 0;
    //string filevdw = "/home/zhangwy/program/dock6/parameters/vdw_AMBER_parm99.defn";
    //string fileflex = "/home/zhangwy/program/dock6/parameters/flex.defn";
    //string fileflex_drive_tbl = "/home/zhangwy/program/dock6/parameters/flex_drive.tbl";
    
    string filevdw = "../../parameters/vdw_AMBER_parm99.defn";
    string fileflex = "../../parameters/flex.defn";
    string fileflex_drive_tbl = "../../parameters/flex_drive.tbl";
    //$$$$$$$$$$$$$$$$$$$$read the BioLip_ligand_torsion_lib$$$$$$$$$$$$$$$$$
    //ifstream biolip_file("/home/zhangwy/program/EDock2019/edock_flex/BioLip_ligand_torsion_lib");
    ifstream biolip_file("../../parameters/BioLip_ligand_torsion_lib");

    /*******************************************************************************************/
    //read the ligand informaiton about the title, mol_info_line, comment1, comment2, comment3,     
    //score_text_data,energy, simplex_text, mol_data, num_atoms, num_bonds, num_residues, and so on.
    //and compulate the ring information
    
    
    cout<<"********read ligand part*****************"<<endl;
    readlig(mol,ligand_in);
    readlig(receptor, receptor_mol);

    //theDensityMap.diffmap(receptor);
    //cout << theDensityMap.cellDimensions[0] << endl;
    //theDensityMap.writeMRC("diffmap.mrc");
    //cout << "finish diffmap" << endl;
    //float init_CC = 10000 * (1.0 - theDensityMap.matchpose(receptor,mol));
    /*******************************************************************************************/
    //this part read some template ligand dat file, and calculte the uij and stdij for template ligands 
    //add by wenyizhang 04252019
    
    cout<<"********read binding sites*****************"<<endl;
    
    vector<vector <int> > binding_sites;
    read_binding_sites(bindingsitesfile.c_str(),binding_sites);
    cout<<"********read binding sites finish*****************"<<endl;

    int ligand_non_H= 0;
    for(int i=0;i<mol.num_atoms;i++)
    {
        if(mol.atom_types[i] != "H")
            ligand_non_H ++; 
    }
    cout<<"ligand_non_H"<<ligand_non_H<<endl;
    int protein_atom_num=0;
    vector<int> protein_atom_index;
    vector<float> cartbinding(3, 0.0);
    vector<float> fartbinding(3, 0.0);
    for (int i=0;i<binding_sites[0].size();i++)
    {
        for(int j=0;j<receptor.num_atoms;j++)
        {
            if(atoi(receptor.atom_residue_numbers[j].c_str()) == binding_sites[0][i])
            {
                protein_atom_num++;
                protein_atom_index.push_back(j);

                cartbinding[0] += receptor.x[j];
                cartbinding[1] += receptor.y[j];
                cartbinding[2] += receptor.z[j];
            }   
        }
    }
    cout<<"protein_atom_num"<<protein_atom_num<<endl;
    cout<<protein_atom_index[0]<<endl;
    //save corresponding all atoms in bindingsite**********************
    
    
    cartbinding[0] = cartbinding[0]/ protein_atom_num;
    cartbinding[1] = cartbinding[0]/ protein_atom_num;
    cartbinding[2] = cartbinding[0]/ protein_atom_num;
    float car_x = 0, car_y = 0, car_z = 0;
    float max = 0;
    for (int i = 0; i < mol.num_atoms; i++)
    {
        for (int j = i; j < mol.num_atoms; j++)
        {
            car_x = (mol.x[i] - mol.x[j]) * (mol.x[i] - mol.x[j]);
            car_y = (mol.y[i] - mol.y[j]) * (mol.y[i] - mol.y[j]);
            car_z = (mol.z[i] - mol.z[j]) * (mol.z[i] - mol.z[j]);
            if (sqrt(car_x + car_x + car_x) > max)
                max = sqrt(car_x + car_x + car_x);
        }
    }
    cout << "max:  " << max << endl;
    theDensityMap.cart(cartbinding, fartbinding,max);


    
    cout<<"********amber initialize****************"<<endl;
    amber.initialize(filevdw.c_str(), fileflex.c_str(),fileflex_drive_tbl.c_str());
    cout<<"*******amber.prepare_molecule***********"<<endl;
    // read the vdw_AMBER_parm99.defn, save in latom
    
    amber.prepare_molecule(mol);// prepare the ligand mol for the vdw
    amber.prepare_molecule(receptor); //prepare the receptor mol for the vdw
    mol.prepare_molecule();//neighbors informations add by wenyi
    //receptor.prepare_molecule();//

    //flexible part: bond checking
    vector<TORSION> torsions;
    FLOATVec  vertex;//this vector I don't use     
    id_torsions(mol, vertex, torsions); 
    cout<<"torsions.size is"<<torsions.size()<<endl;
    if (torsions.size() == 0) flexible_flag = 0;
    cout<<"flexible_flag:"<<flexible_flag<<endl;
    cout<<"ligand number "<<mol.num_atoms<<endl;
    
    if(mol.num_atoms>40 && torsions.size()>10)
    {
        initialnum=20; //REMC number
        swapnumber = 200;
        MC_steps = 201; 
    }
    
    c_nrg.grid_file_name = grid_file.c_str();
    c_nrg.initialize(amber);
    //native energy
  
    /*******************************************************************************************/
    //this part read the BioLip_ligand_torsion_lib, and calculte the bond energy 
    
    int torsiontypes[8]={1,2,3,32,6,62,63,122};
    int row=0;
    int tempnumber=0;
    //int biolip_matrix[8][361]={0};
    vector<vector<int>> biolip_matrix(8,vector<int>(361));
    string pattern = "\t";
    for (string s; getline(biolip_file, s);)
    {
      vector< string> tortimes = split(s, pattern);
      for(int i=0;i<tortimes.size();i++)
      {
        tempnumber = atoi(tortimes[i].c_str());
        biolip_matrix[row][i]=tempnumber;
      }
	    row=row+1;  
		
    }
    //cout<<"row number "<<row<<endl;
    //cout<<biolip_matrix[0][1]*1.0<<"^^^"<<biolip_matrix[0][0]*1.0<<"&&&&"<<-(log((biolip_matrix[0][1]*1.0)/(biolip_matrix[0][0]*1.0)))<<endl;
    //**************************************************   
    if(flexible_flag == 1)
    {
      if (MC_flexible_ligand(amber, mol, 1.0,gen,biolip_matrix))
          {
            cout<<"before matching, flexible ligand finished"<<endl;
            //string mol2file="sampling_" + outfile + ".mol2";
            //ofstream mol2(mol2file.c_str());
            //stringstream ss;
            //ss << 1;
            //string str=ss.str();
            //Write_Mol2(mol, mol2, str);
            //mol2.close();
          }
    }
    /*******************************************************************************************/
    vector<float> centers(3, 0);
    double t[3];
    for(int i=0;i<mol.num_atoms;i++)
    {
      centers[0]+=mol.x[i];
      centers[1]+=mol.y[i];
      centers[2]+=mol.z[i];
    }
    centers[0]=centers[0]/mol.num_atoms;
    centers[1]=centers[1]/mol.num_atoms;
    centers[2]=centers[2]/mol.num_atoms;
  
    /*******************************************************************************************/
    // temp definition add by wenyizng 12/04/2019
    //vector<int> protein_atom_index;
    vector<vector <float> > ave;
    vector<vector <float> > std;
    vector<vector< vector<float> > > distance;
    cout<<"*******match part***********************"<<endl;
    c_orient.match(mol,receptor_in);    
    cout<<"orient fininshed!!!!!!!"<<endl;
    //cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<endl;
    //cout<<"spheres.size()"<<c_orient.spheres.size()<<endl;
    //finding docking center
    float spherecenterx,spherecentery,spherecenterz;
    spherecenterx=spherecentery=spherecenterz=0;
    vector<float> sphx;
    vector<float> sphy;
    vector<float> sphz; 
    
    for(int i=0;i<c_orient.spheres.size();i++)
    {
        //cout<<c_orient.spheres[i].crds.x<<"\t"<<c_orient.spheres[i].crds.y<<"\t"<<c_orient.spheres[i].crds.z<<endl;
      spherecenterx += c_orient.spheres[i].crds.x;
      spherecentery += c_orient.spheres[i].crds.y;
      spherecenterz += c_orient.spheres[i].crds.z;
      sphx.push_back(c_orient.spheres[i].crds.x);
      sphy.push_back(c_orient.spheres[i].crds.y);
      sphz.push_back(c_orient.spheres[i].crds.z);
    }
    spherecenterx=spherecenterx/c_orient.spheres.size();
    spherecentery=spherecentery/c_orient.spheres.size();
    spherecenterz=spherecenterz/c_orient.spheres.size();
    //cout<<spherecenterx<<"\t"<<spherecentery<<"\t"<<spherecenterz<<endl;
    
    vector <DOCKMol> molvec;
    vector <DOCKMol> matchvec;
    /*c_orient.cliques saves the match mols */
    //cout<<c_orient.cliques.size()<<endl;
    int matching_count=0;
    int allmatchcount=0;
    
    cout<<"c_orient.cliques.size() is:"<<c_orient.cliques.size()<<endl;
    //
    int orient_iterations=0;
    if(c_orient.cliques.size() < 5000 && c_orient.cliques.size()>0) //it means that the mol can be oriented
	     orient_iterations=c_orient.cliques.size();
    else if(c_orient.cliques.size() >= 5000)
	     orient_iterations=5000;
    for(c_orient.current_clique = 0; c_orient.current_clique<orient_iterations;c_orient.current_clique++) 
    {
	   //cout<<"c_orient.current_clique is:"<<c_orient.current_clique<<endl;
      if(c_orient.new_next_orientation(mol))	
	  {
          //cout << "aaa" << endl;
        if((energy1(amber, mol, receptor,c_nrg,filevdw.c_str(), fileflex.c_str(),fileflex_drive_tbl.c_str(),grid_file,cutoff,biolip_matrix,sphx,sphy,sphz,bindsite_weight,protein_atom_index,ave,std,theDensityMap)))
        {
		      //cout<<"score is:"<<mol.intral_energy<<endl;
		      //if (mol.current_score < 1e+6)
			  if (mol.current_score < 1e+6)
		      //if (mol.intral_energy <1e+6)
		      {
		        cout <<"check me 000000000   "<<mol.current_score<<"   "<<mol.current_EM_score<<endl;
		        molvec.push_back(mol); 
		        if(molvec.size()>1000)
			      break;
		        // molvec save the matching conformations
		        //cout<<molvec[molvec.size()-1].current_score<<endl;
		      }    
        }
      }      
    }
    cout<<"checking666"<<molvec.size()<<endl;  
    bool REMC_flag;
    if(molvec.size()>=initialnum)
      REMC_flag=1;
    else if(molvec.size()==0)
    {
        int num_molvec=molvec.size();
        cout<<"num_molvec "<<num_molvec<<endl;
        vector<vector<float> > ini_x(initialnum-num_molvec, vector<float> (mol.num_atoms,0.0));
        vector<vector<float> > ini_y(initialnum-num_molvec, vector<float> (mol.num_atoms,0.0));
        vector<vector<float> > ini_z(initialnum-num_molvec, vector<float> (mol.num_atoms,0.0));
        REMC_flag = intialmol2(ini_x,ini_y,ini_z,num_molvec, mol, spherecenterx, spherecentery, spherecenterz, initialnum, gen);
        //cout<<"ini_x len"<<ini_x.size()<<endl;
        for(int i=0;i<initialnum-num_molvec;i++)
        {
            DOCKMol moltemp;
            copy_molecule(moltemp,mol);
            for(int atom=0;atom<mol.num_atoms;atom++)
            {
                moltemp.x[atom]=ini_x[i][atom];
                moltemp.y[atom]=ini_y[i][atom];
                moltemp.z[atom]=ini_z[i][atom];
                
            }
            //cout<<ini_x[i][0]<<"\t"<<moltemp.x[0]<<"\t";
            //cout<<molvec[molvec.size()-1].x[0]<<endl;

            if (!(energy(amber, moltemp, receptor, c_nrg,filevdw.c_str(), fileflex.c_str(),fileflex_drive_tbl.c_str(),grid_file,cutoff,biolip_matrix,sphx,sphy,sphz,bindsite_weight,protein_atom_index,ave,std, theDensityMap)))
            {
                moltemp.current_score=1e+6;
                moltemp.internal_energy=1e+6;
                moltemp.intral_energy=1e+6;
            }
            molvec.push_back(moltemp);
            //cout<<molvec[molvec.size()-1].x[0]<<endl;
       }    
   }
   else if(molvec.size()<initialnum && molvec.size()>0)
   {
       int num_molvec=molvec.size();
       cout<<"num_molvec "<<num_molvec<<endl;
       vector<float> axis1(3, 0);
       vector<float> axis2(3, 0);
       double a,b,c;
       for (int i=0;i<initialnum-num_molvec;i++)
       {   
	        int index=i%num_molvec;
          //cout<<i<<"\t"<<index<<endl;
          for(int m=0;m<mol.num_atoms;m++)
          {
            mol.x[m]=molvec[index].x[m];
            mol.y[m]=molvec[index].y[m];
            mol.z[m]=molvec[index].z[m];
          }
          DOCKMol moltemp;
          copy_molecule(moltemp,mol);
          molvec.push_back(moltemp);
          //cout<<"now molvec.size()"<<molvec.size()<<endl;
          //cout<<"before mol"<<molvec[0].x[0]<<"\t"<<molvec[0].y[0]<<"\t"<<molvec[0].z[0]<<endl;

          {
            float ligcenterx,ligcentery,ligcenterz;
            ligcenterx=ligcentery=ligcenterz=0;
            for(int m=0;m<mol.num_atoms;m++)
            {
             ligcenterx +=mol.x[m];
             ligcentery +=mol.y[m];
             ligcenterz +=mol.z[m];
            }
            ligcenterx=ligcenterx/mol.num_atoms;
            ligcentery=ligcentery/mol.num_atoms;
            ligcenterz=ligcenterz/mol.num_atoms;
            //cout<<ligcenterx<<"\t"<<ligcentery<<"\t"<<ligcenterz<<endl;
             
             axis1[0]=ligcenterx; axis1[1]=ligcentery; axis1[2]=ligcenterz;
             a = rand0to1(gen);
             b = rand0to1(gen);
             //cout<<"a"<<a<<" "<<"b"<<b<<" "<<"c"<<c<<endl;
             double theta = 2 * 3.1415926 * a;
             double phi = acos(2*b-1.0);
             axis2[0] = sin(phi) * cos(theta) + ligcenterx;
             axis2[1] = sin(phi) * sin(theta) + ligcentery;
             axis2[2]= cos(phi) + ligcenterz;
             c=randgeneration(gen);
             float angle = c*180.0;
             
             GroupRotation(axis2,axis1,angle,mol); // random rotate the mol 
             //cout<<"checking5555"<<mol.x[0]<<endl;
             //translate the center, some situation is can not calculate the energy
             double t[3];
             t[0]=spherecenterx-ligcenterx;
             t[1]=spherecentery-ligcentery;
             t[2]=spherecenterz-ligcenterz;
             for(int atom=0;atom<mol.num_atoms;atom++)
             {
                mol.x[atom] +=t[0];
                mol.y[atom] +=t[1];
                mol.z[atom] +=t[2];
             }
             //cout<<"after mol"<<mol.x[0]<<"\t"<<mol.y[0]<<"\t"<<mol.z[0]<<endl;
             if (!(energy(amber, mol, receptor,c_nrg,filevdw.c_str(), fileflex.c_str(),fileflex_drive_tbl.c_str(),grid_file,cutoff,biolip_matrix,sphx,sphy,sphz,bindsite_weight,protein_atom_index,ave,std, theDensityMap)))
             {
                 mol.current_score=1e+6;
                 mol.internal_energy=1e+6;
                 mol.intral_energy=1e+6;
             }
             for(int m=0;m<mol.num_atoms;m++)
             {
                molvec[molvec.size()-1].x[m] = mol.x[m];
                molvec[molvec.size()-1].y[m] = mol.y[m];
                molvec[molvec.size()-1].z[m] = mol.z[m];

             }
             molvec[molvec.size()-1].current_score=mol.current_score;
             molvec[molvec.size()-1].internal_energy=mol.internal_energy;
             molvec[molvec.size()-1].intral_energy = mol.intral_energy;
             molvec[molvec.size()-1].energy_tor_total = mol.energy_tor_total;
             //cout<<"after molvec"<<molvec[molvec.size()-1].x[0]<<endl;
             //vector <float>(axis1).swap(axis1);  
             //vector <float>(axis2).swap(axis2);
          }
       }
       REMC_flag=1;
   }
    //cout<<"REMC_flag is:"<<REMC_flag<<endl; 
    if(initialnum > molvec.size()) 
        initialnum=molvec.size();
    //cout<<"final checking initialnum:"<<initialnum<<endl;
    vector<int> molindex, molindex1;
    vector<DOCKMol> tempmachmols;
    for (int i = 0;i<molvec.size();i++)
    {

      tempmachmols.push_back(molvec[i]);
      //cout<<"checking 9999"<<molvec[i].x[0]<<"\t"<<tempmachmols[tempmachmols.size()-1].x[0]<<endl;
      molindex.push_back(i);
      molindex1.push_back(i);
    }
    sortmol(tempmachmols,molindex);
  
    for(int i=0;i<initialnum;i++)
      cout<<i<<tempmachmols[i].x[0]<<"   " << tempmachmols[i].current_score<<endl;


    vector<float> centers1(3, 0);
    for (int i = 0; i < tempmachmols[molindex[0]].num_atoms; i++)
    {
        centers1[0] += tempmachmols[molindex[0]].x[i];
        centers1[1] += tempmachmols[molindex[0]].y[i];
        centers1[2] += tempmachmols[molindex[0]].z[i];
    }
    centers1[0] = centers1[0] / tempmachmols[molindex[0]].num_atoms;
    centers1[1] = centers1[1] / tempmachmols[molindex[0]].num_atoms;
    centers1[2] = centers1[2] / tempmachmols[molindex[0]].num_atoms;

    for (int i = 0; i < tempmachmols[molindex[0]].num_atoms; i++)
    {
        for (int j = i; j < tempmachmols[molindex[0]].num_atoms; j++)
        {
            car_x = (tempmachmols[molindex[0]].x[i] - tempmachmols[molindex[0]].x[j]) * (tempmachmols[molindex[0]].x[i] - tempmachmols[molindex[0]].x[j]);
            car_y = (tempmachmols[molindex[0]].y[i] - tempmachmols[molindex[0]].y[j]) * (tempmachmols[molindex[0]].y[i] - tempmachmols[molindex[0]].y[j]);
            car_z = (tempmachmols[molindex[0]].z[i] - tempmachmols[molindex[0]].z[j]) * (tempmachmols[molindex[0]].z[i] - tempmachmols[molindex[0]].z[j]);
            if (sqrt(car_x + car_x + car_x) > max)
                max = sqrt(car_x + car_x + car_x);
        }
    }
    cout << "max:  " << max << endl;
    theDensityMap.cart(centers1, fartbinding, max);
    cout << "centers1:  " << centers1[0] << "  " << centers1[1] << "  " << centers1[2] << endl;


    ObjexxFCL::FArray3D< float > rhoC, rhoMask, rhoOmask, rhoO2mask;
    ObjexxFCL::FArray3D< std::complex<float> > FrhoC, FrhoMask, FrhoCmask, FrhoOmask, FrhoO2mask;
    ObjexxFCL::FArray3D< std::complex<float> > FrhoO, FrhoO2;

    // read Density map

    int Denx = theDensityMap.density.u1();
    int Deny = theDensityMap.density.u2();
    int Denz = theDensityMap.density.u3();
    int Denxyz = Denx * Deny * Denz;
    vector<int> grid(3, 0);
    grid[0] = theDensityMap.grid[0];
    grid[1] = theDensityMap.grid[1];
    grid[2] = theDensityMap.grid[2];
    vector<float> origin(3, 0.0);
    origin[0] = theDensityMap.origin[0];
    origin[1] = theDensityMap.origin[1];
    origin[2] = theDensityMap.origin[2];
    float padding = theDensityMap.ATOM_MASK_PADDING;
    float ca_m = theDensityMap.CA_MASK;
    float atom_m = theDensityMap.ATOM_MASK;
    //  theDensityMap.calcRhoC(posey,4.0,rhoC,rhoMask);
    //    theDensityMap.calcRhoCx(posex,4.0,rhoC,rhoMask);
    //  theDensityMap.calcRhoCy(pose,4.0,rhoC,rhoMask);
    //    float CC;
    //  CC=theDensityMap.getRSCC(rhoC,rhoMask);
     //   CC=theDensityMap.getRSCCX(rhoC,rhoMask);
     //   cout<<"CC: "<<CC<<endl; 
    cout << "origin: " << origin[0] << " " << origin[1] << " " << origin[2] << endl;

    DOCKMol fin_pose;
    float KT = 0.001;


    vector<float> axyz;
    axyz = vector<float>(3, 0.0);
    float ang;
    float angle_rotate;
    float old_dE1 = 0.0;
    float new_dE1 = 0.0;
    float old_CC = 0.0;
    float new_CC = 0.0;
    float old_dE2 = 0.0;
    float new_dE2 = 0.0;
    float new_dE = 0.0;
    float old_dE = 0.0;
    float dE = 0.0;
    float old_vwd = 0.0;
    float new_vwd = 0.0;
    float old_clash = 0.0;
    float new_clash = 0.0;

    vector<float> coord_change(3, 0.0);
    DOCKMol best_model;
    float best_E = 10.0;

    //vector<vector<float>> atm_idx(pnum, vector<float>(3, 0.0));
    DOCKMol pointsBx = mol;

    float effReso = std::max(2.4 + 0.8 * resolution, double(resolution));
    float k = (M_PI / effReso) * (M_PI / effReso);
    float a1 = 33.0;  // treat everything as ALA
    float C = a1 * pow(k / 3.1415926, 1.5);
    int tmp_i = 0;
    ObjexxFCL::FArray3D< float > rhoC0;
    //    ObjexxFCL::FArray3D< float > rhoC01;
    ObjexxFCL::FArray3D< float > inv_rho_mask0;
    //    ObjexxFCL::FArray3D< float > inv_rho_mask1;
    rhoC0.dimension(Denx, Deny, Denz, 0.0);
    //    rhoC01.dimension(Denx , Deny , Denz);
    inv_rho_mask0.dimension(Denx, Deny, Denz, 1.0);
    //    inv_rho_mask1.dimension(Denx , Deny , Denz);
    for (int t = 0; t < Denxyz; ++t) {
        rhoC0[t] = 0.0;
        inv_rho_mask0[t] = 1.0;
    }
    vector<float> del_ijx(3, 0.0);
    vector<float> atm_jx(3, 0.0);
    string elt_i;
    float bond0max = 0.0, bond1max = 0.0, bond2max = 0.0;
    float bond0min = 1000.0, bond1min = 1000.0, bond2min = 1000.0;
    int pp_r = pointsBx.num_atoms;
    for (int i = 0; i < pp_r; i++)
    {
        vector<float> cartX1(3, 0.0);
        vector<float> fracX1;
        elt_i = pointsBx.atom_names[i];
        elt_i = elt_i[0];
        OneGaussianScattering sig_j = get_A(elt_i);
        k = sig_j.k(theDensityMap.effectiveB);
        C = sig_j.C(k);
        if (C < 1e-6) continue;

        cartX1[0] = pointsBx.x[i];
        cartX1[1] = pointsBx.y[i];
        cartX1[2] = pointsBx.z[i];
        MatrixTimesTransVector(theDensityMap.c2f, cartX1, fracX1);
        // the location of atom in grid ?
        vector<float> atm_idxt(3, 0.0);
        atm_idxt[0] = (double(fracX1[0] * grid[0] - origin[0] + 1));
        atm_idxt[1] = (double(fracX1[1] * grid[1] - origin[1] + 1));
        atm_idxt[2] = (double(fracX1[2] * grid[2] - origin[2] + 1));
        for (int z = 1; z <= Denz; z++)
        {
            if (z > theDensityMap.z_min && z < theDensityMap.z_max)
            {
                atm_jx[2] = z;
                del_ijx[2] = (atm_idxt[2] - atm_jx[2]) / grid[2];
                del_ijx[0] = del_ijx[1] = 0.0;
                vector<float> frac_tmpz;
                MatrixTimesTransVector(theDensityMap.f2c, del_ijx, frac_tmpz);
                if (square_len(frac_tmpz) > (padding + ca_m) * (padding + ca_m)) continue;
                for (int y = 1; y <= Deny; y++)
                {
                    if (y > theDensityMap.y_min && y < theDensityMap.y_max)
                    {
                        atm_jx[1] = y;
                        del_ijx[1] = (atm_idxt[1] - atm_jx[1]) / grid[1];
                        // wrap-around??                      
                        del_ijx[0] = 0.0;
                        vector<float> frac_tmpy;
                        MatrixTimesTransVector(theDensityMap.f2c, del_ijx, frac_tmpy);
                        if (square_len(frac_tmpy) > (padding + ca_m) * (padding + ca_m)) continue;
                        for (int x = 1; x <= Denx; x++)
                        {
                            if (x > theDensityMap.x_min && x < theDensityMap.x_max)
                            {
                                atm_jx[0] = x;
                                del_ijx[0] = (atm_idxt[0] - atm_jx[0]) / grid[0];
                                // wrap-around??                        
                                vector<float> cart_del_ij2;
                                MatrixTimesTransVector(theDensityMap.f2c, del_ijx, cart_del_ij2);
                                float d2 = square_len(cart_del_ij2);
                                if (d2 > (padding + ca_m) * (padding + ca_m)) continue;

                                float atm = C * exp(-k * d2);
                                float sigmoid_msk = exp(d2 - (atom_m) * (atom_m));

                                float inv_msk = 1 / (1 + sigmoid_msk);
                                rhoC0(x, y, z) += atm;
                                inv_rho_mask0(x, y, z) *= (1 - inv_msk);
                            }
                            else continue;
                        }
                    }
                    else continue;
                }
            }
            else continue;
        }
        tmp_i = tmp_i + 1;
    }
    vector<float> coor_pdb(3, 0.0);
    float sumC_i2 = 0.0, sumO_i2 = 0.0, sumCO_i2 = 0.0, vol_i2 = 0.0, CC_i2 = 0.0;
    float sumO2_i2 = 0.0, sumC2_i2 = 0.0, varC_i2 = 0.0, varO_i2 = 0.0;
    float clc_x2 = 0.0, obs_x2 = 0.0, eps_x2 = 0.0;
    for (int x = 1; x <= Denx; ++x)
    {
        if (x > theDensityMap.x_min && x < theDensityMap.x_max)
        {
            for (int y = 1; y <= Deny; ++y)
            {
                if (y > theDensityMap.y_min && y < theDensityMap.y_max)
                {
                    for (int z = 1; z <= Denz; ++z)
                    {
                        if (z > theDensityMap.z_min && z < theDensityMap.z_max)
                        {
                            clc_x2 = rhoC0(x, y, z);
                            //                clc_x2 = tmp_rhc;
                            obs_x2 = theDensityMap.density(x, y, z);
                            //                obs_x2 = tmp_den;
                            eps_x2 = 1.0 - inv_rho_mask0(x, y, z); //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
                            sumCO_i2 += eps_x2 * clc_x2 * obs_x2;
                            sumO_i2 += eps_x2 * obs_x2;
                            sumO2_i2 += eps_x2 * obs_x2 * obs_x2;
                            sumC_i2 += eps_x2 * clc_x2;
                            sumC2_i2 += eps_x2 * clc_x2 * clc_x2;
                            vol_i2 += eps_x2;
                        }
                        else continue;
                    }
                }
                else continue;
            }
        }
        else continue;
    }
    varC_i2 = (sumC2_i2 - sumC_i2 * sumC_i2 / vol_i2);
    varO_i2 = (sumO2_i2 - sumO_i2 * sumO_i2 / vol_i2);
    if (varC_i2 == 0 || varO_i2 == 0 || vol_i2 == 0) {
        CC_i2 = 0;
    }
    else {
        CC_i2 = (sumCO_i2 - sumC_i2 * sumO_i2 / vol_i2) / sqrt(varC_i2 * varO_i2);
    }
    best_E = 1.0 - CC_i2;
    cout << "best_E: " << best_E << endl;
    best_model = pointsBx;
    coor_pdb = vector<float>(3, 0.0);
    for (int j = 0; j < pointsBx.num_atoms; j++)
    {
        coor_pdb[0] = coor_pdb[0] + pointsBx.x[j];
        coor_pdb[1] = coor_pdb[1] + pointsBx.y[j];
        coor_pdb[2] = coor_pdb[2] + pointsBx.z[j];
    }
    coor_pdb[0] = coor_pdb[0] / (pointsBx.num_atoms);
    coor_pdb[1] = coor_pdb[1] / (pointsBx.num_atoms);
    coor_pdb[2] = coor_pdb[2] / (pointsBx.num_atoms);
    cout << "coor_pdb: " << coor_pdb[0] << " " << coor_pdb[1] << " " << coor_pdb[2] << endl;

    vector<float> frax_pdb;
    MatrixTimesTransVector(theDensityMap.c2f, coor_pdb, frax_pdb);
    vector<float> frac_atm(3, 0.0);
    frac_atm[0] = (double(frax_pdb[0] * grid[0] - origin[0] + 1));
    frac_atm[1] = (double(frax_pdb[1] * grid[1] - origin[1] + 1));
    frac_atm[2] = (double(frax_pdb[2] * grid[2] - origin[2] + 1));

    vector<float> del_tran(3, 0.0);
    del_tran[2] = (frac_atm[2] - float(grid[2]) / 2.0) / float(grid[2]);
    del_tran[1] = (frac_atm[1] - float(grid[1]) / 2.0) / float(grid[1]);
    del_tran[0] = (frac_atm[0] - float(grid[0]) / 2.0) / float(grid[0]);
    vector<float> tran0(3, 0.0);
    tran0[0] = coor_pdb[0] - centers1[0];
    tran0[1] = coor_pdb[1] - centers1[1];
    tran0[2] = coor_pdb[2] - centers1[2];
    //MatrixTimesTransVector(theDensityMap.f2c, del_tran, tran0);

    cout << "tran0: " << tran0[0] << " " << tran0[1] << " " << tran0[2] << endl;
    DOCKMol pointsBxj;
    pointsBxj = pointsBx;
    int p_atm_size = pointsBx.num_atoms;
    // -tran0 平移量把PDB结构中心变成密度图中心
    for (int j = 0; j < p_atm_size; j++) // not last N
    {
        pointsBxj.x[j] = -tran0[0] + pointsBxj.x[j];
        pointsBxj.y[j] = -tran0[1] + pointsBxj.y[j];
        pointsBxj.z[j] = -tran0[2] + pointsBxj.z[j];
    }

    rhoC0.dimension(Denx, Deny, Denz, 0.0);
    inv_rho_mask0.dimension(Denx, Deny, Denz, 1.0);

    for (int t = 0; t < Denxyz; ++t) {
        rhoC0[t] = 0.0;
        inv_rho_mask0[t] = 1.0;
    }
    del_ijx = vector<float>(3, 0.0);
    atm_jx = vector<float>(3, 0.0);

    tmp_i = 0;
    bond0max = -1000.0; bond1max = -1000.0; bond2max = -1000.0;
    bond0min = 1000.0; bond1min = 1000.0; bond2min = 1000.0;
    //    cout<<"FF"<<endl;      
    for (int i = 0; i < p_atm_size; i++)
    {
        vector<float> cartX1(3, 0.0);
        vector<float> fracX1;
        elt_i = pointsBxj.atom_names[i];
        elt_i = elt_i[0];
        OneGaussianScattering sig_j = get_A(elt_i);
        k = sig_j.k(theDensityMap.effectiveB);
        C = sig_j.C(k);
        if (C < 1e-6) continue;

        cartX1[0] = pointsBxj.x[i];
        cartX1[1] = pointsBxj.y[i];
        cartX1[2] = pointsBxj.z[i];
        MatrixTimesTransVector(theDensityMap.c2f, cartX1, fracX1);

        // the location of atom in grid ?
        vector<float> atm_idxt(3, 0.0);
        atm_idxt[0] = (double(fracX1[0] * grid[0] - origin[0] + 1));
        atm_idxt[1] = (double(fracX1[1] * grid[1] - origin[1] + 1));
        atm_idxt[2] = (double(fracX1[2] * grid[2] - origin[2] + 1));

        for (int z = 1; z <= Denz; z++)
        {
            //if (z > theDensityMap.z_min && z < theDensityMap.z_max)
            {
                atm_jx[2] = z;
                del_ijx[2] = (atm_idxt[2] - atm_jx[2]) / grid[2];
                del_ijx[0] = del_ijx[1] = 0.0;
                vector<float> frac_tmpz;
                MatrixTimesTransVector(theDensityMap.f2c, del_ijx, frac_tmpz);
                if (square_len(frac_tmpz) > (padding + ca_m) * (padding + ca_m)) continue;
                for (int y = 1; y <= Deny; y++)
                {
                    //if (y > theDensityMap.y_min && y < theDensityMap.y_max)
                    {
                        atm_jx[1] = y;
                        del_ijx[1] = (atm_idxt[1] - atm_jx[1]) / grid[1];
                        // wrap-around??                    
                        del_ijx[0] = 0.0;
                        vector<float> frac_tmpy;
                        MatrixTimesTransVector(theDensityMap.f2c, del_ijx, frac_tmpy);
                        if (square_len(frac_tmpy) > (padding + ca_m) * (padding + ca_m)) continue;
                        for (int x = 1; x <= Denx; x++)
                        {
                            //if (x > theDensityMap.x_min && x < theDensityMap.x_max)
                            {
                                atm_jx[0] = x;
                                del_ijx[0] = (atm_idxt[0] - atm_jx[0]) / grid[0];
                                // wrap-around??                          
                                vector<float> cart_del_ij2;
                                MatrixTimesTransVector(theDensityMap.f2c, del_ijx, cart_del_ij2);
                                float d2 = square_len(cart_del_ij2);
                                if (d2 > (padding + ca_m) * (padding + ca_m)) continue;

                                float atm = C * exp(-k * d2);
                                float sigmoid_msk = exp(d2 - (atom_m) * (atom_m));

                                float inv_msk = 1 / (1 + sigmoid_msk);
                                rhoC0(x, y, z) += atm;
                                inv_rho_mask0(x, y, z) *= (1 - inv_msk);
                            }
                            //else continue;
                        }
                    }
                    //else continue;
                }
            }
            //else continue;
        }
        tmp_i = tmp_i + 1;
    }
    sumC_i2 = 0.0, sumO_i2 = 0.0, sumCO_i2 = 0.0, vol_i2 = 0.0, CC_i2 = 0.0;
    sumO2_i2 = 0.0, sumC2_i2 = 0.0, varC_i2 = 0.0, varO_i2 = 0.0;
    clc_x2 = 0.0, obs_x2 = 0.0, eps_x2 = 0.0;
    for (int x = 1; x <= Denx; x++)
    {
        //if (x > theDensityMap.x_min && x < theDensityMap.x_max)
        {
            for (int y = 1; y <= Deny; y++)
            {
                //if (y > theDensityMap.y_min && y < theDensityMap.y_max)
                {
                    for (int z = 1; z <= Denz; z++)
                    {
                        //if (z > theDensityMap.z_min && z < theDensityMap.z_max)
                        {
                            clc_x2 = rhoC0(x, y, z);
                            //                clc_x2 = tmp_rhc;
                            obs_x2 = theDensityMap.density(x, y, z);
                            //                obs_x2 = tmp_den;
                            eps_x2 = 1.0 - inv_rho_mask0(x, y, z); //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
                            sumCO_i2 += eps_x2 * clc_x2 * obs_x2;
                            sumO_i2 += eps_x2 * obs_x2;
                            sumO2_i2 += eps_x2 * obs_x2 * obs_x2;
                            sumC_i2 += eps_x2 * clc_x2;
                            sumC2_i2 += eps_x2 * clc_x2 * clc_x2;
                            vol_i2 += eps_x2;
                        }
                        //else continue;
                    }
                }
                //else continue;
            }
        }
        //else continue;
    }
    varC_i2 = (sumC2_i2 - sumC_i2 * sumC_i2 / vol_i2);
    varO_i2 = (sumO2_i2 - sumO_i2 * sumO_i2 / vol_i2);
    if (varC_i2 == 0 || varO_i2 == 0 || vol_i2 == 0) {
        CC_i2 = 0;
    }
    else {
        CC_i2 = (sumCO_i2 - sumC_i2 * sumO_i2 / vol_i2) / sqrt(varC_i2 * varO_i2);
    }
    float E_2 = 1.0 - CC_i2;
    cout << "E_2: " << E_2 << endl;
    if (E_2 <= best_E)
    {
        pointsBx = pointsBxj;
    }
    cout << "CCCCC" << endl;

    string molfile = "test.mol2";
    //ofstream clusterpdb(pdbfile.c_str());
    ofstream clustermol(molfile.c_str());
    string str = string(argv[1]);
    clustermol << "########## Name:" << endl;
    Write_Mol2(pointsBx, clustermol, str);

    ElectronDensity theDensityMapx;
    theDensityMapx.preso = resolution;
    theDensityMapx.pATOM_MASK = 7.0;
    theDensityMapx.pCA_MASK = 7.0;
    theDensityMapx.pforce_apix_on_map_load_ = 0.0;
    theDensityMapx.pWINDOW_ = 1;
    theDensityMapx.pscore_window_context_ = false;
    theDensityMapx.premap_symm_ = false;
    theDensityMapx.pde_edensity = false;
    theDensityMapx.pnkbins_ = 0;
    float mapspx = 0.0;
    theDensityMapx.readMRCandResize(inputMRC, resolution, mapspx);
    int pDenx = theDensityMapx.density.u1();
    int pDeny = theDensityMapx.density.u2();
    int pDenz = theDensityMapx.density.u3();
    int pDenxyz = pDenx * pDeny * pDenz;
    vector<int> pgrid(3, 0);
    pgrid[0] = theDensityMapx.grid[0];
    pgrid[1] = theDensityMapx.grid[1];
    pgrid[2] = theDensityMapx.grid[2];
    vector<float> porigin(3, 0.0);
    porigin[0] = theDensityMapx.origin[0];
    porigin[1] = theDensityMapx.origin[1];
    porigin[2] = theDensityMapx.origin[2];
    float ppadding = theDensityMapx.ATOM_MASK_PADDING;
    float pca_m = theDensityMapx.CA_MASK;
    float patom_m = theDensityMapx.ATOM_MASK;

    ObjexxFCL::FArray3D< float > rhoC0x;
    ObjexxFCL::FArray3D< float > inv_rho_mask0x;
    rhoC0x.dimension(theDensityMapx.density.u1(), theDensityMapx.density.u2(), theDensityMapx.density.u3());
    inv_rho_mask0x.dimension(theDensityMapx.density.u1(), theDensityMapx.density.u2(), theDensityMapx.density.u3());

    cout << "density,x y z: " << pDenx << " " << pDeny << " " << pDenz << endl;
    float new_CC1 = 0.0;
    float old_CC1 = 0.0;
    float new_CC2 = 0.0;

    CC_i2 = theDensityMap.clCC(pointsBx);
    cout << "CC_i2: " << CC_i2 << endl;
    new_CC2 = 1.0 - CC_i2;
    cout << "new_CC2: " << new_CC2 << endl;
    best_E = new_CC2;
    //vector<poseCoord>().swap(best_model);
    best_model = pointsBx;
    cout << "rotate 0 and 180" << endl;
    // rotate 60 and 180
//    vector<poseCoord>().swap(pointsBx);
//    pointsBx = best_model;
    int tmp3 = 0;
    float init_E0 = 0.0;
    float init_E1 = 0.0;
    //vector<poseCoord> init_mat0;
    //vector<poseCoord> init_mat1;
    //int seq_num = pnum;
    float out_sup_all_CC = 2.0;

    float supp_KT0s = 1.5;//;0.0001;
    float supp_KT0e = 0.001;//;0.00001;
    int supp_num1 = 20;
    int supp_num2 = 12;
    int supp_num3 = 500;
    int supp_num4 = supp_num1;
    vector<DOCKMol> decstr_tmp;
    vector<vector<float>> E_Tx(supp_num4);

    vector<float> Eng_T(supp_num1, 0.0);
    int nnx = 0;
    int nny = 0;

    vector<float> best_Ex(supp_num1, 100000.0);
    vector<float> best_Ex_4(4, 100000.0);
    vector<DOCKMol> best_modelx(supp_num1);
    vector<DOCKMol> best_modelx_4(4);
    DOCKMol best_model_u;
    float best_E_u = 1000;
    vector<float> supp_REMC(supp_num1, 0.0);
    for (int i = 0; i < supp_num1; i++)
    {
        float supp_KT0x = pow(float(supp_KT0e / supp_KT0s), float(float(i) / float(supp_num1)));
        supp_REMC[i] = supp_KT0s * float(supp_KT0x);
    }
    vector<DOCKMol > decstrp;
    //    float ang_cus[] = {0.0,0.0,0.0,60.0,60.0,60.0,100.0,100.0,100.0,150.0,150.0,150.0,200.0,200.0,200.0,270.0,270.0,270.0,330.0,330.0,330.0};
    for (int kkk = 0; kkk < supp_num2; kkk++) // 2
    {
        vector<DOCKMol > decstrz;
        vector<DOCKMol >().swap(decstrz);
        vector<float> E_REMC(supp_num1, 0.0);
        for (int jjj = 0; jjj < supp_num1; jjj++)
        {
            DOCKMol point_mat;
            if (kkk > 0)
            {
                point_mat = decstrp[jjj];
                coor_pdb = vector<float>(3, 0.0);
                for (int j = 0; j < point_mat.num_atoms; j++)
                {
                    coor_pdb[0] = coor_pdb[0] + point_mat.x[j];
                    coor_pdb[1] = coor_pdb[1] + point_mat.y[j];
                    coor_pdb[2] = coor_pdb[2] + point_mat.z[j];
                }
                coor_pdb[0] = coor_pdb[0] / (point_mat.num_atoms);
                coor_pdb[1] = coor_pdb[1] / (point_mat.num_atoms);
                coor_pdb[2] = coor_pdb[2] / (point_mat.num_atoms);
            }
            else
            {
                coor_pdb = vector<float>(3, 0.0);
                point_mat = pointsBx;
                for (int j = 0; j < point_mat.num_atoms; j++)
                {
                    coor_pdb[0] = coor_pdb[0] + point_mat.x[j];
                    coor_pdb[1] = coor_pdb[1] + point_mat.y[j];
                    coor_pdb[2] = coor_pdb[2] + point_mat.z[j];
                }
                coor_pdb[0] = coor_pdb[0] / (point_mat.num_atoms);
                coor_pdb[1] = coor_pdb[1] / (point_mat.num_atoms);
                coor_pdb[2] = coor_pdb[2] / (point_mat.num_atoms);
            }
            //            float zang=80.0;
            float zang = float(randIntCustom(30, 180));
            //            float zang=ang_cus[jjj];
            float zang_rd = zang * M_PI / 180.0;

            if (jjj % 3 == 0 && kkk == 0) // rotate around Z axis
            {
                float asin_theta = 2.0 * randf0and1() - 1.0;
                float acos_theta = sqrt(1.0 - asin_theta * asin_theta);
                float apha = 2.0 * PI * randf0and1();
                float awx = 0.0;
                float awy = 0.0;
                float awz = 1.0;
                // Translation Vector
                float t0 = 0.0;
                float angle_rotategg = zang_rd; // rotate angle  
                float t1 = (randf0and1() * 2.0 - 1.0) * t0 + 0.0;
                float t2 = (randf0and1() * 2.0 - 1.0) * t0 + 0.0;
                float t3 = (randf0and1() * 2.0 - 1.0) * t0 + 0.0;
                //    cout<<"ang,t0: "<<angle_rotategg<<" ("<<t1<<","<<t2<<","<<t3<<")"<<endl;          

                float asin = sin(angle_rotategg);
                float acos = cos(angle_rotategg);
                float u[3][3];
                u[0][0] = acos + awx * awx * (1.0 - acos);
                u[0][1] = awx * awy * (1.0 - acos) - awz * asin;
                u[0][2] = awx * awz * (1.0 - acos) + awy * asin;
                u[1][0] = awx * awy * (1.0 - acos) + awz * asin;
                u[1][1] = acos + awy * awy * (1.0 - acos);
                u[1][2] = awy * awz * (1.0 - acos) - awx * asin;
                u[2][0] = awx * awz * (1.0 - acos) - awy * asin;
                u[2][1] = awy * awz * (1.0 - acos) + awx * asin;
                u[2][2] = acos + awz * awz * (1.0 - acos);
                // Rotation points
                axyz[0] = coor_pdb[0];
                axyz[1] = coor_pdb[1];
                axyz[2] = coor_pdb[2];

                DOCKMol fin_matx;
                fin_matx = point_mat;
                tmp3 = 0;
                int pp_q = point_mat.num_atoms;
                for (int j = 0; j < pp_q; j++) // not last N
                {
                    point_mat.x[j] = t1 + axyz[0] + (fin_matx.x[j] - axyz[0]) * u[0][0] + (fin_matx.y[j] - axyz[1]) * u[0][1] + (fin_matx.z[j] - axyz[2]) * u[0][2];
                    point_mat.y[j] = t2 + axyz[1] + (fin_matx.x[j] - axyz[0]) * u[1][0] + (fin_matx.y[j] - axyz[1]) * u[1][1] + (fin_matx.z[j] - axyz[2]) * u[1][2];
                    point_mat.z[j] = t3 + axyz[2] + (fin_matx.x[j] - axyz[0]) * u[2][0] + (fin_matx.y[j] - axyz[1]) * u[2][1] + (fin_matx.z[j] - axyz[2]) * u[2][2];
                    tmp3 = tmp3 + 1;
                }
                //vector<poseCoord>().swap(fin_matx);
            }
            if (jjj % 3 == 1 && kkk == 0)  //rotate around X axis
            {
                float asin_theta = 2.0 * randf0and1() - 1.0;
                float acos_theta = sqrt(1.0 - asin_theta * asin_theta);
                float apha = 2.0 * PI * randf0and1();
                float awx = 1.0;
                float awy = 0.0;
                float awz = 0.0;
                // Translation Vector
                float t0 = 0.0;
                float angle_rotategg = zang_rd; // rotate angle  
                float t1 = (randf0and1() * 2.0 - 1.0) * t0 + 0.0;
                float t2 = (randf0and1() * 2.0 - 1.0) * t0 + 0.0;
                float t3 = (randf0and1() * 2.0 - 1.0) * t0 + 0.0;
                //    cout<<"ang,t0: "<<angle_rotategg<<" ("<<t1<<","<<t2<<","<<t3<<")"<<endl;          

                float asin = sin(angle_rotategg);
                float acos = cos(angle_rotategg);
                float u[3][3];
                u[0][0] = acos + awx * awx * (1.0 - acos);
                u[0][1] = awx * awy * (1.0 - acos) - awz * asin;
                u[0][2] = awx * awz * (1.0 - acos) + awy * asin;
                u[1][0] = awx * awy * (1.0 - acos) + awz * asin;
                u[1][1] = acos + awy * awy * (1.0 - acos);
                u[1][2] = awy * awz * (1.0 - acos) - awx * asin;
                u[2][0] = awx * awz * (1.0 - acos) - awy * asin;
                u[2][1] = awy * awz * (1.0 - acos) + awx * asin;
                u[2][2] = acos + awz * awz * (1.0 - acos);
                // Rotation points
                axyz[0] = coor_pdb[0];
                axyz[1] = coor_pdb[1];
                axyz[2] = coor_pdb[2];

                DOCKMol fin_matx;
                fin_matx = point_mat;
                tmp3 = 0;
                int pp_q = point_mat.num_atoms;
                for (int j = 0; j < pp_q; j++) // not last N
                {
                    point_mat.x[j] = t1 + axyz[0] + (fin_matx.x[j] - axyz[0]) * u[0][0] + (fin_matx.y[j] - axyz[1]) * u[0][1] + (fin_matx.z[j] - axyz[2]) * u[0][2];
                    point_mat.y[j] = t2 + axyz[1] + (fin_matx.x[j] - axyz[0]) * u[1][0] + (fin_matx.y[j] - axyz[1]) * u[1][1] + (fin_matx.z[j] - axyz[2]) * u[1][2];
                    point_mat.z[j] = t3 + axyz[2] + (fin_matx.x[j] - axyz[0]) * u[2][0] + (fin_matx.y[j] - axyz[1]) * u[2][1] + (fin_matx.z[j] - axyz[2]) * u[2][2];
                    tmp3 = tmp3 + 1;
                }
                //vector<poseCoord>().swap(fin_matx);
            }
            if (jjj % 3 == 2 && kkk == 0)  //rotate around Y axis
            {
                float asin_theta = 2.0 * randf0and1() - 1.0;
                float acos_theta = sqrt(1.0 - asin_theta * asin_theta);
                float apha = 2.0 * PI * randf0and1();
                float awx = 0.0;
                float awy = 1.0;
                float awz = 0.0;
                // Translation Vector
                float t0 = 0.0;
                float angle_rotategg = zang_rd; // rotate angle  
                float t1 = (randf0and1() * 2.0 - 1.0) * t0 + 0.0;
                float t2 = (randf0and1() * 2.0 - 1.0) * t0 + 0.0;
                float t3 = (randf0and1() * 2.0 - 1.0) * t0 + 0.0;
                //    cout<<"ang,t0: "<<angle_rotategg<<" ("<<t1<<","<<t2<<","<<t3<<")"<<endl;          

                float asin = sin(angle_rotategg);
                float acos = cos(angle_rotategg);
                float u[3][3];
                u[0][0] = acos + awx * awx * (1.0 - acos);
                u[0][1] = awx * awy * (1.0 - acos) - awz * asin;
                u[0][2] = awx * awz * (1.0 - acos) + awy * asin;
                u[1][0] = awx * awy * (1.0 - acos) + awz * asin;
                u[1][1] = acos + awy * awy * (1.0 - acos);
                u[1][2] = awy * awz * (1.0 - acos) - awx * asin;
                u[2][0] = awx * awz * (1.0 - acos) - awy * asin;
                u[2][1] = awy * awz * (1.0 - acos) + awx * asin;
                u[2][2] = acos + awz * awz * (1.0 - acos);
                // Rotation points
                axyz[0] = coor_pdb[0];
                axyz[1] = coor_pdb[1];
                axyz[2] = coor_pdb[2];

                DOCKMol fin_matx;
                fin_matx = point_mat;
                tmp3 = 0;
                int pp_q = point_mat.num_atoms;
                for (int j = 0; j < pp_q; j++) // not last N
                {
                    point_mat.x[j] = t1 + axyz[0] + (fin_matx.x[j] - axyz[0]) * u[0][0] + (fin_matx.y[j] - axyz[1]) * u[0][1] + (fin_matx.z[j] - axyz[2]) * u[0][2];
                    point_mat.y[j] = t2 + axyz[1] + (fin_matx.x[j] - axyz[0]) * u[1][0] + (fin_matx.y[j] - axyz[1]) * u[1][1] + (fin_matx.z[j] - axyz[2]) * u[1][2];
                    point_mat.z[j] = t3 + axyz[2] + (fin_matx.x[j] - axyz[0]) * u[2][0] + (fin_matx.y[j] - axyz[1]) * u[2][1] + (fin_matx.z[j] - axyz[2]) * u[2][2];
                    tmp3 = tmp3 + 1;
                }
                //vector<poseCoord>().swap(fin_matx);
            }
            //        best_E=10.0;
            coor_pdb = vector<float>(3, 0.0);
            for (int j = 0; j < point_mat.num_atoms; j++)
            {
                coor_pdb[0] = coor_pdb[0] + point_mat.x[j];
                coor_pdb[1] = coor_pdb[1] + point_mat.y[j];
                coor_pdb[2] = coor_pdb[2] + point_mat.z[j];
            }
            coor_pdb[0] = coor_pdb[0] / (point_mat.num_atoms);
            coor_pdb[1] = coor_pdb[1] / (point_mat.num_atoms);
            coor_pdb[2] = coor_pdb[2] / (point_mat.num_atoms);

            if (kkk == 0)
            {
                best_modelx[jjj] = point_mat;
            }

            float E_500 = 0.0;
            int E_500_int = 0;
            KT = supp_REMC[jjj];

            CC_i2 = theDensityMap.clCC(point_mat);

            new_CC2 = 1.0 - CC_i2;
            old_dE = new_CC2;
            //    cout<<"new_CC2: "<<new_CC2<<endl;
            for (int iii = 0; iii < supp_num3; iii++)
            {
                float asin_theta = 2.0 * randf0and1() - 1.0;
                float acos_theta = sqrt(1.0 - asin_theta * asin_theta);
                float apha = 2.0 * PI * randf0and1();
                float awx = acos_theta * cos(apha);
                float awy = acos_theta * sin(apha);
                float awz = asin_theta;
                // Translation Vector
                float t0 = 1.0;
                //    float t0=1.0;

                //    cout<<"t0: "<<t1<<" "<<t2<<" "<<t3<<endl;

                    // Rotation matrix
                float anggg = 30.0;
                float angle_rotategg = (2.0 * randf0and1() - 1.0) * anggg; // rotate angle

                float t1 = (randf0and1() * 2.0 - 1.0) * t0 + 0.0;
                float t2 = (randf0and1() * 2.0 - 1.0) * t0 + 0.0;
                float t3 = (randf0and1() * 2.0 - 1.0) * t0 + 0.0;
                //    cout<<"ang,t0: "<<angle_rotategg<<" ("<<t1<<","<<t2<<","<<t3<<")"<<endl;          

                float asin = sin(angle_rotategg);
                float acos = cos(angle_rotategg);
                float u[3][3];
                u[0][0] = acos + awx * awx * (1.0 - acos);
                u[0][1] = awx * awy * (1.0 - acos) - awz * asin;
                u[0][2] = awx * awz * (1.0 - acos) + awy * asin;
                u[1][0] = awx * awy * (1.0 - acos) + awz * asin;
                u[1][1] = acos + awy * awy * (1.0 - acos);
                u[1][2] = awy * awz * (1.0 - acos) - awx * asin;
                u[2][0] = awx * awz * (1.0 - acos) - awy * asin;
                u[2][1] = awy * awz * (1.0 - acos) + awx * asin;
                u[2][2] = acos + awz * awz * (1.0 - acos);

                axyz[0] = coor_pdb[0];
                axyz[1] = coor_pdb[1];
                axyz[2] = coor_pdb[2];

                //    beg_p=rand_point; 
            //    cout<<"XXX2"<<endl;  

                int tmp3 = 0;
                DOCKMol fin_mat;
                //    fin_mat = tmp_mat;
                //    if(beg_p==0)
                //    {
                tmp3 = 0;
                //vector<poseCoord>().swap(fin_mat);
                fin_mat = point_mat;
                DOCKMol fin_matx;
                fin_matx = fin_mat;
                tmp3 = 0;
                int pp_k = fin_mat.num_atoms;
                for (int j = 0; j < pp_k; j++) // not last N
                {
                    fin_mat.x[j] = t1 + axyz[0] + (fin_matx.x[j] - axyz[0]) * u[0][0] + (fin_matx.y[j] - axyz[1]) * u[0][1] + (fin_matx.z[j] - axyz[2]) * u[0][2];
                    fin_mat.y[j] = t2 + axyz[1] + (fin_matx.x[j] - axyz[0]) * u[1][0] + (fin_matx.y[j] - axyz[1]) * u[1][1] + (fin_matx.z[j] - axyz[2]) * u[1][2];
                    fin_mat.z[j] = t3 + axyz[2] + (fin_matx.x[j] - axyz[0]) * u[2][0] + (fin_matx.y[j] - axyz[1]) * u[2][1] + (fin_matx.z[j] - axyz[2]) * u[2][2];
                    tmp3 = tmp3 + 1;
                }
                //vector<poseCoord>().swap(fin_matx);

                CC_i2 = theDensityMap.clCC(fin_mat);

                new_dE = 1.0 - CC_i2;
                //        new_dE = new_CC1 ;                                   
            //            cout<<"old_CC1,new_CC1:RRRRRRRRRRRRR "<<old_CC1<<" "<<new_CC1<<endl;
                dE = new_dE - old_dE;
                if (new_dE < old_dE)
                {
                    //tmp3 = 0;
                    coor_pdb = vector<float>(3, 0.0);
                    point_mat = fin_mat;
                    for (int j = 0; j < point_mat.num_atoms; j++)
                    {
                        //point_mat[3 * j + 0] = fin_mat[3 * tmp3 + 0];
                        //point_mat[3 * j + 1] = fin_mat[3 * tmp3 + 1];
                        //point_mat[3 * j + 2] = fin_mat[3 * tmp3 + 2];

                        coor_pdb[0] = coor_pdb[0] + fin_mat.x[j];
                        coor_pdb[1] = coor_pdb[1] + fin_mat.y[j];
                        coor_pdb[2] = coor_pdb[2] + fin_mat.z[j];
                        ///tmp3 = tmp3 + 1;
                    }
                    coor_pdb[0] = coor_pdb[0] / (fin_mat.num_atoms);
                    coor_pdb[1] = coor_pdb[1] / (fin_mat.num_atoms);
                    coor_pdb[2] = coor_pdb[2] / (fin_mat.num_atoms);

                    if (new_dE < best_Ex[jjj])
                    {
                        best_Ex[jjj] = new_dE;
                        //vector<poseCoord>().swap(best_modelx[jjj]);
                        best_modelx[jjj] = fin_mat;

                    }
                    if (new_dE < best_E_u)
                    {
                        best_model_u = fin_mat;
                    }
                    if (new_dE < out_sup_all_CC)
                    {
                        out_sup_all_CC = new_dE;
                    }
                    //    cout<<"coor: "<<coor_pdb[0]<<" "<<coor_pdb[1]<<" "<<coor_pdb[2]<<endl;                  
                //        new_CC2 = new_CC1;
                //        cout<<"old_CC1,new_CC1: "<<old_dE<<" "<<new_dE<<endl;
                    old_dE = new_dE;
                    E_500 = new_dE;
                    E_500_int = 1;

                }
                else
                {
                    //    float tmpx=rand()/double(RAND_MAX);
                    float tmpx = randf0and1();
                    float mc_v = exp(-dE / (KT));
                    if (tmpx < mc_v) // CC must be >0
                    {
                        tmp3 = 0;
                        coor_pdb = vector<float>(3, 0.0);
                        point_mat = fin_mat;
                        for (int j = 0; j < fin_mat.num_atoms; j++)
                        {
                            coor_pdb[0] = coor_pdb[0] + fin_mat.x[j];
                            coor_pdb[1] = coor_pdb[1] + fin_mat.y[j];
                            coor_pdb[2] = coor_pdb[2] + fin_mat.z[j];
                            tmp3 = tmp3 + 1;
                        }
                        coor_pdb[0] = coor_pdb[0] / (fin_mat.num_atoms);
                        coor_pdb[1] = coor_pdb[1] / (fin_mat.num_atoms);
                        coor_pdb[2] = coor_pdb[2] / (fin_mat.num_atoms);
                        //    cout<<"coor: "<<coor_pdb[0]<<" "<<coor_pdb[1]<<" "<<coor_pdb[2]<<endl;
                        //    cout<<"old_CC1,new_CC1:XX "<<old_dE<<" "<<new_dE<<endl;                     
                        old_dE = new_dE;
                        E_500 = new_dE;
                        E_500_int = 1;
                        //    new_CC2 = new_CC1;

                    }
                }
                nnx = nnx + 1;
            }
            if (E_500_int == 1)
            {
                E_REMC[jjj] = E_500;
            }
            decstrz.push_back(point_mat);  // change            
        }
        vector<DOCKMol >().swap(decstrp);
        decstrp = vector<DOCKMol >(supp_num1);

        int js = kkk % 2;
        if (js == 0)
        {
            for (int i = 0; i < supp_num1 - 1; i = i + 2)
            {
                int j = i + 1;
                double CH_REMC = exp((1.0 / supp_REMC[i] - 1.0 / supp_REMC[j]) * (E_REMC[i] - E_REMC[j]));
                float Pglobal = min(1.0, CH_REMC);
                float tmpx = randf0and1();
                if (tmpx < Pglobal)
                {
                    decstrp[i] = decstrz[j];
                    decstrp[j] = decstrz[i];
                }
                else
                {
                    decstrp[i] = decstrz[i];
                    decstrp[j] = decstrz[j];
                }
            }
        }
        if (js == 1)
        {
            decstrp[0] = decstrz[0];
            for (int i = 1; i < supp_num1 - 1; i = i + 2)
            {
                int j = i + 1;
                double CH_REMC = exp((1.0 / supp_REMC[i] - 1.0 / supp_REMC[j]) * (E_REMC[i] - E_REMC[j]));
                float Pglobal = min(1.0, CH_REMC);
                float tmpx = randf0and1();
                if (tmpx < Pglobal)
                {
                    decstrp[i] = decstrz[j];
                    decstrp[j] = decstrz[i];
                }
                else
                {
                    decstrp[i] = decstrz[i];
                    decstrp[j] = decstrz[j];
                }
            }
            decstrp[supp_num1 - 1] = decstrz[supp_num1 - 1];
        }
    }

    pointsBx = best_model_u;
    string molfile1 = "test1.mol2";
    //ofstream clusterpdb(pdbfile.c_str());
    ofstream clustermol1(molfile1.c_str());
    string str1 = string(argv[1]);
    clustermol1 << "########## Name:" << endl;
    Write_Mol2(pointsBx, clustermol1, str1);

    for(int i=0;i<initialnum;i++)
    {
      matchvec.push_back(tempmachmols[molindex[i]]);
      DOCKMol moltemp;
      copy_molecule(moltemp,mol);
      matchvec.push_back(moltemp);//just stand one vector, not real conformation.    
    }

    cout<<"matchvec.size()"<<matchvec.size()<<endl;
    for(int i=0;i<matchvec.size();i++)
    {
      cout<<i<<" "<<matchvec[i].x[0]<<"   " << matchvec[i].current_score<<"  " << matchvec[i].current_EM_score<< endl;
    }
    /******************************************/
    //vector <DOCKMol> matchvec2;
    vector<float> axis1(3, 0);
    vector<float> axis2(3, 0);
    double a,b,c;
    for (int i=0;i<matchvec.size()-1;i=i+2)
    {   
       for(int m=0;m<mol.num_atoms;m++)
       {
         mol.x[m]=matchvec[i].x[m];
         mol.y[m]=matchvec[i].y[m];
         mol.z[m]=matchvec[i].z[m];

       }

       float ligcenterx,ligcentery,ligcenterz;
       ligcenterx=ligcentery=ligcenterz=0;
       for(int m=0;m<mol.num_atoms;m++)
       {
           ligcenterx +=mol.x[m];
           ligcentery +=mol.y[m];
           ligcenterz +=mol.z[m];
       }
       ligcenterx=ligcenterx/mol.num_atoms;
       ligcentery=ligcentery/mol.num_atoms;
       ligcenterz=ligcenterz/mol.num_atoms;
       //cout<<ligcenterx<<"\t"<<ligcentery<<"\t"<<ligcenterz<<endl;
       axis1[0]=ligcenterx; axis1[1]=ligcentery; axis1[2]=ligcenterz;
       a = rand0to1(gen);
       b = rand0to1(gen);
       //cout<<"a"<<a<<" "<<"b"<<b<<" "<<"c"<<c<<endl;
       double theta = 2 * 3.1415926 * a;
       double phi = acos(2*b-1.0);
       //cout<<"theta "<<theta<<" phi "<<phi<<endl;
       axis2[0] = sin(phi) * cos(theta) + ligcenterx;
       axis2[1] = sin(phi) * sin(theta) + ligcentery;
       axis2[2]= cos(phi) + ligcenterz;
       c=randgeneration(gen);
       float angle=0.0;
       if(c>=0)
          angle = 180.0;
       else
          angle = -180.0;
       //cout<<"rotation angle "<<angle<<endl;
       //cout<<mol.x[0]<<"****"<<endl;
       GroupRotation(axis2,axis1,angle,mol); // random rotate the mol 
       //cout<<mol.x[0]<<"****"<<endl;
       //%%%%%%%%%%%%%%%%add the large rigid movement
       //cout<<i<<"initial conformation "<<mol.x[0]<<" before one "<<matchvec[i].x[0]<<endl; // different is correct.
       //translate the center to calculate the energy

       //double t[3];
       //t[0]=spherecenterx-ligcenterx;
       //t[1]=spherecentery-ligcentery;
       //t[2]=spherecenterz-ligcenterz;
       //for(int atom=0;atom<mol.num_atoms;atom++)
       //{
       //   mol.x[atom] +=t[0];
       //   mol.y[atom] +=t[1];
       //   mol.z[atom] +=t[2];
       //}    
       MC_rigid_ligand(amber, mol,receptor, c_nrg, filevdw.c_str(), fileflex.c_str(),fileflex_drive_tbl.c_str(),grid_file,1.00, 1.00,gen,biolip_matrix,sphx,sphy,sphz,bindsite_weight,protein_atom_index,ave,std, theDensityMap);
       //cout<<i<<"after MC "<<mol.x[0]<<" "<<mol.current_score<<" before "<<matchvec[i].x[0]<<" "<<matchvec[i].current_score<<endl;
       if (!(energy(amber, mol, receptor, c_nrg,filevdw.c_str(), fileflex.c_str(),fileflex_drive_tbl.c_str(),grid_file,cutoff,biolip_matrix,sphx,sphy,sphz,bindsite_weight,protein_atom_index,ave,std, theDensityMap)))
       {
          mol.current_score=1e+6;
          mol.internal_energy=1e+6;
          mol.intral_energy=1e+6;
       }
       //cout<<"before"<<matchvec[i].x[0]<<"\t"<<matchvec[i+1].x[0]<<endl;
       for(int m=0;m<mol.num_atoms;m++)
       {
          
          matchvec[i+1].x[m] = mol.x[m];
          //cout<<matchvec[i].x[m]<<endl;
          matchvec[i+1].y[m] = mol.y[m];
          matchvec[i+1].z[m] = mol.z[m];
       }
       
       matchvec[i+1].current_score=mol.current_score;
       matchvec[i+1].internal_energy=mol.internal_energy;
       matchvec[i+1].intral_energy = mol.intral_energy;
       matchvec[i+1].energy_tor_total = mol.energy_tor_total;
    }


/*****************************REMC simulation****************************************************/
//REMC_flag = false; //not do remc
  if(REMC_flag==true)
  {
      vector<DOCKMol> clustermols;
      vector<vector<float> > x,y,z;
      vector<float> accenergy;
      if(REMC(amber, matchvec, receptor,c_nrg,filevdw.c_str(), fileflex.c_str(),fileflex_drive_tbl.c_str(),grid_file,1,x,y,z,accenergy,Tmin,Tmax,outfile,gen,flexible_flag, biolip_matrix,swapnumber,MC_steps,sphx,sphy,sphz,bindsite_weight,protein_atom_index,ave,std, theDensityMap))
      {
	       //string pdbfile = string(argv[1])+".pdb";
	       string molfile = string(argv[1])+".mol2";
	       //ofstream clusterpdb(pdbfile.c_str());
	       ofstream clustermol(molfile.c_str());

        for(int i=0;i<x.size();i++)
        {
	         //cout<<"i="<<i<<endl;
	
	        clustermols.push_back(mol);
	        for(int j=0;j<mol.num_atoms;j++)
	        {
	           //cout<<i<<"\t"<<j<<"\t"<<clustermols[i].atom_names[j]<<endl;
	           //cout<< "ATOM  " << setw(5) << j+1 << " " << setw(4)<<mol.atom_names[j] << " " << setw(3)<<"ALA" <<"  "<<setw(4) << j+1 <<"   " ;
	           //cout<<setiosflags(ios::fixed)<<setprecision(3)<<setw(8)<<x[i][j]<<setw(8)<<y[i][j]<<setw(8)<<z[i][j]<<"\n";
	           clustermols[clustermols.size()-1].x[j]=x[i][j];
	           clustermols[clustermols.size()-1].y[j]=y[i][j];
	           clustermols[clustermols.size()-1].z[j]=z[i][j];
	        }
	        //cout<<"\n";
	        clustermols[clustermols.size()-1].current_score=accenergy[i];
	        //outpdb(clustermols[i],clusterpdb,i);
          stringstream ss;
          ss << i;
          string str=string(argv[1])+ss.str();
	        clustermol<<"########## Name:"<<"\t"<<i<<"\t"<<clustermols[i].current_score<<endl;
	        Write_Mol2(clustermols[i],clustermol,str);  
        }
        cout<<"total conformations "<<x.size()<<endl;
        //no conformations can be generated by REMC,using the initial conformations as the final results
        if(x.size()==0)
        {
	   
          for(int i=0;i<matchvec.size();i++)
	        {
	          //outpdb(matchvec[i],clusterpdb,i);
	          stringstream ss;
            ss << i;
            string str=string(argv[1])+ss.str();
	          clustermol<<"########## Name:"<<"\t"<<i<<"\t"<<matchvec[i].current_score<<endl;
	          Write_Mol2(matchvec[i],clustermol,str);  
	        } 
        }
        cout<<"REMC simulation is successful!"<<endl;
        //clusterpdb.close();
        //clustermol.close();	
      }
      else
      {
        cout<<"MC simulation is not successful! Don't cry! fighting!"<<endl;
        REMC_flag=false;
      }
  }
	
  if(REMC_flag==false)
  {
    //string pdbfile = string(argv[1])+".pdb";
	  string molfile = string(argv[1])+".mol2";
	  //ofstream clusterpdb(pdbfile.c_str());
	  ofstream clustermol(molfile.c_str());
    //outpdb(mol,clusterpdb,1);
    for(int i=0;i<matchvec.size();i++)
    {
      stringstream ss;
      ss << i;
      string str=string(argv[1])+ss.str();
	    clustermol<<"########## Name:"<<"\t"<<i<<"\t"<<matchvec[i].current_score<<endl;
	    Write_Mol2(matchvec[i],clustermol,str);
    }
        
  } 
  return 1;
 }



