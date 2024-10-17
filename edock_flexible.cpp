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

// ####################################################################
  /*
    // template part for EDock ############################################
    //read dat file**************************************
    string cofactor_dat = "template_lig_pre.dat";
    
    vector <vector <float> > cofactor_x;
    vector <vector <float> > cofactor_y;
    vector <vector <float> > cofactor_z;
    vector <vector <int> > cofactor_atomindex;
    
    read_dat(cofactor_dat.c_str(), cofactor_x,cofactor_y,cofactor_z,cofactor_atomindex);
    cout<<"read_dat finish"<<endl;
    
    // test read dat 
    //cout<<"x.size()"<<cofactor_x.size()<<" "<<"y.size()"<<cofactor_y.size()<<" "
    //  <<"z.size()"<<cofactor_z.size()<<" "<<"atomindex.size()"<<cofactor_atomindex.size()<<endl; 
    //for(int i=0;i<cofactor_atomindex.size();i++)
    //{
        //for(int j=0;j<cofactor_atomindex[i].size();j++)
        //{
            //cout<<cofactor_atomindex[i][j]<<"**"<<cofactor_x[i][j]<<"**"<<cofactor_y[i][j]<<"**"<<cofactor_z[i][j]<<endl;
        //}
        //cout<<"TRE"<<endl;
    //}
    //**************calculate dij**********************
    
    vector<vector< vector<float> > > distance; 
    vector<vector <float> > ave;
    vector<vector <float> > std;
    
    distance.resize(protein_atom_index.size());
    ave.resize(protein_atom_index.size(),vector <float> (ligand_non_H, 0.0));
    std.resize(protein_atom_index.size(),vector <float> (ligand_non_H, 0.0));
    for (int i=0;i<protein_atom_index.size();i++)
    {
        distance[i].resize(ligand_non_H);
    }
    //cout<<distance[0].size()<<endl;   
       
    for (int protein_i=0;protein_i<protein_atom_index.size();protein_i++)
    {
        for(int template_j=0;template_j<cofactor_atomindex.size();template_j++)
        {
            for(int ligand_m= 0;ligand_m<cofactor_atomindex[template_j].size();ligand_m++)
            {
                int p_index = protein_atom_index[protein_i];
                float distance2= 
                    (receptor.x[p_index] - cofactor_x[template_j][ligand_m])*(receptor.x[p_index] - cofactor_x[template_j][ligand_m])\
                    +(receptor.y[p_index] - cofactor_y[template_j][ligand_m])*(receptor.y[p_index] - cofactor_y[template_j][ligand_m])\
                    +(receptor.z[p_index] - cofactor_z[template_j][ligand_m])*(receptor.z[p_index] - cofactor_z[template_j][ligand_m]);
                float distance_temp = sqrt(distance2);
                //cout<<receptor.x[p_index]<<" "<<receptor.y[p_index]<<" "<<receptor.z[p_index]<<endl;
                //cout<<cofactor_x[template_j][ligand_m]<<" "<<cofactor_y[template_j][ligand_m]<<" "<<cofactor_z[template_j][ligand_m]<<endl;
                //cout<<"protein "<<protein_i<<"ligand "<<cofactor_atomindex[template_j][ligand_m]-1<<endl;
                distance[protein_i][cofactor_atomindex[template_j][ligand_m]-1].push_back(distance_temp);
                //cout<<distance_temp<<" ";
                //cout<<cofactor_atomindex[template_j][ligand_m]<<endl;
                //break;
            }
            //break;
        }   
        //cout<<"TER"<<endl;
        //break;
    }
    for(int i=0;i<distance.size();i++)
    { 
        for (int j=0;j<distance[i].size();j++)
        {
            if(distance[i][j].size()<=1)
            {
                ave[i][j] = -1.0;
                //std[i][j] = -1.0;
                continue;
            }

            for(int m=0; m<distance[i][j].size();m++)
            {
                //cout<<"protein "<<i<<"ligand "<<j<<" "<<distance[i][j].size()<<" ";
                //cout<<distance[i][j][m]<<" ";
                ave[i][j] = ave[i][j] + distance[i][j][m];
            }
            //cout<<"\n";
            ave[i][j] = ave[i][j]/distance[i][j].size();
            //cout<<i<<" "<<j<<" ave "<<ave[i][j]<<endl;
            //break;
            
            //cout<<"TER"<<endl;
            //cout<<distance[i][j]<<endl;
        }
        //break;
    }
    for (int i = 0;i<distance.size();i++)
    {
        for (int j=0; j<distance[i].size();j++)
        {
            if(distance[i][j].size()<=1)
            {
                std[i][j] = -1.0;
                continue;
            }

            for (int m=0;m<distance[i][j].size();m++)
            {
                std[i][j] = std[i][j] + (distance[i][j][m] -ave[i][j])*(distance[i][j][m] -ave[i][j]);
            }
            std[i][j] = std[i][j]/(distance[i][j].size()-1);
            //std[i][j] = sqrt(std[i][j]);
            //cout<<i<<" "<<j<<" std "<<std[i][j]<<endl;
        }
    }

*/
/*******************************************************************************************/
    
    /*******************************************************************************************/
    //this part read some parameters, if input is native ligand, can output the energy, 
    //this part is used    
    //for train the energy 
    //add by wenyizhang 07052017
    
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
    //for (int tt=0;tt<torsions.size();tt++)
    //{
    //    cout<<tt<<"tpye:"<<amber.bond_typer.types[amber.bond_typer.flex_ids[torsions[tt].bond_num]].drive_id<<" ihedral"<<torsions[tt].ihedral<<endl;
    //}
    
    /******idea 13*******/
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
	   /*    	    
      {
	    allmatchcount++;
	    //cout<<(energy(amber, mol, c_nrg,filevdw.c_str(), fileflex.c_str(),fileflex_drive_tbl.c_str(),grid_file,cutoff))<<endl;
	    cout<<"allmatchcout "<<allmatchcount<<endl;
	    outpdb(mol,allmatchmol,allmatchcount); 
      }
	   */        
    }
    //allmatchmol.close(); 
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

            if (!(energy1(amber, moltemp, receptor, c_nrg,filevdw.c_str(), fileflex.c_str(),fileflex_drive_tbl.c_str(),grid_file,cutoff,biolip_matrix,sphx,sphy,sphz,bindsite_weight,protein_atom_index,ave,std, theDensityMap)))
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
             if (!(energy1(amber, mol, receptor,c_nrg,filevdw.c_str(), fileflex.c_str(),fileflex_drive_tbl.c_str(),grid_file,cutoff,biolip_matrix,sphx,sphy,sphz,bindsite_weight,protein_atom_index,ave,std, theDensityMap)))
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
    cout << "CCCCCCCCCC" << theDensityMap.matchcentroidposex(mol) << endl;
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
    //theDensityMap.cart(centers1, fartbinding, max);
    cout << "centers1:  " << centers1[0] << "  " << centers1[1] << "  " << centers1[2] << endl;
    
    for (int i = 0; i < tempmachmols.size(); i++)
    {
        if ((energy(amber, tempmachmols[i], receptor, c_nrg, filevdw.c_str(), fileflex.c_str(), fileflex_drive_tbl.c_str(), grid_file, cutoff, biolip_matrix, sphx, sphy, sphz, bindsite_weight, protein_atom_index, ave, std, theDensityMap)))
        {
            cout << "check me 000000000   " << tempmachmols[i].current_score << "   " << tempmachmols[i].current_EM_score << endl;
        }
    }
    sortmol_em(tempmachmols, molindex1);

    for(int i=0;i<initialnum;i++)
    {
      matchvec.push_back(tempmachmols[molindex1[i]]);
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
       //cout<<"after"<<i+1<<" "<<matchvec[i+1].x[0]<<"\t"<<i<<" "<<matchvec[i].x[0]<<endl;

       //vector <float>(axis1).swap(axis1);  
       //vector <float>(axis2).swap(axis2); 
    }
   
  ///////////////////////////////////////////////////////
  ////////////////////////////////////////
  /*
   string matchfile = "match_" + outfile+".pdb";
   ofstream matchmol(matchfile.c_str());
   for(int i=0;i<matchvec.size();i++)
   {
      outpdb(matchvec[i],matchmol,i); 
   }
   matchmol.close();
  */
  //string mol2file="match_" + outfile + ".mol2";
  //ofstream mol2(mol2file.c_str());
  //for(int i=0;i<matchvec.size();i++)
  //{
  //  stringstream ss;
  //  ss << i;
  //  string str=ss.str();
  //  Write_Mol2(matchvec[i], mol2, str);
  //}
  //mol2.close();
  //############################hbond testing################################
  //
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



