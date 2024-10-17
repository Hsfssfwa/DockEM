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
#include "precon.h"
#include "initialmol.h"
#include "flex_dock.h"
#include "generaterandnum.h"
using namespace std;


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
    
    float cutoff = 100.0;
    bool flexible_flag =1; //whethe do flexible docking
    //bool flexible_flag = 0;

    char * ligand_in = argv[1]; //ligand mol2 file
    char * receptor_mol = argv[2];// receptor mol2 file
    string grid_file = argv[3]; // receptor grid file
    
    //const char * ligand_in = "lig_charge.mol2";
    //char * ligand_in = argv[8];
    
    //const char * receptor_mol = "rec_ITASSER.mol2";
    //char * receptor_mol = argv[10];
    //const char * receptor_mol = "rec_native.mol2";
    //char * receptor_in = argv[5]; // is a .sph file
    
    string filevdw = "/home/zhangwy/program/dock6/parameters/vdw_AMBER_parm99.defn";
    string fileflex = "/home/zhangwy/program/dock6/parameters/flex.defn";
    string fileflex_drive_tbl = "/home/zhangwy/program/dock6/parameters/flex_drive.tbl";
    //string filevdw = "../../parameters/vdw_AMBER_parm99.defn";
    //string fileflex = "../../parameters/flex.defn";
    //string fileflex_drive_tbl = "../../parameters/flex_drive.tbl";
    //$$$$$$$$$$$$$$$$$$$$read the BioLip_ligand_torsion_lib$$$$$$$$$$$$$$$$$
    //ifstream biolip_file("/home/zhangwy/program/EDock2019/edock_flex/BioLip_ligand_torsion_lib");
    //ifstream biolip_file("../../torsion/BioLip_ligand_torsion_lib");
    

    int initialnum=40; //REMC number
    int swapnumber = 100;
    int MC_steps = 201;

    //char * vdw_weight_c = argv[6]; //user input vdw weight
    //char * vdw_weight_c = "1";
    //float vdw_weight = atof(vdw_weight_c);
    float vdw_weight = 1.0;
    //char * bindsite_weight_c = argv[7];//user input binding sites weight
    //float  bindsite_weight = atof(bindsite_weight_c);
    c_nrg.vdw_scale = vdw_weight;



    /*******************************************************************************************/
    //read the ligand informaiton about the title, mol_info_line, comment1, comment2, comment3,     
    //score_text_data,energy, simplex_text, mol_data, num_atoms, num_bonds, num_residues, and so on.
    //and compulate the ring information
    
    
    cout<<"********read ligand part*****************"<<endl;
    readlig(mol,ligand_in);
    readlig(receptor, receptor_mol);

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
    cout <<"flexible_flag:"<<flexible_flag<<endl;
    for (int tt=0;tt<torsions.size();tt++)
    {
        cout<<tt<<"tpye:"<<amber.bond_typer.types[amber.bond_typer.flex_ids[torsions[tt].bond_num]].drive_id<<" ihedral"<<torsions[tt].ihedral<<endl;
    }
    
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
    
    if((energy(amber, mol, receptor,c_nrg,filevdw.c_str(), fileflex.c_str(),fileflex_drive_tbl.c_str(),grid_file,cutoff)))
    {     
      cout<<"native ligand energy is:"<<mol.current_score<<"\t"<<"internal_energy is:" <<mol.internal_energy<<"\t"<<"intral_energy is:"<<mol.intral_energy<<endl;   
    }    
    
   
  return 1;
}



