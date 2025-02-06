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
    string inputMRC = argv[5];
    float resolution = atof(argv[7]);
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
    //cout << theDensityMap.density.u3() << "  " << theDensityMap.density.u2() << "  " << theDensityMap.density.u1()<< endl;
    /*for (int i = 1; i <= theDensityMap.density.u1(); i++)
    {
        for (int j = 1; j <= theDensityMap.density.u2(); j++)
        {
            for (int k = 1; k <= theDensityMap.density.u3(); k++)
            {
                cout<< theDensityMap.density(i, j, k)<<"  ";
            }
            cout << endl;
        }
        cout << endl <<i<< endl;
    }*/

    //char * smin=argv[4]; // min tempreture for remc
    //char * smax=argv[5]; // max tempreture for remc
    //char * vdw_weight_c = argv[6]; //user input vdw weight
    //char * bindsite_weight_c = argv[7];//user input binding sites weight
    //float Tmin=atof(smin);
    //float Tmax=atof(smax);
    float Tmin=1.0;
    float Tmax=50.0;//60
    float cutoff = 100.0; // total energy cutoff
    //float vdw_weight = atof(vdw_weight_c);
    float vdw_weight = 0.001;
    //float  bindsite_weight = atof(bindsite_weight_c);
    float  bindsite_weight = 1;
    c_nrg.vdw_scale = vdw_weight;
    int initialnum=40; //REMC number
    int swapnumber = 25;//100
    int MC_steps = 500;//201
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

    {

        vector<int> gridp(3, 0);
        gridp[0] = theDensityMap.grid[0];
        gridp[1] = theDensityMap.grid[1];
        gridp[2] = theDensityMap.grid[2];
        vector<float> originp(3, 0.0);
        originp[0] = theDensityMap.origin[0];
        originp[1] = theDensityMap.origin[1];
        originp[2] = theDensityMap.origin[2];

        int tmp3 = 0;
        float init_E0 = 0.0;
        float init_E1 = 0.0;
        float KTp = 0.001;

        float old_CCp = 0.0;
        float new_CCp = 0.0;
        float out_sup_all_CC = 2.0;
        //CC 0.05 0.001
        float supp_KT0sp = 0.05;//5  //;0.0001;1.5
        float supp_KT0ep = 0.001;//;0.00001;0.001
        //1 0.01
        int supp_num1p = 6;
        int supp_num3p = 250;
        int supp_num4p = supp_num1p;
        int supp_num2p = 5;//12 20
        // 300 50
        vector<DOCKMol> decstr_tmp;
        vector<vector<float>> E_Tx(supp_num4p);

        vector<float> Eng_T(supp_num1p, 0.0);
        int nnx = 0;
        int nny = 0;

        vector<float> best_Ex(supp_num1p, 100000.0);
        vector<float> best_Ex_4(4, 100000.0);
        vector<DOCKMol> best_modelx;//best_modelx(supp_num1);
        vector<DOCKMol> best_modelx_4(4);
        //vector<DOCKMol> best_model_u;
        DOCKMol best_model_up;
        copy_molecule(best_model_up, receptor);
        vector<float> coor_pdbp(3, 0.0);

        vector<float> coor_pro(3, 0.0);
        coor_pro = vector<float>(3, 0.0);
        for (int j = 0; j < receptor.num_atoms; j++)
        {
            coor_pro[0] = coor_pro[0] + receptor.x[j];
            coor_pro[1] = coor_pro[1] + receptor.y[j];
            coor_pro[2] = coor_pro[2] + receptor.z[j];
        }
        coor_pro[0] = coor_pro[0] / (receptor.num_atoms);
        coor_pro[1] = coor_pro[1] / (receptor.num_atoms);
        coor_pro[2] = coor_pro[2] / (receptor.num_atoms);
        cout << "coor_pdb: " << coor_pro[0] << " " << coor_pro[1] << " " << coor_pro[2] << endl;

        vector<float> frax_pdbp;
        MatrixTimesTransVector(theDensityMap.c2f, coor_pro, frax_pdbp);
        vector<float> frac_atmp(3, 0.0);
        frac_atmp[0] = (double(frax_pdbp[0] * gridp[0] - originp[0] + 1));
        frac_atmp[1] = (double(frax_pdbp[1] * gridp[1] - originp[1] + 1));
        frac_atmp[2] = (double(frax_pdbp[2] * gridp[2] - originp[2] + 1));
        vector<float> del_tranp(3, 0.0);
        del_tranp[2] = (frac_atmp[2] - float(gridp[2]) / 2.0) / float(gridp[2]);
        del_tranp[1] = (frac_atmp[1] - float(gridp[1]) / 2.0) / float(gridp[1]);
        del_tranp[0] = (frac_atmp[0] - float(gridp[0]) / 2.0) / float(gridp[0]);
        vector<float> tran0p;
        MatrixTimesTransVector(theDensityMap.f2c, del_tranp, tran0p);
        cout << "tran0: " << tran0p[0] << " " << tran0p[1] << " " << tran0p[2] << endl;
        for (int j = 0; j < receptor.num_atoms; j++) // not last N
        {
            receptor.x[j] = -tran0p[0] + receptor.x[j];
            receptor.y[j] = -tran0p[1] + receptor.y[j];
            receptor.z[j] = -tran0p[2] + receptor.z[j];
        }

        DOCKMol pointsBxp;
        copy_molecule(pointsBxp, receptor);
        float best_E_up = 10;
        vector<float> supp_REMCp(supp_num1p, 0.0);
        for (int i = 0; i < supp_num1p; i++)
        {
            float supp_KT0xp = pow(float(supp_KT0ep / supp_KT0sp), float(float(i) / float(supp_num1p)));
            supp_REMCp[i] = supp_KT0sp * float(supp_KT0xp);
        }
        vector<DOCKMol > decstrpp;
        //    float ang_cus[] = {0.0,0.0,0.0,60.0,60.0,60.0,100.0,100.0,100.0,150.0,150.0,150.0,200.0,200.0,200.0,270.0,270.0,270.0,330.0,330.0,330.0};
        for (int kkk = 0; kkk < supp_num2p; kkk++) // 2
        {
            vector<DOCKMol > decstrz;
            vector<DOCKMol >().swap(decstrz);
            vector<float> E_REMC(supp_num1p, 0.0);
            for (int jjj = 0; jjj < supp_num1p; jjj++)
            {
                //vector<DOCKMol> point_mat;
                DOCKMol point_mat;
                copy_molecule(point_mat, receptor);
                if (kkk > 0)
                {
                    //point_mat.push_back(decstrp[jjj]);
                    cout << point_mat.num_atoms <<"      "<< decstrpp[jjj].num_atoms << endl;
                    swapmol(point_mat, decstrpp[jjj]);
                    //point_mat = decstrp[jjj];
                    coor_pdbp = vector<float>(3, 0.0);
                    for (int j = 0; j < point_mat.num_atoms; j++)
                    {
                        coor_pdbp[0] = coor_pdbp[0] + point_mat.x[j];
                        coor_pdbp[1] = coor_pdbp[1] + point_mat.y[j];
                        coor_pdbp[2] = coor_pdbp[2] + point_mat.z[j];
                    }
                    coor_pdbp[0] = coor_pdbp[0] / (point_mat.num_atoms);
                    coor_pdbp[1] = coor_pdbp[1] / (point_mat.num_atoms);
                    coor_pdbp[2] = coor_pdbp[2] / (point_mat.num_atoms);
                }
                else
                {
                    coor_pdbp = vector<float>(3, 0.0);
                    //point_mat.push_back(pointsBx[0]);
                    swapmol(point_mat, pointsBxp);
                    for (int j = 0; j < point_mat.num_atoms; j++)
                    {
                        coor_pdbp[0] = coor_pdbp[0] + point_mat.x[j];
                        coor_pdbp[1] = coor_pdbp[1] + point_mat.y[j];
                        coor_pdbp[2] = coor_pdbp[2] + point_mat.z[j];
                    }
                    coor_pdbp[0] = coor_pdbp[0] / (point_mat.num_atoms);
                    coor_pdbp[1] = coor_pdbp[1] / (point_mat.num_atoms);
                    coor_pdbp[2] = coor_pdbp[2] / (point_mat.num_atoms);
                }

                DOCKMol fin_matx;
                copy_molecule(fin_matx, receptor);
                swapmol(fin_matx, point_mat);
                generatationmol(fin_matx, point_mat, 2, gen);

                float zang = float(randIntCustom(30, 180));
                //            float zang=ang_cus[jjj];
                float zang_rd = zang * M_PI / 180.0;
                coor_pdbp = vector<float>(3, 0.0);
                for (int j = 0; j < point_mat.num_atoms; j++)
                {
                    coor_pdbp[0] = coor_pdbp[0] + point_mat.x[j];
                    coor_pdbp[1] = coor_pdbp[1] + point_mat.y[j];
                    coor_pdbp[2] = coor_pdbp[2] + point_mat.z[j];
                }
                coor_pdbp[0] = coor_pdbp[0] / (point_mat.num_atoms);
                coor_pdbp[1] = coor_pdbp[1] / (point_mat.num_atoms);
                coor_pdbp[2] = coor_pdbp[2] / (point_mat.num_atoms);

                if (kkk == 0)
                {
                    //best_modelx[jjj] = point_mat[0];
                    DOCKMol best;
                    copy_molecule(best, receptor);
                    swapmol(best, point_mat);
                    best_modelx.push_back(best);
                }
                float E_500 = 0.0;
                int E_500_int = 0;
                KTp = supp_REMCp[jjj];


                vector<float> centers_R1(3, 0);
                for (int i = 0; i < point_mat.num_atoms; i++)
                {
                    centers_R1[0] += point_mat.x[i];
                    centers_R1[1] += point_mat.y[i];
                    centers_R1[2] += point_mat.z[i];
                }
                centers_R1[0] = centers_R1[0] / point_mat.num_atoms;
                centers_R1[1] = centers_R1[1] / point_mat.num_atoms;
                centers_R1[2] = centers_R1[2] / point_mat.num_atoms;

                //CC_i2 = theDensityMap.clCC_2(point_mat, receptor_mapsit);
                float CC_i2p = theDensityMap.matchpose(point_mat);
                new_CCp = 1 - CC_i2p;
                //new_CC2 = 1.0 - CC_i2;
                old_CCp = new_CCp;
                //cout<<"new_CC2: "<<new_CC2<<"  " << kkk <<"    " << jjj << endl;
                int acc_rate1 = 0;
                for (int iii = 0; iii < supp_num3p; iii++)
                {
                    
                    DOCKMol fin_mat;
                    copy_molecule(fin_mat, receptor);
                    swapmol(fin_mat, point_mat);

                    DOCKMol fin_matx;
                    copy_molecule(fin_matx, receptor);
                    swapmol(fin_matx, fin_mat);
                    tmp3 = 0;
                    int pp_k = fin_mat.num_atoms;

                    generatationmol(fin_matx, fin_mat, 2, gen);
                    //CC_i2 = theDensityMap.clCC_2(fin_mat, receptor_mapsit);
                    float CC_i2p = theDensityMap.matchpose(fin_mat);



                    vector<float> centers_R(3, 0);
                    for (int i = 0; i < fin_mat.num_atoms; i++)
                    {
                        centers_R[0] += fin_mat.x[i];
                        centers_R[1] += fin_mat.y[i];
                        centers_R[2] += fin_mat.z[i];
                    }
                    centers_R[0] = centers_R[0] / fin_mat.num_atoms;
                    centers_R[1] = centers_R[1] / fin_mat.num_atoms;
                    centers_R[2] = centers_R[2] / fin_mat.num_atoms;


                    new_CCp = 1 - CC_i2p;

                    float dEp = new_CCp - old_CCp;
                    if (new_CCp < old_CCp)// && new_CC<old_CC)
                    {
                        cout << "old_CC1,new_CC1:RRRRRRRRRRRRR "  << new_CCp << "  " << old_CCp << endl;
                        //tmp3 = 0;
                        coor_pdbp = vector<float>(3, 0.0);

                        //point_mat.clear();
                        //point_mat.push_back(fin_mat[0]);
                        swapmol(point_mat, fin_mat);
                        //point_mat = fin_mat;
                        for (int j = 0; j < fin_mat.num_atoms; j++)
                        {
                            coor_pdbp[0] = coor_pdbp[0] + fin_mat.x[j];
                            coor_pdbp[1] = coor_pdbp[1] + fin_mat.y[j];
                            coor_pdbp[2] = coor_pdbp[2] + fin_mat.z[j];
                            ///tmp3 = tmp3 + 1;
                        }
                        coor_pdbp[0] = coor_pdbp[0] / (fin_mat.num_atoms);
                        coor_pdbp[1] = coor_pdbp[1] / (fin_mat.num_atoms);
                        coor_pdbp[2] = coor_pdbp[2] / (fin_mat.num_atoms);

                        if (new_CCp < best_Ex[jjj])
                        {
                            best_Ex[jjj] = new_CCp;
                            //vector<poseCoord>().swap(best_modelx[jjj]);
                            //best_modelx[jjj] = fin_mat;
                            swapmol(best_modelx[jjj], fin_mat);
                        }
                        if (new_CCp < best_E_up)
                        {
     
                            swapmol(best_model_up, fin_mat);

                            best_E_up = new_CCp;
                            cout << "best" << endl;
                            cout << "AAAAAA" << best_model_up.x[0] << "  " << best_model_up.y[0] << "  " << best_model_up.z[0] << endl;
                        }


                        old_CCp = new_CCp;
                        E_500 = new_CCp;
                        E_500_int = 1;

                    }
                    else if (new_CCp == old_CCp) continue;
                    else
                    {
                        //    float tmpx=rand()/double(RAND_MAX);
                        float tmpx = randf0and1();
                        float mc_v = exp(-dEp / (KTp));
                        if (tmpx < mc_v) // CC must be >0
                        {
                            cout << "RRRRRRRRRRRRR old_CC1,new_CC1: " << new_CCp << " "<< old_CCp << endl;
                            //cout << "1301" << endl;
                            acc_rate1 = acc_rate1 + 1;
                            //cout << acc_rate1 << endl;
                            tmp3 = 0;
                            coor_pdbp = vector<float>(3, 0.0);
                            swapmol(point_mat, fin_mat);
                            //point_mat.clear();
                            //point_mat.push_back(fin_mat[0]);
                            //point_mat = fin_mat;
                            for (int j = 0; j < fin_mat.num_atoms; j++)
                            {
                                coor_pdbp[0] = coor_pdbp[0] + fin_mat.x[j];
                                coor_pdbp[1] = coor_pdbp[1] + fin_mat.y[j];
                                coor_pdbp[2] = coor_pdbp[2] + fin_mat.z[j];
                                tmp3 = tmp3 + 1;
                            }
                            coor_pdbp[0] = coor_pdbp[0] / (fin_mat.num_atoms);
                            coor_pdbp[1] = coor_pdbp[1] / (fin_mat.num_atoms);
                            coor_pdbp[2] = coor_pdbp[2] / (fin_mat.num_atoms);
                            //    cout<<"coor: "<<coor_pdb[0]<<" "<<coor_pdb[1]<<" "<<coor_pdb[2]<<endl;
                            //    cout<<"old_CC1,new_CC1:XX "<<old_dE<<" "<<new_dE<<endl;                     

                            old_CCp = new_CCp;
                            E_500 = new_CCp;
                            E_500_int = 1;
                            //    new_CC2 = new_CC1;

                        }
                    }
                    nnx = nnx + 1;
                    //cout << "BBBBB" << best_model_u.x[0] << "  " << best_model_u.y[0] << "  " << best_model_u.z[0] << endl;
                }
                cout << "acc_rate:  " << acc_rate1 << endl;
                cout << best_E_up << endl;
                if (E_500_int == 1)
                {
                    E_REMC[jjj] = E_500;
                }
                decstrz.push_back(point_mat);  // change    
                cout << "end" << endl;
            }
            if (best_E_up < 0.15)
            {
                cout << "find min CC" << endl;
                break;
            }
            vector<DOCKMol >().swap(decstrpp);
            decstrpp = vector<DOCKMol >(supp_num1p);

            int js = kkk % 2;
            if (js == 0)
            {
                for (int i = 0; i < supp_num1p - 1; i = i + 2)
                {
                    int j = i + 1;
                    double CH_REMC = exp((1.0 / supp_REMCp[i] - 1.0 / supp_REMCp[j]) * (E_REMC[i] - E_REMC[j]));
                    float Pglobal = min(1.0, CH_REMC);
                    float tmpx = randf0and1();
                    if (tmpx < Pglobal)
                    {
                        decstrpp[i] = decstrz[j];
                        decstrpp[j] = decstrz[i];
                        cout << "aaaaaaa" << endl;
                    }
                    else
                    {
                        decstrpp[i] = decstrz[i];
                        decstrpp[j] = decstrz[j];
                    }
                }
            }
            if (js == 1)
            {
                decstrpp[0] = decstrz[0];
                for (int i = 1; i < supp_num1p - 1; i = i + 2)
                {
                    int j = i + 1;
                    double CH_REMC = exp((1.0 / supp_REMCp[i] - 1.0 / supp_REMCp[j]) * (E_REMC[i] - E_REMC[j]));
                    float Pglobal = min(1.0, CH_REMC);
                    float tmpx = randf0and1();
                    if (tmpx < Pglobal)
                    {
                        decstrpp[i] = decstrz[j];
                        decstrpp[j] = decstrz[i];
                        cout << "aaaaaaa" << endl;
                    }
                    else
                    {
                        decstrpp[i] = decstrz[i];
                        decstrpp[j] = decstrz[j];
                    }
                }
                decstrpp[supp_num1p - 1] = decstrz[supp_num1p - 1];
            }
        }
        swapmol(receptor, best_model_up);
    }
    string molfile11 = "test_protein.mol2";
    //ofstream clusterpdb(pdbfile.c_str());
    ofstream clustermol11(molfile11.c_str());
    string str11 = string(argv[1]);
    clustermol11 << "########## Name:" << endl;
    Write_Mol2(receptor, clustermol11, str11);

    //float X, Y, Z;
    //vector<float> x, y, z;
    ////cout << argv[12] << endl;
    //string keyfile = string(argv[12]) + ".pdb";
    //ifstream inputFile(keyfile.c_str());
    //if (inputFile.is_open())
    //{
    //    string PDBline;
    //    while (getline(inputFile, PDBline))
    //    {
    //        if (PDBline.substr(0, 4) == "ATOM")
    //        {
    //            X = (strtof(PDBline.substr(30, 8).c_str(), NULL));
    //            Y = (strtof(PDBline.substr(38, 8).c_str(), NULL));
    //            Z = (strtof(PDBline.substr(46, 8).c_str(), NULL));
    //            x.push_back(X);
    //            y.push_back(Y);
    //            z.push_back(Z);
    //        }
    //        else continue;
    //    }
    //    inputFile.close();
    //}
    //else std::cout << "文件打开失败" << endl;
    //cout << "centers7" << x.size() << endl;
    //vector<float> centers7(3, 0);
    //for (int i = 0; i < x.size(); i++)
    //{
    //    centers7[0] += x[i];
    //    centers7[1] += y[i];
    //    centers7[2] += z[i];
    //}
    //centers7[0] = centers7[0] / x.size();
    //centers7[1] = centers7[1] / y.size();
    //centers7[2] = centers7[2] / z.size();
    //cout << "ccccccc" << centers7[0] << "  " << centers7[1] << "  " << centers7[2] << endl;

    float X_1, Y_1, Z_1;
    vector<float> x_cen, y_cen, z_cen;
    string keyfile1 = string(argv[5]) + ".pdb";
    ifstream inputFile1(keyfile1.c_str());
    if (inputFile1.is_open())
    {
        string PDBline;
        while (getline(inputFile1, PDBline))
        {
            if (PDBline.substr(0, 4) == "ATOM")
            {
                X_1 = (strtof(PDBline.substr(30, 8).c_str(), NULL));
                Y_1 = (strtof(PDBline.substr(38, 8).c_str(), NULL));
                Z_1 = (strtof(PDBline.substr(46, 8).c_str(), NULL));
                x_cen.push_back(X_1);
                y_cen.push_back(Y_1);
                z_cen.push_back(Z_1);
            }
            else continue;
        }
        inputFile1.close();

    }
    else std::cout << "文件打开失败" << endl;

    //cout << "centers7" << x.size() << endl;
    vector<float> centers7(3, 0);
    for (int i = 0; i < x_cen.size(); i++)
    {
        centers7[0] += x_cen[i];
        centers7[1] += y_cen[i];
        centers7[2] += z_cen[i];
    }
    centers7[0] = centers7[0] / x_cen.size();
    centers7[1] = centers7[1] / y_cen.size();
    centers7[2] = centers7[2] / z_cen.size();

    cout << "ccccccc" << centers7[0] << "  " << centers7[1] << "  " << centers7[2] << endl;

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
    vector<poseCoord> receptor_bindingsit;

    int num = binding_sites[0].size();
    for (int i = 0; i < num; i++)
    {
        int idx = 0;
        bool q = 1;
        idx = binding_sites[0][i] - 1;
        int num1 = binding_sites[0].size();
        for (int j = 0; j < num1; j++)
        {
            if (idx == binding_sites[0][j])
            {
                q = 0;
                break;
            }
        }
        if (q)
        {
            binding_sites[0].push_back(idx);
        }

        q = 1;
        int num2 = binding_sites[0].size();
        idx = binding_sites[0][i] + 1;
        for (int j = 0; j < num2; j++)
        {
            if (idx == binding_sites[0][j])
            {
                q = 0;
                break;
            }
        }
        if (q)
        {
            binding_sites[0].push_back(idx);
        }
    }
    for (int j = 0; j < binding_sites[0].size(); j++)
    {
        cout << binding_sites[0][j] << ":";
    }
    cout << endl;

    if(binding_sites[1][0] != 0 && binding_sites[1][1] != 0 && binding_sites[1][2] != 0)
    {
        for (int j = 0; j < binding_sites[1].size(); j++)
        {
            cout << binding_sites[1][j] << ":";
        }
        cout << endl;
    }

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

                poseCoord atom;
                atom.atom_name = receptor.atom_names[j];
                atom.residue_numbers = receptor.atom_residue_numbers[j];
                atom.x = receptor.x[j];
                atom.y = receptor.y[j];
                atom.z = receptor.z[j];
                receptor_bindingsit.push_back(atom);
            }   
        }
    }
    //for (int i = 0; i < receptor_bindingsit.size(); i++)
    //{
    //    cout << receptor_bindingsit[i].atom_name << "  " << receptor_bindingsit[i].residue_numbers << "  " << receptor_bindingsit[i].x << "  " << receptor_bindingsit[i].y << "  " << receptor_bindingsit[i].z << endl;
    //}

    cout<<"protein_atom_num"<<protein_atom_num<<endl;
    cout<<protein_atom_index[0]<<endl;
    //save corresponding all atoms in bindingsite**********************
     
    cartbinding[0] = cartbinding[0]/ protein_atom_num;
    cartbinding[1] = cartbinding[1]/ protein_atom_num;
    cartbinding[2] = cartbinding[2]/ protein_atom_num;

    if (binding_sites[1][0] != 0 && binding_sites[1][1] != 0 && binding_sites[1][2] != 0)
    {
        cartbinding[0] = binding_sites[1][0];
        cartbinding[1] = binding_sites[1][1];
        cartbinding[2] = binding_sites[1][2];
    }

    cout << "QQQQQ" << cartbinding[0] << "  " << cartbinding[1] << "  " << cartbinding[2] << endl;
    float car_x = 0, car_y = 0, car_z = 0;
    float max = 0;
    for (int i = 0; i < receptor.num_atoms; i++)
    {
        for (int j = i; j < receptor.num_atoms; j++)
        {
            car_x = (receptor.x[i] - receptor.x[j]);
            car_y = (receptor.y[i] - receptor.y[j]);
            car_z = (receptor.z[i] - receptor.z[j]);
            if (sqrt(car_x* car_x + car_y* car_y + car_z * car_z) > max)
                max = sqrt(car_x * car_x + car_y * car_y + car_z * car_z);
        }
    }
    cout << "max:  " << max << endl;
    //theDensityMap.cart(cartbinding, fartbinding,max);

    
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
        swapnumber = 50;
        MC_steps = 501; 
    }
    
    //c_nrg.grid_file_name = grid_file.c_str();
    //c_nrg.initialize(amber);
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
  
    cout << "ccccccc" << centers[0] << "  " << centers[1] << "  " << centers[2] << endl;

    max = 0;
    for (int i = 0; i < mol.num_atoms; i++)
    {
        for (int j = i; j < mol.num_atoms; j++)
        {
            car_x = (mol.x[i] - mol.x[j]);
            car_y = (mol.y[i] - mol.y[j]);
            car_z = (mol.z[i] - mol.z[j]);
            if (sqrt(car_x * car_x + car_y * car_y + car_z * car_z) > max)
                max = sqrt(car_x * car_x + car_y * car_y + car_z * car_z);
        }
    }
    cout << "max:  " << max << endl;
    //theDensityMap.x_min = 1, theDensityMap.y_min = 1, theDensityMap.z_min = 1;
    //theDensityMap.x_max = theDensityMap.density.u1(), theDensityMap.y_max = theDensityMap.density.u2(), theDensityMap.z_max = theDensityMap.density.u3();
    //
    //cout << "clCC_1clCC_1clCC_1" << theDensityMap.clCC_1(receptor) << endl;
    //cout << "clCC_2clCC_2clCC_2" << theDensityMap.clCC_2(mol, receptor_bindingsit) << endl;
    //cout << "clCC_3clCC_3clCC_3" << theDensityMap.clCC_3(mol, receptor) << endl;
    //cout << endl;
    //theDensityMap.x_min = 17, theDensityMap.y_min = 15, theDensityMap.z_min = 19;
    //theDensityMap.x_max = 22, theDensityMap.y_max = 20, theDensityMap.z_max = 24;
    ////theDensityMap.cart(centers, fartbinding, max);
    //cout << "clCC_1clCC_1clCC_1" << theDensityMap.clCC_1(mol) << endl;
    //cout << "clCC_2clCC_2clCC_2" << theDensityMap.clCC_2(mol, receptor_bindingsit) << endl;
    //cout << "clCC_3clCC_3clCC_3" << theDensityMap.clCC_3(mol, receptor) << endl;
    //cout << "matchcentroidposex" << theDensityMap.matchcentroidposex(mol) << endl;
    //cout << "matchcentroidzzzzz" << theDensityMap.matchcentroidposexz(mol) << endl;


    /*******************************************************************************************/
    // temp definition add by wenyizng 12/04/2019
    //vector<int> protein_atom_index;
    vector<vector <float> > ave;
    vector<vector <float> > std;
    vector<vector< vector<float> > > distance;
    cout<<"*******match part***********************"<<endl;
    //c_orient.match(mol,receptor_in);    
    cout<<"orient fininshed!!!!!!!"<<endl;
    float spherecenterx, spherecentery, spherecenterz;
    spherecenterx = spherecentery = spherecenterz = 0;
    vector<float> sphx;
    vector<float> sphy;
    vector<float> sphz;

    //for (int i = 0; i < c_orient.spheres.size(); i++)
    //{
    //    //cout<<c_orient.spheres[i].crds.x<<"\t"<<c_orient.spheres[i].crds.y<<"\t"<<c_orient.spheres[i].crds.z<<endl;
    //    spherecenterx += c_orient.spheres[i].crds.x;
    //    spherecentery += c_orient.spheres[i].crds.y;
    //    spherecenterz += c_orient.spheres[i].crds.z;
    //    sphx.push_back(c_orient.spheres[i].crds.x);
    //    sphy.push_back(c_orient.spheres[i].crds.y);
    //    sphz.push_back(c_orient.spheres[i].crds.z);
    //}
    //spherecenterx = spherecenterx / c_orient.spheres.size();
    //spherecentery = spherecentery / c_orient.spheres.size();
    //spherecenterz = spherecenterz / c_orient.spheres.size();

    //cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<endl;
    //vector<float> centersa(3, 0);
    ////double t[3];
    //for (int i = 0; i < mol.num_atoms; i++)
    //{
    //    centersa[0] += mol.x[i];
    //    centersa[1] += mol.y[i];
    //    centersa[2] += mol.z[i];
    //}
    //centersa[0] = centersa[0] / mol.num_atoms;
    //centersa[1] = centersa[1] / mol.num_atoms;
    //centersa[2] = centersa[2] / mol.num_atoms;
    vector <DOCKMol> molvec;
    vector <DOCKMol> matchvec;
    /*c_orient.cliques saves the match mols */
    //cout<<c_orient.cliques.size()<<endl;
    int matching_count = 0;
    int allmatchcount = 0;

    cout << "c_orient.cliques.size() is:" << c_orient.cliques.size() << endl;
    //
    int orient_iterations = 0;
    if (c_orient.cliques.size() < 5000 && c_orient.cliques.size() > 0) //it means that the mol can be oriented
        orient_iterations = c_orient.cliques.size();
    else if (c_orient.cliques.size() >= 5000)
        orient_iterations = 5000;
    for (c_orient.current_clique = 0; c_orient.current_clique < orient_iterations; c_orient.current_clique++)
    {
        //cout<<"c_orient.current_clique is:"<<c_orient.current_clique<<endl;
        if (c_orient.new_next_orientation(mol))
        {
            //cout << "aaa" << endl;
            if ((energy1(amber, mol, receptor, c_nrg, filevdw.c_str(), fileflex.c_str(), fileflex_drive_tbl.c_str(), grid_file, cutoff, biolip_matrix, sphx, sphy, sphz, bindsite_weight, protein_atom_index, ave, std, theDensityMap)))
            {
                //cout<<"score is:"<<mol.intral_energy<<endl;
                //if (mol.current_score < 1e+6)
                if (mol.current_score < 1e+6)
                    //if (mol.intral_energy <1e+6)
                {
                    //cout << "check me 000000000   " << mol.current_score << "   " << mol.current_EM_score << endl;
                    molvec.push_back(mol);
                    if (molvec.size() > 1000)
                        break;
                    // molvec save the matching conformations
                    //cout<<molvec[molvec.size()-1].current_score<<endl;
                }
            }
        }
    }
    cout << "checking666" << molvec.size() << endl;
    vector<float> centersa(3, 0);
    //double t[3];
    for (int i = 0; i < mol.num_atoms; i++)
    {
        centersa[0] += mol.x[i];
        centersa[1] += mol.y[i];
        centersa[2] += mol.z[i];
    }
    centersa[0] = centersa[0] / mol.num_atoms;
    centersa[1] = centersa[1] / mol.num_atoms;
    centersa[2] = centersa[2] / mol.num_atoms;

    cout << "dddddddd" << centersa[0] << "  " << centersa[1] << "  " << centersa[2] << endl;
    bool REMC_flag = 1;
    //if (molvec.size() >= initialnum)
    //    REMC_flag = 1;
    //else if (molvec.size() == 0)
    //{
    //    int num_molvec = molvec.size();
    //    cout << "num_molvec " << num_molvec << endl;
    //    vector<vector<float> > ini_x(initialnum - num_molvec, vector<float>(mol.num_atoms, 0.0));
    //    vector<vector<float> > ini_y(initialnum - num_molvec, vector<float>(mol.num_atoms, 0.0));
    //    vector<vector<float> > ini_z(initialnum - num_molvec, vector<float>(mol.num_atoms, 0.0));
    //    REMC_flag = intialmol2(ini_x, ini_y, ini_z, num_molvec, mol, spherecenterx, spherecentery, spherecenterz, initialnum, gen);
    //    //cout<<"ini_x len"<<ini_x.size()<<endl;
    //    for (int i = 0; i < initialnum - num_molvec; i++)
    //    {
    //        DOCKMol moltemp;
    //        copy_molecule(moltemp, mol);
    //        for (int atom = 0; atom < mol.num_atoms; atom++)
    //        {
    //            moltemp.x[atom] = ini_x[i][atom];
    //            moltemp.y[atom] = ini_y[i][atom];
    //            moltemp.z[atom] = ini_z[i][atom];
    //        }
    //        //cout<<ini_x[i][0]<<"\t"<<moltemp.x[0]<<"\t";
    //        //cout<<molvec[molvec.size()-1].x[0]<<endl;
    //        if (!(energy1(amber, moltemp, receptor, c_nrg, filevdw.c_str(), fileflex.c_str(), fileflex_drive_tbl.c_str(), grid_file, cutoff, biolip_matrix, sphx, sphy, sphz, bindsite_weight, protein_atom_index, ave, std, theDensityMap)))
    //        {
    //            moltemp.current_score = 1e+6;
    //            moltemp.internal_energy = 1e+6;
    //            moltemp.intral_energy = 1e+6;
    //        }
    //        molvec.push_back(moltemp);
    //        //cout<<molvec[molvec.size()-1].x[0]<<endl;
    //    }
    //}
    //else if (molvec.size() < initialnum && molvec.size() > 0)
    //{
    //    int num_molvec = molvec.size();
    //    cout << "num_molvec " << num_molvec << endl;
    //    //cout << initialnum << "  " << mol.num_atoms << endl;
    //    vector<float> axis1(3, 0);
    //    vector<float> axis2(3, 0);
    //    double a, b, c;
    //    for (int i = 0; i < initialnum - num_molvec; i++)
    //    {
    //        int index = i % num_molvec;
    //        //cout<<i<<"\t"<<index<<endl;
    //        for (int m = 0; m < mol.num_atoms; m++)
    //        {
    //            mol.x[m] = molvec[index].x[m];
    //            mol.y[m] = molvec[index].y[m];
    //            mol.z[m] = molvec[index].z[m];
    //        }
    //        DOCKMol moltemp;
    //        copy_molecule(moltemp, mol);
    //        molvec.push_back(moltemp);
    //        cout<<"now molvec.size()"<<molvec.size()<<endl;
    //        //cout<<"before mol"<<molvec[0].x[0]<<"\t"<<molvec[0].y[0]<<"\t"<<molvec[0].z[0]<<endl;
    //        {
    //            float ligcenterx, ligcentery, ligcenterz;
    //            ligcenterx = ligcentery = ligcenterz = 0;
    //            for (int m = 0; m < mol.num_atoms; m++)
    //            {
    //                ligcenterx += mol.x[m];
    //                ligcentery += mol.y[m];
    //                ligcenterz += mol.z[m];
    //            }
    //            ligcenterx = ligcenterx / mol.num_atoms;
    //            ligcentery = ligcentery / mol.num_atoms;
    //            ligcenterz = ligcenterz / mol.num_atoms;
    //            //cout<<ligcenterx<<"\t"<<ligcentery<<"\t"<<ligcenterz<<endl;
    //            axis1[0] = ligcenterx; axis1[1] = ligcentery; axis1[2] = ligcenterz;
    //            a = rand0to1(gen);
    //            b = rand0to1(gen);
    //            //cout<<"a"<<a<<" "<<"b"<<b<<" "<<"c"<<c<<endl;
    //            double theta = 2 * 3.1415926 * a;
    //            double phi = acos(2 * b - 1.0);
    //            axis2[0] = sin(phi) * cos(theta) + ligcenterx;
    //            axis2[1] = sin(phi) * sin(theta) + ligcentery;
    //            axis2[2] = cos(phi) + ligcenterz;
    //            c = randgeneration(gen);
    //            float angle = c * 180.0;
    //            GroupRotation(axis2, axis1, angle, mol); // random rotate the mol 
    //            //cout<<"checking5555"<<mol.x[0]<<endl;
    //            //translate the center, some situation is can not calculate the energy
    //            double t[3];
    //            t[0] = spherecenterx - ligcenterx;
    //            t[1] = spherecentery - ligcentery;
    //            t[2] = spherecenterz - ligcenterz;
    //            for (int atom = 0; atom < mol.num_atoms; atom++)
    //            {
    //                mol.x[atom] += t[0];
    //                mol.y[atom] += t[1];
    //                mol.z[atom] += t[2];
    //            }
    //            //cout<<"after mol"<<mol.x[0]<<"\t"<<mol.y[0]<<"\t"<<mol.z[0]<<endl;
    //            if (!(energy1(amber, mol, receptor, c_nrg, filevdw.c_str(), fileflex.c_str(), fileflex_drive_tbl.c_str(), grid_file, cutoff, biolip_matrix, sphx, sphy, sphz, bindsite_weight, protein_atom_index, ave, std, theDensityMap)))              
    //            {
    //                mol.current_score = 1e+6;
    //                mol.internal_energy = 1e+6;
    //                mol.intral_energy = 1e+6;
    //            }
    //            for (int m = 0; m < mol.num_atoms; m++)
    //            {
    //                molvec[molvec.size() - 1].x[m] = mol.x[m];
    //                molvec[molvec.size() - 1].y[m] = mol.y[m];
    //                molvec[molvec.size() - 1].z[m] = mol.z[m];
    //            }
    //            molvec[molvec.size() - 1].current_score = mol.current_score;
    //            molvec[molvec.size() - 1].internal_energy = mol.internal_energy;
    //            molvec[molvec.size() - 1].intral_energy = mol.intral_energy;
    //            molvec[molvec.size() - 1].energy_tor_total = mol.energy_tor_total;
    //            //cout<<"after molvec"<<molvec[molvec.size()-1].x[0]<<endl;
    //            //vector <float>(axis1).swap(axis1);  
    //            //vector <float>(axis2).swap(axis2);
    //        }
    //    }
    //    REMC_flag = 1;
    //}
    //cout << "655" << endl;
    ////cout<<"REMC_flag is:"<<REMC_flag<<endl; 
    //if (initialnum > molvec.size())
    //    initialnum = molvec.size();
    //cout<<"final checking initialnum:"<<initialnum<<endl;
    cout << "660" << endl;
    //vector<int> molindex, molindex1;
    //vector<DOCKMol> tempmachmols;
    //for (int i = 0; i < molvec.size(); i++)
    //{
    //    tempmachmols.push_back(molvec[i]);
    //    //cout<<"checking 9999"<<molvec[i].x[0]<<"\t"<<tempmachmols[tempmachmols.size()-1].x[0]<<endl;
    //    molindex.push_back(i);
    //    molindex1.push_back(i);
    //}
    //sortmol(tempmachmols, molindex);
    ////for (int i = 0; i < initialnum; i++)
    //    //cout << i << tempmachmols[i].x[0] << "   " << tempmachmols[i].current_score << endl;
    //vector<float> centers1(3, 0);
    //for (int i = 0; i < tempmachmols[molindex[0]].num_atoms; i++)
    //{
    //    centers1[0] += tempmachmols[molindex[0]].x[i];
    //    centers1[1] += tempmachmols[molindex[0]].y[i];
    //    centers1[2] += tempmachmols[molindex[0]].z[i];
    //}
    //centers1[0] = centers1[0] / tempmachmols[molindex[0]].num_atoms;
    //centers1[1] = centers1[1] / tempmachmols[molindex[0]].num_atoms;
    //centers1[2] = centers1[2] / tempmachmols[molindex[0]].num_atoms;
    //cout << "aaa" << centers1[0] << "  " << centers1[1] << "  " << centers1[2] << endl;

    vector<float> molcen(3, 0);
    for (int i = 0; i < mol.num_atoms; i++)
    {
        molcen[0] += mol.x[i];
        molcen[1] += mol.y[i];
        molcen[2] += mol.z[i];
    }
    molcen[0] = molcen[0] / mol.num_atoms;
    molcen[1] = molcen[1] / mol.num_atoms;
    molcen[2] = molcen[2] / mol.num_atoms;
    cout << "aaa" << molcen[0] << "  " << molcen[1] << "  " << molcen[2] << endl;

    //theDensityMap.cart(centers1, fartbinding, max);
    //theDensityMap.cart(centers7, fartbinding, max);
    theDensityMap.cart(cartbinding, fartbinding, 0.8*max);

    cout << "  " << endl;
    //theDensityMap.cart(cartbinding, fartbinding, 0.5 * max);
    
    vector<poseCoord> receptor_mapsit, receptor_mapsit_1;
    vector<float> fra_max(3, 0.0),fra_min(3, 0.0);
    vector<float> car_max(3, 0.0), car_min(3, 0.0);
    //fra_max[0] = theDensityMap.x_max, fra_max[1] = theDensityMap.y_max, fra_max[2] = theDensityMap.z_max;
    //fra_min[0] = theDensityMap.x_min, fra_min[1] = theDensityMap.y_min, fra_min[2] = theDensityMap.z_min;
    //MatrixTimesTransVector(theDensityMap.f2c, fra_max, car_max);
    //MatrixTimesTransVector(theDensityMap.f2c, fra_min, car_min);

    //for (int j = 0; j < receptor.num_atoms; j++)
    //{
    //    car_max[0] = receptor.x[j], car_max[1] = receptor.x[j], car_max[2] = receptor.x[j];
    //    MatrixTimesTransVector(theDensityMap.c2f, car_max, fra_max);
    //    // the location of atom in grid ?
    //    vector<float> atm_idxt(3, 0.0);
    //    atm_idxt[0] = (double(fra_max[0] * theDensityMap.grid[0] - theDensityMap.origin[0] + 1));
    //    atm_idxt[1] = (double(fra_max[1] * theDensityMap.grid[1] - theDensityMap.origin[1] + 1));
    //    atm_idxt[2] = (double(fra_max[2] * theDensityMap.grid[2] - theDensityMap.origin[2] + 1));
    //    cout << j << endl;
    //    cout << atm_idxt[0] << "  " << theDensityMap.x_min << "  " << theDensityMap.x_max << endl;
    //    cout << atm_idxt[1] << "  " << theDensityMap.y_min << "  " << theDensityMap.y_max << endl;
    //    cout << atm_idxt[2] << "  " << theDensityMap.z_min << "  " << theDensityMap.z_max << endl;
    //    if (atm_idxt[0]< theDensityMap.x_max + 0.5 && atm_idxt[0] > theDensityMap.x_min - 0.5 &&
    //        atm_idxt[1]< theDensityMap.y_max + 0.5 && atm_idxt[1] > theDensityMap.y_min - 0.5 &&
    //        atm_idxt[2]< theDensityMap.z_max + 0.5 && atm_idxt[2] > theDensityMap.z_min - 0.5)
    //    //if (receptor.x[j]< car_max[0] && receptor.x[j] > car_min[0] &&
    //    //    receptor.y[j]< car_max[1] && receptor.y[j] > car_min[1] &&
    //    //    receptor.z[j]< car_max[2] && receptor.z[j] > car_min[2])
    //    {
    //        poseCoord atom;
    //        atom.atom_name = receptor.atom_names[j];
    //        atom.residue_numbers = receptor.atom_residue_numbers[j];
    //        atom.x = receptor.x[j];
    //        atom.y = receptor.y[j];
    //        atom.z = receptor.z[j];
    //        receptor_mapsit.push_back(atom);
    //    }
    //}
    //cout << car_max[0] << "  " << car_max[1] << "  " << car_max[2] << endl;
    //cout << car_min[0] << "  " << car_min[1] << "  " << car_min[2] << endl;
    for (int j = 0; j < receptor.num_atoms; j++)
    {
        car_max[0] = receptor.x[j], car_max[1] = receptor.x[j], car_max[2] = receptor.x[j];
        //double distance;
        //distance = sqrt((receptor.x[j] - centers7[0]) * (receptor.x[j] - centers7[0]) +
        //    (receptor.y[j] - centers7[1]) * (receptor.y[j] - centers7[1]) +
        //    (receptor.z[j] - centers7[2]) * (receptor.z[j] - centers7[2]));

        //if (abs(receptor.x[j] - centers7[0])<0.6*max && abs(receptor.y[j] - centers7[1])<0.6*max&&
        //    abs(receptor.z[j] - centers7[2])<0.6*max)
        if (abs(receptor.x[j] - cartbinding[0]) < 0.8 * max && abs(receptor.y[j] - cartbinding[1]) < 0.8 * max &&
            abs(receptor.z[j] - cartbinding[2]) < 0.8 * max)
        {
            poseCoord atom;
            atom.atom_name = receptor.atom_names[j];
            atom.residue_numbers = receptor.atom_residue_numbers[j];
            atom.x = receptor.x[j];
            atom.y = receptor.y[j];
            atom.z = receptor.z[j];
            receptor_mapsit.push_back(atom);
        }
    }
    cout << receptor_mapsit.size() << endl;
    string old_res, new_res = "aa";
    for (int i = 0; i < receptor_mapsit.size(); i++)
    {
        old_res = receptor_mapsit[i].residue_numbers;
        if (old_res != new_res)
        {
            cout << ":" << old_res ;
            new_res = old_res;
        }
    }
    cout << endl;

    for (int j = 0; j < receptor.num_atoms; j++)
    {
        car_max[0] = receptor.x[j], car_max[1] = receptor.x[j], car_max[2] = receptor.x[j];
        //double distance;
        //distance = sqrt((receptor.x[j] - centers7[0]) * (receptor.x[j] - centers7[0]) +
        //    (receptor.y[j] - centers7[1]) * (receptor.y[j] - centers7[1]) +receptor_mapsit_1
        //    (receptor.z[j] - centers7[2]) * (receptor.z[j] - centers7[2]));
        if (abs(receptor.x[j] - centers7[0]) < 0.35 * max && abs(receptor.y[j] - centers7[1]) < 0.35 * max &&
            abs(receptor.z[j] - centers7[2]) < 0.35 * max)
        {
            poseCoord atom;
            atom.atom_name = receptor.atom_names[j];
            atom.residue_numbers = receptor.atom_residue_numbers[j];
            atom.x = receptor.x[j];
            atom.y = receptor.y[j];
            atom.z = receptor.z[j];
            receptor_mapsit_1.push_back(atom);
        }
    }
    cout << receptor_mapsit_1.size() << endl;
    string old_res1, new_res1 = "aa";
    for (int i = 0; i < receptor_mapsit_1.size(); i++)
    {
        old_res1 = receptor_mapsit_1[i].residue_numbers;
        if (old_res1 != new_res1)
        {
            cout << ":" << old_res1;
            new_res1 = old_res1;
        }
    }
    cout << endl;
    //float xmin = 1000, xmax = -1000, ymin = 1000, ymax = -1000, zmin = 1000, zmax = -1000;
    //for (int i = 0; i < tempmachmols[molindex[0]].num_atoms; i++)
    //{
    //    if (tempmachmols[molindex[0]].x[i] > xmax) xmax = tempmachmols[molindex[0]].x[i];
    //    if (tempmachmols[molindex[0]].y[i] > ymax) ymax = tempmachmols[molindex[0]].y[i];
    //    if (tempmachmols[molindex[0]].z[i] > zmax) zmax = tempmachmols[molindex[0]].z[i];
    //    if (tempmachmols[molindex[0]].x[i] < xmin) xmin = tempmachmols[molindex[0]].x[i];
    //    if (tempmachmols[molindex[0]].y[i] < ymin) ymin = tempmachmols[molindex[0]].y[i];
    //    if (tempmachmols[molindex[0]].z[i] < zmin) zmin = tempmachmols[molindex[0]].z[i];
    //}
    //vector<float> c(3, 0.0),f(3,0.0),c1(3, 0.0), f1(3, 0.0),F(3,0.0),F1(3,0.0);
    //c[0] = xmin, c[1] = ymin, c[2] = zmin, c1[0] = xmax, c1[1] = ymax, c1[2] = zmax;
    //MatrixTimesTransVector(theDensityMap.c2f, c, f);
    //MatrixTimesTransVector(theDensityMap.c2f, c1, f1);
    //// the location of atom in grid ?
    //F[0] = (int)(double(f[0] * theDensityMap.grid[0] - theDensityMap.origin[0] + 1));
    //F[1] = (int)(double(f[1] * theDensityMap.grid[1] - theDensityMap.origin[1] + 1));
    //F[2] = (int)(double(f[2] * theDensityMap.grid[2] - theDensityMap.origin[2] + 1));
    //F1[0] = (int)(double(f1[0] * theDensityMap.grid[0] - theDensityMap.origin[0] + 1));
    //F1[1] = (int)(double(f1[1] * theDensityMap.grid[1] - theDensityMap.origin[1] + 1));
    //F1[2] = (int)(double(f1[2] * theDensityMap.grid[2] - theDensityMap.origin[2] + 1));
    //cout << "jubujubujubujubu" << F[0] << "  " << F1[0] << "  " << F[1] << "  " << F1[1] << "  " << F[2] << "  " << F1[2] << endl;

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
    //vector<DOCKMol> best_model;
    DOCKMol best_model;
    copy_molecule(best_model, mol);
    float best_E = 10.0;

    //vector<vector<float>> atm_idx(pnum, vector<float>(3, 0.0));
    //vector<DOCKMol> pointsBx;
    //pointsBx.push_back(mol);
    DOCKMol pointsBx;
    copy_molecule(pointsBx, mol);
    //pointsBx.push_back(mol);

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
        for (int z = theDensityMap.z_min; z <= theDensityMap.z_max; z++)
        //for (int z = 1; z <= Denz; z++)
        {
            atm_jx[2] = z;
            del_ijx[2] = (atm_idxt[2] - atm_jx[2]) / grid[2];
            del_ijx[0] = del_ijx[1] = 0.0;
            vector<float> frac_tmpz;
            MatrixTimesTransVector(theDensityMap.f2c, del_ijx, frac_tmpz);
            if (square_len(frac_tmpz) > (padding + ca_m) * (padding + ca_m)) continue;
            for (int y = theDensityMap.y_min; y <= theDensityMap.y_max; y++)
            //for (int y = 1; y <= Deny; y++)
            {

                atm_jx[1] = y;
                del_ijx[1] = (atm_idxt[1] - atm_jx[1]) / grid[1];
                // wrap-around??                      
                del_ijx[0] = 0.0;
                vector<float> frac_tmpy;
                MatrixTimesTransVector(theDensityMap.f2c, del_ijx, frac_tmpy);
                if (square_len(frac_tmpy) > (padding + ca_m) * (padding + ca_m)) continue;
                for (int x = theDensityMap.x_min; x <= theDensityMap.x_max; x++)
                //for (int x = 1; x <= Denx; x++)
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

            }

        }
        tmp_i = tmp_i + 1;
    }
    vector<float> coor_pdb(3, 0.0);
    float sumC_i2 = 0.0, sumO_i2 = 0.0, sumCO_i2 = 0.0, vol_i2 = 0.0, CC_i2 = 0.0;
    float sumO2_i2 = 0.0, sumC2_i2 = 0.0, varC_i2 = 0.0, varO_i2 = 0.0;
    float clc_x2 = 0.0, obs_x2 = 0.0, eps_x2 = 0.0;
    for (int z = theDensityMap.z_min; z <= theDensityMap.z_max; z++)
    //for (int z = 1; z <= Denz; z++)
    {
        for (int y = theDensityMap.y_min; y <= theDensityMap.y_max; y++)
        //for (int y = 1; y <= Deny; y++)
        {
            for (int x = theDensityMap.x_min; x <= theDensityMap.x_max; x++)
            //for (int x = 1; x <= Denx; x++)
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
    best_E = 1.0 - CC_i2;

    float dis_p = 0.0;
    for (int i = 0; i < receptor_mapsit.size(); i++)
    {
        for (int j = 0; j < pointsBx.num_atoms; j++)
        {
            string VR1 = receptor_mapsit[i].atom_name;//fin_maty[jj].elt_;
            string VR2 = pointsBx.atom_names[j];// fin_maty[t].elt_;
            float dis = 0.0;
            dis = dis + pow((receptor_mapsit[i].x - pointsBx.x[j]),2);
            dis = dis + pow((receptor_mapsit[i].y - pointsBx.y[j]),2);
            dis = dis + pow((receptor_mapsit[i].z - pointsBx.z[j]),2);
            dis_p = sqrt(dis);
            //dis_p = Distance_point(fin_maty[jj].x_, fin_maty[t].x_);
            new_clash = new_clash + GetVdwEgCG(VR1, VR2, dis_p);
        }
    }
    cout << "clash" << new_clash << endl;
    cout << "best_E: " << best_E << endl;
    swapmol(best_model, pointsBx);
    //best_model.push_back(pointsBx);
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
    tran0[0] = coor_pdb[0] - cartbinding[0];
    tran0[1] = coor_pdb[1] - cartbinding[1];
    tran0[2] = coor_pdb[2] - cartbinding[2];
    //MatrixTimesTransVector(theDensityMap.f2c, del_tran, tran0);

    cout << "tran0: " << tran0[0] << " " << tran0[1] << " " << tran0[2] << endl;
    //vector<DOCKMol> pointsBxj;
    //pointsBxj.push_back(pointsBx);
    DOCKMol pointsBxj;
    copy_molecule(pointsBxj, mol);
    int p_atm_size = pointsBx.num_atoms;
    // -tran0 平移量把PDB结构中心变成密度图中心
    for (int j = 0; j < p_atm_size; j++) // not last N
    {
        pointsBxj.x[j] = -tran0[0] + pointsBxj.x[j];
        pointsBxj.y[j] = -tran0[1] + pointsBxj.y[j];
        pointsBxj.z[j] = -tran0[2] + pointsBxj.z[j];
    }

    for (int j = 0; j < pointsBx.num_atoms; j++)
    {
        coor_pdb[0] = coor_pdb[0] + pointsBxj.x[j];
        coor_pdb[1] = coor_pdb[1] + pointsBxj.y[j];
        coor_pdb[2] = coor_pdb[2] + pointsBxj.z[j];
    }
    coor_pdb[0] = coor_pdb[0] / (pointsBxj.num_atoms);
    coor_pdb[1] = coor_pdb[1] / (pointsBxj.num_atoms);
    coor_pdb[2] = coor_pdb[2] / (pointsBxj.num_atoms);
    cout << "coor: " << coor_pdb[0] << " " << coor_pdb[1] << " " << coor_pdb[2] << endl;

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

        for (int z = theDensityMap.z_min; z <= theDensityMap.z_max; z++)
        //for (int z = 1; z <= Denz; z++)
        {
            atm_jx[2] = z;
            del_ijx[2] = (atm_idxt[2] - atm_jx[2]) / grid[2];
            del_ijx[0] = del_ijx[1] = 0.0;
            vector<float> frac_tmpz;
            MatrixTimesTransVector(theDensityMap.f2c, del_ijx, frac_tmpz);
            if (square_len(frac_tmpz) > (padding + ca_m) * (padding + ca_m)) continue;
            for (int y = theDensityMap.y_min; y <= theDensityMap.y_max; y++)
            //for (int y = 1; y <= Deny; y++)
            {
                atm_jx[1] = y;
                del_ijx[1] = (atm_idxt[1] - atm_jx[1]) / grid[1];
                // wrap-around??                    
                del_ijx[0] = 0.0;
                vector<float> frac_tmpy;
                MatrixTimesTransVector(theDensityMap.f2c, del_ijx, frac_tmpy);
                if (square_len(frac_tmpy) > (padding + ca_m) * (padding + ca_m)) continue;
                for (int x = theDensityMap.x_min; x <= theDensityMap.x_max; x++)
                //for (int x = 1; x <= Denx; x++)
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

            }
        }
        tmp_i = tmp_i + 1;
    }
    sumC_i2 = 0.0, sumO_i2 = 0.0, sumCO_i2 = 0.0, vol_i2 = 0.0, CC_i2 = 0.0;
    sumO2_i2 = 0.0, sumC2_i2 = 0.0, varC_i2 = 0.0, varO_i2 = 0.0;
    clc_x2 = 0.0, obs_x2 = 0.0, eps_x2 = 0.0;
    for (int z = theDensityMap.z_min; z <= theDensityMap.z_max; z++)
    //for (int z = 1; z <= Denz; z++)
    {
        for (int y = theDensityMap.y_min; y <= theDensityMap.y_max; y++)
        //for (int y = 1; y <= Deny; y++)
        {
            for (int x = theDensityMap.x_min; x <= theDensityMap.x_max; x++)
            //for (int x = 1; x <= Denx; x++)
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

    new_clash = 0.0;
    dis_p = 0.0;
    for (int i = 0; i < receptor_mapsit.size(); i++)
    {
        for (int j = 0; j < pointsBxj.num_atoms; j++)
        {
            string VR1 = receptor_mapsit[i].atom_name;//fin_maty[jj].elt_;
            string VR2 = pointsBxj.atom_names[j];// fin_maty[t].elt_;
            float dis = 0.0;
            dis = dis + pow((receptor_mapsit[i].x - pointsBxj.x[j]) ,2);
            dis = dis + pow((receptor_mapsit[i].y - pointsBxj.y[j]) ,2);
            dis = dis + pow((receptor_mapsit[i].z - pointsBxj.z[j]) ,2);
            dis_p = sqrt(dis);
            //dis_p = Distance_point(fin_maty[jj].x_, fin_maty[t].x_);
            new_clash = new_clash + GetVdwEgCG(VR1, VR2, dis_p);
        }
    }
    cout << "clash" << new_clash << endl;

    float E_2 = 1.0 - CC_i2;
    cout << "E_2: " << E_2 << endl;
    if (E_2 <= best_E)
    {
        //pointsBx.clear();
        //pointsBx.push_back(pointsBxj[0]);
        swapmol(pointsBx, pointsBxj);
        //pointsBx = pointsBxj;
    }
    swapmol(pointsBx, pointsBxj);
    cout << "CCCCC" << endl;
    cout << pointsBxj.x[0] << "  " << pointsBxj.y[0] << "  " << pointsBxj.z[0] << endl;
    cout << pointsBx.x[0] << "  " << pointsBx.y[0] << "  " << pointsBx.z[0] << endl;
    string molfile = "test.mol2";
    //ofstream clusterpdb(pdbfile.c_str());
    ofstream clustermol(molfile.c_str());
    string str = string(argv[1]);
    clustermol << "########## Name:" << endl;
    Write_Mol2(pointsBxj, clustermol, str);

    //cout << "density,x y z: " << pDenx << " " << pDeny << " " << pDenz << endl;
    float new_CC1 = 0.0;
    float old_CC1 = 0.0;
    float new_CC2 = 0.0;

    cout << "clCC_1clCC_1clCC_1" << theDensityMap.clCC_1(pointsBx) << endl;
    cout << "dandu" << theDensityMap.clCC_4(receptor_mapsit) << endl;
    cout << "dandu1" << theDensityMap.clCC_4(receptor_mapsit_1) << endl;
    cout << "bindingsite" << theDensityMap.clCC_4(receptor_bindingsit) << endl;
    //cout << "matchcentroidposex" << theDensityMap.matchcentroidposex(pointsBx) << endl;
    cout << "clCC_2clCC_2clCC_2" << theDensityMap.clCC_2(pointsBx, receptor_bindingsit) << endl;
    cout << "clCC_3clCC_3clCC_3" << theDensityMap.clCC_2(pointsBx, receptor_mapsit) << endl;
    cout << "clCC_4clCC_4clCC_4" << theDensityMap.clCC_2(pointsBx, receptor_mapsit_1) << endl;
    CC_i2 = theDensityMap.clCC_1(pointsBx);
    cout << "CC_i2: " << CC_i2 << endl;
    new_CC2 = 1.0 - CC_i2;
    cout << "new_CC2: " << new_CC2 << endl;
    best_E = new_CC2;
    //vector<poseCoord>().swap(best_model);
    //best_model.clear();
    //best_model.push_back(pointsBx);
    swapmol(best_model, pointsBx);
    //best_model = pointsBx[0];
    cout << "rotate 0 and 180" << endl;

    //generatationmol(pointsBx, pointsBx, 1, gen);
    //generatationmol(pointsBx, pointsBx, 2, gen);
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
    //CC 0.05 0.001
    float supp_KT0s = 10;//5  //;0.0001;1.5
    float supp_KT0e = 0.1;//;0.00001;0.001
    //1 0.01
    int supp_num1 = 20;
    int supp_num3 = 200;
    int supp_num4 = supp_num1;
    int supp_num2 = 10;//12 20
    // 300 50
    vector<DOCKMol> decstr_tmp;
    vector<vector<float>> E_Tx(supp_num4);

    vector<float> Eng_T(supp_num1, 0.0);
    int nnx = 0;
    int nny = 0;

    vector<float> best_Ex(supp_num1, 100000.0);
    vector<float> best_Ex_4(4, 100000.0);
    vector<DOCKMol> best_modelx;//best_modelx(supp_num1);
    vector<DOCKMol> best_modelx_4(4);
    //vector<DOCKMol> best_model_u;
    DOCKMol best_model_u;
    copy_molecule(best_model_u, mol);

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
            //vector<DOCKMol> point_mat;
            DOCKMol point_mat;
            copy_molecule(point_mat, mol);
            if (kkk > 0)
            {
                //point_mat.push_back(decstrp[jjj]);
                swapmol(point_mat, decstrp[jjj]);
                //point_mat = decstrp[jjj];
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
                //point_mat.push_back(pointsBx[0]);
                swapmol(point_mat, pointsBx);
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

            DOCKMol fin_matx;
            copy_molecule(fin_matx, mol);
            swapmol(fin_matx, point_mat);
            //vector<DOCKMol> fin_matx;
            //fin_matx.push_back(point_mat[0]);
            //fin_matx = point_mat[0];
            /*if (kkk == 0) generatationmol(fin_matx, point_mat, 2, gen);
            else generatationmol(fin_matx, point_mat, 1, gen);*/
            generatationmol(fin_matx, point_mat, 2, gen);

            float zang = float(randIntCustom(30, 180));
            //            float zang=ang_cus[jjj];
            float zang_rd = zang * M_PI / 180.0;
            //if (jjj % 3 == 0 && kkk == 0) // rotate around Z axis
            //{
            //    float asin_theta = 2.0 * randf0and1() - 1.0;
            //    float acos_theta = sqrt(1.0 - asin_theta * asin_theta);
            //    float apha = 2.0 * PI * randf0and1();
            //    float awx = 0.0;
            //    float awy = 0.0;
            //    float awz = 1.0;
            //    // Translation Vector
            //    float t0 = 0.0;
            //    float angle_rotategg = zang_rd; // rotate angle  
            //    float t1 = (randf0and1() * 2.0 - 1.0) * t0 + 0.0;
            //    float t2 = (randf0and1() * 2.0 - 1.0) * t0 + 0.0;
            //    float t3 = (randf0and1() * 2.0 - 1.0) * t0 + 0.0;
            //    //    cout<<"ang,t0: "<<angle_rotategg<<" ("<<t1<<","<<t2<<","<<t3<<")"<<endl;          
            //    float asin = sin(angle_rotategg);
            //    float acos = cos(angle_rotategg);
            //    float u[3][3];
            //    u[0][0] = acos + awx * awx * (1.0 - acos);
            //    u[0][1] = awx * awy * (1.0 - acos) - awz * asin;
            //    u[0][2] = awx * awz * (1.0 - acos) + awy * asin;
            //    u[1][0] = awx * awy * (1.0 - acos) + awz * asin;
            //    u[1][1] = acos + awy * awy * (1.0 - acos);
            //    u[1][2] = awy * awz * (1.0 - acos) - awx * asin;
            //    u[2][0] = awx * awz * (1.0 - acos) - awy * asin;
            //    u[2][1] = awy * awz * (1.0 - acos) + awx * asin;
            //    u[2][2] = acos + awz * awz * (1.0 - acos);
            //    // Rotation points
            //    axyz[0] = coor_pdb[0];
            //    axyz[1] = coor_pdb[1];
            //    axyz[2] = coor_pdb[2];
            //    vector<DOCKMol> fin_matx;
            //    fin_matx.push_back(point_mat[0]);
            //    //fin_matx = point_mat[0];
            //    tmp3 = 0;
            //    int pp_q = point_mat[0].num_atoms;
            //    //for (int j = 0; j < pp_q; j++) // not last N
            //    //{
            //    //    point_mat.x[j] = t1 + axyz[0] + (fin_matx.x[j] - axyz[0]) * u[0][0] + (fin_matx.y[j] - axyz[1]) * u[0][1] + (fin_matx.z[j] - axyz[2]) * u[0][2];
            //    //    point_mat.y[j] = t2 + axyz[1] + (fin_matx.x[j] - axyz[0]) * u[1][0] + (fin_matx.y[j] - axyz[1]) * u[1][1] + (fin_matx.z[j] - axyz[2]) * u[1][2];
            //    //    point_mat.z[j] = t3 + axyz[2] + (fin_matx.x[j] - axyz[0]) * u[2][0] + (fin_matx.y[j] - axyz[1]) * u[2][1] + (fin_matx.z[j] - axyz[2]) * u[2][2];
            //    //    tmp3 = tmp3 + 1;
            //    //}
            //    generatationmol(fin_matx[0], point_mat[0], 1, gen);
            //    //vector<poseCoord>().swap(fin_matx);
            //}
            //if (jjj % 3 == 1 && kkk == 0)  //rotate around X axis
            //{
            //    float asin_theta = 2.0 * randf0and1() - 1.0;
            //    float acos_theta = sqrt(1.0 - asin_theta * asin_theta);
            //    float apha = 2.0 * PI * randf0and1();
            //    float awx = 1.0;
            //    float awy = 0.0;
            //    float awz = 0.0;
            //    // Translation Vector
            //    float t0 = 0.0;
            //    float angle_rotategg = zang_rd; // rotate angle  
            //    float t1 = (randf0and1() * 2.0 - 1.0) * t0 + 0.0;
            //    float t2 = (randf0and1() * 2.0 - 1.0) * t0 + 0.0;
            //    float t3 = (randf0and1() * 2.0 - 1.0) * t0 + 0.0;
            //    //    cout<<"ang,t0: "<<angle_rotategg<<" ("<<t1<<","<<t2<<","<<t3<<")"<<endl;          
            //    float asin = sin(angle_rotategg);
            //    float acos = cos(angle_rotategg);
            //    float u[3][3];
            //    u[0][0] = acos + awx * awx * (1.0 - acos);
            //    u[0][1] = awx * awy * (1.0 - acos) - awz * asin;
            //    u[0][2] = awx * awz * (1.0 - acos) + awy * asin;
            //    u[1][0] = awx * awy * (1.0 - acos) + awz * asin;
            //    u[1][1] = acos + awy * awy * (1.0 - acos);
            //    u[1][2] = awy * awz * (1.0 - acos) - awx * asin;
            //    u[2][0] = awx * awz * (1.0 - acos) - awy * asin;
            //    u[2][1] = awy * awz * (1.0 - acos) + awx * asin;
            //    u[2][2] = acos + awz * awz * (1.0 - acos);
            //    // Rotation points
            //    axyz[0] = coor_pdb[0];
            //    axyz[1] = coor_pdb[1];
            //    axyz[2] = coor_pdb[2];
            //    DOCKMol fin_matx;
            //    fin_matx = point_mat;
            //    tmp3 = 0;
            //    int pp_q = point_mat.num_atoms;
            //    //for (int j = 0; j < pp_q; j++) // not last N
            //    //{
            //    //    point_mat.x[j] = t1 + axyz[0] + (fin_matx.x[j] - axyz[0]) * u[0][0] + (fin_matx.y[j] - axyz[1]) * u[0][1] + (fin_matx.z[j] - axyz[2]) * u[0][2];
            //    //    point_mat.y[j] = t2 + axyz[1] + (fin_matx.x[j] - axyz[0]) * u[1][0] + (fin_matx.y[j] - axyz[1]) * u[1][1] + (fin_matx.z[j] - axyz[2]) * u[1][2];
            //    //    point_mat.z[j] = t3 + axyz[2] + (fin_matx.x[j] - axyz[0]) * u[2][0] + (fin_matx.y[j] - axyz[1]) * u[2][1] + (fin_matx.z[j] - axyz[2]) * u[2][2];
            //    //    tmp3 = tmp3 + 1;
            //    //}
            //    generatationmol(fin_matx, point_mat, 1, gen);
            //    //vector<poseCoord>().swap(fin_matx);
            //}
            //if (jjj % 3 == 2 && kkk == 0)  //rotate around Y axis
            //{
            //    float asin_theta = 2.0 * randf0and1() - 1.0;
            //    float acos_theta = sqrt(1.0 - asin_theta * asin_theta);
            //    float apha = 2.0 * PI * randf0and1();
            //    float awx = 0.0;
            //    float awy = 1.0;
            //    float awz = 0.0;
            //    // Translation Vector
            //    float t0 = 0.0;
            //    float angle_rotategg = zang_rd; // rotate angle  
            //    float t1 = (randf0and1() * 2.0 - 1.0) * t0 + 0.0;
            //    float t2 = (randf0and1() * 2.0 - 1.0) * t0 + 0.0;
            //    float t3 = (randf0and1() * 2.0 - 1.0) * t0 + 0.0;
            //    //    cout<<"ang,t0: "<<angle_rotategg<<" ("<<t1<<","<<t2<<","<<t3<<")"<<endl;          
            //    float asin = sin(angle_rotategg);
            //    float acos = cos(angle_rotategg);
            //    float u[3][3];
            //    u[0][0] = acos + awx * awx * (1.0 - acos);
            //    u[0][1] = awx * awy * (1.0 - acos) - awz * asin;
            //    u[0][2] = awx * awz * (1.0 - acos) + awy * asin;
            //    u[1][0] = awx * awy * (1.0 - acos) + awz * asin;
            //    u[1][1] = acos + awy * awy * (1.0 - acos);
            //    u[1][2] = awy * awz * (1.0 - acos) - awx * asin;
            //    u[2][0] = awx * awz * (1.0 - acos) - awy * asin;
            //    u[2][1] = awy * awz * (1.0 - acos) + awx * asin;
            //    u[2][2] = acos + awz * awz * (1.0 - acos);
            //    // Rotation points
            //    axyz[0] = coor_pdb[0];
            //    axyz[1] = coor_pdb[1];
            //    axyz[2] = coor_pdb[2];
            //    DOCKMol fin_matx;
            //    fin_matx = point_mat;
            //    tmp3 = 0;
            //    int pp_q = point_mat.num_atoms;
            //    //for (int j = 0; j < pp_q; j++) // not last N
            //    //{
            //    //    point_mat.x[j] = t1 + axyz[0] + (fin_matx.x[j] - axyz[0]) * u[0][0] + (fin_matx.y[j] - axyz[1]) * u[0][1] + (fin_matx.z[j] - axyz[2]) * u[0][2];
            //    //    point_mat.y[j] = t2 + axyz[1] + (fin_matx.x[j] - axyz[0]) * u[1][0] + (fin_matx.y[j] - axyz[1]) * u[1][1] + (fin_matx.z[j] - axyz[2]) * u[1][2];
            //    //    point_mat.z[j] = t3 + axyz[2] + (fin_matx.x[j] - axyz[0]) * u[2][0] + (fin_matx.y[j] - axyz[1]) * u[2][1] + (fin_matx.z[j] - axyz[2]) * u[2][2];
            //    //    tmp3 = tmp3 + 1;
            //    //}
            //    generatationmol(fin_matx, point_mat, 1, gen);
            //    //vector<poseCoord>().swap(fin_matx);
            //}
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
                //best_modelx[jjj] = point_mat[0];
                DOCKMol best;
                copy_molecule(best, mol);
                swapmol(best, point_mat);
                best_modelx.push_back(best);
            }

            float E_500 = 0.0;
            int E_500_int = 0;
            KT = supp_REMC[jjj];

            new_clash = 0.0;
            dis_p = 0.0;
            for (int i = 0; i < receptor_mapsit.size(); i++)
            {
                for (int j = 0; j < point_mat.num_atoms; j++)
                {
                    string VR1 = receptor_mapsit[i].atom_name;//fin_maty[jj].elt_;
                    string VR2 = point_mat.atom_names[j];// fin_maty[t].elt_;
                    float dis = 0.0;
                    dis = dis + pow((receptor_mapsit[i].x - point_mat.x[j]),2);
                    dis = dis + pow((receptor_mapsit[i].y - point_mat.y[j]) ,2);
                    dis = dis + pow((receptor_mapsit[i].z - point_mat.z[j]) ,2);
                    dis_p = sqrt(dis);
                    //dis_p = Distance_point(fin_maty[jj].x_, fin_maty[t].x_);
                    new_clash = new_clash + GetVdwEgCG(VR1, VR2, dis_p);
                }
            }

            float old_dist = 0.0;
            vector<float> centers_R1(3, 0);
            for (int i = 0; i < point_mat.num_atoms; i++)
            {
                centers_R1[0] += point_mat.x[i];
                centers_R1[1] += point_mat.y[i];
                centers_R1[2] += point_mat.z[i];
            }
            centers_R1[0] = centers_R1[0] / point_mat.num_atoms;
            centers_R1[1] = centers_R1[1] / point_mat.num_atoms;
            centers_R1[2] = centers_R1[2] / point_mat.num_atoms;

            old_dist = sqrt(pow(cartbinding[0] - centers_R1[0], 2) + pow(cartbinding[1] - centers_R1[1], 2) +
                pow(cartbinding[2] - centers_R1[2], 2));
            if (old_dist >  max) old_dist = 500* old_dist;

            //CC_i2 = theDensityMap.clCC_2(point_mat, receptor_mapsit);
            CC_i2 = theDensityMap.clCC_1(point_mat);
            new_CC2 = 0.5 - CC_i2;
            //new_CC2 = 1.0 - CC_i2;
            old_CC = new_CC2;
            old_clash = new_clash;
            if (point_mat.num_atoms > 40) old_dE = 2000 * new_CC2 + old_clash + 10 * old_dist;
            else old_dE = 500 * new_CC2 + old_clash + 10 * old_dist;
            //cout<<"new_CC2: "<<new_CC2<<"  " << kkk <<"    " << jjj << endl;
            int acc_rate1 = 0;
            for (int iii = 0; iii < supp_num3; iii++)
            {
            //    float asin_theta = 2.0 * randf0and1() - 1.0;
            //    float acos_theta = sqrt(1.0 - asin_theta * asin_theta);
            //    float apha = 2.0 * PI * randf0and1();
            //    float awx = acos_theta * cos(apha);
            //    float awy = acos_theta * sin(apha);
            //    float awz = asin_theta;
            //    // Translation Vector
            //    float t0 = 1.0;
            //    //    float t0=1.0;
            //    //    cout<<"t0: "<<t1<<" "<<t2<<" "<<t3<<endl;
            //    // Rotation matrix
            //    float anggg = 30.0;
            //    float angle_rotategg = (2.0 * randf0and1() - 1.0) * anggg; // rotate angle
            //    float t1 = (randf0and1() * 2.0 - 1.0) * t0 + 0.0;
            //    float t2 = (randf0and1() * 2.0 - 1.0) * t0 + 0.0;
            //    float t3 = (randf0and1() * 2.0 - 1.0) * t0 + 0.0;
            //    //    cout<<"ang,t0: "<<angle_rotategg<<" ("<<t1<<","<<t2<<","<<t3<<")"<<endl;          
            //    float asin = sin(angle_rotategg);
            //    float acos = cos(angle_rotategg);
            //    float u[3][3];
            //    u[0][0] = acos + awx * awx * (1.0 - acos);
            //    u[0][1] = awx * awy * (1.0 - acos) - awz * asin;
            //    u[0][2] = awx * awz * (1.0 - acos) + awy * asin;
            //    u[1][0] = awx * awy * (1.0 - acos) + awz * asin;
            //    u[1][1] = acos + awy * awy * (1.0 - acos);
            //    u[1][2] = awy * awz * (1.0 - acos) - awx * asin;
            //    u[2][0] = awx * awz * (1.0 - acos) - awy * asin;
            //    u[2][1] = awy * awz * (1.0 - acos) + awx * asin;
            //    u[2][2] = acos + awz * awz * (1.0 - acos);
            //    axyz[0] = coor_pdb[0];
            //    axyz[1] = coor_pdb[1];
            //    axyz[2] = coor_pdb[2];
            //    //    beg_p=rand_point; 
            ////    cout<<"XXX2"<<endl;  
            //    int tmp3 = 0;
            
                DOCKMol fin_mat;
                copy_molecule(fin_mat, mol);
                swapmol(fin_mat, point_mat);
                //vector<DOCKMol> fin_mat;
                //fin_mat.push_back(point_mat[0]);
                //fin_mat = point_mat;
                //DOCKMol fin_matx;
                //fin_matx = fin_mat;
                DOCKMol fin_matx;
                copy_molecule(fin_matx, mol);
                swapmol(fin_matx, fin_mat);
                tmp3 = 0;
                int pp_k = fin_mat.num_atoms;
                //for (int j = 0; j < pp_k; j++) // not last N
                //{
                //    fin_mat.x[j] = t1 + axyz[0] + (fin_matx.x[j] - axyz[0]) * u[0][0] + (fin_matx.y[j] - axyz[1]) * u[0][1] + (fin_matx.z[j] - axyz[2]) * u[0][2];
                //    fin_mat.y[j] = t2 + axyz[1] + (fin_matx.x[j] - axyz[0]) * u[1][0] + (fin_matx.y[j] - axyz[1]) * u[1][1] + (fin_matx.z[j] - axyz[2]) * u[1][2];
                //    fin_mat.z[j] = t3 + axyz[2] + (fin_matx.x[j] - axyz[0]) * u[2][0] + (fin_matx.y[j] - axyz[1]) * u[2][1] + (fin_matx.z[j] - axyz[2]) * u[2][2];
                //    tmp3 = tmp3 + 1;
                //}
                generatationmol(fin_matx, fin_mat, 2, gen);
                //CC_i2 = theDensityMap.clCC_2(fin_mat, receptor_mapsit);
                CC_i2 = theDensityMap.clCC_1(fin_mat);
                new_clash = 0.0;
                dis_p = 0.0;
                for (int i = 0; i < receptor_mapsit.size(); i++)
                {
                    for (int j = 0; j < fin_mat.num_atoms; j++)
                    {
                        string VR1 = receptor_mapsit[i].atom_name;//fin_maty[jj].elt_;
                        string VR2 = fin_mat.atom_names[j];// fin_maty[t].elt_;
                        float dis = 0.0;
                        dis = dis + pow((receptor_mapsit[i].x - fin_mat.x[j]) ,2);
                        dis = dis + pow((receptor_mapsit[i].y - fin_mat.y[j]) ,2);
                        dis = dis + pow((receptor_mapsit[i].z - fin_mat.z[j]) ,2);
                        dis_p = sqrt(dis);
                        //dis_p = Distance_point(fin_maty[jj].x_, fin_maty[t].x_);
                        new_clash = new_clash + GetVdwEgCG(VR1, VR2, dis_p);
                    }
                }

/* for (int i = 0; i < fin_mat.num_atoms; i++)
                 {
                     cout << fin_mat.x[i] << "  " << fin_mat.y[i] << "  " << fin_mat.z[i] << endl;
                 }*/
                float new_dist = 0.0;
                vector<float> centers_R(3, 0);
                for (int i = 0; i < fin_mat.num_atoms; i++)
                {
                    centers_R[0] += fin_mat.x[i];
                    centers_R[1] += fin_mat.y[i];
                    centers_R[2] += fin_mat.z[i];
                }
                centers_R[0] = centers_R[0] / fin_mat.num_atoms;
                centers_R[1] = centers_R[1] / fin_mat.num_atoms;
                centers_R[2] = centers_R[2] / fin_mat.num_atoms;

                new_dist = sqrt(pow(cartbinding[0] - centers_R[0], 2) + pow(cartbinding[1] - centers_R[1], 2) +
                    pow(cartbinding[2] - centers_R[2], 2));
                if (new_dist >  max) new_dist = 500* new_dist;

                new_CC = 0.5 - CC_i2;
                //new_CC = 1.0 - CC_i2;
                if(fin_mat.num_atoms>40) new_dE = 2000 * new_CC + new_clash + 10 * new_dist;
                else new_dE = 500 * new_CC + new_clash + 10 * new_dist;
                //new_dE = new_CC1 ;                                   
                //cout<<"old_CC1,new_CC1:RRRRRRRRRRRRR "<< new_dE <<" "<< old_dE <<endl;

                //for (int j = 0; j < point_mat.num_atoms; j++)
                //{
                //    cout << "CCCC" << point_mat.x[0] << "  " << point_mat.y[0] << "  " << point_mat.z[0] << endl;
                //}

                dE = new_dE - old_dE;
                if (new_dE < old_dE)// && new_CC<old_CC)
                {
                    cout << "old_CC1,new_CC1:RRRRRRRRRRRRR " << new_dE << " " <<new_CC<<"  " << old_dE <<"  "<<old_CC<< endl;
                    //tmp3 = 0;
                    coor_pdb = vector<float>(3, 0.0);

                    //point_mat.clear();
                    //point_mat.push_back(fin_mat[0]);
                    swapmol(point_mat, fin_mat);
                    //point_mat = fin_mat;
                    for (int j = 0; j < fin_mat.num_atoms; j++)
                    {
                        //point_mat[3 * j + 0] = fin_mat[3 * tmp3 + 0];
                        //point_mat[3 * j + 1] = fin_mat[3 * tmp3 + 1];
                        //point_mat[3 * j + 2] = fin_mat[3 * tmp3 + 2];
                        //cout << "CCCC" << point_mat.x[j] << "  " << point_mat.y[j] << "  " << point_mat.z[j] << endl;
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
                        //best_modelx[jjj] = fin_mat;
                        swapmol(best_modelx[jjj], fin_mat);
                    }
                    if (new_dE < best_E_u)
                    {
                        //for (int i = 0; i < fin_mat[0].num_atoms; i++)
                        //{
                        //    best_model_u.x[i] = fin_mat[0].x[i];
                        //    best_model_u.y[i] = fin_mat[0].y[i];
                        //    best_model_u.z[i] = fin_mat[0].z[i];
                        //}
                        //best_model_u = fin_mat;
                        swapmol(best_model_u, fin_mat);
                        //best_model_u.clear();
                        //best_model_u.push_back(fin_mat[0]);
                        best_E_u = new_dE;
                        cout << "best" << endl;
                        cout << "AAAAAA" << best_model_u.x[0] << "  " << best_model_u.y[0] << "  " << best_model_u.z[0] << endl;
                    }
                    if (new_dE < out_sup_all_CC)
                    {
                        out_sup_all_CC = new_dE;
                    }
                    //    cout<<"coor: "<<coor_pdb[0]<<" "<<coor_pdb[1]<<" "<<coor_pdb[2]<<endl;                  
                //        new_CC2 = new_CC1;
                //        cout<<"old_CC1,new_CC1: "<<old_dE<<" "<<new_dE<<endl;
                    old_dE = new_dE;
                    old_CC = new_CC;
                    E_500 = new_dE;
                    E_500_int = 1;

                }
                else if (new_dE == old_dE) continue;
                else
                {
                    //    float tmpx=rand()/double(RAND_MAX);
                    float tmpx = randf0and1();
                    float mc_v = exp(-dE / (KT));
                    if (tmpx < mc_v) // CC must be >0
                    {
                        cout << "RRRRRRRRRRRRR old_CC1,new_CC1: " << new_dE << " " << new_CC << "  " << old_dE << "  " << old_CC << endl;
                        //cout << "1301" << endl;
                        acc_rate1 = acc_rate1 + 1;
                        //cout << acc_rate1 << endl;
                        tmp3 = 0;
                        coor_pdb = vector<float>(3, 0.0);
                        swapmol(point_mat, fin_mat);
                        //point_mat.clear();
                        //point_mat.push_back(fin_mat[0]);
                        //point_mat = fin_mat;
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
                        old_CC = new_CC;
                        E_500 = new_dE;
                        E_500_int = 1;
                        //    new_CC2 = new_CC1;

                    }
                }
                nnx = nnx + 1;
                //cout << "BBBBB" << best_model_u.x[0] << "  " << best_model_u.y[0] << "  " << best_model_u.z[0] << endl;
            }
            cout << "acc_rate:  " << acc_rate1 << endl;
            cout << best_E_u << endl;
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
                    cout << "aaaaaaa" << endl;
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
                    cout << "aaaaaaa" << endl;
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

    vector<int> ccindex;
    for (int i = 0; i < best_Ex.size(); i++)
    {
        ccindex.push_back(i);
    }

    sortmol_em(best_modelx, ccindex);
    //sortmol_CC(best_Ex, ccindex);

    vector<float> centers_Fin(3, 0);
    for (int i = 0; i < best_modelx[ccindex[0]].num_atoms; i++)
    {
        centers_Fin[0] += best_modelx[ccindex[0]].x[i];
        centers_Fin[1] += best_modelx[ccindex[0]].y[i];
        centers_Fin[2] += best_modelx[ccindex[0]].z[i];
    }
    centers_Fin[0] = centers_Fin[0] / best_modelx[ccindex[0]].num_atoms;
    centers_Fin[1] = centers_Fin[1] / best_modelx[ccindex[0]].num_atoms;
    centers_Fin[2] = centers_Fin[2] / best_modelx[ccindex[0]].num_atoms;
    cout << "FIN center  " << centers_Fin[0] << "  " << centers_Fin[1] << "  " << centers_Fin[2] << endl;
    cout << "Center distance  " << abs(centers_Fin[0] - cartbinding[0]) << "  " << abs(centers_Fin[1] - cartbinding[1]) << "  " << abs(centers_Fin[2] - cartbinding[2]) << endl;
    cout << "Center distance Native " << abs(centers_Fin[0] - centers7[0]) << "  " << abs(centers_Fin[1] - centers7[1]) << "  " << abs(centers_Fin[2] - centers7[2]) << endl;

    float max_idx = 0;
    if (mol.num_atoms > 40) max_idx = 0.55 * max;
    else max_idx = 0.5 * max;
    theDensityMap.cart(centers_Fin, fartbinding, max_idx);

    if (mol.num_atoms > 40) max_idx = 0.7 * max;
    else max_idx = 0.6 * max;
    receptor_mapsit.clear();
    for (int j = 0; j < receptor.num_atoms; j++)
    {
        if (abs(receptor.x[j] - centers_Fin[0]) < max_idx && abs(receptor.y[j] - centers_Fin[1]) < max_idx &&
            abs(receptor.z[j] - centers_Fin[2]) < max_idx)
        {
            poseCoord atom;
            atom.atom_name = receptor.atom_names[j];
            atom.residue_numbers = receptor.atom_residue_numbers[j];
            atom.x = receptor.x[j];
            atom.y = receptor.y[j];
            atom.z = receptor.z[j];
            receptor_mapsit.push_back(atom);
        }
    }
    cout << receptor_mapsit.size() << endl;
    old_res, new_res = "aa";
    for (int i = 0; i < receptor_mapsit.size(); i++)
    {
        old_res = receptor_mapsit[i].residue_numbers;
        if (old_res != new_res)
        {
            cout << ":" << old_res;
            new_res = old_res;
        }
    }
    cout << endl;

    //decstrp.clear();
    //decstr_tmp.clear();
    //best_modelx.clear();
    //cout << best_Ex.size() << endl;
    //for (int i = 0; i < best_Ex.size(); i++)
    //{
    //    best_Ex[i] = 1000000;
    //}
    //for (int i = 0; i < supp_num1; i++)
    //{
    //    float supp_KT0x = pow(float(supp_KT0e / supp_KT0s), float(float(i) / float(supp_num1)));
    //    supp_REMC[i] = supp_KT0s * float(supp_KT0x);
    //}
    //for (int kkk = 0; kkk < supp_num2; kkk++) // 2
    //{
    //    vector<DOCKMol > decstrz;
    //    vector<DOCKMol >().swap(decstrz);
    //    vector<float> E_REMC(supp_num1, 0.0);
    //    for (int jjj = 0; jjj < supp_num1; jjj++)
    //    {
    //        //vector<DOCKMol> point_mat;
    //        DOCKMol point_mat;
    //        copy_molecule(point_mat, mol);
    //        if (kkk > 0)
    //        {
    //            //point_mat.push_back(decstrp[jjj]);
    //            swapmol(point_mat, decstrp[jjj]);
    //            //point_mat = decstrp[jjj];
    //            coor_pdb = vector<float>(3, 0.0);
    //            for (int j = 0; j < point_mat.num_atoms; j++)
    //            {
    //                coor_pdb[0] = coor_pdb[0] + point_mat.x[j];
    //                coor_pdb[1] = coor_pdb[1] + point_mat.y[j];
    //                coor_pdb[2] = coor_pdb[2] + point_mat.z[j];
    //            }
    //            coor_pdb[0] = coor_pdb[0] / (point_mat.num_atoms);
    //            coor_pdb[1] = coor_pdb[1] / (point_mat.num_atoms);
    //            coor_pdb[2] = coor_pdb[2] / (point_mat.num_atoms);
    //        }
    //        else
    //        {
    //            coor_pdb = vector<float>(3, 0.0);
    //            //point_mat.push_back(pointsBx[0]);
    //            swapmol(point_mat, pointsBx);
    //            for (int j = 0; j < point_mat.num_atoms; j++)
    //            {
    //                coor_pdb[0] = coor_pdb[0] + point_mat.x[j];
    //                coor_pdb[1] = coor_pdb[1] + point_mat.y[j];
    //                coor_pdb[2] = coor_pdb[2] + point_mat.z[j];
    //            }
    //            coor_pdb[0] = coor_pdb[0] / (point_mat.num_atoms);
    //            coor_pdb[1] = coor_pdb[1] / (point_mat.num_atoms);
    //            coor_pdb[2] = coor_pdb[2] / (point_mat.num_atoms);
    //        }
    //        DOCKMol fin_matx;
    //        copy_molecule(fin_matx, mol);
    //        swapmol(fin_matx, point_mat);
    //        //vector<DOCKMol> fin_matx;
    //        //fin_matx.push_back(point_mat[0]);
    //        //fin_matx = point_mat[0];
    //        /*if (kkk == 0) generatationmol(fin_matx, point_mat, 2, gen);
    //        else generatationmol(fin_matx, point_mat, 1, gen);*/
    //        generatationmol(fin_matx, point_mat, 2, gen);
    //        float zang = float(randIntCustom(30, 180));
    //        //            float zang=ang_cus[jjj];
    //        float zang_rd = zang * M_PI / 180.0;
    //        
    //        coor_pdb = vector<float>(3, 0.0);
    //        for (int j = 0; j < point_mat.num_atoms; j++)
    //        {
    //            coor_pdb[0] = coor_pdb[0] + point_mat.x[j];
    //            coor_pdb[1] = coor_pdb[1] + point_mat.y[j];
    //            coor_pdb[2] = coor_pdb[2] + point_mat.z[j];
    //        }
    //        coor_pdb[0] = coor_pdb[0] / (point_mat.num_atoms);
    //        coor_pdb[1] = coor_pdb[1] / (point_mat.num_atoms);
    //        coor_pdb[2] = coor_pdb[2] / (point_mat.num_atoms);
    //        if (kkk == 0)
    //        {
    //            //best_modelx[jjj] = point_mat[0];
    //            DOCKMol best;
    //            copy_molecule(best, mol);
    //            swapmol(best, point_mat);
    //            best_modelx.push_back(best);
    //        }
    //        float E_500 = 0.0;
    //        int E_500_int = 0;
    //        KT = supp_REMC[jjj];
    //        new_clash = 0.0;
    //        dis_p = 0.0;
    //        for (int i = 0; i < receptor_mapsit.size(); i++)
    //        {
    //            for (int j = 0; j < point_mat.num_atoms; j++)
    //            {
    //                string VR1 = receptor_mapsit[i].atom_name;//fin_maty[jj].elt_;
    //                string VR2 = point_mat.atom_names[j];// fin_maty[t].elt_;
    //                float dis = 0.0;
    //                dis = dis + pow((receptor_mapsit[i].x - point_mat.x[j]), 2);
    //                dis = dis + pow((receptor_mapsit[i].y - point_mat.y[j]), 2);
    //                dis = dis + pow((receptor_mapsit[i].z - point_mat.z[j]), 2);
    //                dis_p = sqrt(dis);
    //                //dis_p = Distance_point(fin_maty[jj].x_, fin_maty[t].x_);
    //                new_clash = new_clash + GetVdwEgCG(VR1, VR2, dis_p);
    //            }
    //        }
    //        float old_dist = 0.0;
    //        vector<float> centers_R1(3, 0);
    //        for (int i = 0; i < point_mat.num_atoms; i++)
    //        {
    //            centers_R1[0] += point_mat.x[i];
    //            centers_R1[1] += point_mat.y[i];
    //            centers_R1[2] += point_mat.z[i];
    //        }
    //        centers_R1[0] = centers_R1[0] / point_mat.num_atoms;
    //        centers_R1[1] = centers_R1[1] / point_mat.num_atoms;
    //        centers_R1[2] = centers_R1[2] / point_mat.num_atoms;
    //        old_dist = sqrt(pow(centers_Fin[0] - centers_R1[0], 2) + pow(centers_Fin[1] - centers_R1[1], 2) +
    //            pow(centers_Fin[2] - centers_R1[2], 2));
    //        if (old_dist < 0.5) old_dist = 0;
    //        //CC_i2 = theDensityMap.clCC_2(point_mat, receptor_mapsit);
    //        CC_i2 = theDensityMap.clCC_1(point_mat);
    //        new_CC2 = 0.5 - CC_i2;
    //        //new_CC2 = 1.0 - CC_i2;
    //        old_CC = new_CC2;
    //        old_clash = new_clash;
    //        old_dE = 500 * new_CC2 + old_clash + 500 * old_dist;
    //        //cout<<"new_CC2: "<<new_CC2<<"  " << kkk <<"    " << jjj << endl;
    //        int acc_rate1 = 0;
    //        for (int iii = 0; iii < supp_num3; iii++)
    //        {
    //            
    //            DOCKMol fin_mat;
    //            copy_molecule(fin_mat, mol);
    //            swapmol(fin_mat, point_mat);
    //            //vector<DOCKMol> fin_mat;
    //            //fin_mat.push_back(point_mat[0]);
    //            //fin_mat = point_mat;
    //            //DOCKMol fin_matx;
    //            //fin_matx = fin_mat;
    //            DOCKMol fin_matx;
    //            copy_molecule(fin_matx, mol);
    //            swapmol(fin_matx, fin_mat);
    //            tmp3 = 0;
    //            int pp_k = fin_mat.num_atoms;
    //            //for (int j = 0; j < pp_k; j++) // not last N
    //            //{
    //            //    fin_mat.x[j] = t1 + axyz[0] + (fin_matx.x[j] - axyz[0]) * u[0][0] + (fin_matx.y[j] - axyz[1]) * u[0][1] + (fin_matx.z[j] - axyz[2]) * u[0][2];
    //            //    fin_mat.y[j] = t2 + axyz[1] + (fin_matx.x[j] - axyz[0]) * u[1][0] + (fin_matx.y[j] - axyz[1]) * u[1][1] + (fin_matx.z[j] - axyz[2]) * u[1][2];
    //            //    fin_mat.z[j] = t3 + axyz[2] + (fin_matx.x[j] - axyz[0]) * u[2][0] + (fin_matx.y[j] - axyz[1]) * u[2][1] + (fin_matx.z[j] - axyz[2]) * u[2][2];
    //            //    tmp3 = tmp3 + 1;
    //            //}
    //            generatationmol(fin_matx, fin_mat, 2, gen);
    //            //CC_i2 = theDensityMap.clCC_2(fin_mat, receptor_mapsit);
    //            CC_i2 = theDensityMap.clCC_1(fin_mat);
    //            new_clash = 0.0;
    //            dis_p = 0.0;
    //            for (int i = 0; i < receptor_mapsit.size(); i++)
    //            {
    //                for (int j = 0; j < fin_mat.num_atoms; j++)
    //                {
    //                    string VR1 = receptor_mapsit[i].atom_name;//fin_maty[jj].elt_;
    //                    string VR2 = fin_mat.atom_names[j];// fin_maty[t].elt_;
    //                    float dis = 0.0;
    //                    dis = dis + pow((receptor_mapsit[i].x - fin_mat.x[j]), 2);
    //                    dis = dis + pow((receptor_mapsit[i].y - fin_mat.y[j]), 2);
    //                    dis = dis + pow((receptor_mapsit[i].z - fin_mat.z[j]), 2);
    //                    dis_p = sqrt(dis);
    //                    //dis_p = Distance_point(fin_maty[jj].x_, fin_maty[t].x_);
    //                    new_clash = new_clash + GetVdwEgCG(VR1, VR2, dis_p);
    //                }
    //            }
    //            /* for (int i = 0; i < fin_mat.num_atoms; i++)
    //                             {
    //                                 cout << fin_mat.x[i] << "  " << fin_mat.y[i] << "  " << fin_mat.z[i] << endl;
    //                             }*/
    //            float new_dist = 0.0;
    //            vector<float> centers_R(3, 0);
    //            for (int i = 0; i < fin_mat.num_atoms; i++)
    //            {
    //                centers_R[0] += fin_mat.x[i];
    //                centers_R[1] += fin_mat.y[i];
    //                centers_R[2] += fin_mat.z[i];
    //            }
    //            centers_R[0] = centers_R[0] / fin_mat.num_atoms;
    //            centers_R[1] = centers_R[1] / fin_mat.num_atoms;
    //            centers_R[2] = centers_R[2] / fin_mat.num_atoms;
    //            new_dist = sqrt(pow(centers_Fin[0] - centers_R[0], 2) + pow(centers_Fin[1] - centers_R[1], 2) +
    //                pow(centers_Fin[2] - centers_R[2], 2));
    //            if (new_dist < 0.5) new_dist = 0;
    //            new_CC = 0.5 - CC_i2;
    //            //new_CC = 1.0 - CC_i2;
    //            new_dE = 500 * new_CC + new_clash + 500 * new_dist;
    //            //new_dE = new_CC1 ;                                   
    //            //cout<<"old_CC1,new_CC1:RRRRRRRRRRRRR "<< new_dE <<" "<< old_dE <<endl;
    //            //for (int j = 0; j < point_mat.num_atoms; j++)
    //            //{
    //            //    cout << "CCCC" << point_mat.x[0] << "  " << point_mat.y[0] << "  " << point_mat.z[0] << endl;
    //            //}
    //            dE = new_dE - old_dE;
    //            if (new_dE < old_dE)// && new_CC<old_CC)
    //            {
    //                cout << "old_CC1,new_CC1:RRRRRRRRRRRRR " << new_dE << " " << new_CC << "  " << old_dE << "  " << old_CC << endl;
    //                //tmp3 = 0;
    //                coor_pdb = vector<float>(3, 0.0);
    //                //point_mat.clear();
    //                //point_mat.push_back(fin_mat[0]);
    //                swapmol(point_mat, fin_mat);
    //                //point_mat = fin_mat;
    //                for (int j = 0; j < fin_mat.num_atoms; j++)
    //                {
    //                    //point_mat[3 * j + 0] = fin_mat[3 * tmp3 + 0];
    //                    //point_mat[3 * j + 1] = fin_mat[3 * tmp3 + 1];
    //                    //point_mat[3 * j + 2] = fin_mat[3 * tmp3 + 2];
    //                    //cout << "CCCC" << point_mat.x[j] << "  " << point_mat.y[j] << "  " << point_mat.z[j] << endl;
    //                    coor_pdb[0] = coor_pdb[0] + fin_mat.x[j];
    //                    coor_pdb[1] = coor_pdb[1] + fin_mat.y[j];
    //                    coor_pdb[2] = coor_pdb[2] + fin_mat.z[j];
    //                    ///tmp3 = tmp3 + 1;
    //                }
    //                coor_pdb[0] = coor_pdb[0] / (fin_mat.num_atoms);
    //                coor_pdb[1] = coor_pdb[1] / (fin_mat.num_atoms);
    //                coor_pdb[2] = coor_pdb[2] / (fin_mat.num_atoms);
    //                if (new_dE < best_Ex[jjj])
    //                {
    //                    best_Ex[jjj] = new_dE;
    //                    //vector<poseCoord>().swap(best_modelx[jjj]);
    //                    //best_modelx[jjj] = fin_mat;
    //                    swapmol(best_modelx[jjj], fin_mat);
    //                }
    //                if (new_dE < best_E_u)
    //                {
    //                    //for (int i = 0; i < fin_mat[0].num_atoms; i++)
    //                    //{
    //                    //    best_model_u.x[i] = fin_mat[0].x[i];
    //                    //    best_model_u.y[i] = fin_mat[0].y[i];
    //                    //    best_model_u.z[i] = fin_mat[0].z[i];
    //                    //}
    //                    //best_model_u = fin_mat;
    //                    swapmol(best_model_u, fin_mat);
    //                    //best_model_u.clear();
    //                    //best_model_u.push_back(fin_mat[0]);
    //                    best_E_u = new_dE;
    //                    cout << "best" << endl;
    //                    cout << "AAAAAA" << best_model_u.x[0] << "  " << best_model_u.y[0] << "  " << best_model_u.z[0] << endl;
    //                }
    //                if (new_dE < out_sup_all_CC)
    //                {
    //                    out_sup_all_CC = new_dE;
    //                }
    //                //    cout<<"coor: "<<coor_pdb[0]<<" "<<coor_pdb[1]<<" "<<coor_pdb[2]<<endl;                  
    //            //        new_CC2 = new_CC1;
    //            //        cout<<"old_CC1,new_CC1: "<<old_dE<<" "<<new_dE<<endl;
    //                old_dE = new_dE;
    //                old_CC = new_CC;
    //                E_500 = new_dE;
    //                E_500_int = 1;
    //            }
    //            else if (new_dE == old_dE) continue;
    //            else
    //            {
    //                //    float tmpx=rand()/double(RAND_MAX);
    //                float tmpx = randf0and1();
    //                float mc_v = exp(-dE / (KT));
    //                if (tmpx < mc_v) // CC must be >0
    //                {
    //                    cout << "RRRRRRRRRRRRR old_CC1,new_CC1: " << new_dE << " " << new_CC << "  " << old_dE << "  " << old_CC << endl;
    //                    //cout << "1301" << endl;
    //                    acc_rate1 = acc_rate1 + 1;
    //                    //cout << acc_rate1 << endl;
    //                    tmp3 = 0;
    //                    coor_pdb = vector<float>(3, 0.0);
    //                    swapmol(point_mat, fin_mat);
    //                    //point_mat.clear();
    //                    //point_mat.push_back(fin_mat[0]);
    //                    //point_mat = fin_mat;
    //                    for (int j = 0; j < fin_mat.num_atoms; j++)
    //                    {
    //                        coor_pdb[0] = coor_pdb[0] + fin_mat.x[j];
    //                        coor_pdb[1] = coor_pdb[1] + fin_mat.y[j];
    //                        coor_pdb[2] = coor_pdb[2] + fin_mat.z[j];
    //                        tmp3 = tmp3 + 1;
    //                    }
    //                    coor_pdb[0] = coor_pdb[0] / (fin_mat.num_atoms);
    //                    coor_pdb[1] = coor_pdb[1] / (fin_mat.num_atoms);
    //                    coor_pdb[2] = coor_pdb[2] / (fin_mat.num_atoms);
    //                    //    cout<<"coor: "<<coor_pdb[0]<<" "<<coor_pdb[1]<<" "<<coor_pdb[2]<<endl;
    //                    //    cout<<"old_CC1,new_CC1:XX "<<old_dE<<" "<<new_dE<<endl;                     
    //                    old_dE = new_dE;
    //                    old_CC = new_CC;
    //                    E_500 = new_dE;
    //                    E_500_int = 1;
    //                    //    new_CC2 = new_CC1;
    //                }
    //            }
    //            nnx = nnx + 1;
    //            //cout << "BBBBB" << best_model_u.x[0] << "  " << best_model_u.y[0] << "  " << best_model_u.z[0] << endl;
    //        }
    //        cout << "acc_rate:  " << acc_rate1 << endl;
    //        cout << best_E_u << endl;
    //        if (E_500_int == 1)
    //        {
    //            E_REMC[jjj] = E_500;
    //        }
    //        decstrz.push_back(point_mat);  // change            
    //    }
    //    vector<DOCKMol >().swap(decstrp);
    //    decstrp = vector<DOCKMol >(supp_num1);
    //    int js = kkk % 2;
    //    if (js == 0)
    //    {
    //        for (int i = 0; i < supp_num1 - 1; i = i + 2)
    //        {
    //            int j = i + 1;
    //            double CH_REMC = exp((1.0 / supp_REMC[i] - 1.0 / supp_REMC[j]) * (E_REMC[i] - E_REMC[j]));
    //            float Pglobal = min(1.0, CH_REMC);
    //            float tmpx = randf0and1();
    //            if (tmpx < Pglobal)
    //            {
    //                decstrp[i] = decstrz[j];
    //                decstrp[j] = decstrz[i];
    //                cout << "aaaaaaa" << endl;
    //            }
    //            else
    //            {
    //                decstrp[i] = decstrz[i];
    //                decstrp[j] = decstrz[j];
    //            }
    //        }
    //    }
    //    if (js == 1)
    //    {
    //        decstrp[0] = decstrz[0];
    //        for (int i = 1; i < supp_num1 - 1; i = i + 2)
    //        {
    //            int j = i + 1;
    //            double CH_REMC = exp((1.0 / supp_REMC[i] - 1.0 / supp_REMC[j]) * (E_REMC[i] - E_REMC[j]));
    //            float Pglobal = min(1.0, CH_REMC);
    //            float tmpx = randf0and1();
    //            if (tmpx < Pglobal)
    //            {
    //                decstrp[i] = decstrz[j];
    //                decstrp[j] = decstrz[i];
    //                cout << "aaaaaaa" << endl;
    //            }
    //            else
    //            {
    //                decstrp[i] = decstrz[i];
    //                decstrp[j] = decstrz[j];
    //            }
    //        }
    //        decstrp[supp_num1 - 1] = decstrz[supp_num1 - 1];
    //    }
    //}
    //ccindex.clear();
    //for (int i = 0; i < best_Ex.size(); i++)
    //{
    //    ccindex.push_back(i);
    //}
    //sortmol_CC(best_Ex, ccindex);
    ////vector<float> centers_Fin(3, 0);
    //for (int i = 0; i < 3; i++)
    //{
    //    centers_Fin[i] = 0;
    //}
    //for (int i = 0; i < best_modelx[0].num_atoms; i++)
    //{
    //    centers_Fin[0] += best_modelx[0].x[i];
    //    centers_Fin[1] += best_modelx[0].y[i];
    //    centers_Fin[2] += best_modelx[0].z[i];
    //}
    //centers_Fin[0] = centers_Fin[0] / best_modelx[0].num_atoms;
    //centers_Fin[1] = centers_Fin[1] / best_modelx[0].num_atoms;
    //centers_Fin[2] = centers_Fin[2] / best_modelx[0].num_atoms;
    //cout << "FIN center  " << centers_Fin[0] << "  " << centers_Fin[1] << "  " << centers_Fin[2] << endl;
    //max_idx = 0;
    //if (mol.num_atoms > 40) max_idx = 0.75 * max;
    //else max_idx = 0.5 * max;
    //theDensityMap.cart(centers_Fin, fartbinding, max_idx);

    ccindex.clear();
    for (int i = 0; i < best_Ex.size(); i++)
    {
        ccindex.push_back(i);
    }
    sortmol_CC(best_Ex, ccindex);

    int j = 0;
    for (int i = 0; i < best_Ex.size(); i++)
    {
        if (i <= 4 || i >= 15)
        //if(i<=1 || i>=18 || i == 10)
        {
            cout << "eeeee" << best_Ex[ccindex[i]] << endl;
            cout << "clCC_1clCC_1clCC_1" << theDensityMap.clCC_1(best_modelx[ccindex[i]]) << endl;
            string molfile1 = "test_" + std::to_string(j) + ".mol2";
            //ofstream clusterpdb(pdbfile.c_str());
            ofstream clustermol1(molfile1.c_str());
            string str1 = string(argv[1]);
            clustermol1 << "########## Name:" << endl;
            Write_Mol2(best_modelx[ccindex[i]], clustermol1, str1);
            j++;
        }
        //for (int j = 0; j < best_modelx[i].num_atoms; j++)
        //{
        //    //cout << best_modelx[i].x[j] << "  " << best_modelx[i].y[j] << "  " << best_modelx[i].z[j] << endl;
        //}
    }

    //cout << "AAAAAA" << best_model_u.x[0] << "  " << best_model_u.y[0] << "  " << best_model_u.z[0] << endl;
    //pointsBx.clear();
    //pointsBx.push_back(best_model_u);
    //swapmol(pointsBx, best_model_u);
    //pointsBx = best_model_u;
    //cout << "AAAAAA" << pointsBx.x[0] << "  " << pointsBx.y[0] << "  " << pointsBx.z[0] << endl;
    //string molfile1 = "test1.mol2";
    ////ofstream clusterpdb(pdbfile.c_str());
    //ofstream clustermol1(molfile1.c_str());
    //string str1 = string(argv[1]);
    //clustermol1 << "########## Name:" << endl;
    //Write_Mol2(pointsBx, clustermol1, str1);


    vector <DOCKMol> copy_mol2;
    for (int i = 0; i < best_Ex.size(); i++)
    {
        DOCKMol trans_mol;
        copy_molecule(trans_mol, best_modelx[ccindex[i]]);
        generatationmol(best_modelx[ccindex[i]], trans_mol, 5, gen);
        copy_mol2.push_back(trans_mol);
    }

    //for (int i = 0; i < initialnum; i++)
    //{
    //    matchvec.push_back(pointsBx);
    //    DOCKMol moltemp;
    //    copy_molecule(moltemp, mol);
    //    matchvec.push_back(moltemp);//just stand one vector, not real conformation.    
    //}
    for (int i = 0; i < best_modelx.size(); i++)
    {
        cout << "2339" << endl;
        energy(amber, best_modelx[i], receptor, c_nrg, filevdw.c_str(), fileflex.c_str(), fileflex_drive_tbl.c_str(), grid_file, cutoff, biolip_matrix, sphx, sphy, sphz, bindsite_weight, protein_atom_index, ave, std, theDensityMap, receptor_mapsit, centers_Fin);
        cout << i << "  current_score: " << best_modelx[i].current_score << "  intral_energy: " << best_modelx[i].intral_energy << "  internal_energy: " << best_modelx[i].internal_energy << "  contact_energy: " << best_modelx[i].contact_energy << "  hb_energy: " << best_modelx[i].hb_energy << "  CC: " << best_modelx[i].current_EM_score << endl;
        if (i <= 4 || i >= 15)
        //if (i <= 1 || i >= 18 || i == 10)
        {
            //DOCKMol trans_mol;
            //copy_molecule(trans_mol, best_modelx[ccindex[i]]);
            //generatationmol(best_modelx[ccindex[i]], trans_mol,5,gen);
            //copy_mol2.push_back(trans_mol);
            matchvec.push_back(best_modelx[ccindex[i]]);
            //string molfile1 = "test_" + std::to_string(j) + ".mol2";
            //j++;
            ////ofstream clusterpdb(pdbfile.c_str());
            //ofstream clustermol1(molfile1.c_str());
            //string str1 = string(argv[1]);
            //clustermol1 << "########## Name:" << endl;
            //Write_Mol2(trans_mol, clustermol1, str1);
            //if(i%2 == 0)
            //{
            //    DOCKMol trans_mol;
            //    copy_molecule(trans_mol, best_modelx[ccindex[i]]);
            //    generatationmol(best_modelx[ccindex[i]], trans_mol,5,gen);
            //    matchvec.push_back(trans_mol);
            //}
            //else matchvec.push_back(best_modelx[ccindex[i]]);
        }
    }
    for (int i = 0; i < copy_mol2.size(); i++)
    {
        energy(amber, copy_mol2[i], receptor, c_nrg, filevdw.c_str(), fileflex.c_str(), fileflex_drive_tbl.c_str(), grid_file, cutoff, biolip_matrix, sphx, sphy, sphz, bindsite_weight, protein_atom_index, ave, std, theDensityMap, receptor_mapsit, centers_Fin);
        if (i <= 4 || i >= 15)
        {
            matchvec.push_back(copy_mol2[i]);
            string molfile1 = "test_" + std::to_string(j) + ".mol2";
            j++;
            //ofstream clusterpdb(pdbfile.c_str());
            ofstream clustermol1(molfile1.c_str());
            string str1 = string(argv[1]);
            clustermol1 << "########## Name:" << endl;
            Write_Mol2(copy_mol2[i], clustermol1, str1);
        }
    }
    cout << copy_mol2.size() << endl;
    cout << matchvec.size() << endl;
        //matchvec.push_back(best_modelx[ccindex[i]]);
        //DOCKMol moltemp;
        //copy_molecule(moltemp, mol);
        //matchvec.push_back(moltemp);

    cout << "matchvec.size()" << matchvec.size() << endl;
    for (int i = 0; i < matchvec.size(); i++)
    {
        cout << i << " " << matchvec[i].x[0] << "   " << matchvec[i].current_score << "  " << matchvec[i].current_EM_score << endl;
    }
    /******************************************/
    //vector <DOCKMol> matchvec2;
    //vector<float> axis1(3, 0);
    //vector<float> axis2(3, 0);
    //double a, b, c;
    //for (int i = 0; i < matchvec.size() - 1; i = i + 2)
    //{
    //    for (int m = 0; m < mol.num_atoms; m++)
    //    {
    //        mol.x[m] = matchvec[i].x[m];
    //        mol.y[m] = matchvec[i].y[m];
    //        mol.z[m] = matchvec[i].z[m];
    //    }
    //    float ligcenterx, ligcentery, ligcenterz;
    //    ligcenterx = ligcentery = ligcenterz = 0;
    //    for (int m = 0; m < mol.num_atoms; m++)
    //    {
    //        ligcenterx += mol.x[m];
    //        ligcentery += mol.y[m];
    //        ligcenterz += mol.z[m];
    //    }
    //    ligcenterx = ligcenterx / mol.num_atoms;
    //    ligcentery = ligcentery / mol.num_atoms;
    //    ligcenterz = ligcenterz / mol.num_atoms;
    //    //cout<<ligcenterx<<"\t"<<ligcentery<<"\t"<<ligcenterz<<endl;
    //    axis1[0] = ligcenterx; axis1[1] = ligcentery; axis1[2] = ligcenterz;
    //    a = rand0to1(gen);
    //    b = rand0to1(gen);
    //    //cout<<"a"<<a<<" "<<"b"<<b<<" "<<"c"<<c<<endl;
    //    double theta = 2 * 3.1415926 * a;
    //    double phi = acos(2 * b - 1.0);
    //    //cout<<"theta "<<theta<<" phi "<<phi<<endl;
    //    axis2[0] = sin(phi) * cos(theta) + ligcenterx;
    //    axis2[1] = sin(phi) * sin(theta) + ligcentery;
    //    axis2[2] = cos(phi) + ligcenterz;
    //    c = randgeneration(gen);
    //    float angle = 0.0;
    //    if (c >= 0)
    //        angle = 180.0;
    //    else
    //        angle = -180.0;
    //    //cout<<"rotation angle "<<angle<<endl;
    //    //cout<<mol.x[0]<<"****"<<endl;
    //    GroupRotation(axis2, axis1, angle, mol); // random rotate the mol 
    //    MC_rigid_ligand(amber, mol, receptor, c_nrg, filevdw.c_str(), fileflex.c_str(), fileflex_drive_tbl.c_str(), grid_file, 1.00, 1.00, gen, biolip_matrix, sphx, sphy, sphz, bindsite_weight, protein_atom_index, ave, std, theDensityMap,receptor_bindingsit);
    //    //cout<<i<<"after MC "<<mol.x[0]<<" "<<mol.current_score<<" before "<<matchvec[i].x[0]<<" "<<matchvec[i].current_score<<endl;
    //    if (!(energy(amber, mol, receptor, c_nrg, filevdw.c_str(), fileflex.c_str(), fileflex_drive_tbl.c_str(), grid_file, cutoff, biolip_matrix, sphx, sphy, sphz, bindsite_weight, protein_atom_index, ave, std, theDensityMap, receptor_bindingsit)))
    //    {
    //        mol.current_score = 1e+6;
    //        mol.internal_energy = 1e+6;
    //        mol.intral_energy = 1e+6;
    //    }
    //    //cout<<"before"<<matchvec[i].x[0]<<"\t"<<matchvec[i+1].x[0]<<endl;
    //    for (int m = 0; m < mol.num_atoms; m++)
    //    {
    //        matchvec[i + 1].x[m] = mol.x[m];
    //        //cout<<matchvec[i].x[m]<<endl;
    //        matchvec[i + 1].y[m] = mol.y[m];
    //        matchvec[i + 1].z[m] = mol.z[m];
    //    }
    //    matchvec[i + 1].current_score = mol.current_score;
    //    matchvec[i + 1].internal_energy = mol.internal_energy;
    //    matchvec[i + 1].intral_energy = mol.intral_energy;
    //    matchvec[i + 1].energy_tor_total = mol.energy_tor_total;

    //}

  /*****************************REMC simulation****************************************************/
  //REMC_flag = false; //not do remc
    //bool REMC_flag = 1;
    vector <DOCKMol> REMC_best_model;
    for(int i =0;i< matchvec.size();i++)
    {
        DOCKMol best;
        copy_molecule(best, mol);
        REMC_best_model.push_back(best);
        //REMC_best_model.push_back(best);
    }
    vector<float> REMC_best_E(matchvec.size(), 1000000000.0);

    if (REMC_flag == true)
    {
        vector<DOCKMol> clustermols, clustermols1, clustermols2;
        vector<vector<float> > x, y, z;
        vector<float> accenergy;
        if (REMC(amber, matchvec, receptor, c_nrg, filevdw.c_str(), fileflex.c_str(), fileflex_drive_tbl.c_str(), grid_file, 1, x, y, z, accenergy, Tmin, Tmax, outfile, gen, flexible_flag, biolip_matrix, swapnumber, MC_steps, sphx, sphy, sphz, bindsite_weight, protein_atom_index, ave, std, theDensityMap, REMC_best_model, REMC_best_E, receptor_mapsit, centers_Fin))
        {
            //string pdbfile = string(argv[1])+".pdb";
            string molfile = string(argv[1]) + "_13.mol2";
            //ofstream clusterpdb(pdbfile.c_str());
            ofstream clustermol(molfile.c_str());
            for (int i = 0; i < REMC_best_model.size(); i++)
            {
                //cout<<"i="<<i<<endl;
                clustermols.push_back(mol);
                swapmol(clustermols[i], REMC_best_model[i]);
                //for (int j = 0; j < mol.num_atoms; j++)
                //{
                //    //cout<<i<<"\t"<<j<<"\t"<<clustermols[i].atom_names[j]<<endl;
                //    //cout<< "ATOM  " << setw(5) << j+1 << " " << setw(4)<<mol.atom_names[j] << " " << setw(3)<<"ALA" <<"  "<<setw(4) << j+1 <<"   " ;
                //    //cout<<setiosflags(ios::fixed)<<setprecision(3)<<setw(8)<<x[i][j]<<setw(8)<<y[i][j]<<setw(8)<<z[i][j]<<"\n";
                //    clustermols[clustermols.size() - 1].x[j] = x[i][j];
                //    clustermols[clustermols.size() - 1].y[j] = y[i][j];
                //    clustermols[clustermols.size() - 1].z[j] = z[i][j];               
                //}
                //cout<<"\n";
                //clustermols[clustermols.size() - 1].current_score = accenergy[i];
                //outpdb(clustermols[i],clusterpdb,i);
                stringstream ss;
                ss << i;
                string str = string(argv[1]) + ss.str();
                clustermol << "########## Name:" << "\t" << i << "\t" << clustermols[i].current_score <<"  "<< clustermols[i].current_EM_score << endl;
                Write_Mol2(clustermols[i], clustermol, str);
            }
            for (int i = 0; i < best_modelx.size(); i++)
            {
                clustermols1.push_back(mol);
                swapmol(clustermols1[i], best_modelx[i]);
                stringstream ss;
                ss << i;
                string str = string(argv[1]) + ss.str();
                clustermol << "########## Name:" << "\t" << i << "\t" << clustermols1[i].current_score << "  " << clustermols[i].current_EM_score << endl;
                Write_Mol2(clustermols1[i], clustermol, str);               
            }
            for (int i = 0; i < copy_mol2.size(); i++)
            {
                clustermols2.push_back(mol);
                swapmol(clustermols2[i], copy_mol2[i]);
                stringstream ss;
                ss << i;
                string str = string(argv[1]) + ss.str();
                clustermol << "########## Name:" << "\t" << i << "\t" << clustermols2[i].current_score << "  " << clustermols[i].current_EM_score << endl;
                Write_Mol2(clustermols2[i], clustermol, str);
            }

            //for (int i = 0; i < x.size(); i++)
            //{
            //    //cout<<"i="<<i<<endl;
            //    clustermols.push_back(mol);
            //    for (int j = 0; j < mol.num_atoms; j++)
            //    {
            //        //cout<<i<<"\t"<<j<<"\t"<<clustermols[i].atom_names[j]<<endl;
            //        //cout<< "ATOM  " << setw(5) << j+1 << " " << setw(4)<<mol.atom_names[j] << " " << setw(3)<<"ALA" <<"  "<<setw(4) << j+1 <<"   " ;
            //        //cout<<setiosflags(ios::fixed)<<setprecision(3)<<setw(8)<<x[i][j]<<setw(8)<<y[i][j]<<setw(8)<<z[i][j]<<"\n";
            //        clustermols[clustermols.size() - 1].x[j] = x[i][j];
            //        clustermols[clustermols.size() - 1].y[j] = y[i][j];
            //        clustermols[clustermols.size() - 1].z[j] = z[i][j];
            //    }
            //    //cout<<"\n";
            //    clustermols[clustermols.size() - 1].current_score = accenergy[i];
            //    //outpdb(clustermols[i],clusterpdb,i);
            //    stringstream ss;
            //    ss << i;
            //    string str = string(argv[1]) + ss.str();
            //    clustermol << "########## Name:" << "\t" << i << "\t" << clustermols[i].current_score << endl;
            //    Write_Mol2(clustermols[i], clustermol, str);
            //}
            cout << "total conformations " << x.size() << endl;

            for (int i = 0; i < REMC_best_model.size(); i++)
            {
                cout << REMC_best_model[i].current_score <<"  "<< REMC_best_model[i].current_EM_score<< endl;
                string molfile1 = "test_" + std::to_string(i) + "_F.mol2";
                //ofstream clusterpdb(pdbfile.c_str());
                ofstream clustermol1(molfile1.c_str());
                string str1 = string(argv[1]);
                clustermol1 << "########## Name:" << endl;
                Write_Mol2(REMC_best_model[i], clustermol1, str1);

                //int j = i;
                //DOCKMol trans_mol;
                //copy_molecule(trans_mol, REMC_best_model[i]);
                //generatationmol(REMC_best_model[i], trans_mol,5,gen);
                //string molfile2 = "test_" + std::to_string(j+10) + "_F.mol2";
                ////ofstream clusterpdb(pdbfile.c_str());
                //ofstream clustermol2(molfile2.c_str());
                //string str2 = string(argv[1]);
                //clustermol2 << "########## Name:" << endl;
                //Write_Mol2(trans_mol, clustermol2, str2);
            }

            //no conformations can be generated by REMC,using the initial conformations as the final results
            if (x.size() == 0)
            {
                for (int i = 0; i < matchvec.size(); i++)
                {
                    //outpdb(matchvec[i],clusterpdb,i);
                    stringstream ss;
                    ss << i;
                    string str = string(argv[1]) + ss.str();
                    clustermol << "########## Name:" << "\t" << i << "\t" << matchvec[i].current_score << endl;
                    Write_Mol2(matchvec[i], clustermol, str);
                }
            }
            cout << "REMC simulation is successful!" << endl;
            //clusterpdb.close();
            //clustermol.close();	
        }
        else
        {
            cout << "MC simulation is not successful! Don't cry! fighting!" << endl;
            REMC_flag = false;
        }
    }

    if (REMC_flag == false)
    {
        //string pdbfile = string(argv[1])+".pdb";
        string molfile = string(argv[1]) + ".mol2";
        //ofstream clusterpdb(pdbfile.c_str());
        ofstream clustermol(molfile.c_str());
        //outpdb(mol,clusterpdb,1);
        for (int i = 0; i < matchvec.size(); i++)
        {
            stringstream ss;
            ss << i;
            string str = string(argv[1]) + ss.str();
            clustermol << "########## Name:" << "\t" << i << "\t" << matchvec[i].current_score << endl;
            Write_Mol2(matchvec[i], clustermol, str);
        }
    }
  return 1;
 }



