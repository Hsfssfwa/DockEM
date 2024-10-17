/*
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <cstring>
#include <iomanip>
#include <vector>
#include <math.h>  */
#include "Rotbuilder.h"

// using namespace std;


vector<float> calculateCoordinates(vector<float>& refA,vector<float>& refB,vector<float>& refC,float L,float ang,float di)
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

    float F = (sqrt(BX*BX + BY*BY + BZ*BZ)) * L * cos(ang * (3.1415926/180.0));
//  float F = sqrt(BX*BX + BY*BY + BZ*BZ) * L * cos(ang);

    float constx = sqrt((pow(B*BZ-BY*G,2))*(-(F*F)*(A*A+B*B+G*G)+(B*B*(BX*BX+BZ*BZ) + A*A*(BY*BY+BZ*BZ)- (2*A*BX*BZ*G) + (BX*BX+ BY*BY)*G*G - (2*B*BY)*(A*BX+BZ*G))*L*L));
    float denom = (B*B)*(BX*BX+BZ*BZ)+ (A*A)*(BY*BY+BZ*BZ) - (2*A*BX*BZ*G) + (BX*BX+BY*BY)*(G*G) - (2*B*BY)*(A*BX+BZ*G);

    float X= ((B*B*BX*F)-(A*B*BY*F)+(F*G)*(-A*BZ+BX*G)+constx)/denom;

    float Y=0.0,Z=0.0;
    if((B==0.0 || BZ==0.0) && (BY==0.0 || G==0.0))
    {
        float const1=sqrt( G*G*(-A*A*X*X+(B*B+G*G)*(L-X)*(L+X)));
        Y= ((-A*B*X)+const1)/(B*B+G*G);
        Z= -(A*G*G*X+B*const1)/(G*(B*B+G*G));
    }
    else
    {
        Y= ((A*A*BY*F)*(B*BZ-BY*G)+ G*( -F*pow(B*BZ-BY*G,2.0) + BX*constx) - A*( B*B*BX*BZ*F- B*BX*BY*F*G + BZ*constx)) / ((B*BZ-BY*G)*denom);
        Z= ((A*A*BZ*F)*(B*BZ-BY*G) + (B*F)*pow(B*BZ-BY*G,2.0) + (A*BX*F*G)*(-B*BZ+BY*G) - B*BX*constx + A*BY*constx) / ((B*BZ-BY*G)*denom);
    }
    vector<float> D(3,0.0);
    D[0] = X + CV[0];
    D[1] = Y + CV[1];
    D[2] = Z + CV[2];

    float temp=0.0;

    temp=Points2Dihedral(AV,BV,CV,D) * (180.0/3.1415926);
//    cout<<"temp: "<<temp<<endl;

    di = di - temp;
    vector<float> Dx(3,0.0);
//    Dx[0] = D[0] - BV[0];
//    Dx[1] = D[1] - BV[1];
//    Dx[2] = D[2] - BV[2];
    Dx[0] = D[0];
    Dx[1] = D[1];
    Dx[2] = D[2];    

//    di = (di * (3.1415926/180.0));
    GroupRotationt(CV,BV,di,Dx);
//    D[0] = Dx[0] + BV[0];
//    D[1] = Dx[1] + BV[1];
//    D[2] = Dx[2] + BV[2];
    D[0] = Dx[0] ;
    D[1] = Dx[1] ;
    D[2] = Dx[2] ;    

    return D;
}

Rotm makeGly(vector<float > N, vector<float > CA, vector<float > C)
{
 ///   '''Creates a Glycine residue'''
 ///   ##Create Residue Data Structure
    Rotm resx;

    float C_O_length=1.23;
    float CA_C_O_angle=120.5117;
    float N_CA_C_O_diangle=180.0;

    resx.Ocd=calculateCoordinates(N, CA, C, C_O_length, CA_C_O_angle, N_CA_C_O_diangle);
    return resx;    
}

Rotm makeAla(vector<float > N, vector<float > CA, vector<float > C)
{
//    '''Creates an Alanine residue'''
//    ##R-Group

    float C_O_length=1.23;
    float CA_C_O_angle=120.5;
    float N_CA_C_O_diangle=-60.5;

    float CA_CB_length=1.52;
    float C_CA_CB_angle=109.5;
    float N_C_CA_CB_diangle=122.6860;

    vector<float> carbon_b= calculateCoordinates(N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle);
    Rotm resx;

    resx.Rotcd.push_back(carbon_b);
    resx.Rotnm.push_back("CB");

    resx.Ocd=calculateCoordinates(N, CA, C, C_O_length, CA_C_O_angle, N_CA_C_O_diangle);

    return resx;    
}


Rotm makeSer(vector<float > N, vector<float > CA, vector<float > C)
{
//    '''Creates a Serine residue'''
 //   ##R-Group
    float C_O_length=1.23;
    float CA_C_O_angle=120.5;
    float N_CA_C_O_diangle=-60.0;

    float CA_CB_length=1.52;
    float C_CA_CB_angle=109.5;
    float N_C_CA_CB_diangle=122.6618;

    float CB_OG_length=1.417;
    float CA_CB_OG_angle=110.773;
    float N_CA_CB_OG_diangle=-63.3;

    Rotm resx;
    
    vector<float > carbon_b= calculateCoordinates(N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle);
//    CB= Atom("CB", carbon_b, 0.0 , 1.0, " "," CB", 0,"C")
    resx.Rotcd.push_back(carbon_b);
    resx.Rotnm.push_back("CB");
    vector<float> CB = carbon_b;  
    vector<float > oxygen_g= calculateCoordinates(N, CA, CB, CB_OG_length, CA_CB_OG_angle, N_CA_CB_OG_diangle);
//    OG= Atom("OG", oxygen_g, 0.0, 1.0, " ", " OG", 0, "O")
    resx.Rotcd.push_back(oxygen_g);
    resx.Rotnm.push_back("OG");

    resx.Ocd=calculateCoordinates(N, CA, C, C_O_length, CA_C_O_angle, N_CA_C_O_diangle);

//    ##print res
    return resx;    
}


Rotm makeCys(vector<float > N, vector<float > CA,vector<float > C)
{
//    '''Creates a Cysteine residue'''
//    ##R-Group
    float C_O_length=1.23;
    float CA_C_O_angle=120.5;
    float N_CA_C_O_diangle=-60.0;

    float CA_CB_length=1.52;
    float C_CA_CB_angle=109.5;
    float N_C_CA_CB_diangle=122.5037;
    
    float CB_SG_length= 1.808;
    float CA_CB_SG_angle= 113.8169;
    float N_CA_CB_SG_diangle= -62.2;

    Rotm resx;

    vector<float > carbon_b= calculateCoordinates(N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle);
//    CB= Atom("CB", carbon_b, 0.0 , 1.0, " "," CB", 0,"C")
    resx.Rotcd.push_back(carbon_b);
    resx.Rotnm.push_back("CB");
    vector<float> CB = carbon_b;      
    vector<float > sulfur_g= calculateCoordinates(N, CA, CB, CB_SG_length, CA_CB_SG_angle, N_CA_CB_SG_diangle);
//    SG= Atom("SG", sulfur_g, 0.0, 1.0, " ", " SG", 0, "S")
    resx.Rotcd.push_back(sulfur_g);
    resx.Rotnm.push_back("SG");  

    resx.Ocd=calculateCoordinates(N, CA, C, C_O_length, CA_C_O_angle, N_CA_C_O_diangle);  

    return resx;
}

Rotm makeVal(vector<float > N, vector<float > CA, vector<float > C)
{
///    '''Creates a Valine residue'''
//    ##R-Group
    float C_O_length=1.23;
    float CA_C_O_angle=120.5686;
    float N_CA_C_O_diangle=-60.0;

    float CA_CB_length=1.52;
    float C_CA_CB_angle=109.5;
    float N_C_CA_CB_diangle=123.2347;
    
    float CB_CG1_length=1.527;
    float CA_CB_CG1_angle=110.7;
    float N_CA_CB_CG1_diangle=177.2;
    
    float CB_CG2_length=1.527;
    float CA_CB_CG2_angle=110.4;
    float N_CA_CB_CG2_diangle=-63.3;

    Rotm resx;

    vector<float > carbon_b= calculateCoordinates(N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle);
//    CB= Atom("CB", carbon_b, 0.0 , 1.0, " "," CB", 0,"C")
    resx.Rotcd.push_back(carbon_b);
    resx.Rotnm.push_back("CB");
    vector<float> CB = carbon_b;        
    vector<float > carbon_g1= calculateCoordinates(N, CA, CB, CB_CG1_length, CA_CB_CG1_angle, N_CA_CB_CG1_diangle);
//    CG1= Atom("CG1", carbon_g1, 0.0, 1.0, " ", " CG1", 0, "C")
    resx.Rotcd.push_back(carbon_g1);
    resx.Rotnm.push_back("CG1");      
    vector<float > carbon_g2= calculateCoordinates(N, CA, CB, CB_CG2_length, CA_CB_CG2_angle, N_CA_CB_CG2_diangle);
//    CG2= Atom("CG2", carbon_g2, 0.0, 1.0, " ", " CG2", 0, "C")
    resx.Rotcd.push_back(carbon_g2);
    resx.Rotnm.push_back("CG2");     

    resx.Ocd=calculateCoordinates(N, CA, C, C_O_length, CA_C_O_angle, N_CA_C_O_diangle); 

    return resx; 
}


Rotm makeIle(vector<float > N, vector<float > CA, vector<float > C)
{
//    '''Creates an Isoleucine residue'''
//    ##R-group
    float C_O_length=1.23;
    float CA_C_O_angle=120.5403;
    float N_CA_C_O_diangle= -60.0;

    float CA_CB_length=1.52;
    float C_CA_CB_angle=109.5;
    float N_C_CA_CB_diangle=123.2347;
    
    float CB_CG1_length=1.527;
    float CA_CB_CG1_angle=110.7;
    float N_CA_CB_CG1_diangle=-61.6;
        
    float CB_CG2_length=1.527;
    float CA_CB_CG2_angle=110.4;
    float N_CA_CB_CG2_diangle= 59.7;

    float CG1_CD1_length= 1.52;
    float CB_CG1_CD1_angle= 113.97;
    float CA_CB_CG1_CD1_diangle= 169.8;

    Rotm resx;
    
    vector<float > carbon_b= calculateCoordinates(N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle);
//    CB= Atom("CB", carbon_b, 0.0 , 1.0, " "," CB", 0,"C")
    resx.Rotcd.push_back(carbon_b);
    resx.Rotnm.push_back("CB");
    vector<float> CB = carbon_b;       
    vector<float > carbon_g1= calculateCoordinates(N, CA, CB, CB_CG1_length, CA_CB_CG1_angle, N_CA_CB_CG1_diangle);
//    CG1= Atom("CG1", carbon_g1, 0.0, 1.0, " ", " CG1", 0, "C")
    resx.Rotcd.push_back(carbon_g1);
    resx.Rotnm.push_back("CG1"); 
    vector<float > carbon_g2= calculateCoordinates(N, CA, CB, CB_CG2_length, CA_CB_CG2_angle, N_CA_CB_CG2_diangle);
//    CG2= Atom("CG2", carbon_g2, 0.0, 1.0, " ", " CG2", 0, "C")
    resx.Rotcd.push_back(carbon_g2);
    resx.Rotnm.push_back("CG2");
    vector<float> CG1 = carbon_g1;    
    vector<float > carbon_d1= calculateCoordinates(CA, CB, CG1, CG1_CD1_length, CB_CG1_CD1_angle, CA_CB_CG1_CD1_diangle);
 //   CD1= Atom("CD1", carbon_d1, 0.0, 1.0, " ", " CD1", 0, "C")
    resx.Rotcd.push_back(carbon_d1);
    resx.Rotnm.push_back("CD1");  

    resx.Ocd=calculateCoordinates(N, CA, C, C_O_length, CA_C_O_angle, N_CA_C_O_diangle);

    return resx;
}

Rotm makeLeu(vector<float > N, vector<float > CA, vector<float > C)
{
//    '''Creates a Leucine residue'''
//    ##R-Group
    float C_O_length=1.23;
    float CA_C_O_angle=120.4647;
    float N_CA_C_O_diangle=120.0;

    float CA_CB_length=1.52;
    float C_CA_CB_angle=109.5;
    float N_C_CA_CB_diangle=122.4948;

    float CB_CG_length=1.53;
    float CA_CB_CG_angle= 116.10;
    float N_CA_CB_CG_diangle=-60.1;
    
    float CG_CD1_length=1.524;
    float CB_CG_CD1_angle=110.27;
    float CA_CB_CG_CD1_diangle=174.9;

    float CG_CD2_length=1.525;
    float CB_CG_CD2_angle=110.58;
    float CA_CB_CG_CD2_diangle=66.7;

    Rotm resx;

    vector<float > carbon_b= calculateCoordinates(N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle);
//    CB= Atom("CB", carbon_b, 0.0 , 1.0, " "," CB", 0,"C")
    resx.Rotcd.push_back(carbon_b);
    resx.Rotnm.push_back("CB"); 
    vector<float> CB = carbon_b;  
    vector<float > carbon_g1= calculateCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle);
//    CG= Atom("CG", carbon_g1, 0.0, 1.0, " ", " CG", 0, "C")
    resx.Rotcd.push_back(carbon_g1);
    resx.Rotnm.push_back("CG"); 
    vector<float> CG = carbon_g1;
    vector<float > carbon_d1= calculateCoordinates(CA, CB, CG, CG_CD1_length, CB_CG_CD1_angle, CA_CB_CG_CD1_diangle);
//    CD1= Atom("CD1", carbon_d1, 0.0, 1.0, " ", " CD1", 0, "C")
    resx.Rotcd.push_back(carbon_d1);
    resx.Rotnm.push_back("CD1"); 
    vector<float > carbon_d2= calculateCoordinates(CA, CB, CG, CG_CD2_length, CB_CG_CD2_angle, CA_CB_CG_CD2_diangle);
//    CD2= Atom("CD2", carbon_d2, 0.0, 1.0, " ", " CD2", 0, "C")
    resx.Rotcd.push_back(carbon_d2);
    resx.Rotnm.push_back("CD2"); 

    resx.Ocd=calculateCoordinates(N, CA, C, C_O_length, CA_C_O_angle, N_CA_C_O_diangle);

    return resx;
}
    
Rotm makeThr(vector<float > N, vector<float > CA, vector<float > C)
{
//    '''Creates a Threonine residue'''
//    ##R-Group
    float C_O_length=1.23;
    float CA_C_O_angle=120.5359;
    float N_CA_C_O_diangle=120.0;

    float CA_CB_length=1.52;
    float C_CA_CB_angle=109.5;
    float N_C_CA_CB_diangle=123.0953;
    
    float CB_OG1_length=1.43;
    float CA_CB_OG1_angle=109.18;
    float N_CA_CB_OG1_diangle=-60.3;
        
    float CB_CG2_length=1.53;
    float CA_CB_CG2_angle=111.13;
    float N_CA_CB_CG2_diangle= 60.0;

    Rotm resx;

    vector<float > carbon_b= calculateCoordinates(N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle);
//    CB= Atom("CB", carbon_b, 0.0 , 1.0, " "," CB", 0,"C")
    resx.Rotcd.push_back(carbon_b);
    resx.Rotnm.push_back("CB"); 
    vector<float> CB = carbon_b;   
    vector<float > oxygen_g1= calculateCoordinates(N, CA, CB, CB_OG1_length, CA_CB_OG1_angle, N_CA_CB_OG1_diangle);
//    OG1= Atom("OG1", oxygen_g1, 0.0, 1.0, " ", " OG1", 0, "O")
    resx.Rotcd.push_back(oxygen_g1);
    resx.Rotnm.push_back("OG1");    
    vector<float > carbon_g2= calculateCoordinates(N, CA, CB, CB_CG2_length, CA_CB_CG2_angle, N_CA_CB_CG2_diangle);
//    CG2= Atom("CG2", carbon_g2, 0.0, 1.0, " ", " CG2", 0, "C")
    resx.Rotcd.push_back(carbon_g2);
    resx.Rotnm.push_back("CG2"); 

    resx.Ocd=calculateCoordinates(N, CA, C, C_O_length, CA_C_O_angle, N_CA_C_O_diangle);

    return resx;
}

Rotm makeArg(vector<float > N, vector<float > CA, vector<float > C)
{
//    '''Creates an Arginie residue'''
//    ##R-Group
    float C_O_length=1.23;
    float CA_C_O_angle=120.54;
    float N_CA_C_O_diangle=120.0;

    float CA_CB_length=1.52;
    float C_CA_CB_angle=109.5;
    float N_C_CA_CB_diangle=122.76;

    float CB_CG_length=1.52;
    float CA_CB_CG_angle= 113.83;
    float N_CA_CB_CG_diangle=-65.2;
    
    float CG_CD_length=1.52;
    float CB_CG_CD_angle=111.79;
    float CA_CB_CG_CD_diangle=-179.2;
    
    float CD_NE_length=1.46;
    float CG_CD_NE_angle=111.68;
    float CB_CG_CD_NE_diangle=-179.3;

    float NE_CZ_length=1.33;
    float CD_NE_CZ_angle=124.79;
    float CG_CD_NE_CZ_diangle=-178.7;

    float CZ_NH1_length=1.33;
    float NE_CZ_NH1_angle=120.64;
    float CD_NE_CZ_NH1_diangle=0.0;

    float CZ_NH2_length=1.33;
    float NE_CZ_NH2_angle=119.63;
    float CD_NE_CZ_NH2_diangle=180.0;

    Rotm resx;
    
    vector<float > carbon_b= calculateCoordinates(N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle);
//    CB= Atom("CB", carbon_b, 0.0 , 1.0, " "," CB", 0,"C")
    resx.Rotcd.push_back(carbon_b);
    resx.Rotnm.push_back("CB"); 
    vector<float> CB = carbon_b;    
    vector<float > carbon_g= calculateCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle);
//    CG= Atom("CG", carbon_g, 0.0, 1.0, " ", " CG", 0, "C")
    resx.Rotcd.push_back(carbon_g);
    resx.Rotnm.push_back("CG");
    vector<float> CG = carbon_g; 
    vector<float > carbon_d= calculateCoordinates(CA, CB, CG, CG_CD_length, CB_CG_CD_angle, CA_CB_CG_CD_diangle);
//    CD= Atom("CD", carbon_d, 0.0, 1.0, " ", " CD", 0, "C")
    resx.Rotcd.push_back(carbon_d);
    resx.Rotnm.push_back("CD");
    vector<float> CD = carbon_d; 
    vector<float > nitrogen_e= calculateCoordinates(CB, CG, CD, CD_NE_length, CG_CD_NE_angle, CB_CG_CD_NE_diangle);
//    NE= Atom("NE", nitrogen_e, 0.0, 1.0, " ", " NE", 0, "N")
    resx.Rotcd.push_back(nitrogen_e);
    resx.Rotnm.push_back("NE"); 
    vector<float> NE = nitrogen_e;    
    vector<float > carbon_z= calculateCoordinates(CG, CD, NE, NE_CZ_length, CD_NE_CZ_angle, CG_CD_NE_CZ_diangle);
//    CZ= Atom("CZ", carbon_z, 0.0, 1.0, " ", " CZ", 0, "C")
    resx.Rotcd.push_back(carbon_z);
    resx.Rotnm.push_back("CZ"); 
    vector<float> CZ = carbon_z;     
    vector<float > nitrogen_h1= calculateCoordinates(CD, NE, CZ, CZ_NH1_length, NE_CZ_NH1_angle, CD_NE_CZ_NH1_diangle);
//    NH1= Atom("NH1", nitrogen_h1, 0.0, 1.0, " ", " NH1", 0, "N")
    resx.Rotcd.push_back(nitrogen_h1);
    resx.Rotnm.push_back("NH1");      
    vector<float > nitrogen_h2= calculateCoordinates(CD, NE, CZ, CZ_NH2_length, NE_CZ_NH2_angle, CD_NE_CZ_NH2_diangle);
//    NH2= Atom("NH2", nitrogen_h2, 0.0, 1.0, " ", " NH2", 0, "N")
    resx.Rotcd.push_back(nitrogen_h2);
    resx.Rotnm.push_back("NH2");      

    resx.Ocd=calculateCoordinates(N, CA, C, C_O_length, CA_C_O_angle, N_CA_C_O_diangle);

    return resx;
}

Rotm makeLys(vector<float > N, vector<float > CA, vector<float > C)
{
//    '''Creates a Lysine residue'''
 //   ##R-Group
    float C_O_length=1.23;
    float CA_C_O_angle=120.54;
    float N_CA_C_O_diangle=120.0;

    float CA_CB_length=1.52;
    float C_CA_CB_angle=109.5;
    float N_C_CA_CB_diangle=122.76;

    float CB_CG_length=1.52;
    float CA_CB_CG_angle=113.83;
    float N_CA_CB_CG_diangle=-64.5;

    float CG_CD_length=1.52;
    float CB_CG_CD_angle=111.79;
    float CA_CB_CG_CD_diangle=-178.1;

    float CD_CE_length=1.46;
    float CG_CD_CE_angle=111.68;
    float CB_CG_CD_CE_diangle=-179.6;

    float CE_NZ_length=1.33;
    float CD_CE_NZ_angle=124.79;
    float CG_CD_CE_NZ_diangle=179.6;

    Rotm resx;
    
    vector<float > carbon_b= calculateCoordinates(N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle);
//    CB= Atom("CB", carbon_b, 0.0 , 1.0, " "," CB", 0,"C")
    resx.Rotcd.push_back(carbon_b);
    resx.Rotnm.push_back("CB"); 
    vector<float> CB = carbon_b;  
    vector<float > carbon_g= calculateCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle);
//    CG= Atom("CG", carbon_g, 0.0, 1.0, " ", " CG", 0, "C")
    resx.Rotcd.push_back(carbon_g);
    resx.Rotnm.push_back("CG");
    vector<float> CG = carbon_g;   
    vector<float > carbon_d= calculateCoordinates(CA, CB, CG, CG_CD_length, CB_CG_CD_angle, CA_CB_CG_CD_diangle);
//    CD= Atom("CD", carbon_d, 0.0, 1.0, " ", " CD", 0, "C")
    resx.Rotcd.push_back(carbon_d);
    resx.Rotnm.push_back("CD");
    vector<float> CD = carbon_d;       
    vector<float > carbon_e= calculateCoordinates(CB, CG, CD, CD_CE_length, CG_CD_CE_angle, CB_CG_CD_CE_diangle);
//    CE= Atom("CE", carbon_e, 0.0, 1.0, " ", " CE", 0, "C")
    resx.Rotcd.push_back(carbon_e);
    resx.Rotnm.push_back("CE"); 
    vector<float> CE = carbon_e;
    vector<float > nitrogen_z= calculateCoordinates(CG, CD, CE, CE_NZ_length, CD_CE_NZ_angle, CG_CD_CE_NZ_diangle);
 //   NZ= Atom("NZ", nitrogen_z, 0.0, 1.0, " ", " NZ", 0, "N")
    resx.Rotcd.push_back(nitrogen_z);
    resx.Rotnm.push_back("NZ");   

    resx.Ocd=calculateCoordinates(N, CA, C, C_O_length, CA_C_O_angle, N_CA_C_O_diangle);

    return resx;    
}


Rotm makeAsp(vector<float > N, vector<float > CA, vector<float > C)
{
//    '''Creates an Aspartic Acid residue'''
//    ##R-Group
    float C_O_length=1.23;
    float CA_C_O_angle=120.51;
    float N_CA_C_O_diangle=120.0;

    float CA_CB_length=1.52;
    float C_CA_CB_angle=109.5;
    float N_C_CA_CB_diangle=122.82;

    float CB_CG_length=1.52;
    float CA_CB_CG_angle=113.06;
    float N_CA_CB_CG_diangle=-66.4;

    float CG_OD1_length=1.25;
    float CB_CG_OD1_angle=119.22;
    float CA_CB_CG_OD1_diangle=-46.7;

    float CG_OD2_length=1.25;
    float CB_CG_OD2_angle=118.208;
    float CA_CB_CG_OD2_diangle=133.3;

    Rotm resx;

    vector<float > carbon_b= calculateCoordinates(N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle);
//    CB= Atom("CB", carbon_b, 0.0 , 1.0, " "," CB", 0,"C")
    resx.Rotcd.push_back(carbon_b);
    resx.Rotnm.push_back("CB"); 
    vector<float> CB = carbon_b; 
    vector<float > carbon_g= calculateCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle);
//    CG= Atom("CG", carbon_g, 0.0, 1.0, " ", " CG", 0, "C")
    resx.Rotcd.push_back(carbon_g);
    resx.Rotnm.push_back("CG");
    vector<float> CG = carbon_g;  
    vector<float > oxygen_d1= calculateCoordinates(CA, CB, CG, CG_OD1_length, CB_CG_OD1_angle, CA_CB_CG_OD1_diangle);
//    OD1= Atom("OD1", oxygen_d1, 0.0, 1.0, " ", " OD1", 0, "O")
    resx.Rotcd.push_back(oxygen_d1);
    resx.Rotnm.push_back("OD1");  
    vector<float > oxygen_d2= calculateCoordinates(CA, CB, CG, CG_OD2_length, CB_CG_OD2_angle, CA_CB_CG_OD2_diangle);
//    OD2= Atom("OD2", oxygen_d2, 0.0, 1.0, " ", " OD2", 0, "O")
    resx.Rotcd.push_back(oxygen_d2);
    resx.Rotnm.push_back("OD2");  

    resx.Ocd=calculateCoordinates(N, CA, C, C_O_length, CA_C_O_angle, N_CA_C_O_diangle);

    return resx;    
}


Rotm makeAsn(vector<float > N, vector<float > CA, vector<float > C)
{
 //   '''Creates an Asparagine residue'''
//    ##R-Group
    float C_O_length=1.23;
    float CA_C_O_angle=120.4826;
    float N_CA_C_O_diangle=-60.0;

    float CA_CB_length=1.52;
    float C_CA_CB_angle=109.5;
    float N_C_CA_CB_diangle=123.2254;
    
    float CB_CG_length=1.52;
    float CA_CB_CG_angle=112.62;
    float N_CA_CB_CG_diangle=-65.5;
    
    float CG_OD1_length=1.23;
    float CB_CG_OD1_angle=120.85;
    float CA_CB_CG_OD1_diangle=-58.3;
    
    float CG_ND2_length=1.33;
    float CB_CG_ND2_angle=116.48;
    float CA_CB_CG_ND2_diangle=121.7;

    Rotm resx;

    vector<float > carbon_b= calculateCoordinates(N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle);
//    CB= Atom("CB", carbon_b, 0.0 , 1.0, " "," CB", 0,"C")
    resx.Rotcd.push_back(carbon_b);
    resx.Rotnm.push_back("CB");  
    vector<float> CB = carbon_b;
    vector<float > carbon_g= calculateCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle);
//    CG= Atom("CG", carbon_g, 0.0, 1.0, " ", " CG", 0, "C")
    resx.Rotcd.push_back(carbon_g);
    resx.Rotnm.push_back("CG");
    vector<float> CG = carbon_g;  
    vector<float > oxygen_d1= calculateCoordinates(CA, CB, CG, CG_OD1_length, CB_CG_OD1_angle, CA_CB_CG_OD1_diangle);
//    OD1= Atom("OD1", oxygen_d1, 0.0, 1.0, " ", " OD1", 0, "O")
    resx.Rotcd.push_back(oxygen_d1);
    resx.Rotnm.push_back("OD1");  
    vector<float > nitrogen_d2= calculateCoordinates(CA, CB, CG, CG_ND2_length, CB_CG_ND2_angle, CA_CB_CG_ND2_diangle);
//    ND2= Atom("ND2", nitrogen_d2, 0.0, 1.0, " ", " ND2", 0, "N")
    resx.Rotcd.push_back(nitrogen_d2);
    resx.Rotnm.push_back("ND2");  
//    res= Residue((' ', segID, ' '), "ASN", '    ')

    resx.Ocd=calculateCoordinates(N, CA, C, C_O_length, CA_C_O_angle, N_CA_C_O_diangle);

    return resx;    
}


Rotm makeGlu(vector<float > N, vector<float > CA, vector<float > C)
{
//    '''Creates a Glutamic Acid residue'''
//    ##R-Group
    float C_O_length=1.23;
    float CA_C_O_angle=120.511;
    float N_CA_C_O_diangle=120.0;

    float CA_CB_length=1.52;
    float C_CA_CB_angle =109.5;
    float N_C_CA_CB_diangle=122.8702;
    
    float CB_CG_length=1.52;
    float CA_CB_CG_angle=113.82;
    float N_CA_CB_CG_diangle=-63.8;

    float CG_CD_length=1.52;
    float CB_CG_CD_angle=113.31;
    float CA_CB_CG_CD_diangle=-179.8;

    float CD_OE1_length=1.25;
    float CG_CD_OE1_angle=119.02;
    float CB_CG_CD_OE1_diangle=-6.2;

    float CD_OE2_length=1.25;
    float CG_CD_OE2_angle=118.08;
    float CB_CG_CD_OE2_diangle=173.8;

    Rotm resx;

    vector<float > carbon_b= calculateCoordinates(N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle);
//    CB= Atom("CB", carbon_b, 0.0 , 1.0, " "," CB", 0,"C")
    resx.Rotcd.push_back(carbon_b);
    resx.Rotnm.push_back("CB"); 
    vector<float> CB = carbon_b;
    vector<float > carbon_g= calculateCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle);
//    CG= Atom("CG", carbon_g, 0.0, 1.0, " ", " CG", 0, "C")
    resx.Rotcd.push_back(carbon_g);
    resx.Rotnm.push_back("CG");
    vector<float> CG = carbon_g; 
    vector<float > carbon_d= calculateCoordinates(CA, CB, CG, CG_CD_length, CB_CG_CD_angle, CA_CB_CG_CD_diangle);
//    CD= Atom("CD", carbon_d, 0.0, 1.0, " ", " CD", 0, "C")
    resx.Rotcd.push_back(carbon_d);
    resx.Rotnm.push_back("CD"); 
    vector<float> CD = carbon_d;
    vector<float > oxygen_e1= calculateCoordinates(CB, CG, CD, CD_OE1_length, CG_CD_OE1_angle, CB_CG_CD_OE1_diangle);
//    OE1= Atom("OE1", oxygen_e1, 0.0, 1.0, " ", " OE1", 0, "O")
    resx.Rotcd.push_back(oxygen_e1);
    resx.Rotnm.push_back("OE1"); 
    vector<float > oxygen_e2= calculateCoordinates(CB, CG, CD, CD_OE2_length, CG_CD_OE2_angle, CB_CG_CD_OE2_diangle);
//    OE2= Atom("OE2", oxygen_e2, 0.0, 1.0, " ", " OE2", 0, "O")
    resx.Rotcd.push_back(oxygen_e2);
    resx.Rotnm.push_back("OE2"); 

    resx.Ocd=calculateCoordinates(N, CA, C, C_O_length, CA_C_O_angle, N_CA_C_O_diangle);

    return resx;    
}


Rotm makeGln(vector<float > N, vector<float > CA, vector<float > C)
{
//    '''Creates a Glutamine residue'''
//    ##R-Group
    float C_O_length=1.23;
    float CA_C_O_angle=120.5029;
    float N_CA_C_O_diangle=120.0;

    float CA_CB_length=1.52;
    float C_CA_CB_angle=109.5;
    float N_C_CA_CB_diangle=122.8134;
    
    float CB_CG_length=1.52;
    float CA_CB_CG_angle=113.75;
    float N_CA_CB_CG_diangle=-60.2;

    float CG_CD_length=1.52;
    float CB_CG_CD_angle=112.78;
    float CA_CB_CG_CD_diangle=-69.6;
    
    float CD_OE1_length=1.24;
    float CG_CD_OE1_angle=120.86;
    float CB_CG_CD_OE1_diangle=-50.5;
    
    float CD_NE2_length=1.33;
    float CG_CD_NE2_angle=116.50;
    float CB_CG_CD_NE2_diangle=129.5;

    Rotm resx;

    vector<float > carbon_b= calculateCoordinates(N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle);
//    CB= Atom("CB", carbon_b, 0.0 , 1.0, " "," CB", 0,"C")
    resx.Rotcd.push_back(carbon_b);
    resx.Rotnm.push_back("CB"); 
    vector<float> CB = carbon_b;
    vector<float > carbon_g= calculateCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle);
//    CG= Atom("CG", carbon_g, 0.0, 1.0, " ", " CG", 0, "C")
    resx.Rotcd.push_back(carbon_g);
    resx.Rotnm.push_back("CG"); 
    vector<float> CG = carbon_g;
    vector<float > carbon_d= calculateCoordinates(CA, CB, CG, CG_CD_length, CB_CG_CD_angle, CA_CB_CG_CD_diangle);
//    CD= Atom("CD", carbon_d, 0.0, 1.0, " ", " CD", 0, "C")
    resx.Rotcd.push_back(carbon_d);
    resx.Rotnm.push_back("CD"); 
    vector<float> CD = carbon_d;
    vector<float > oxygen_e1= calculateCoordinates(CB, CG, CD, CD_OE1_length, CG_CD_OE1_angle, CB_CG_CD_OE1_diangle);
//    OE1= Atom("OE1", oxygen_e1, 0.0, 1.0, " ", " OE1", 0, "O")
    resx.Rotcd.push_back(oxygen_e1);
    resx.Rotnm.push_back("OE1"); 
    vector<float > nitrogen_e2= calculateCoordinates(CB, CG, CD, CD_NE2_length, CG_CD_NE2_angle, CB_CG_CD_NE2_diangle);
 //   NE2= Atom("NE2", nitrogen_e2, 0.0, 1.0, " ", " NE2", 0, "N")
    resx.Rotcd.push_back(nitrogen_e2);
    resx.Rotnm.push_back("NE2"); 

    resx.Ocd=calculateCoordinates(N, CA, C, C_O_length, CA_C_O_angle, N_CA_C_O_diangle);

    return resx;
}

Rotm makeMet(vector<float > N, vector<float > CA, vector<float > C)
{
//    '''Creates a Methionine residue'''
//    ##R-Group
    float C_O_length=1.23;
    float CA_C_O_angle=120.4816;
    float N_CA_C_O_diangle=120.0;

    float CA_CB_length=1.52;
    float C_CA_CB_angle=109.5;
    float N_C_CA_CB_diangle=122.6733;

    float CB_CG_length=1.52;
    float CA_CB_CG_angle=113.68;
    float N_CA_CB_CG_diangle=-64.4;
    
    float CG_SD_length=1.81;
    float CB_CG_SD_angle=112.69;
    float CA_CB_CG_SD_diangle=-179.6;
    
    float SD_CE_length=1.79;
    float CG_SD_CE_angle=100.61;
    float CB_CG_SD_CE_diangle=70.1;

    Rotm resx;
    
    vector<float > carbon_b= calculateCoordinates(N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle);
//    CB= Atom("CB", carbon_b, 0.0 , 1.0, " "," CB", 0,"C")
    resx.Rotcd.push_back(carbon_b);
    resx.Rotnm.push_back("CB");
    vector<float> CB = carbon_b; 
    vector<float > carbon_g= calculateCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle);
//    CG= Atom("CG", carbon_g, 0.0, 1.0, " ", " CG", 0, "C")
    resx.Rotcd.push_back(carbon_g);
    resx.Rotnm.push_back("CG"); 
    vector<float> CG = carbon_g;
    vector<float > sulfur_d= calculateCoordinates(CA, CB, CG, CG_SD_length, CB_CG_SD_angle, CA_CB_CG_SD_diangle);
//    SD= Atom("SD", sulfur_d, 0.0, 1.0, " ", " SD", 0, "S")
    resx.Rotcd.push_back(sulfur_d);
    resx.Rotnm.push_back("SD"); 
    vector<float> SD = sulfur_d;
    vector<float > carbon_e= calculateCoordinates(CB, CG, SD, SD_CE_length, CG_SD_CE_angle, CB_CG_SD_CE_diangle);
//    CE= Atom("CE", carbon_e, 0.0, 1.0, " ", " CE", 0, "C")
    resx.Rotcd.push_back(carbon_e);
    resx.Rotnm.push_back("CE"); 

    resx.Ocd=calculateCoordinates(N, CA, C, C_O_length, CA_C_O_angle, N_CA_C_O_diangle);

    return resx;
}


Rotm makeHis(vector<float > N, vector<float> CA, vector<float> C)
{
//    '''Creates a Histidine residue'''
//    ##R-Group
    float C_O_length=1.23;
    float CA_C_O_angle=120.4732;
    float N_CA_C_O_diangle=120.0;

    float CA_CB_length=1.52;
    float C_CA_CB_angle=109.5;
    float N_C_CA_CB_diangle=122.6711;

    float CB_CG_length=1.49;
    float CA_CB_CG_angle=113.74;
    float N_CA_CB_CG_diangle=-63.2;
    
    float CG_ND1_length=1.38;
    float CB_CG_ND1_angle=122.85;
    float CA_CB_CG_ND1_diangle=-75.7;
    
    float CG_CD2_length=1.35;
    float CB_CG_CD2_angle=130.61;
    float CA_CB_CG_CD2_diangle=104.3;
    
    float ND1_CE1_length=1.32;
    float CG_ND1_CE1_angle=108.5;
    float CB_CG_ND1_CE1_diangle=180.0;
    
    float CD2_NE2_length=1.35;
    float CG_CD2_NE2_angle=108.5;
    float CB_CG_CD2_NE2_diangle=180.0;

    Rotm resx;
    
    vector<float > carbon_b= calculateCoordinates(N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle);
//    CB= Atom("CB", carbon_b, 0.0 , 1.0, " "," CB", 0,"C")
    resx.Rotcd.push_back(carbon_b);
    resx.Rotnm.push_back("CB"); 
    vector<float> CB = carbon_b;
    vector<float > carbon_g= calculateCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle);
//    CG= Atom("CG", carbon_g, 0.0, 1.0, " ", " CG", 0, "C")
    resx.Rotcd.push_back(carbon_g);
    resx.Rotnm.push_back("CG"); 
    vector<float> CG = carbon_g;
    vector<float > nitrogen_d1= calculateCoordinates(CA, CB, CG, CG_ND1_length, CB_CG_ND1_angle, CA_CB_CG_ND1_diangle);
//    ND1= Atom("ND1", nitrogen_d1, 0.0, 1.0, " ", " ND1", 0, "N")
    resx.Rotcd.push_back(nitrogen_d1);
    resx.Rotnm.push_back("ND1"); 
    vector<float > carbon_d2= calculateCoordinates(CA, CB, CG, CG_CD2_length, CB_CG_CD2_angle, CA_CB_CG_CD2_diangle);
//    CD2= Atom("CD2", carbon_d2, 0.0, 1.0, " ", " CD2", 0, "C")
    resx.Rotcd.push_back(carbon_d2);
    resx.Rotnm.push_back("CD2"); 
    vector<float> ND1 = nitrogen_d1;
    vector<float > carbon_e1= calculateCoordinates(CB, CG, ND1, ND1_CE1_length, CG_ND1_CE1_angle, CB_CG_ND1_CE1_diangle);
 //   CE1= Atom("CE1", carbon_e1, 0.0, 1.0, " ", " CE1", 0, "C")
    resx.Rotcd.push_back(carbon_e1);
    resx.Rotnm.push_back("CE1"); 
    vector<float> CD2 = carbon_d2;
    vector<float > nitrogen_e2= calculateCoordinates(CB, CG, CD2, CD2_NE2_length, CG_CD2_NE2_angle, CB_CG_CD2_NE2_diangle);
//    NE2= Atom("NE2", nitrogen_e2, 0.0, 1.0, " ", " NE2", 0, "N")
    resx.Rotcd.push_back(nitrogen_e2);
    resx.Rotnm.push_back("NE2"); 

    resx.Ocd=calculateCoordinates(N, CA, C, C_O_length, CA_C_O_angle, N_CA_C_O_diangle);

    return resx;
}

Rotm makePro(vector<float> N, vector<float> CA, vector<float> C)
{
//    '''Creates a Proline residue'''
//    ##R-Group
    float C_O_length=1.23;
    float CA_C_O_angle=120.2945;
    float N_CA_C_O_diangle= -45.0;

    float CA_CB_length=1.52;
    float C_CA_CB_angle=109.5;
    float N_C_CA_CB_diangle=115.2975;
    
    float CB_CG_length=1.49;
    float CA_CB_CG_angle=104.21;
    float N_CA_CB_CG_diangle=29.6;
    
    float CG_CD_length=1.50;
    float CB_CG_CD_angle=105.03;
    float CA_CB_CG_CD_diangle=-34.8;

    Rotm resx;
    
    vector<float > carbon_b= calculateCoordinates(N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle);
//    CB= Atom("CB", carbon_b, 0.0 , 1.0, " "," CB", 0,"C")
    resx.Rotcd.push_back(carbon_b);
    resx.Rotnm.push_back("CB"); 
    vector<float> CB = carbon_b;
    vector<float > carbon_g= calculateCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle);
//    CG= Atom("CG", carbon_g, 0.0, 1.0, " ", " CG", 0, "C")
    resx.Rotcd.push_back(carbon_g);
    resx.Rotnm.push_back("CG"); 
    vector<float> CG = carbon_g;
    vector<float > carbon_d= calculateCoordinates(CA, CB, CG, CG_CD_length, CB_CG_CD_angle, CA_CB_CG_CD_diangle);
//    CD= Atom("CD", carbon_d, 0.0, 1.0, " ", " CD", 0, "C")
    resx.Rotcd.push_back(carbon_d);
    resx.Rotnm.push_back("CD"); 

    resx.Ocd=calculateCoordinates(N, CA, C, C_O_length, CA_C_O_angle, N_CA_C_O_diangle);

    return resx;
}

Rotm makePhe(vector<float> N, vector<float> CA, vector<float> C)
{
//    '''Creates a Phenylalanine residue'''
//    ##R-Group
    float C_O_length=1.23;
    float CA_C_O_angle=120.5316;
    float N_CA_C_O_diangle=120.0;

    float CA_CB_length=1.52;
    float C_CA_CB_angle=109.5;
    float N_C_CA_CB_diangle=122.6054;
    
    float CB_CG_length=1.50;
    float CA_CB_CG_angle=113.85;
    float N_CA_CB_CG_diangle=-64.7;

    float CG_CD1_length=1.39;
    float CB_CG_CD1_angle=120.0;
    float CA_CB_CG_CD1_diangle=93.3;

    float CG_CD2_length=1.39;
    float CB_CG_CD2_angle=120.0;
    float CA_CB_CG_CD2_diangle= -86.7;
    
    float CD1_CE1_length=1.39;
    float CG_CD1_CE1_angle=120.0;
    float CB_CG_CD1_CE1_diangle=180.0;

    float CD2_CE2_length=1.39;
    float CG_CD2_CE2_angle=120.0;
    float CB_CG_CD2_CE2_diangle=180.0;

    float CE1_CZ_length=1.39;
    float CD1_CE1_CZ_angle=120.0;
    float CG_CD1_CE1_CZ_diangle=0.0;

    Rotm resx;

    vector<float > carbon_b= calculateCoordinates(N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle);
//    CB= Atom("CB", carbon_b, 0.0 , 1.0, " "," CB", 0,"C")
    resx.Rotcd.push_back(carbon_b);
    resx.Rotnm.push_back("CB"); 
    vector<float> CB = carbon_b;
    vector<float > carbon_g= calculateCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle);
//    CG= Atom("CG", carbon_g, 0.0, 1.0, " ", " CG", 0, "C")
    resx.Rotcd.push_back(carbon_g);
    resx.Rotnm.push_back("CG");
    vector<float> CG = carbon_g; 
    vector<float > carbon_d1= calculateCoordinates(CA, CB, CG, CG_CD1_length, CB_CG_CD1_angle, CA_CB_CG_CD1_diangle);
//    CD1= Atom("CD1", carbon_d1, 0.0, 1.0, " ", " CD1", 0, "C")
    resx.Rotcd.push_back(carbon_d1);
    resx.Rotnm.push_back("CD1"); 
    vector<float > carbon_d2= calculateCoordinates(CA, CB, CG, CG_CD2_length, CB_CG_CD2_angle, CA_CB_CG_CD2_diangle);
//    CD2= Atom("CD2", carbon_d2, 0.0, 1.0, " ", " CD2", 0, "C")
    resx.Rotcd.push_back(carbon_d2);
    resx.Rotnm.push_back("CD2"); 
    vector<float> CD1 = carbon_d1;
    vector<float > carbon_e1= calculateCoordinates(CB, CG, CD1, CD1_CE1_length, CG_CD1_CE1_angle, CB_CG_CD1_CE1_diangle);
//    CE1= Atom("CE1", carbon_e1, 0.0, 1.0, " ", " CE1", 0, "C")
    resx.Rotcd.push_back(carbon_e1);
    resx.Rotnm.push_back("CE1"); 
    vector<float> CD2 = carbon_d2;
    vector<float > carbon_e2= calculateCoordinates(CB, CG, CD2, CD2_CE2_length, CG_CD2_CE2_angle, CB_CG_CD2_CE2_diangle);
//    CE2= Atom("CE2", carbon_e2, 0.0, 1.0, " ", " CE2", 0, "C")
    resx.Rotcd.push_back(carbon_e2);
    resx.Rotnm.push_back("CE2"); 
    vector<float> CE1 = carbon_e1;
    vector<float > carbon_z= calculateCoordinates(CG, CD1, CE1, CE1_CZ_length, CD1_CE1_CZ_angle, CG_CD1_CE1_CZ_diangle);
//    CZ= Atom("CZ", carbon_z, 0.0, 1.0, " ", " CZ", 0, "C")
    resx.Rotcd.push_back(carbon_z);
    resx.Rotnm.push_back("CZ"); 

    resx.Ocd=calculateCoordinates(N, CA, C, C_O_length, CA_C_O_angle, N_CA_C_O_diangle);

    return resx;

}

Rotm makeTyr(vector<float> N, vector<float> CA, vector<float> C)
{
//    '''Creates a Tyrosine residue'''
//    ##R-Group
    float C_O_length=1.23;
    float CA_C_O_angle=120.5434;
    float N_CA_C_O_diangle=120.0;

    float CA_CB_length=1.52;
    float C_CA_CB_angle=109.5;
    float N_C_CA_CB_diangle=122.6023;
    
    float CB_CG_length=1.51;
    float CA_CB_CG_angle=113.8;
    float N_CA_CB_CG_diangle=-64.3;

    float CG_CD1_length=1.39;
    float CB_CG_CD1_angle=120.98;
    float CA_CB_CG_CD1_diangle=93.1;
    
    float CG_CD2_length=1.39;
    float CB_CG_CD2_angle=120.82;
    float CA_CB_CG_CD2_diangle=273.1; // 86.9
    
    float CD1_CE1_length=1.39;
    float CG_CD1_CE1_angle=120.0;
    float CB_CG_CD1_CE1_diangle=180.0;

    float CD2_CE2_length=1.39;
    float CG_CD2_CE2_angle=120.0;
    float CB_CG_CD2_CE2_diangle=180.0;

    float CE1_CZ_length=1.39;
    float CD1_CE1_CZ_angle=120.0;
    float CG_CD1_CE1_CZ_diangle=0.0;

    float CZ_OH_length=1.39;
    float CE1_CZ_OH_angle=119.78;
    float CD1_CE1_CZ_OH_diangle=180.0;

    Rotm resx;

    vector<float > carbon_b= calculateCoordinates(N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle);
//    CB= Atom("CB", carbon_b, 0.0 , 1.0, " "," CB", 0,"C")
    resx.Rotcd.push_back(carbon_b);
    resx.Rotnm.push_back("CB"); 
    vector<float> CB = carbon_b;
    vector<float > carbon_g= calculateCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle);
//    CG= Atom("CG", carbon_g, 0.0, 1.0, " ", " CG", 0, "C")
    resx.Rotcd.push_back(carbon_g);
    resx.Rotnm.push_back("CG"); 
    vector<float> CG = carbon_g;
    vector<float > carbon_d1= calculateCoordinates(CA, CB, CG, CG_CD1_length, CB_CG_CD1_angle, CA_CB_CG_CD1_diangle);
//    CD1= Atom("CD1", carbon_d1, 0.0, 1.0, " ", " CD1", 0, "C")
    resx.Rotcd.push_back(carbon_d1);
    resx.Rotnm.push_back("CD1"); 
    vector<float > carbon_d2= calculateCoordinates(CA, CB, CG, CG_CD2_length, CB_CG_CD2_angle, CA_CB_CG_CD2_diangle);
//    CD2= Atom("CD2", carbon_d2, 0.0, 1.0, " ", " CD2", 0, "C")
    resx.Rotcd.push_back(carbon_d2);
    resx.Rotnm.push_back("CD2"); 
    vector<float> CD1 = carbon_d1;
    vector<float > carbon_e1= calculateCoordinates(CB, CG, CD1, CD1_CE1_length, CG_CD1_CE1_angle, CB_CG_CD1_CE1_diangle);
//    CE1= Atom("CE1", carbon_e1, 0.0, 1.0, " ", " CE1", 0, "C")
    resx.Rotcd.push_back(carbon_e1);
    resx.Rotnm.push_back("CE1"); 
    vector<float> CD2 = carbon_d2;
    vector<float > carbon_e2= calculateCoordinates(CB, CG, CD2, CD2_CE2_length, CG_CD2_CE2_angle, CB_CG_CD2_CE2_diangle);
//    CE2= Atom("CE2", carbon_e2, 0.0, 1.0, " ", " CE2", 0, "C")
    resx.Rotcd.push_back(carbon_e2);
    resx.Rotnm.push_back("CE2"); 
    vector<float> CE1 = carbon_e1;
    vector<float > carbon_z= calculateCoordinates(CG, CD1, CE1, CE1_CZ_length, CD1_CE1_CZ_angle, CG_CD1_CE1_CZ_diangle);
//    CZ= Atom("CZ", carbon_z, 0.0, 1.0, " ", " CZ", 0, "C")
    resx.Rotcd.push_back(carbon_z);
    resx.Rotnm.push_back("CZ"); 
    vector<float> CZ = carbon_z;
    vector<float > oxygen_h= calculateCoordinates(CD1, CE1, CZ, CZ_OH_length, CE1_CZ_OH_angle, CD1_CE1_CZ_OH_diangle);
//    OH= Atom("OH", oxygen_h, 0.0, 1.0, " ", " OH", 0, "O")
    resx.Rotcd.push_back(oxygen_h);
    resx.Rotnm.push_back("OH"); 

    resx.Ocd=calculateCoordinates(N, CA, C, C_O_length, CA_C_O_angle, N_CA_C_O_diangle);

    return resx;
}

Rotm makeTrp(vector<float> N, vector<float> CA, vector<float> C)
{
//    '''Creates a Tryptophan residue'''
 //   ##R-Group
    float C_O_length=1.23;
    float CA_C_O_angle=120.5117;
    float N_CA_C_O_diangle=120.0;

    float CA_CB_length=1.52;
    float C_CA_CB_angle=109.5;
    float N_C_CA_CB_diangle=122.611;

    float CB_CG_length=1.50;
    float CA_CB_CG_angle=114.10;
    float N_CA_CB_CG_diangle=-66.4;

    float CG_CD1_length=1.37;
    float CB_CG_CD1_angle=127.07;
    float CA_CB_CG_CD1_diangle=96.3;

    float CG_CD2_length=1.43;
    float CB_CG_CD2_angle=126.66;
    float CA_CB_CG_CD2_diangle=-83.7;
    
    float CD1_NE1_length=1.38;
    float CG_CD1_NE1_angle=108.5;
    float CB_CG_CD1_NE1_diangle=180.0;

    float CD2_CE2_length=1.40;
    float CG_CD2_CE2_angle=108.5;
    float CB_CG_CD2_CE2_diangle=180.0;

    float CD2_CE3_length=1.40;
    float CG_CD2_CE3_angle=133.83;
    float CB_CG_CD2_CE3_diangle=0.0;

    float CE2_CZ2_length=1.40;
    float CD2_CE2_CZ2_angle=120.0;
    float CG_CD2_CE2_CZ2_diangle=180.0;

    float CE3_CZ3_length=1.40;
    float CD2_CE3_CZ3_angle=120.0;
    float CG_CD2_CE3_CZ3_diangle=180.0;

    float CZ2_CH2_length=1.40;
    float CE2_CZ2_CH2_angle=120.0;
    float CD2_CE2_CZ2_CH2_diangle=0.0;

    Rotm resx;

    vector<float > carbon_b= calculateCoordinates(N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle);
//    CB= Atom("CB", carbon_b, 0.0 , 1.0, " "," CB", 0,"C")
    resx.Rotcd.push_back(carbon_b);
    resx.Rotnm.push_back("CB"); 
    vector<float> CB = carbon_b;
    vector<float > carbon_g= calculateCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle);
//    CG= Atom("CG", carbon_g, 0.0, 1.0, " ", " CG", 0, "C")
    resx.Rotcd.push_back(carbon_g);
    resx.Rotnm.push_back("CG"); 
    vector<float> CG = carbon_g;
    vector<float > carbon_d1= calculateCoordinates(CA, CB, CG, CG_CD1_length, CB_CG_CD1_angle, CA_CB_CG_CD1_diangle);
//    CD1= Atom("CD1", carbon_d1, 0.0, 1.0, " ", " CD1", 0, "C")
    resx.Rotcd.push_back(carbon_d1);
    resx.Rotnm.push_back("CD1"); 
    vector<float > carbon_d2= calculateCoordinates(CA, CB, CG, CG_CD2_length, CB_CG_CD2_angle, CA_CB_CG_CD2_diangle);
//    CD2= Atom("CD2", carbon_d2, 0.0, 1.0, " ", " CD2", 0, "C")
    resx.Rotcd.push_back(carbon_d2);
    resx.Rotnm.push_back("CD2"); 
    vector<float> CD1 = carbon_d1;
    vector<float > nitrogen_e1= calculateCoordinates(CB, CG, CD1, CD1_NE1_length, CG_CD1_NE1_angle, CB_CG_CD1_NE1_diangle);
//    NE1= Atom("NE1", nitrogen_e1, 0.0, 1.0, " ", " NE1", 0, "N")
    resx.Rotcd.push_back(nitrogen_e1);
    resx.Rotnm.push_back("NE1"); 
    vector<float> CD2 = carbon_d2;
    vector<float > carbon_e2= calculateCoordinates(CB, CG, CD2, CD2_CE2_length, CG_CD2_CE2_angle, CB_CG_CD2_CE2_diangle);
//    CE2= Atom("CE2", carbon_e2, 0.0, 1.0, " ", " CE2", 0, "C")
    resx.Rotcd.push_back(carbon_e2);
    resx.Rotnm.push_back("CE2"); 
    vector<float > carbon_e3= calculateCoordinates(CB, CG, CD2, CD2_CE3_length, CG_CD2_CE3_angle, CB_CG_CD2_CE3_diangle);
//    CE3= Atom("CE3", carbon_e3, 0.0, 1.0, " ", " CE3", 0, "C")
    resx.Rotcd.push_back(carbon_e3);
    resx.Rotnm.push_back("CE3"); 
    vector<float> CE2 = carbon_e2;
    vector<float > carbon_z2= calculateCoordinates(CG, CD2, CE2, CE2_CZ2_length, CD2_CE2_CZ2_angle, CG_CD2_CE2_CZ2_diangle);
//    CZ2= Atom("CZ2", carbon_z2, 0.0, 1.0, " ", " CZ2", 0, "C")
    resx.Rotcd.push_back(carbon_z2);
    resx.Rotnm.push_back("CZ2"); 
    vector<float> CE3 = carbon_e3;
    vector<float > carbon_z3= calculateCoordinates(CG, CD2, CE3, CE3_CZ3_length, CD2_CE3_CZ3_angle, CG_CD2_CE3_CZ3_diangle);
//    CZ3= Atom("CZ3", carbon_z3, 0.0, 1.0, " ", " CZ3", 0, "C")
    resx.Rotcd.push_back(carbon_z3);
    resx.Rotnm.push_back("CZ3"); 
    vector<float> CZ2 = carbon_z2;
    vector<float > carbon_h2= calculateCoordinates(CD2, CE2, CZ2, CZ2_CH2_length, CE2_CZ2_CH2_angle, CD2_CE2_CZ2_CH2_diangle);
//    CH2= Atom("CH2", carbon_h2, 0.0, 1.0, " ", " CH2", 0, "C")
    resx.Rotcd.push_back(carbon_h2);
    resx.Rotnm.push_back("CH2"); 

    resx.Ocd=calculateCoordinates(N, CA, C, C_O_length, CA_C_O_angle, N_CA_C_O_diangle);
    
    return resx;
}


Rotm CalculateRotm(string AA,vector<float> N,vector<float> CA,vector<float> C)
{
    Rotm res;
    if(strcmp(AA.c_str(),"G")==0||strcmp(AA.c_str()," GLY")==0)
        res=makeGly(N, CA, C);
    else if(strcmp(AA.c_str(),"A")==0||strcmp(AA.c_str()," ALA")==0)
        res=makeAla(N, CA, C);
    else if(strcmp(AA.c_str(),"S")==0||strcmp(AA.c_str()," SER")==0)
        res=makeSer(N, CA, C);
    else if(strcmp(AA.c_str(),"C")==0||strcmp(AA.c_str()," CYS")==0)
        res=makeCys(N, CA, C);
    else if(strcmp(AA.c_str(),"V")==0||strcmp(AA.c_str()," VAL")==0)
        res=makeVal(N, CA, C);
    else if(strcmp(AA.c_str(),"I")==0||strcmp(AA.c_str()," ILE")==0)
        res=makeIle(N, CA, C);
    else if(strcmp(AA.c_str(),"L")==0||strcmp(AA.c_str()," LEU")==0)
        res=makeLeu(N, CA, C);
    else if(strcmp(AA.c_str(),"T")==0||strcmp(AA.c_str()," THR")==0)
        res=makeThr(N, CA, C);
    else if(strcmp(AA.c_str(),"R")==0||strcmp(AA.c_str()," ARG")==0)
        res=makeArg(N, CA, C);
    else if(strcmp(AA.c_str(),"K")==0||strcmp(AA.c_str()," LYS")==0)
        res=makeLys(N, CA, C);
    else if(strcmp(AA.c_str(),"D")==0||strcmp(AA.c_str()," ASP")==0)
        res=makeAsp(N, CA, C);
    else if(strcmp(AA.c_str(),"E")==0||strcmp(AA.c_str()," GLU")==0)
        res=makeGlu(N, CA, C);
    else if(strcmp(AA.c_str(),"N")==0||strcmp(AA.c_str()," ASN")==0)
        res=makeAsn(N, CA, C);
    else if(strcmp(AA.c_str(),"Q")==0||strcmp(AA.c_str()," GLN")==0)
        res=makeGln(N, CA, C);
    else if(strcmp(AA.c_str(),"M")==0||strcmp(AA.c_str()," MET")==0)
        res=makeMet(N, CA, C);
    else if(strcmp(AA.c_str(),"H")==0||strcmp(AA.c_str()," HIS")==0)
        res=makeHis(N, CA, C);
    else if(strcmp(AA.c_str(),"P")==0||strcmp(AA.c_str()," PRO")==0)
        res=makePro(N, CA, C);
    else if(strcmp(AA.c_str(),"F")==0||strcmp(AA.c_str()," PHE")==0)
        res=makePhe(N, CA, C);
    else if(strcmp(AA.c_str(),"Y")==0||strcmp(AA.c_str()," TYR")==0)
        res=makeTyr(N, CA, C);
    else if(strcmp(AA.c_str(),"W")==0||strcmp(AA.c_str()," TRP")==0)
        res=makeTrp(N, CA, C);
    else
        res=makeGly(N, CA, C);    
/*    if(AA.c_str()=="G"||AA.c_str()==" GLY")
        res=makeGly(N, CA, C);
    else if(AA.c_str()=="A"||AA.c_str()==" ALA")
        res=makeAla(N, CA, C);
    else if(AA.c_str()=="S"||AA.c_str()==" SER")
        res=makeSer(N, CA, C);
    else if(AA.c_str()=="C"||AA.c_str()==" CYS")
        res=makeCys(N, CA, C);
    else if(AA.c_str()=="V"||AA.c_str()==" VAL")
        res=makeVal(N, CA, C);
    else if(AA.c_str()=="I"||AA.c_str()==" ILE")
        res=makeIle(N, CA, C);
    else if(AA.c_str()=="L"||AA.c_str()==" LEU")
        res=makeLeu(N, CA, C);
    else if(AA.c_str()=="T"||AA.c_str()==" THR")
        res=makeThr(N, CA, C);
    else if(AA.c_str()=="R"||AA.c_str()==" ARG")
        res=makeArg(N, CA, C);
    else if(AA.c_str()=="K"||AA.c_str()==" LYS")
        res=makeLys(N, CA, C);
    else if(AA.c_str()=="D"||AA.c_str()==" ASP")
        res=makeAsp(N, CA, C);
    else if(AA.c_str()=="E"||AA.c_str()==" GLU")
        res=makeGlu(N, CA, C);
    else if(AA.c_str()=="N"||AA.c_str()==" ASN")
        res=makeAsn(N, CA, C);
    else if(AA.c_str()=="Q"||AA.c_str()==" GLN")
        res=makeGln(N, CA, C);
    else if(AA.c_str()=="M"||AA.c_str()==" MET")
        res=makeMet(N, CA, C);
    else if(AA.c_str()=="H"||AA.c_str()==" HIS")
        res=makeHis(N, CA, C);
    else if(AA.c_str()=="P"||AA.c_str()==" PRO")
        res=makePro(N, CA, C);
    else if(AA.c_str()=="F"||AA.c_str()==" PHE")
        res=makePhe(N, CA, C);
    else if(AA.c_str()=="Y"||AA.c_str()==" TYR")
        res=makeTyr(N, CA, C);
    else if(AA.c_str()=="W"||AA.c_str()==" TRP")
        res=makeTrp(N, CA, C);
    else
        res=makeGly(N, CA, C);  */

    return res;
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

