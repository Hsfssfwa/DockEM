#include "GetVdwRadius.h"
#include "Kabsch.h"
#include "GeometryTools.h"

using namespace std;

const float PI=3.141592653589793;
const float Extra=1.0e-4;
const float UpMax=1.0e+10;

#define epsilon  1e-7

static float expval[]={
1.0000000,0.9900498,0.9801987,0.9704455,0.9607894,0.9512294,0.9417645,0.9323938,0.9231163,0.9139312,
0.9048374,0.8958341,0.8869204,0.8780954,0.8693582,0.8607080,0.8521438,0.8436648,0.8352702,0.8269591,
0.8187308,0.8105842,0.8025188,0.7945336,0.7866279,0.7788008,0.7710516,0.7633795,0.7557837,0.7482636,
0.7408182,0.7334470,0.7261490,0.7189237,0.7117703,0.7046881,0.6976763,0.6907343,0.6838614,0.6770569,
0.6703200,0.6636503,0.6570468,0.6505091,0.6440364,0.6376282,0.6312836,0.6250023,0.6187834,0.6126264,
0.6065307,0.6004956,0.5945205,0.5886050,0.5827483,0.5769498,0.5712091,0.5655254,0.5598984,0.5543273,
0.5488116,0.5433509,0.5379444,0.5325918,0.5272924,0.5220458,0.5168513,0.5117086,0.5066170,0.5015761,
0.4965853,0.4916442,0.4867523,0.4819090,0.4771139,0.4723666,0.4676664,0.4630131,0.4584060,0.4538448,
0.4493290,0.4448581,0.4404317,0.4360493,0.4317105,0.4274149,0.4231621,0.4189515,0.4147829,0.4106558,
0.4065697,0.4025242,0.3985190,0.3945537,0.3906278,0.3867410,0.3828929,0.3790830,0.3753111,0.3715767,
0.3678794,0.3642190,0.3605949,0.3570070,0.3534547,0.3499377,0.3464558,0.3430085,0.3395955,0.3362165,
0.3328711,0.3295590,0.3262798,0.3230333,0.3198190,0.3166368,0.3134862,0.3103669,0.3072787,0.3042213,
0.3011942,0.2981973,0.2952302,0.2922926,0.2893842,0.2865048,0.2836540,0.2808316,0.2780373,0.2752708,
0.2725318,0.2698201,0.2671353,0.2644773,0.2618457,0.2592403,0.2566608,0.2541070,0.2515786,0.2490753,
0.2465970,0.2441433,0.2417140,0.2393089,0.2369278,0.2345703,0.2322363,0.2299255,0.2276377,0.2253727,
0.2231302,0.2209100,0.2187119,0.2165357,0.2143811,0.2122480,0.2101361,0.2080452,0.2059751,0.2039256,
0.2018965,0.1998876,0.1978987,0.1959296,0.1939800,0.1920499,0.1901390,0.1882471,0.1863740,0.1845195,
0.1826835,0.1808658,0.1790661,0.1772844,0.1755204,0.1737739,0.1720449,0.1703330,0.1686381,0.1669602,
0.1652989,0.1636541,0.1620258,0.1604136,0.1588174,0.1572372,0.1556726,0.1541237,0.1525901,0.1510718,
0.1495686,0.1480804,0.1466070,0.1451482,0.1437039,0.1422741,0.1408584,0.1394569,0.1380692,0.1366954,
0.1353353,0.1339887,0.1326555,0.1313355,0.1300287,0.1287349,0.1274540,0.1261858,0.1249302,0.1236871,
0.1224564,0.1212380,0.1200316,0.1188373,0.1176548,0.1164842,0.1153251,0.1141776,0.1130415,0.1119167,
0.1108032,0.1097006,0.1086091,0.1075284,0.1064585,0.1053992,0.1043505,0.1033122,0.1022842,0.1012665,
0.1002588,0.0992613,0.0982736,0.0972957,0.0963276,0.0953692,0.0944202,0.0934807,0.0925506,0.0916297,
0.0907180,0.0898153,0.0889216,0.0880368,0.0871609,0.0862936,0.0854350,0.0845849,0.0837432,0.0829100,
0.0820850,0.0812682,0.0804596,0.0796590,0.0788664,0.0780817,0.0773047,0.0765355,0.0757740,0.0750200,
0.0742736,0.0735345,0.0728029,0.0720785,0.0713613,0.0706512,0.0699482,0.0692522,0.0685632,0.0678809,
0.0672055,0.0665368,0.0658748,0.0652193,0.0645703,0.0639279,0.0632918,0.0626620,0.0620385,0.0614212,
0.0608101,0.0602050,0.0596059,0.0590129,0.0584257,0.0578443,0.0572688,0.0566989,0.0561348,0.0555762,
0.0550232,0.0544757,0.0539337,0.0533970,0.0528657,0.0523397,0.0518189,0.0513033,0.0507928,0.0502874,
0.0497871,0.0492917,0.0488012,0.0483156,0.0478349,0.0473589,0.0468877,0.0464212,0.0459593,0.0455020,
0.0450492,0.0446010,0.0441572,0.0437178,0.0432828,0.0428521,0.0424257,0.0420036,0.0415857,0.0411719,
0.0407622,0.0403566,0.0399551,0.0395575,0.0391639,0.0387742,0.0383884,0.0380064,0.0376283,0.0372538,
0.0368832,0.0365162,0.0361528,0.0357931,0.0354370,0.0350844,0.0347353,0.0343896,0.0340475,0.0337087,
0.0333733,0.0330412,0.0327124,0.0323869,0.0320647,0.0317456,0.0314298,0.0311170,0.0308074,0.0305009,
0.0301974,0.0298969,0.0295994,0.0293049,0.0290133,0.0287246,0.0284388,0.0281559,0.0278757,0.0275983,
0.0273237,0.0270518,0.0267827,0.0265162,0.0262523,0.0259911,0.0257325,0.0254765,0.0252230,0.0249720,
0.0247235,0.0244775,0.0242340,0.0239928,0.0237541,0.0235177,0.0232837,0.0230521,0.0228227,0.0225956,
0.0223708,0.0221482,0.0219278,0.0217096,0.0214936,0.0212797,0.0210680,0.0208584,0.0206508,0.0204453,
0.0202419,0.0200405,0.0198411,0.0196437,0.0194482,0.0192547,0.0190631,0.0188734,0.0186856,0.0184997,
0.0183156,0.0181334,0.0179530,0.0177743,0.0175975,0.0174224,0.0172490,0.0170774,0.0169075,0.0167392,
0.0165727,0.0164078,0.0162445,0.0160829,0.0159229,0.0157644,0.0156076,0.0154523,0.0152985,0.0151463,
0.0149956,0.0148464,0.0146986,0.0145524,0.0144076,0.0142642,0.0141223,0.0139818,0.0138427,0.0137049,
0.0135686,0.0134335,0.0132999,0.0131675,0.0130365,0.0129068,0.0127784,0.0126512,0.0125254,0.0124007,
0.0122773,0.0121552,0.0120342,0.0119145,0.0117959,0.0116786,0.0115624,0.0114473,0.0113334,0.0112206,
0.0111090,0.0109985,0.0108890,0.0107807,0.0106734,0.0105672,0.0104621,0.0103580,0.0102549,0.0101529,
0.0100518,0.0099518,0.0098528,0.0097548,0.0096577,0.0095616,0.0094665,0.0093723,0.0092790,0.0091867,
0.0090953,0.0090048,0.0089152,0.0088265,0.0087386,0.0086517,0.0085656,0.0084804,0.0083960,0.0083125,
0.0082297,0.0081479,0.0080668,0.0079865,0.0079071,0.0078284,0.0077505,0.0076734,0.0075970,0.0075214,
0.0074466,0.0073725,0.0072991,0.0072265,0.0071546,0.0070834,0.0070129,0.0069431,0.0068741,0.0068057,
0.0067379,0.0066709,0.0066045,0.0065388,0.0064737,0.0064093,0.0063456,0.0062824,0.0062199,0.0061580,
0.0060967,0.0060361,0.0059760,0.0059166,0.0058577,0.0057994,0.0057417,0.0056846,0.0056280,0.0055720,
0.0055166,0.0054617,0.0054073,0.0053535,0.0053003,0.0052475,0.0051953,0.0051436,0.0050924,0.0050418,
0.0049916,0.0049419,0.0048928,0.0048441,0.0047959,0.0047482,0.0047009,0.0046541,0.0046078,0.0045620,
0.0045166,0.0044716,0.0044271,0.0043831,0.0043395,0.0042963,0.0042536,0.0042112,0.0041693,0.0041278,
0.0040868,0.0040461,0.0040058,0.0039660,0.0039265,0.0038875,0.0038488,0.0038105,0.0037726,0.0037350,
0.0036979,0.0036611,0.0036246,0.0035886,0.0035529,0.0035175,0.0034825,0.0034479,0.0034136,0.0033796,
0.0033460,0.0033127,0.0032797,0.0032471,0.0032148,0.0031828,0.0031511,0.0031198,0.0030887,0.0030580,
0.0030276,0.0029974,0.0029676,0.0029381,0.0029088,0.0028799,0.0028512,0.0028229,0.0027948,0.0027670,
0.0027394,0.0027122,0.0026852,0.0026585,0.0026320,0.0026058,0.0025799,0.0025542,0.0025288,0.0025037,
0.0024788,0.0024541,0.0024297,0.0024055,0.0023816,0.0023579,0.0023344,0.0023112,0.0022882,0.0022654,
0.0022429,0.0022206,0.0021985,0.0021766,0.0021549,0.0021335,0.0021123,0.0020912,0.0020704,0.0020498,
0.0020294,0.0020092,0.0019892,0.0019695,0.0019499,0.0019305,0.0019112,0.0018922,0.0018734,0.0018548,
0.0018363,0.0018180,0.0017999,0.0017820,0.0017643,0.0017467,0.0017294,0.0017122,0.0016951,0.0016783,
0.0016616,0.0016450,0.0016287,0.0016125,0.0015964,0.0015805,0.0015648,0.0015492,0.0015338,0.0015185,
0.0015034,0.0014885,0.0014737,0.0014590,0.0014445,0.0014301,0.0014159,0.0014018,0.0013878,0.0013740,
0.0013604,0.0013468,0.0013334,0.0013202,0.0013070,0.0012940,0.0012811,0.0012684,0.0012558,0.0012433,
0.0012309,0.0012187,0.0012065,0.0011945,0.0011826,0.0011709,0.0011592,0.0011477,0.0011363,0.0011250,
0.0011138,0.0011027,0.0010917,0.0010809,0.0010701,0.0010595,0.0010489,0.0010385,0.0010281,0.0010179,
0.0010078,0.0009978,0.0009878,0.0009780,0.0009683,0.0009586,0.0009491,0.0009397,0.0009303,0.0009210,
0.0009119,0.0009028,0.0008938,0.0008849,0.0008761,0.0008674,0.0008588,0.0008502,0.0008418,0.0008334,
0.0008251,0.0008169,0.0008088,0.0008007,0.0007928,0.0007849,0.0007771,0.0007693,0.0007617,0.0007541,
0.0007466,0.0007392,0.0007318,0.0007245,0.0007173,0.0007102,0.0007031,0.0006961,0.0006892,0.0006823,
0.0006755,0.0006688,0.0006622,0.0006556,0.0006491,0.0006426,0.0006362,0.0006299,0.0006236,0.0006174,
0.0006113,0.0006052,0.0005991,0.0005932,0.0005873,0.0005814,0.0005757,0.0005699,0.0005643,0.0005586,
0.0005531,0.0005476,0.0005421,0.0005367,0.0005314,0.0005261,0.0005209,0.0005157,0.0005106,0.0005055,
0.0005005,0.0004955,0.0004905,0.0004857,0.0004808,0.0004760,0.0004713,0.0004666,0.0004620,0.0004574,
0.0004528,0.0004483,0.0004439,0.0004394,0.0004351,0.0004307,0.0004265,0.0004222,0.0004180,0.0004139,
0.0004097,0.0004057,0.0004016,0.0003976,0.0003937,0.0003898,0.0003859,0.0003820,0.0003782,0.0003745,
0.0003707,0.0003671,0.0003634,0.0003598,0.0003562,0.0003527,0.0003492,0.0003457,0.0003422,0.0003388,
0.0003355,0.0003321,0.0003288,0.0003255,0.0003223,0.0003191,0.0003159,0.0003128,0.0003097,0.0003066,
0.0003035,0.0003005,0.0002975,0.0002946,0.0002916,0.0002887,0.0002859,0.0002830,0.0002802,0.0002774,
0.0002747,0.0002719,0.0002692,0.0002665,0.0002639,0.0002613,0.0002587,0.0002561,0.0002535,0.0002510,
0.0002485,0.0002460,0.0002436,0.0002412,0.0002388,0.0002364,0.0002340,0.0002317,0.0002294,0.0002271,
0.0002249,0.0002226,0.0002204,0.0002182,0.0002161,0.0002139,0.0002118,0.0002097,0.0002076,0.0002055,
0.0002035,0.0002014,0.0001994,0.0001975,0.0001955,0.0001935,0.0001916,0.0001897,0.0001878,0.0001860,
0.0001841,0.0001823,0.0001805,0.0001787,0.0001769,0.0001751,0.0001734,0.0001717,0.0001700,0.0001683,
0.0001666,0.0001649,0.0001633,0.0001617,0.0001601,0.0001585,0.0001569,0.0001553,0.0001538,0.0001522,
0.0001507,0.0001492,0.0001477,0.0001463,0.0001448,0.0001434,0.0001420,0.0001405,0.0001391,0.0001378,
0.0001364,0.0001350,0.0001337,0.0001324,0.0001310,0.0001297,0.0001284,0.0001272,0.0001259,0.0001247,
0.0001234,0.0001222,0.0001210,0.0001198,0.0001186,0.0001174,0.0001162,0.0001151,0.0001139,0.0001128,
0.0001117,0.0001106,0.0001095,0.0001084,0.0001073,0.0001062,0.0001052,0.0001041,0.0001031,0.0001021,
0.0001010,0.0001000,0.0000990,0.0000981,0.0000971,0.0000961,0.0000952,0.0000942,0.0000933,0.0000923,
0.0000914,0.0000905,0.0000896,0.0000887,0.0000878,0.0000870,0.0000861,0.0000852,0.0000844,0.0000836,
0.0000827,0.0000819,0.0000811,0.0000803,0.0000795,0.0000787,0.0000779,0.0000771,0.0000764,0.0000756,
0.0000749,0.0000741,0.0000734,0.0000726,0.0000719,0.0000712,0.0000705,0.0000698,0.0000691,0.0000684,
0.0000677,0.0000671,0.0000664,0.0000657,0.0000651,0.0000644,0.0000638,0.0000631,0.0000625,0.0000619,
0.0000613,0.0000607,0.0000601,0.0000595,0.0000589,0.0000583,0.0000577,0.0000571,0.0000566,0.0000560,
0.0000555,0.0000549,0.0000544,0.0000538,0.0000533,0.0000527,0.0000522,0.0000517,0.0000512,0.0000507,
0.0000502,0.0000497,0.0000492,0.0000487,0.0000482,0.0000477,0.0000473,0.0000468,0.0000463,0.0000459
};

float dist(vector<float> x1,vector<float> x2)
{
    float dis=0.0;
    for(int i=0;i<3;i++)
    {
        dis=dis + (x1[i]-x2[i])*(x1[i]-x2[i]);
    }
    dis = (dis);
    return dis;
}

float GetVdwRadius(const string& atm)
{
    if(atm[0] =='H')  return 1.00;
//    if(atm == "CA")   return 1.85;
    if(atm[0] =='C')  return 1.80;
    if(atm[0] =='N')  return 1.65;
    if(atm[0] =='O')  return 1.40;
    if(atm[0] =='S')  return 1.85;
    if(atm[0] =='P')  return 1.90;
    else return 1.70; 
}

// Anbor enger form wenyi
vwdt GetVdwRadiusx(const string& atm)
{
    vwdt atmt;
    if(atm[0] =='H') 
    {
        atmt.radius = 1.00;
        atmt.well_depth = 0.0157;
        return atmt;
    }
//    if(atm == "CA")   return 1.85;
    if(atm[0] =='C')
    {
        atmt.radius = 1.908;
        atmt.well_depth = 0.1094;
        return atmt;
    }
    if(atm[0] =='N') 
    {
        atmt.radius = 1.824;
        atmt.well_depth = 0.170;
        return atmt;
    }    
    if(atm[0] =='O') 
    {
        atmt.radius = 1.6612;
        atmt.well_depth = 0.210;
        return atmt;
    }
    if(atm[0] =='S') 
    {
        atmt.radius = 2.000;
        atmt.well_depth = 0.250;
        return atmt;
    }
    if(atm[0] =='P')
    {
        atmt.radius = 2.100;
        atmt.well_depth = 0.200;
        return atmt;
    }
    else
    {
        atmt.radius = 1.908;
        atmt.well_depth = 0.086;
        return atmt;
    }     
}

// charm enger from charm energy
vwdt GetVdwRadiusC(const string& atm)
{
    vwdt atmt;
    if(atm[0] =='H') 
    {
        atmt.radius = 0.2245;
        atmt.well_depth = -0.046;
        return atmt;
    }
//    if(atm == "CA")   return 1.85;
    if(atm =="CA")
    {
        atmt.radius = 2.275;
        atmt.well_depth = -0.02;
        return atmt;
    }
    if(atm[0] == 'C')
    {
        atmt.radius = 2.00;
        atmt.well_depth = -0.11;
        return atmt;
    }    
    if(atm[0] =='N') 
    {
        atmt.radius = 1.85;
        atmt.well_depth = -0.2;
        return atmt;
    }    
    if(atm[0] =='O') 
    {
        atmt.radius = 1.7;
        atmt.well_depth = -0.12;
        return atmt;
    }
    if(atm[0] =='S') 
    {
        atmt.radius = 2.000;
        atmt.well_depth = -0.45;
        return atmt;
    }
    if(atm[0] =='P')
    {
        atmt.radius = 2.15;
        atmt.well_depth = -0.585;
        return atmt;
    }
    else
    {
        atmt.radius = 1.9924;
        atmt.well_depth = -0.07;
    //    cout<<"error"<<endl;
        return atmt;
    }     
}

/*
float GetVdwEg(string& R1,string& R2,float Dij)
{
	float A = GetVdwRadius(R1)/Dij;
	float B = GetVdwRadius(R2)/Dij;
	float E = 4.0*pow(A,12.0)-4.0*pow(B,6.0);
	return E;
} */

float GetVdwEg(string& R1,string& R2,float Dij)
{
    vwdt A1 = GetVdwRadiusx(R1);
    vwdt B1 = GetVdwRadiusx(R2);
    float vwdA1 = sqrt((A1.well_depth)*pow(2*A1.radius,12));
    float vwdA2 = sqrt((A1.well_depth)*2*pow(2*A1.radius,6));
    float vwdB1 = sqrt((B1.well_depth)*pow(2*B1.radius,12));
    float vwdB2 = sqrt((B1.well_depth)*2*pow(2*B1.radius,6));

    float E = vwdA1*vwdB1/(pow(Dij,12))-vwdA2*vwdB2/(pow(Dij,6));  

    return E;
}

float GetVdwEgC(string& R1,string& R2,float Dij)
{
    vwdt A1 = GetVdwRadiusC(R1);
    vwdt B1 = GetVdwRadiusC(R2);
    float eij=sqrt(A1.well_depth*B1.well_depth);
    float rij=0.9*(A1.radius+B1.radius);
//    float rij=(A1.radius+B1.radius);

    float E = eij*((pow(rij/Dij,12))-2.0*(pow(rij/Dij,6)));  

    return E;
}

float GetVdwEgCT(string& R1,string& R2,float Dij)
{
    vwdt A1 = GetVdwRadiusC(R1);
    vwdt B1 = GetVdwRadiusC(R2);
    float eij=sqrt(A1.well_depth*B1.well_depth);
    float rij=0.92*(A1.radius+B1.radius);
//    float rij=(A1.radius+B1.radius);

    float E = eij*((pow(rij/Dij,12))-2.0*(pow(rij/Dij,6)));  

    return E;
}

float GetVdwEgCG(string& R1,string& R2,float Dij)
{
    float E=0.0;
    vwdt A1 = GetVdwRadiusC(R1);
    vwdt B1 = GetVdwRadiusC(R2);
    float eij=sqrt(A1.well_depth*B1.well_depth);
    float rij=0.95*(A1.radius+B1.radius);
//    float rij=(A1.radius+B1.radius);
    float rdij=rij/Dij;
    if(rdij<=1.1225)
    {
        E = eij*((pow(rij/Dij,12))-2.0*(pow(rij/Dij,6)));  // big
    } else
    {
        E = 50.0-56.125*(Dij/rij);
//        E = 10.0-10.6125*(Dij/rij);  // small
//        E = 10.0*cos(0.56125*(PI)*Dij/rij); // medium
    } 

    return E;
}


float GetVdwEgCH(string& R1,string& R2,float Dij)
{
    vwdt A1 = GetVdwRadiusC(R1);
    vwdt B1 = GetVdwRadiusC(R2);
    float Aii = A1.well_depth*(pow(2.0*A1.radius*0.9,12));
    float Bii = A1.well_depth*(pow(2.0*A1.radius*0.9,6));
    float Ajj = B1.well_depth*(pow(2.0*B1.radius*0.9,12));
    float Bjj = B1.well_depth*(pow(2.0*B1.radius*0.9,6));
    float Aij = sqrt(Aii*Ajj);
    float Bij = sqrt(Bii*Bjj);
    float E = Aij/(pow(Dij,12))-Bij/(pow(Dij,6)); 

    return E;
}  

float GetVdwEgCx(string& R1,string& R2,float Dij)
{
    vwdt A1 = GetVdwRadiusC(R1);
    vwdt B1 = GetVdwRadiusC(R2);
    float eij=sqrt(A1.well_depth*B1.well_depth);
    float rij=0.9*(A1.radius+B1.radius);
//    float rij=(A1.radius+B1.radius);

    float E = eij*((pow(rij/Dij,9))-2.0*(pow(rij/Dij,6)));  

    return E;
}

float GetVdwEg(float& R1,float& R2,float Dij)
{
    float A = R1*R2/Dij;
    float B = R1*R2/Dij;
//    float E = 4.0*pow(A,12.0)-4.0*pow(B,6.0);
    float E = 4.0*pow(A,12.0)-4.0*pow(B,6.0);
    return E;
}

void Spicker(vector<vector<point3f>> vect_decstr,vector<vector<point3f>> &out_model)
{
    int pp_t = vect_decstr.size();
    int res_num = vect_decstr[0].size();
    vector<vector<float>> vect_rmsd(pp_t,vector<float>(pp_t,0.0));
    int num_clu=0;
//    vector<vector<int> > out_model_num;
    float tran[3];
    float rotatx[3][3];
    // calculate initial RMSD matrix
    for(int kk=0;kk<pp_t;kk++)
    {
        vector<point3f> decstrx;
        vector<point3f>().swap(decstrx);
        decstrx=vect_decstr[kk];
        int detcx = decstrx.size();
        vector<vector<float> > x_vect(detcx,vector<float>(3,0.0));
        for(int i=0;i<detcx;i++)
        {
            vector<float> tmk(3,0.0);
            tmk[0] = decstrx[i].x;
            tmk[1] = decstrx[i].y;
            tmk[2] = decstrx[i].z;
            x_vect[i] = tmk;
        }
        for(int tt=kk+1;tt<pp_t;tt++)
        {
            vector<point3f> decstry;
            vector<point3f>().swap(decstry);
            decstry=vect_decstr[tt];
            int detcy = decstry.size();
            float sumx=0.0;
            vector<vector<float> > y_vect(detcx,vector<float>(3,0.0));
            for(int i=0;i<detcy;i++)
            {
                vector<float> tmk(3,0.0);
                tmk[0] = decstry[i].x;
                tmk[1] = decstry[i].y;
                tmk[2] = decstry[i].z;
                y_vect[i] = tmk;
            }            
            Kabsch(x_vect,y_vect,detcx,0,sumx,tran,rotatx);
/*            for(int k=0;k<res_num;k++)
            {
                vector<float> vectA(3,0.0);
                vectA[0] = decstrx[k].x; 
                vectA[1] = decstrx[k].y;
                vectA[2] = decstrx[k].z;
                vector<float> vectB(3,0.0);
                vectB[0] = decstry[k].x;
                vectB[1] = decstry[k].y;
                vectB[2] = decstry[k].z; 
                sumx = sumx +  (vectA[0] - vectB[0])*(vectA[0] - vectB[0]) + (vectA[1] - vectB[1])*(vectA[1] - vectB[1]) +(vectA[2] - vectB[2])*(vectA[2] - vectB[2]);
            }
            sumx = sumx/float(res_num);
            sumx = sqrt(sumx); */
            vect_rmsd[kk][tt] = sumx;
            vect_rmsd[tt][kk] = sumx;
        }
    }
    vector<int> neib_num(pp_t,0);
    cout<<"GGGG"<<endl;

    // calculate the maximum cluster, R_cut
    int max_num = 0;
    int max_index = 0;
    float R_cut = 7.5;
    int N=0;
    vector<int> max_clu;
    while(N<2000)
    {
        vector<vector<int>> neib_index(pp_t);

        neib_num = vector<int>(pp_t,0);
        vector<int>().swap(max_clu);
        max_num = 0;
        max_index = 0;
        for(int kk=0;kk<pp_t;kk++)
        {
            for(int tt=kk+1;tt<pp_t;tt++)
            {
                float tmp_rmsd = vect_rmsd[kk][tt];
                if(tmp_rmsd < R_cut) 
                {
                    neib_index[kk].push_back(tt);
                    neib_index[tt].push_back(kk);
                    neib_num[kk] = neib_num[kk] + 1;
                    neib_num[tt] = neib_num[tt] + 1;
                }
            }
            if(neib_num[kk]> max_num) 
            {
                max_num = neib_num[kk];
                max_index = kk;
            }
        }
        float tmp_cut = float(max_num)/float(pp_t);
        if(tmp_cut>0.7 && R_cut>3.5)
        {
            R_cut = R_cut-0.1;
        } else
        {
            if(tmp_cut < 0.15 && R_cut<12.0)
            {
                R_cut = R_cut + 0.1;
            } else
            {
                max_clu = neib_index[max_index];
                num_clu = num_clu + 1;
//                out_model_num.push_back(max_clu);
                break;
            }
        }
        N = N+1;
        vector<vector<int>>().swap(neib_index);        
    }  // maximum cluster
    out_model.push_back(vect_decstr[max_index]);
    cout<<"hhkk"<<endl;

    vector<int> all_clu;
    // calculate 5 cluster 
    while(pp_t>5 )
    {
        int pp_u = 0;
        int tmp_clu = max_clu.size();
        for(int i=0;i<tmp_clu;i++)
        {
            all_clu.push_back(max_clu[i]);
        }
        int pp_g = all_clu.size();

        vector<vector<int>> neib_indey(pp_t);
        vector<int> neib_numy(pp_t,0);
        max_num = 0;
        int max_indey = 0;
        for(int kk=0;kk<pp_t;kk++)
        {
            bool tmp_FT = false;
            for(int k=0;k<pp_g;k++)
            {
                if(kk == all_clu[k]) 
                {
                    tmp_FT = true;
                    break;
                }
            }
            if(tmp_FT) continue;
            for(int tt=kk+1;tt<pp_t;tt++)
            {
                bool tmp_FTF = false;
                for(int k=0;k<pp_g;k++)
                {
                    if(tt == all_clu[k]) 
                    {
                        tmp_FTF = true;
                        break;
                    }
                }
                if(tmp_FTF) continue;

                float tmp_rmsd = vect_rmsd[kk][tt];
//                vect_rmsdy[xx][yy] = tmp_rmsd;

                if(tmp_rmsd < R_cut) 
                {
                    neib_indey[kk].push_back(tt);
                    neib_indey[tt].push_back(kk);
                    neib_numy[kk] = neib_numy[kk] + 1;
                    neib_numy[tt] = neib_numy[tt] + 1;
                }
            }
            if(neib_numy[kk]> max_num) 
            {
                max_num = neib_numy[kk];
                max_indey = kk; 
                pp_u = kk;             
            }
        } 
        cout<<"KKK"<<endl;
        vector<int>().swap(max_clu);
        max_clu = neib_indey[max_indey];
        num_clu = num_clu + 1;
        out_model.push_back(vect_decstr[pp_u]);
//        out_model_num.push_back(max_clu);  
        if(num_clu == 5) break;

        vector<vector<int>>().swap(neib_indey);
        vector<int>().swap(neib_numy);
    }


/*    while(pp_t>5 )
    {
        int pp_u = 0;
        int tmp_clu = max_clu.size();
        int pp_y = pp_t - tmp_clu;
        if(pp_y<0) break;
        vector<vector<float>> vect_rmsdy(pp_y,vector<float>(pp_y,0.0));
        vector<vector<int>> neib_indey(pp_y);
        vector<int> neib_numy(pp_y,0);
        int xx =0;
        int yy =0;
        max_num = 0;
        int max_indey = 0;
        for(int kk=0;kk<pp_t;kk++)
        {
            bool tmp_FT = false;
            for(int k=0;k<tmp_clu;k++)
            {
                if(kk == max_clu[k]) 
                {
                    tmp_FT = true;
                    break;
                }
            }
            if(tmp_FT) continue;
            for(int tt=kk+1;tt<pp_t;tt++)
            {
                bool tmp_FTF = false;
                for(int k=0;k<tmp_clu;k++)
                {
                    if(tt == max_clu[k]) 
                    {
                        tmp_FTF = true;
                        break;
                    }
                }
                if(tmp_FTF) continue;

                float tmp_rmsd = vect_rmsd[kk][tt];
                vect_rmsdy[xx][yy] = tmp_rmsd;

                if(tmp_rmsd < R_cut) 
                {
                    neib_indey[xx].push_back(xx);
                    neib_indey[yy].push_back(yy);
                    neib_numy[xx] = neib_numy[xx] + 1;
                    neib_numy[yy] = neib_numy[yy] + 1;
                }
                yy = yy + 1;
            }
            if(neib_numy[xx]> max_num) 
            {
                max_num = neib_numy[xx];
                max_indey = xx; 
                pp_u = kk;             
            }
            xx = xx + 1;
        } 
        cout<<"KKK"<<endl;
        vector<int>().swap(max_clu);
        max_clu = neib_indey[max_indey];
        num_clu = num_clu + 1;
        out_model.push_back(vect_decstr[pp_u]);
//        out_model_num.push_back(max_clu);  
        if(num_clu == 5) break;

        vector<vector<float>>().swap(vect_rmsd);   
        vect_rmsd = vect_rmsdy;
        pp_t = vect_rmsd.size();
        vector<vector<int>>().swap(neib_indey);
        vector<int>().swap(neib_numy);
    } */
    return ;

}

void Spickerx(vector<vector<point3f>> vect_decstr,vector<vector<point3f>> &out_model)
{
    int pp_t = vect_decstr.size();
    int res_num = vect_decstr[0].size();
    vector<vector<float>> vect_rmsd(pp_t,vector<float>(pp_t,0.0));
    int num_clu=0;
//    vector<vector<int> > out_model_num;
    float tran[3];
    float rotatx[3][3];
    // calculate initial RMSD matrix
    for(int kk=0;kk<pp_t;kk++)
    {
        vector<point3f> decstrx;
        vector<point3f>().swap(decstrx);
        decstrx=vect_decstr[kk];
        int detcx = decstrx.size();
    /*    vector<vector<float> > x_vect(detcx,vector<float>(3,0.0));
        for(int i=0;i<detcx;i++)
        {
            vector<float> tmk(3,0.0);
            tmk[0] = decstrx[i].x;
            tmk[1] = decstrx[i].y;
            tmk[2] = decstrx[i].z;
            x_vect[i] = tmk;
        }  */
        for(int tt=kk+1;tt<pp_t;tt++)
        {
            vector<point3f> decstry;
            vector<point3f>().swap(decstry);
            decstry=vect_decstr[tt];
            float sumx=0.0;
    /*        int detcy = decstry.size();
            vector<vector<float> > y_vect(detcx,vector<float>(3,0.0));
            for(int i=0;i<detcy;i++)
            {
                vector<float> tmk(3,0.0);
                tmk[0] = decstry[i].x;
                tmk[1] = decstry[i].y;
                tmk[2] = decstry[i].z;
                y_vect[i] = tmk;
            }           */ 
//            Kabsch(x_vect,y_vect,detcx,0,sumx,tran,rotatx);
            for(int k=0;k<res_num;k++)
            {
                vector<float> vectA(3,0.0);
                vectA[0] = decstrx[k].x; 
                vectA[1] = decstrx[k].y;
                vectA[2] = decstrx[k].z;
                vector<float> vectB(3,0.0);
                vectB[0] = decstry[k].x;
                vectB[1] = decstry[k].y;
                vectB[2] = decstry[k].z; 
                sumx = sumx +  (vectA[0] - vectB[0])*(vectA[0] - vectB[0]) + (vectA[1] - vectB[1])*(vectA[1] - vectB[1]) +(vectA[2] - vectB[2])*(vectA[2] - vectB[2]);
            }
            sumx = sumx/float(res_num);
            sumx = sqrt(sumx); 
            vect_rmsd[kk][tt] = sumx;
            vect_rmsd[tt][kk] = sumx;
        }
    }
    vector<int> neib_num(pp_t,0);
    cout<<"GGGG"<<endl;

    // calculate the maximum cluster, R_cut
    int max_num = 0;
    int max_index = 0;
    float R_cut = 7.5;
    int N=0;
    vector<int> max_clu;
    while(N<2000)
    {
        vector<vector<int>> neib_index(pp_t);

        neib_num = vector<int>(pp_t,0);
        vector<int>().swap(max_clu);
        max_num = 0;
        max_index = 0;
        for(int kk=0;kk<pp_t;kk++)
        {
            for(int tt=kk+1;tt<pp_t;tt++)
            {
                float tmp_rmsd = vect_rmsd[kk][tt];
                if(tmp_rmsd < R_cut) 
                {
                    neib_index[kk].push_back(tt);
                    neib_index[tt].push_back(kk);
                    neib_num[kk] = neib_num[kk] + 1;
                    neib_num[tt] = neib_num[tt] + 1;
                }
            }
            if(neib_num[kk]> max_num) 
            {
                max_num = neib_num[kk];
                max_index = kk;
            }
        }
        float tmp_cut = float(max_num)/float(pp_t);
        if(tmp_cut>0.7 && R_cut>3.5)
        {
            R_cut = R_cut-0.1;
        } else
        {
            if(tmp_cut < 0.15 && R_cut<12.0)
            {
                R_cut = R_cut + 0.1;
            } else
            {
                max_clu = neib_index[max_index];
                num_clu = num_clu + 1;
//                out_model_num.push_back(max_clu);
                break;
            }
        }
        N = N+1;
        vector<vector<int>>().swap(neib_index);        
    }  // maximum cluster
    out_model.push_back(vect_decstr[max_index]);
    cout<<"hhkk"<<endl;

    vector<int> all_clu;
    // calculate 5 cluster 
    while(pp_t>5 )
    {
        int pp_u = 0;
        int tmp_clu = max_clu.size();
        for(int i=0;i<tmp_clu;i++)
        {
            all_clu.push_back(max_clu[i]);
        }
        int pp_g = all_clu.size();

        vector<vector<int>> neib_indey(pp_t);
        vector<int> neib_numy(pp_t,0);
        max_num = 0;
        int max_indey = 0;
        for(int kk=0;kk<pp_t;kk++)
        {
            bool tmp_FT = false;
            for(int k=0;k<pp_g;k++)
            {
                if(kk == all_clu[k]) 
                {
                    tmp_FT = true;
                    break;
                }
            }
            if(tmp_FT) continue;
            for(int tt=kk+1;tt<pp_t;tt++)
            {
                bool tmp_FTF = false;
                for(int k=0;k<pp_g;k++)
                {
                    if(tt == all_clu[k]) 
                    {
                        tmp_FTF = true;
                        break;
                    }
                }
                if(tmp_FTF) continue;

                float tmp_rmsd = vect_rmsd[kk][tt];
//                vect_rmsdy[xx][yy] = tmp_rmsd;

                if(tmp_rmsd < R_cut) 
                {
                    neib_indey[kk].push_back(tt);
                    neib_indey[tt].push_back(kk);
                    neib_numy[kk] = neib_numy[kk] + 1;
                    neib_numy[tt] = neib_numy[tt] + 1;
                }
            }
            if(neib_numy[kk]> max_num) 
            {
                max_num = neib_numy[kk];
                max_indey = kk; 
                pp_u = kk;             
            }
        } 
        cout<<"KKK"<<endl;
        vector<int>().swap(max_clu);
        max_clu = neib_indey[max_indey];
        num_clu = num_clu + 1;
        out_model.push_back(vect_decstr[pp_u]);
//        out_model_num.push_back(max_clu);  
        if(num_clu == 5) break;

        vector<vector<int>>().swap(neib_indey);
        vector<int>().swap(neib_numy);
    }


/*    while(pp_t>5 )
    {
        int pp_u = 0;
        int tmp_clu = max_clu.size();
        int pp_y = pp_t - tmp_clu;
        if(pp_y<0) break;
        vector<vector<float>> vect_rmsdy(pp_y,vector<float>(pp_y,0.0));
        vector<vector<int>> neib_indey(pp_y);
        vector<int> neib_numy(pp_y,0);
        int xx =0;
        int yy =0;
        max_num = 0;
        int max_indey = 0;
        for(int kk=0;kk<pp_t;kk++)
        {
            bool tmp_FT = false;
            for(int k=0;k<tmp_clu;k++)
            {
                if(kk == max_clu[k]) 
                {
                    tmp_FT = true;
                    break;
                }
            }
            if(tmp_FT) continue;
            for(int tt=kk+1;tt<pp_t;tt++)
            {
                bool tmp_FTF = false;
                for(int k=0;k<tmp_clu;k++)
                {
                    if(tt == max_clu[k]) 
                    {
                        tmp_FTF = true;
                        break;
                    }
                }
                if(tmp_FTF) continue;

                float tmp_rmsd = vect_rmsd[kk][tt];
                vect_rmsdy[xx][yy] = tmp_rmsd;

                if(tmp_rmsd < R_cut) 
                {
                    neib_indey[xx].push_back(xx);
                    neib_indey[yy].push_back(yy);
                    neib_numy[xx] = neib_numy[xx] + 1;
                    neib_numy[yy] = neib_numy[yy] + 1;
                }
                yy = yy + 1;
            }
            if(neib_numy[xx]> max_num) 
            {
                max_num = neib_numy[xx];
                max_indey = xx; 
                pp_u = kk;             
            }
            xx = xx + 1;
        } 
        cout<<"KKK"<<endl;
        vector<int>().swap(max_clu);
        max_clu = neib_indey[max_indey];
        num_clu = num_clu + 1;
        out_model.push_back(vect_decstr[pp_u]);
//        out_model_num.push_back(max_clu);  
        if(num_clu == 5) break;

        vector<vector<float>>().swap(vect_rmsd);   
        vect_rmsd = vect_rmsdy;
        pp_t = vect_rmsd.size();
        vector<vector<int>>().swap(neib_indey);
        vector<int>().swap(neib_numy);
    } */
    return ;

}

int energyclash(vector<point3f> &decstr,int numseq,float &fene)
{
    int i,j,ii,jj;
    lableproblematic lp1;
    point3d pin[8],pjn[8],tp;
    double tdist;
    double totene=0;
    bool *flagpos=new bool[numseq]; 
    int ncla=0;
    int nbro=0;
    for(i=0;i<numseq;i++)
    {
        flagpos[i]=false;
    }
    for(i=0;i<numseq;i++)
    {
        pin[0]=setv(decstr[i].x,decstr[i].y,decstr[i].z);
        pin[1]=setv(decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z);
        pin[2]=setv(decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z);
        pin[3]=setv(decstr[i].pto.x,decstr[i].pto.y,decstr[i].pto.z);
        pin[4]=setv(decstr[i].ptb.x,decstr[i].ptb.y,decstr[i].ptb.z);
        pin[5]=setv(decstr[i].pth.x,decstr[i].pth.y,decstr[i].pth.z);
        pin[6]=setv(decstr[i].ptg.x,decstr[i].ptg.y,decstr[i].ptg.z);
        pin[7]=setv(decstr[i].ptsg.x,decstr[i].ptsg.y,decstr[i].ptsg.z);
        for(j=i+1;j<numseq;j++)
        {
            pjn[0]=setv(decstr[j].x,decstr[j].y,decstr[j].z);
            pjn[1]=setv(decstr[j].ptn.x,decstr[j].ptn.y,decstr[j].ptn.z);
            pjn[2]=setv(decstr[j].ptc.x,decstr[j].ptc.y,decstr[j].ptc.z);
            pjn[3]=setv(decstr[j].pto.x,decstr[j].pto.y,decstr[j].pto.z);
            pjn[4]=setv(decstr[j].ptb.x,decstr[j].ptb.y,decstr[j].ptb.z);
            pjn[5]=setv(decstr[j].pth.x,decstr[j].pth.y,decstr[j].pth.z);
            pjn[6]=setv(decstr[j].ptg.x,decstr[j].ptg.y,decstr[j].ptg.z);
            pjn[7]=setv(decstr[j].ptsg.x,decstr[j].ptsg.y,decstr[j].ptsg.z);
            for(ii=0;ii<8;ii++)
            {
                if(ii==4 || ii==5  || ii==6) continue;
                for(jj=0;jj<8;jj++)
                {
                    if(jj==4 || jj==5  || jj==6) continue;
                    tp=minu(pin[ii],pjn[jj]);
                    tdist=tp.x*tp.x+tp.y*tp.y+tp.z*tp.z;
                    tdist-=vdwds[ii][jj];
                    if(tdist<0)
                    {
                         flagpos[i]=true;
                         flagpos[j]=true;
                         totene-=tdist;
                    } 
                    if(ii==0 && jj==0)
                    {
                        tdist=norm(tp);
                        if(j==i+1)
                        {
                            if(tdist>4.2) 
                            {
                                nbro++;
                            }
                        }
                        if(tdist<3.6) 
                        {
                            ncla++;
                        }
                    }
                }
            }
        }
    }
    lp1.nn[8]=0;
    lp1.indn[8] = new int [numseq];
    lp1.indn[28] = new int [numseq];  
    memset(lp1.indn[8],0,numseq*sizeof(int));  
    memset(lp1.indn[28],0,numseq*sizeof(int));   
    for(i=0;i<numseq;i++)
    {
        if(flagpos[i])
        {
            lp1.indn[8][lp1.nn[8]]=i;
            lp1.nn[8]++;
            lp1.indn[28][i]=1;
        }
    }
    int tmp_nn= lp1.nn[8];
    delete[]flagpos;
    delete [] lp1.indn[8];
    delete [] lp1.indn[28];    
    fene=totene;
    return tmp_nn+nbro+5*ncla;
//    return lp1.nn[8]+nbro+5*ncla;
}
double energyexcludedvolume(vector<point3f> decstr,int numseq)
{
    int i,j,ii,jj;
    point3d pin[8],pjn[8],tp;
//  BasicFunc bf;
    int numcon =8;
/*    if(numcon>8)
    {
        printf("should less than 8 %d\n",numcon);
        return 0;
    } */
    //////////////////  ca    n    c    o   cb   h   ha   sg
//  double dthresh[8]={2.45,2.75,1.35,2.25,2.90,1.55,1.65,2.50};
//  double dthresh[8]={6.125,7.5625,1.8225,5.0625,8.41,2.4025,2.7225,6.25};
//  double rvdw[8]={1.70,1.55,1.70,1.52,1.70,1.20,1.20,2.50};
    double ttot=0;
    double tdist;
//  double th;
//  int totnum=0;
    for(i=0;i<numseq-1;i++)
    {
        pin[0]=setv(decstr[i].x,decstr[i].y,decstr[i].z);
        pin[1]=setv(decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z);
        pin[2]=setv(decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z);
        pin[3]=setv(decstr[i].pto.x,decstr[i].pto.y,decstr[i].pto.z);
        pin[4]=setv(decstr[i].ptb.x,decstr[i].ptb.y,decstr[i].ptb.z);
//      pin[5]=bf.setv(decstr[i].pth.x,decstr[i].pth.y,decstr[i].pth.z);
//      pin[6]=bf.setv(decstr[i].ptha.x,decstr[i].ptha.y,decstr[i].ptha.z);
        pin[7]=setv(decstr[i].ptsg.x,decstr[i].ptsg.y,decstr[i].ptsg.z);
        for(j=i+1;j<numseq;j++)
        {
            pjn[0]=setv(decstr[j].x,decstr[j].y,decstr[j].z);
            pjn[1]=setv(decstr[j].ptn.x,decstr[j].ptn.y,decstr[j].ptn.z);
            pjn[2]=setv(decstr[j].ptc.x,decstr[j].ptc.y,decstr[j].ptc.z);
            pjn[3]=setv(decstr[j].pto.x,decstr[j].pto.y,decstr[j].pto.z);
            pjn[4]=setv(decstr[j].ptb.x,decstr[j].ptb.y,decstr[j].ptb.z);
//          pjn[5]=bf.setv(decstr[j].pth.x,decstr[j].pth.y,decstr[j].pth.z);
//          pjn[6]=bf.setv(decstr[j].ptha.x,decstr[j].ptha.y,decstr[j].ptha.z);
            pjn[7]=setv(decstr[j].ptsg.x,decstr[j].ptsg.y,decstr[j].ptsg.z);
            for(ii=0;ii<numcon;ii++)
            {
                if(ii>4 && ii<7) continue;
                for(jj=0;jj<numcon;jj++)
                {
                    if(jj>4 && jj<7) continue;
//                  tdist=bf.norm(bf.minu(pin[ii],pjn[jj]));
//                  if(tdist<dthresh[ii])
//                  {
//                      th=(dthresh[ii]-tdist)/2.0;
//                      tdist=th*th*(3*rvdw[ii]-th);
//                      ttot+=tdist;
//                      tdist=th*th*(3*rvdw[jj]-th);
//                      ttot+=tdist;
//                      totnum++;
//                  }//act volume
                    tp=minu(pin[ii],pjn[jj]);
                    tdist=tp.x*tp.x+tp.y*tp.y+tp.z*tp.z;
                    tdist-=vdwds[ii][jj];
                    if(tdist<0)
                    {
                        ttot-=tdist;
//                      totnum++;
                    } 
                }
            }
        }
    }
//  if(totnum>0)
//  ttot*=100*PI/3.0/double(totnum);
    return ttot;
}

int sec_str(float dis13, float dis14, float dis15, float dis24, float dis25, float dis35)
{
    int s=1;
    
    float delta=2.1;
    if(fabs(dis15-6.37)<delta)
    {
        if(fabs(dis14-5.18)<delta)
        {
            if(fabs(dis25-5.18)<delta)
            {
                if(fabs(dis13-5.45)<delta)
                {
                    if(fabs(dis24-5.45)<delta)
                    {
                        if(fabs(dis35-5.45)<delta)
                        {
                            s=2; //helix                        
                            return s;
                        }
                    }
                }
            }
        }
    }

    delta=1.42;
    if(fabs(dis15-13)<delta)
    {
        if(fabs(dis14-10.4)<delta)
        {
            if(fabs(dis25-10.4)<delta)
            {
                if(fabs(dis13-6.1)<delta)
                {
                    if(fabs(dis24-6.1)<delta)
                    {
                        if(fabs(dis35-6.1)<delta)
                        {
                            s=4; //strand
                            return s;
                        }
                    }
                }
            }
        }
    }

/*    if(dis15 < 8)
    {
        s=3; //turn
    }     */


    return s;
}

string sec_strx(float dis13, float dis14, float dis15, float dis24, float dis25, float dis35)
{
//    int s=1;
    string ssx="C";
    
    float delta=2.1;
    if(fabs(dis15-6.37)<delta)
    {
        if(fabs(dis14-5.18)<delta)
        {
            if(fabs(dis25-5.18)<delta)
            {
                if(fabs(dis13-5.45)<delta)
                {
                    if(fabs(dis24-5.45)<delta)
                    {
                        if(fabs(dis35-5.45)<delta)
                        {
                            ssx="H"; //helix                        
                            return ssx;
                        }
                    }
                }
            }
        }
    }

    delta=1.42;
    if(fabs(dis15-13)<delta)
    {
        if(fabs(dis14-10.4)<delta)
        {
            if(fabs(dis25-10.4)<delta)
            {
                if(fabs(dis13-6.1)<delta)
                {
                    if(fabs(dis24-6.1)<delta)
                    {
                        if(fabs(dis35-6.1)<delta)
                        {
                            ssx="E"; //strand
                            return ssx;
                        }
                    }
                }
            }
        }
    }

    if(dis15 < 8)
    {
        ssx="T"; //turn
    }     


    return ssx;
}

//1->coil, 2->helix, 3->turn, 4->strand
void make_sec(vector<vector<float> > x, int len, vector<int> &sec)
{
    int j1, j2, j3, j4, j5;
    float d13, d14, d15, d24, d25, d35;
    for(int i=0; i<len; i++)
    {   
        sec[i]=1;
        j1=i-2;
        j2=i-1;
        j3=i;
        j4=i+1;
        j5=i+2;     
        
        if(j1>=0 && j5<len)
        {
            d13=sqrt(dist(x[j1], x[j3]));
            d14=sqrt(dist(x[j1], x[j4]));
            d15=sqrt(dist(x[j1], x[j5]));
            d24=sqrt(dist(x[j2], x[j4]));
            d25=sqrt(dist(x[j2], x[j5]));
            d35=sqrt(dist(x[j3], x[j5]));
            sec[i]=sec_str(d13, d14, d15, d24, d25, d35);           
        }    
    } 
}

void make_secx(vector<vector<float> > x, int len, vector<string> &sec)
{
    int j1, j2, j3, j4, j5;
    float d13, d14, d15, d24, d25, d35;
    for(int i=0; i<len; i++)
    {   
        sec[i]="C";
        j1=i-2;
        j2=i-1;
        j3=i;
        j4=i+1;
        j5=i+2;     
        
        if(j1>=0 && j5<len)
        {
            d13=sqrt(dist(x[j1], x[j3]));
            d14=sqrt(dist(x[j1], x[j4]));
            d15=sqrt(dist(x[j1], x[j5]));
            d24=sqrt(dist(x[j2], x[j4]));
            d25=sqrt(dist(x[j2], x[j5]));
            d35=sqrt(dist(x[j3], x[j5]));
            sec[i]=sec_strx(d13, d14, d15, d24, d25, d35);           
        }    
    } 
}

point3d setv(float tx,float ty,float tz)
{
    point3d temp;
    temp.x=tx;
    temp.y=ty;
    temp.z=tz;
    return temp;
}

point3d minu(point3d p1,point3d p2)
{
    point3d temp;
    temp.x=p1.x-p2.x;
    temp.y=p1.y-p2.y;
    temp.z=p1.z-p2.z;
    return temp;
}

float norm(point3d p1)
{
    float t=sqrt(p1.x*p1.x+p1.y*p1.y+p1.z*p1.z);
    return t;
}

point3d prod(point3d p1,point3d p2)//from p1 to p2
{
    point3d temp;
    temp.x=p1.y*p2.z-p1.z*p2.y;
    temp.y=p1.z*p2.x-p1.x*p2.z;
    temp.z=p1.x*p2.y-p1.y*p2.x;
    return temp;
}

point3d unit(point3d p1)
{
    point3d temp;
    float t=sqrt(p1.x*p1.x+p1.y*p1.y+p1.z*p1.z);
    if(t<1e-40) 
    {
        temp.x=0;
        temp.y=0;
        temp.z=0;
        return temp;
    }
    temp.x=p1.x/t;
    temp.y=p1.y/t;
    temp.z=p1.z/t;
    return temp;
}

int getmovetype(vector<float> tmov,int totmov,float trandnum)
{
    int i;
    for(i=0;i<totmov;i++)
    {
        if(trandnum>=tmov[i] && trandnum<tmov[i+1])
        {
            return i;
        }
    }
    return 0;
}

double phi(double xi,double yi,double zi,double xj,double yj,double zj,double xk,
           double yk,double zk,double xl,double yl,double zl)
{
    double xij,yij,zij,
        xkj,ykj,zkj,
        xkl,ykl,zkl,
        dxi,dyi,dzi,
        gxi,gyi,gzi,
        bi,bk,ct,
        boi2,boj2,
        z1,z2,ap,s,
        bioj,bjoi;
    
    
    /* Calculate the vectors C,B,C                                       */
    xij = xi - xj;
    yij = yi - yj;
    zij = zi - zj;
    xkj = xk - xj;
    ykj = yk - yj;
    zkj = zk - zj;
    xkl = xk - xl;
    ykl = yk - yl;
    zkl = zk - zl;
    
    /* Calculate the normals to the two planes n1 and n2
    this is given as the cross products:
    AB x BC
    --------- = n1
    |AB x BC|
    
      BC x CD
      --------- = n2
      |BC x CD|
    */
    dxi = yij * zkj - zij * ykj;     /* Normal to plane 1                */
    dyi = zij * xkj - xij * zkj;
    dzi = xij * ykj - yij * xkj;
    gxi = zkj * ykl - ykj * zkl;     /* Mormal to plane 2                */
    gyi = xkj * zkl - zkj * xkl;
    gzi = ykj * xkl - xkj * ykl;
    
    /* Calculate the length of the two normals                           */
    bi = dxi * dxi + dyi * dyi + dzi * dzi;
    bk = gxi * gxi + gyi * gyi + gzi * gzi;
    ct = dxi * gxi + dyi * gyi + dzi * gzi;
    
    boi2 = 1./bi;
    boj2 = 1./bk;
    bi   = (double)sqrt((double)bi);
    bk   = (double)sqrt((double)bk);
    if(bi<epsilon*0.01 || bk<epsilon*0.01) return 180;
    z1   = 1./bi;
    z2   = 1./bk;
    bioj = bi * z2;
    bjoi = bk * z1;
    ct   = ct * z1 * z2;
    if (ct >  1.0)   ct = 1.0;
    if (ct < (-1.0)) ct = -1.0;
    ap   = acos(ct);
    
    s = xkj * (dzi * gyi - dyi * gzi)
        + ykj * (dxi * gzi - dzi * gxi)
        + zkj * (dyi * gxi - dxi * gyi);
    
    if (s < 0.0) ap = -ap;
    
    ap = (ap > 0.0) ? PI-ap : -(PI+ap);

    // angle
    ap *= (double)180.0/PI;
    if(ap<0) ap+=360;
    return(ap);
}

float phi(float xi,float yi,float zi,float xj,float yj,float zj,float xk,
           float yk,float zk,float xl,float yl,float zl)
{
    double xij,yij,zij,
        xkj,ykj,zkj,
        xkl,ykl,zkl,
        dxi,dyi,dzi,
        gxi,gyi,gzi,
        bi,bk,ct,
        boi2,boj2,
        z1,z2,s,
        bioj,bjoi;
    float ap;
    
    /* Calculate the vectors C,B,C                                       */
    xij = xi - xj;
    yij = yi - yj;
    zij = zi - zj;
    xkj = xk - xj;
    ykj = yk - yj;
    zkj = zk - zj;
    xkl = xk - xl;
    ykl = yk - yl;
    zkl = zk - zl;
    
    /* Calculate the normals to the two planes n1 and n2
    this is given as the cross products:
    AB x BC
    --------- = n1
    |AB x BC|
    
      BC x CD
      --------- = n2
      |BC x CD|
    */
    dxi = yij * zkj - zij * ykj;     /* Normal to plane 1                */
    dyi = zij * xkj - xij * zkj;
    dzi = xij * ykj - yij * xkj;
    gxi = zkj * ykl - ykj * zkl;     /* Mormal to plane 2                */
    gyi = xkj * zkl - zkj * xkl;
    gzi = ykj * xkl - xkj * ykl;
    
    /* Calculate the length of the two normals                           */
    bi = dxi * dxi + dyi * dyi + dzi * dzi;
    bk = gxi * gxi + gyi * gyi + gzi * gzi;
    ct = dxi * gxi + dyi * gyi + dzi * gzi;
    
    boi2 = 1./bi;
    boj2 = 1./bk;
    bi   = (double)sqrt((double)bi);
    bk   = (double)sqrt((double)bk);
    if(bi<epsilon*0.01 || bk<epsilon*0.01) return 180;
    z1   = 1./bi;
    z2   = 1./bk;
    bioj = bi * z2;
    bjoi = bk * z1;
    ct   = ct * z1 * z2;
    if (ct >  1.0)   ct = 1.0;
    if (ct < (-1.0)) ct = -1.0;
    ap   = acos(ct);
    
    s = xkj * (dzi * gyi - dyi * gzi)
        + ykj * (dxi * gzi - dzi * gxi)
        + zkj * (dyi * gxi - dxi * gyi);
    
    if (s < 0.0) ap = -ap;
    
    ap = (ap > 0.0) ? PI-ap : -(PI+ap);

    // angle
    ap *= (float)180.0/PI;
    if(ap<0) ap+=360;
    return(ap);
}

float squgaussian(float x,float sigma,float miu)
{
    float tup=-(x-miu)*(x-miu)/(2.0*sigma*sigma);
    return tup;
}

float expgaussian2(float x,float sigma,float miu)
{
    float tdist;
    float tup=(x-miu)*(x-miu)/(2.0*sigma*sigma);
    int tind=(int)(tup*100.0);
    if(tind>=1000) tind=999;
    tdist=expval[tind];
    return tdist;
}

float expgaussian(float x,float sigma,float miu)
{
    float tdist;
    float tup=-(x-miu)*(x-miu)/(2.0*sigma*sigma);
    tdist=exp(tup);
    return tdist;
}
float expgaussianq(float xx,float sigsig)
{
    float tdist;
    float tup=xx/sigsig;
    tdist=exp(tup);
    return tdist;
}

float dotv(point3d p1,point3d p2)
{
    float temp;
    temp=p1.x*p2.x+p1.y*p2.y+p1.z*p2.z;
    return temp;
}

float angv(point3d p1,point3d p2)
{
    float t1=norm(p1);
    float t2=norm(p2);
    if(t1<0.0000001 || t2<0.0000001)
    {
        return 0;
    }
    float t3=dotv(p1,p2)/t1/t2;
    if(t3<-1.0)
    {
        t3=-1.0;
    }
    else if(t3>1.0)
    {
        t3=1.0;
    }
    return(acos(t3));
}

point3d rotv(point3d p1,double rot[])
{
    point3d tv;
    tv.x=rot[0]*p1.x+rot[1]*p1.y+rot[2]*p1.z;
    tv.y=rot[3]*p1.x+rot[4]*p1.y+rot[5]*p1.z;
    tv.z=rot[6]*p1.x+rot[7]*p1.y+rot[8]*p1.z;
    return tv;
}

point3d rana(double theta)
{
    point3d tp;
    double dphi;
//  dphi=(2.0*rand()/double(RAND_MAX+1.0))*PI;
    dphi=(2.0*randf0and1())*PI;
    tp.x=sin(theta)*cos(dphi);
    tp.y=sin(theta)*sin(dphi);
    tp.z=cos(theta);
    return(tp);
}

void v2rot(point3d p1,double rot[])
{
    //the rotmatrix for [0,0,1] to p1
    double pr=sqrt(p1.x*p1.x+p1.y*p1.y+p1.z*p1.z);
    double theta=acos(p1.z/pr);
    double phi=atan2(p1.y,p1.x)+PI/2.0;
    double cosphi=cos(phi);
    double sinphi=sin(phi);
    double costheta=cos(theta);
    double sintheta=sin(theta);
    rot[0]=cosphi;rot[1]=-sinphi*costheta;rot[2]=sinphi*sintheta;
    rot[3]=sinphi;rot[4]=cosphi*costheta; rot[5]=-cosphi*sintheta;
    rot[6]=0;     rot[7]=sintheta;        rot[8]=costheta;
}

float energyhbondcanc(vector<point3f> &decstr,int numseq)
{
    int i,j,k,l,m;
//  BasicFunc bf;
    point3d tp[30],ap[30],np[60];
    float lamda;
    float diss[30];
    float totenergy=0;
    float totenergy2=0;
//  double threshval[5]={8.51719,9.21034,8.51719,8.11173,8.51719};
//  double threshval[5]={11.5129255,12.7168983,10.8197783,9.7211660,8.51719};
//  double threshval[5]={12.2060726,14.5086577,11.5129255,10.4143132,8.51719};
//  double threshval[5]={5.5214609,6.2146081,5.2983174,5.8091430,8.51719};
    float threshval[5]={17.0343864,17.7275336,12.4292162,13.8155106,8.51719};
    for(i=0;i<numseq;i++)
    {
        decstr[i].vpos=0;
        decstr[i].tpos=0;
        decstr[i].ssm='C';
        decstr[i].indl=-1;
        decstr[i].indr=-1;
        decstr[i].tpr=-1;
        decstr[i].tpl=-1;
    }

    float bmean[ ][30]={
    5.005258,4.812388,5.021389,8.549857,7.358543,0.313033,0.892573,0.434546,0.10,1.031090,1.005027,1.56,1.56,195.8602,188.2889,2.1134,2.1414,0.12,0.6143,2.3610,1.6780,0.6083,5.8643,5.8880,6.2619,6.1894,4.8166,5.3186,5.3258,4.8050,
    4.900784,4.863861,4.999239,7.266375,8.554238,0.348184,2.134800,0.443322,3.03,2.149802,2.227194,1.59,1.59,200.7937,194.2951,2.1542,2.1237,3.06,0.5957,1.6530,2.4022,0.6273,6.2469,6.1463,5.8328,5.9151,4.8085,5.2970,5.3377,4.8108,
    8.365075,4.491176,7.735837,5.339629,5.320091,2.691990,0.536371,2.615259,0.18,0.373315,2.583255,1.68,1.45,195.7795,205.6715,2.1604,2.1455,2.97,0.7491,2.4668,2.4778,1.5352,5.7154,6.2214,5.7354,6.2418,5.7361,5.5213,5.5205,6.0509,
    6.914873,5.235785,8.768489,4.648967,4.582170,2.336092,2.427004,2.741880,2.86,2.717415,0.766849,1.70,1.44,212.0750,185.1968,2.2088,2.2202,0.27,1.3055,2.5039,2.4991,0.7340,5.6729,6.3980,5.6743,6.3585,4.4500,4.0924,4.0764,4.9571};
    float bsigma[ ][30]={
    0.506992,0.277813,0.463843,0.386039,0.769766,0.329624,0.103274,0.229149,0.210566,0.097765,0.157291,0.091175,0.110380, 42.3417, 62.3685,0.1873,0.2068,0.3450,0.3303,0.2746,0.3560,0.3367,0.4202,0.4111,0.3630,0.3539,0.2653,0.2306,0.3050,0.2174,
    0.495829,0.281285,0.489831,0.769988,0.366858,0.342107,0.098529,0.242618,0.241717,0.203498,0.160810,0.109177,0.088545, 35.8454, 45.1044,0.1974,0.1815,0.3165,0.3169,0.3472,0.2630,0.3452,0.3145,0.3637,0.4117,0.4182,0.2478,0.3001,0.2290,0.2179,
    0.391946,0.316595,1.051334,0.412897,0.388383,0.403889,0.250276,0.278586,0.238441,0.252056,0.364956,0.108726,0.110183, 33.6450, 45.5233,0.1975,0.1930,0.3910,0.3124,0.3525,0.3375,0.4186,0.3643,0.5255,0.3653,0.5009,0.2996,0.2934,0.2997,0.4543,
    0.686441,0.297945,0.421034,0.543513,0.527419,0.505260,0.395371,0.216206,0.262468,0.276552,0.503292,0.112229,0.105024, 48.3067, 51.6512,0.2010,0.1934,0.4407,0.3769,0.3261,0.3278,0.3086,0.3525,0.4469,0.3832,0.4593,0.3052,0.2829,0.2854,0.2721};

    float hmean[2][18]={5.911418,5.866403,5.701616,10.037283,3.812445,247.498333,249.593577,242.269149,3.847155,4.219875,3.346957,3.102981,1.830908,2.660360,2.846917,88.5,110.0,199.007496,
    5.177662,5.150089,5.181840,8.646340,3.806791,230.448930,229.593231,229.379739,3.847155,4.219875,3.346957,3.102981,1.830908,2.660360,2.846917,88.5,110.0,199.007496};
    //static double tmean[2][8]={5.911418,5.866403,5.701616,10.037283,3.812445,247.498333,249.593577,242.269149,
    float hsigma[2][18]={0.654601,0.668008,0.711406,0.600463,0.024618,28.657094,26.636027,27.211977,0.231673,0.408075,0.315905,0.487618,0.137187,0.177197,0.353712,8.279363,11.562845,12.839635,
    0.297962,0.276835,0.305942,0.352136,0.025402,8.135139,7.403763,8.548291,0.231673,0.408075,0.315905,0.487618,0.137187,0.177197,0.353712,8.279363,11.562845,12.839635};
/*    
static  double bsigma[ ][36]={
            0.506992,0.277813,0.463843,0.386039,0.769766,0.329624,0.103274,0.229149,0.210566,0.097765,0.157291,0.091175,0.110380, 42.3417, 62.3685,0.1873,0.2068,0.3450,0.3303,0.2746,0.3560,0.3367,0.4202,0.4111,0.3630,0.3539,0.2653,0.2306,0.3050,0.2174,0.3114,0.3133,0.4834,0.4522,0.3182,0.2812,
            0.495829,0.281285,0.489831,0.769988,0.366858,0.342107,0.098529,0.242618,0.241717,0.203498,0.160810,0.109177,0.088545, 35.8454, 45.1044,0.1974,0.1815,0.3165,0.3169,0.3472,0.2630,0.3452,0.3145,0.3637,0.4117,0.4182,0.2478,0.3001,0.2290,0.2179,0.2931,0.3037,0.4205,0.4342,0.2688,0.2779,
            0.391946,0.316595,1.051334,0.412897,0.388383,0.403889,0.250276,0.278586,0.238441,0.252056,0.364956,0.108726,0.110183, 33.6450, 45.5233,0.1975,0.1930,0.3910,0.3124,0.3525,0.3375,0.4186,0.3643,0.5255,0.3653,0.5009,0.2996,0.2934,0.2997,0.4543,0.3369,0.3389,0.5017,0.4669,0.2991,0.2611,
            0.686441,0.297945,0.421034,0.543513,0.527419,0.505260,0.395371,0.216206,0.262468,0.276552,0.503292,0.112229,0.105024, 48.3067, 51.6512,0.2010,0.1934,0.4407,0.3769,0.3261,0.3278,0.3086,0.3525,0.4469,0.3832,0.4593,0.3052,0.2829,0.2854,0.2721,0.3308,0.3516,0.4713,0.4764,0.2618,0.2613};
static  double bmean[ ][36]={
            5.005258,4.812388,5.021389,8.549857,7.358543,0.313033,0.892573,0.434546,0.10,1.031090,1.005027,1.56,1.56,195.8602,188.2889,2.1134,2.1414,0.12,0.6143,2.3610,1.6780,0.6083,5.8643,5.8880,6.2619,6.1894,4.8166,5.3186,5.3258,4.8050,2.7088,2.7445,3.0439,2.9924,4.1640,4.1254,
            4.900784,4.863861,4.999239,7.266375,8.554238,0.348184,2.134800,0.443322,3.03,2.149802,2.227194,1.59,1.59,200.7937,194.2951,2.1542,2.1237,3.06,0.5957,1.6530,2.4022,0.6273,6.2469,6.1463,5.8328,5.9151,4.8085,5.2970,5.3377,4.8108,2.7326,2.7429,2.9931,2.9738,4.1251,4.1081,
            8.365075,4.491176,7.735837,5.339629,5.320091,2.691990,0.536371,2.615259,0.18,0.373315,2.583255,1.68,1.45,195.7795,205.6715,2.1604,2.1455,2.97,0.7491,2.4668,2.4778,1.5352,5.7154,6.2214,5.7354,6.2418,5.7361,5.5213,5.5205,6.0509,2.6494,2.6524,2.9774,2.9686,4.0726,4.0656,
            6.914873,5.235785,8.768489,4.648967,4.582170,2.336092,2.427004,2.741880,2.86,2.717415,0.766849,1.70,1.44,212.0750,185.1968,2.2088,2.2202,0.27,1.3055,2.5039,2.4991,0.7340,5.6729,6.3980,5.6743,6.3585,4.4500,4.0924,4.0764,4.9571,2.6294,2.5994,3.0052,3.0058,4.0919,4.0737};

    double bmean[ ][18]={
    5.005258,4.812388,5.021389,8.549857,7.358543,0.313033,0.892573,0.434546,0.10,1.031090,1.005027,1.56,1.56,195.8602,188.2889,2.1134,2.1414,0.12,
    4.900784,4.863861,4.999239,7.266375,8.554238,0.348184,2.134800,0.443322,3.03,2.149802,2.227194,1.59,1.59,200.7937,194.2951,2.1542,2.1237,3.06,
    8.365075,4.491176,7.735837,5.339629,5.320091,2.691990,0.536371,2.615259,0.18,0.373315,2.583255,1.68,1.45,195.7795,205.6715,2.1604,2.1455,2.97,
    6.914873,5.235785,8.768489,4.648967,4.582170,2.336092,2.427004,2.741880,2.86,2.717415,0.766849,1.70,1.44,212.0750,185.1968,2.2088,2.2202,0.27};

    double bsigma[ ][18]={
    0.506992,0.277813,0.463843,0.386039,0.769766,0.329624,0.103274,0.229149,0.210566,0.097765,0.157291,0.091175,0.110380,42.3417,62.3685,0.1873,0.2068,0.3450,
    0.495829,0.281285,0.489831,0.769988,0.366858,0.342107,0.098529,0.242618,0.241717,0.203498,0.160810,0.109177,0.088545,35.8454,45.1044,0.1974,0.1815,0.3165,
    0.391946,0.316595,1.051334,0.412897,0.388383,0.403889,0.250276,0.278586,0.238441,0.252056,0.364956,0.108726,0.110183,33.6450,45.5233,0.1975,0.1930,0.3910,
    0.686441,0.297945,0.421034,0.543513,0.527419,0.505260,0.395371,0.216206,0.262468,0.276552,0.503292,0.112229,0.105024,48.3067,51.6512,0.2010,0.1934,0.4407}; */            

    for(i=1;i<numseq-4;i++)
    {
        j=i+3;  
//        if((decstr[i].ss2=='E' &&  decstr[i].stype!='H') || (decstr[j].ss2=='E' &&  decstr[j].stype!='H'))
//            continue;
        if((decstr[i].ss2=='E')||(decstr[j].ss2=='E'))
            continue;        
        l=0;
        for(k=0;k<4;k++)
        {
            if(decstr[i+k].ss2=='H') 
                l++;
        }
        if(l==3 && (decstr[i+1].ss2!='H' || decstr[i+2].ss2!='H'))
            continue;
        else if(l==2 && 
            ((decstr[i+1].ss2=='H' && decstr[i+2].ss2=='H') || (decstr[i+1].ss2!='H' && decstr[i+2].ss2!='H') ))
            continue;

        if(l==0 ) continue;
        else if(l==1 ) continue;
        else if(l==2 ) lamda=0.05;
        else if(l==3) lamda=0.4;
        else if(l==4) lamda=1.5;
//        else if(m==3) lamda=0.6;
//        else if(m==4) lamda=1.0;
        else lamda=0.3;
        
        tp[1]=setv(decstr[i].x,decstr[i].y,decstr[i].z);
        ap[1]=setv(decstr[j].x,decstr[j].y,decstr[j].z);
        tp[7]=minu(tp[1],ap[1]);
        diss[1]=norm(tp[7]); 
        if(diss[1]>6.6) continue;

        tp[0]=setv(decstr[i-1].x,decstr[i-1].y,decstr[i-1].z);
        tp[2]=setv(decstr[i+1].x,decstr[i+1].y,decstr[i+1].z);
        ap[0]=setv(decstr[j-1].x,decstr[j-1].y,decstr[j-1].z);
        ap[2]=setv(decstr[j+1].x,decstr[j+1].y,decstr[j+1].z);
        tp[10]=minu(tp[2],ap[0]);
        diss[4]=norm(tp[10]);
        if(diss[4]<3.4 || diss[4]>4.2) continue;
        
        tp[3]=minu(tp[0],tp[1]);
        tp[4]=minu(tp[2],tp[1]);
        tp[5]=prod(tp[3],tp[4]);
        tp[5]=unit(tp[5]);
        tp[6]=minu(tp[0],ap[0]);
        tp[8]=minu(tp[2],ap[2]);
        tp[9]=minu(tp[0],ap[2]);
            
        diss[0]=norm(tp[6]);         
        diss[2]=norm(tp[8]); 
        diss[3]=norm(tp[9]); 
            
        diss[5]=phi(tp[0].x,tp[0].y,tp[0].z,tp[1].x,tp[1].y,tp[1].z,tp[2].x,tp[2].y,tp[2].z,
            ap[0].x,ap[0].y,ap[0].z);//phi 
        if(diss[5]>180) diss[5]-=180;
        else diss[5]+=180;
        diss[6]=phi(tp[1].x,tp[1].y,tp[1].z,tp[2].x,tp[2].y,tp[2].z,ap[0].x,ap[0].y,ap[0].z,
            ap[1].x,ap[1].y,ap[1].z);//phi 
        if(diss[6]>180) diss[6]-=180;
        else diss[6]+=180;
        diss[7]=phi(tp[2].x,tp[2].y,tp[2].z,ap[0].x,ap[0].y,ap[0].z,
            ap[1].x,ap[1].y,ap[1].z,ap[2].x,ap[2].y,ap[2].z);//phi 
        if(diss[7]>180) diss[7]-=180;
        else diss[7]+=180;
        
    //  double tval=1;
        double tval=0;
        for(k=0;k<8;k++)
        {
            if(k==1 || k==3 || k==5 || k==6 || k==7 || k==0 || k==2)
            {
            //  tval*=bf.expgaussian(diss[k],hsigma[1][k],hmean[1][k]);
                tval+=squgaussian(diss[k],hsigma[1][k],hmean[1][k]);
            }
        }
        tval=exp(tval);
        if(tval>2e-4)
        {
            if(tval>1e-2) tval=1e-2;
            tval=lamda*(threshval[4]+log(tval));
            totenergy+=tval;
            for(k=0;k<=3;k++)
            {
                decstr[i+k].ssm='H';
            }
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////beta
    for(j=1;j<numseq-4;j++) 
    {
        for(k=j+3;k<numseq-1;k++)
        {   
            if( (decstr[j].ss2=='H') || (decstr[k].ss2=='H'))
                continue;
            lamda=0.2;
            if((decstr[j].ss2=='E') && (decstr[k].ss2=='E'))
                lamda+=0.1;
            else if((decstr[j].ss2=='E') || (decstr[k].ss2=='E'))
                lamda+=0.2;
            else if(decstr[j].ss2=='E' && decstr[k].ss2=='E') 
                lamda+=1.6;
            else if(decstr[k].ss2=='E' || decstr[j].ss2=='E')
                lamda+=0.6;
            tp[1]=setv(decstr[j].x,decstr[j].y,decstr[j].z);
            ap[1]=setv(decstr[k].x,decstr[k].y,decstr[k].z);
            tp[7]=minu(tp[1],ap[1]);
            if(norm(tp[7])>7.8) continue;
            tp[0]=setv(decstr[j-1].x,decstr[j-1].y,decstr[j-1].z);   
            tp[2]=setv(decstr[j+1].x,decstr[j+1].y,decstr[j+1].z);
            tp[3]=minu(tp[0],tp[1]);
            tp[4]=minu(tp[2],tp[1]);
            tp[5]=prod(tp[3],tp[4]);
            tp[5]=unit(tp[5]);

            ap[0]=setv(decstr[k-1].x,decstr[k-1].y,decstr[k-1].z);           
            ap[2]=setv(decstr[k+1].x,decstr[k+1].y,decstr[k+1].z);
            ap[3]=minu(ap[0],ap[1]);
            ap[4]=minu(ap[2],ap[1]);
            ap[5]=prod(ap[3],ap[4]);
            ap[5]=unit(ap[5]);           
            tp[6]=minu(tp[0],ap[0]);
            
            tp[8]=minu(tp[2],ap[2]);
            tp[9]=minu(tp[0],ap[2]);
            tp[10]=minu(tp[2],ap[0]);
            tp[11]=minu(tp[0],tp[2]);
            tp[12]=minu(ap[0],ap[2]);
            
//          np[0]=bf.minu(tp[1],ap[0]);
//          np[1]=bf.minu(tp[1],ap[2]);
//          np[2]=bf.minu(ap[1],tp[0]);
//          np[3]=bf.minu(ap[1],tp[2]);
            np[4]=setv(decstr[j-1].ptn.x,decstr[j-1].ptn.y,decstr[j-1].ptn.z);
            np[5]=setv(decstr[j].ptn.x,decstr[j].ptn.y,decstr[j].ptn.z);
            np[6]=setv(decstr[j+1].ptn.x,decstr[j+1].ptn.y,decstr[j+1].ptn.z);
            np[7]=setv(decstr[j-1].ptc.x,decstr[j-1].ptc.y,decstr[j-1].ptc.z);
            np[8]=setv(decstr[j].ptc.x,decstr[j].ptc.y,decstr[j].ptc.z);
            np[9]=setv(decstr[j+1].ptc.x,decstr[j+1].ptc.y,decstr[j+1].ptc.z);
            np[10]=setv(decstr[k-1].ptn.x,decstr[k-1].ptn.y,decstr[k-1].ptn.z);
            np[11]=setv(decstr[k].ptn.x,decstr[k].ptn.y,decstr[k].ptn.z);
            np[12]=setv(decstr[k+1].ptn.x,decstr[k+1].ptn.y,decstr[k+1].ptn.z);
            np[13]=setv(decstr[k-1].ptc.x,decstr[k-1].ptc.y,decstr[k-1].ptc.z);
            np[14]=setv(decstr[k].ptc.x,decstr[k].ptc.y,decstr[k].ptc.z);
            np[15]=setv(decstr[k+1].ptc.x,decstr[k+1].ptc.y,decstr[k+1].ptc.z);
    
//          np[16]=bf.minu(tp[0],np[7]);
//          np[17]=bf.minu(np[5],np[7]);
//          np[18]=bf.minu(ap[0],np[13]);
//          np[19]=bf.minu(np[11],np[13]);
//          np[20]=bf.prod(np[16],np[17]);
//          np[20]=bf.unit(np[20]);
//          np[21]=bf.prod(np[18],np[19]);
//          np[21]=bf.unit(np[21]); 
//          np[26]=bf.minu(tp[1],np[8]);
//          np[27]=bf.minu(np[6],np[8]);
//          np[28]=bf.minu(ap[1],np[14]);
//          np[29]=bf.minu(np[12],np[14]);
//          np[30]=bf.prod(np[26],np[27]);
//          np[30]=bf.unit(np[30]);
//          np[31]=bf.prod(np[28],np[29]);
//          np[31]=bf.unit(np[31]);
            
            np[32]=minu(np[5],np[11]);//nn
            np[33]=minu(np[5],np[14]);//nc
            np[34]=minu(np[8],np[11]);//cn
            np[35]=minu(np[8],np[14]);//cc

            diss[0]=norm(tp[6]);//dist caca-1
            diss[1]=norm(tp[7]);//dist caca 
            diss[2]=norm(tp[8]);//dist caca+1            
            diss[3]=norm(tp[9]);//dist caca-1 +1
            diss[4]=norm(tp[10]);//dist caca+1 -1
            diss[5]=angv(tp[5],ap[5]);//plane plane good
            diss[6]=angv(tp[5],tp[10]);//plane1 cn good 
            diss[7]=angv(tp[11],tp[12]);//caca caca
            diss[8]=angv(tp[5],tp[7]);//plane1 caca      
            diss[9]=angv(tp[5],tp[9]);//dist nc
            diss[10]=angv(ap[5],tp[9]);//plane2 nc
            diss[11]=angv(tp[11],tp[7]);//ca1 caca  
            diss[12]=angv(tp[12],tp[7]);//ca2 caca  

            ap[11]=setv(decstr[j+2].x,decstr[j+2].y,decstr[j+2].z);
            ap[12]=setv(decstr[k-2].x,decstr[k-2].y,decstr[k-2].z);
            diss[13]=phi(tp[0].x,tp[0].y,tp[0].z,tp[1].x,tp[1].y,tp[1].z,tp[2].x,tp[2].y,tp[2].z,ap[11].x,ap[11].y,ap[11].z);//phi 
            diss[14]=phi(ap[12].x,ap[12].y,ap[12].z,ap[0].x,ap[0].y,ap[0].z,ap[1].x,ap[1].y,ap[1].z,ap[2].x,ap[2].y,ap[2].z);//phi 

            diss[15]=angv(tp[3],tp[4]);
            diss[16]=angv(ap[3],ap[4]);
            diss[17]=angv(ap[5],tp[7]);//plane2 nc
            
//          diss[18]=bf.angv(np[20],np[21]);// 4 pepetide plane
//          diss[19]=bf.angv(np[20],np[31]);// 4 pepetide plane
//          diss[20]=bf.angv(np[30],np[21]);// 4 pepetide plane
//          diss[21]=bf.angv(np[30],np[31]);// 4 pepetide plane
//          diss[22]=bf.norm(np[0]);//dist
//          diss[23]=bf.norm(np[1]);//dist
//          diss[24]=bf.norm(np[2]);//dist
//          diss[25]=bf.norm(np[3]);//dist      
            diss[26]=norm(np[32]);//dist
            diss[27]=norm(np[33]);//dist
            diss[28]=norm(np[34]);//dist
            diss[29]=norm(np[35]);//dist


            
            int delseg=-1;
            if(diss[7]<PI/2.0 && diss[8]<PI/2.0)//parallel1 left
            {
                delseg=0;
            }
            else if(diss[7]<PI/2.0 && diss[8]>=PI/2.0)//parallel2 right
            {
                delseg=1;
            }
            else if(diss[7]>=PI/2.0 && diss[8]<PI/2.0)//anti1 left
            {
                delseg=2;
            }
            else if(diss[7]>=PI/2.0 && diss[8]>=PI/2.0)//anti2 right
            {
                delseg=3;
            }
            
            if(delseg==0 && (diss[15]<1.57 || diss[15]>2.70)) continue;
            else if(delseg==1 && (diss[15]<1.57 || diss[15]>2.70)) continue;
            else if(delseg==2 && (diss[15]<1.57 || diss[15]>2.70)) continue;
            else if(delseg==3 && (diss[15]<1.57 || diss[15]>2.70)) continue;
            
            if(delseg==0 && (diss[16]<1.57 || diss[16]>2.70)) continue;
            else if(delseg==1 && (diss[16]<1.57 || diss[16]>2.70)) continue;
            else if(delseg==2 && (diss[16]<1.57 || diss[16]>2.70)) continue;
            else if(delseg==3 && (diss[16]<1.57 || diss[16]>2.70)) continue;//exclude false

            if(delseg==0 && (diss[1]<3.8 || diss[1]>6.2)) continue;
            else if(delseg==1 && (diss[1]<3.8 || diss[1]>6.2)) continue;
            else if(delseg==2 && (diss[1]<3.2 || diss[1]>6.0)) continue;
            else if(delseg==3 && (diss[1]<3.4 || diss[1]>6.4)) continue;

            if(delseg==0 && (diss[3]<6.6 || diss[3]>10.0)) continue;
            else if(delseg==1 && (diss[4]<6.6 || diss[4]>9.8)) continue;
            else if(delseg==2 && (diss[0]<6.4 || diss[0]>9.7)) continue;
            else if(delseg==3 && (diss[2]<6.4 || diss[2]>9.8)) continue;

            if(delseg==0 && ((diss[13]>30 &&diss[13]<100) || diss[13]>360)) continue;
            else if(delseg==1 && ((diss[13]>50 &&diss[13]<100) || diss[13]>310)) continue;
            else if(delseg==2 && ((diss[13]>0  &&diss[13]<80)  || diss[13]>360)) continue;
            else if(delseg==3 && ((diss[13]>50 &&diss[13]<100) || diss[13]>360)) continue;
            
            if(delseg==0 && (diss[5]>1.0)) continue;
            else if(delseg==1 && (diss[5]>1.2)) continue;
            else if(delseg==2 && (diss[5]<1.6)) continue;
            else if(delseg==3 && (diss[5]<1.0)) continue;


        //  double tval=1;
            float tval=0;
            if(delseg==0)
            {
                for(l=0;l<30;l++) //if(l!=5  && l!=7)
                //  if( l==0 || l==1 || l==2 || l==3 || l==4 || l==8 || l==17 || l==26 || l==29 || l==27 || l==28)
                if(l==11 || l==12 || l==0 || l==1 || l==2 || l==3  || l==26 || l==29 || l==27 || l==28 || l==8 || l==17)
                {
                //  tval*=bf.expgaussian2(diss[l],bsigma[delseg][l],bmean[delseg][l]);
                    tval+=squgaussian(diss[l],bsigma[delseg][l],bmean[delseg][l]);
                }
            }
            else if(delseg==1)
            {
                for(l=0;l<30;l++)  //if(l!=5  && l!=7)
                //  if( l==0 || l==1 || l==2 || l==3 || l==4  || l==8 || l==17 || l==26 || l==29 || l==27 || l==28)
                if(l==11 || l==12 || l==0 || l==1 || l==2   || l==4  || l==26 || l==29 || l==27 || l==28 || l==8 || l==17)
                {
                //  tval*=bf.expgaussian2(diss[l],bsigma[delseg][l],bmean[delseg][l]);
                    tval+=squgaussian(diss[l],bsigma[delseg][l],bmean[delseg][l]);
                }
            }
            else if(delseg==2)
            {
                for(l=0;l<30;l++) //if( l!=5 && l!=7)
                //  if( l==3 || l==1 || l==4 || l==0 || l==2 || l==8 || l==17 || l==27 || l==28 || l==26 || l==29)
                if(l==11 || l==12 || l==3 || l==1 || l==4 || l==0  || l==26 || l==29 || l==27 || l==28 || l==8 || l==17)
                {
                //  tval*=bf.expgaussian2(diss[l],bsigma[delseg][l],bmean[delseg][l]);
                    tval+=squgaussian(diss[l],bsigma[delseg][l],bmean[delseg][l]);
                }
            }
            else if(delseg==3)
            {
                for(l=0;l<30;l++) //if(l!=5 && l!=6 && l!=7 && l!=10)
                //  if( l==3 || l==1 || l==4 || l==0 || l==2 || l==8 || l==17 || l==27 || l==28 || l==26 || l==29)
                if(l==11 || l==12 || l==3 || l==1 || l==4  || l==2  || l==26 || l==29 || l==27 || l==28 || l==8 || l==17)
                {
                //  tval*=bf.expgaussian2(diss[l],bsigma[delseg][l],bmean[delseg][l]);
                    tval+=squgaussian(diss[l],bsigma[delseg][l],bmean[delseg][l]);
                }
            }
            tval=exp(tval);
            //24 23 76 75
//          if( (delseg==0 && tval>2e-4)|| (delseg==1 && tval>1e-4)||(delseg==2 && tval>2e-4)|| (delseg==3 && tval>3e-4))
//          if( (delseg==0 && tval>1e-5)|| (delseg==1 && tval>3e-6)||(delseg==2 && tval>2e-5)|| (delseg==3 && tval>6e-5))
//          if( (delseg==0 && tval>4e-3)|| (delseg==1 && tval>2e-3)||(delseg==2 && tval>5e-3)|| (delseg==3 && tval>3e-3))
            if( (delseg==0 && tval>4e-8)|| (delseg==1 && tval>2e-8)||(delseg==2 && tval>4e-6)|| (delseg==3 && tval>1e-6))
            {
                if(tval>1e-3) tval=1e-3;
                tval=lamda*(threshval[delseg]+log(tval))/threshval[delseg]*threshval[4];
//              if(tval>decstr[j].vpos)
//              {
//                  decstr[j].tpos=decstr[j].vpos;
//                  decstr[j].vpos=tval;        
//              }
//              else if(tval>decstr[j].tpos)
//              {
//                  decstr[j].tpos=tval;
//              }
//              if(tval>decstr[k].vpos)
//              {
//                  decstr[k].tpos=decstr[k].vpos;
//                  decstr[k].vpos=tval;        
//              }
//              else if(tval>decstr[k].tpos)
//              {
//                  decstr[k].tpos=tval;
//              }
//              decstr[j].ssm='E';
//              decstr[k].ssm='E';
                if(delseg==0 || delseg==2)//left of j
                {
                    if(tval>decstr[j].vpos)
                    {
                        decstr[j].vpos=tval;
                        if(decstr[j].indl>0) 
                        {       
                            if(decstr[decstr[j].indl].indr==j)//set to zero
                            {
                                decstr[decstr[j].indl].tpos=0;
                                decstr[decstr[j].indl].indr=-1;
                                if(decstr[decstr[j].indl].indl==-1)
                                    decstr[decstr[j].indl].ssm='C';
                            }
                        }
                        decstr[j].indl=k;
                        decstr[j].tpl=delseg;
                        decstr[k].ssm='E';
                        decstr[j].ssm='E';
                    }
                }
                else if(delseg==1 || delseg==3)//right of j
                {
                    if(tval>decstr[j].tpos)
                    {
                        decstr[j].tpos=tval;
                        if(decstr[j].indr>0) 
                        {
                            if(decstr[decstr[j].indr].indl==j)//set to zero
                            {
                                decstr[decstr[j].indr].vpos=0;
                                decstr[decstr[j].indr].indl=-1;
                                if(decstr[decstr[j].indr].indr==-1)
                                    decstr[decstr[j].indr].ssm='C';
                            }
                        }
                        decstr[j].indr=k;
                        decstr[j].tpr=delseg;
                        decstr[k].ssm='E';
                        decstr[j].ssm='E';
                    }
                }
                if(delseg==1 || delseg==2)//left of k
                {
                    if(tval>decstr[k].vpos)
                    {
                        decstr[k].vpos=tval;
                        if(decstr[k].indl>0) 
                        {
                            if(decstr[decstr[k].indl].indr==k)//set to zero
                            {
                                decstr[decstr[k].indl].tpos=0;
                                decstr[decstr[k].indl].indr=-1;
                                if(decstr[decstr[k].indl].indl==-1)
                                    decstr[decstr[k].indl].ssm='C';
                            }
                        }
                        decstr[k].indl=j;
                        decstr[k].tpl=delseg;
                        decstr[j].ssm='E';
                        decstr[k].ssm='E';
                    }
                }
                else if(delseg==0 || delseg==3)//right of k
                {
                    if(tval>decstr[k].tpos)
                    {
                        decstr[k].tpos=tval;
                        if(decstr[k].indr>0) 
                        {
                            if(decstr[decstr[k].indr].indl==k)//set to zero
                            {
                                decstr[decstr[k].indr].vpos=0;
                                decstr[decstr[k].indr].indl=-1;
                                if(decstr[decstr[k].indr].indr==-1)
                                    decstr[decstr[k].indr].ssm='C';
                            }
                        }
                        decstr[k].indr=j;
                        decstr[k].tpr=delseg;
                        decstr[j].ssm='E';
                        decstr[k].ssm='E';
                    }
                }   
            }
        }
    }
    for(i=0;i<numseq;i++)
    {
        totenergy2+=decstr[i].vpos;
        totenergy2+=decstr[i].tpos;
        if(decstr[i].indl!=-1 && decstr[i].indr!=-1)
        {
            totenergy2+=2.0;
            if((decstr[i].tpl==0 || decstr[i].tpl==1) && (decstr[i].tpr==0 || decstr[i].tpr==1))
                totenergy2+=1.0;//+decstr[i].vpos+decstr[i].tpos;
            else if((decstr[i].tpl==2 || decstr[i].tpl==3) && (decstr[i].tpr==2 || decstr[i].tpr==3))
                totenergy2+=1.5;//+decstr[i].vpos+decstr[i].tpos;
            else if((decstr[i].tpl==0 || decstr[i].tpl==1) && (decstr[i].tpr==2 || decstr[i].tpr==3))
                totenergy2-=0.0;
            else if((decstr[i].tpl==2 || decstr[i].tpl==3) && (decstr[i].tpr==0 || decstr[i].tpr==1))
                totenergy2-=0.0;
        }
        else if(decstr[i].indl!=-1 || decstr[i].indr!=-1)
        {
            totenergy2+=1.5;
        //  if(decstr[i].indl!=-1 && (decstr[i].tpl==2 || decstr[i].tpl==3))
        //      totenergy2+=0.5;
        //  else if(decstr[i].indr!=-1 && (decstr[i].tpr==2 || decstr[i].tpr==3))
        //      totenergy2+=0.5;
        }
    }
//  totenergy2/=2.0;
    return -(1.0*totenergy+3.0*totenergy2);
}
/*
void calcsse2x(vector<boneinfo> &bb, int numbb,vector<vector<float>> &proseq)
{
    int i,j,k,l;
    point3d tp[30],ap[30],np[60];
    double diss[30];
double tmean[2][8]={5.911418,5.866403,5.701616,10.037283,3.812445,247.498333,249.593577,242.269149,
5.177662,5.150089,5.181840,8.646340,3.806791,230.448930,229.593231,229.379739};
double tsigma[2][8]={0.654601,0.668008,0.711406,0.600463,0.024618,28.657094,26.636027,27.211977,
0.297962,0.276835,0.305942,0.352136,0.025402,8.135139,7.403763,8.548291};
//double tmean[2][8]={5.911418,5.866403,5.701616,10.037283,3.812445,247.498333,249.593577,242.269149,
//5.1,5.1,5.1,8.7,3.9,231,230,230};

    for(i=0;i<numbb;i++)
    {
        bb[i].sst=0;
    }
    for(i=1;i<numbb-4;i++)
    {
        j=i+3;  
        if(bb[i-1].indca==-1 || bb[i].indca==-1 || bb[i+1].indca==-1 || bb[i+2].indca==-1 
            || bb[i+3].indca==-1 || bb[i+4].indca==-1) continue;
        tp[1]=setv(proseq[bb[i].indca].x,proseq[bb[i].indca].y,proseq[bb[i].indca].z);
        ap[1]=setv(proseq[bb[j].indca].x,proseq[bb[j].indca].y,proseq[bb[j].indca].z);
        tp[7]=minu(tp[1],ap[1]);
        diss[1]=norm(tp[7]); 
        if(diss[1]>6.6) continue;

        tp[0]=setv(proseq[bb[i-1].indca].x,proseq[bb[i-1].indca].y,proseq[bb[i-1].indca].z);
        tp[2]=setv(proseq[bb[i+1].indca].x,proseq[bb[i+1].indca].y,proseq[bb[i+1].indca].z);
        ap[0]=setv(proseq[bb[j-1].indca].x,proseq[bb[j-1].indca].y,proseq[bb[j-1].indca].z);
        ap[2]=setv(proseq[bb[j+1].indca].x,proseq[bb[j+1].indca].y,proseq[bb[j+1].indca].z);
        tp[10]=minu(tp[2],ap[0]);
        diss[4]=norm(tp[10]);
        if(diss[4]<3.4 || diss[4]>4.2) continue;
        
        tp[3]=minu(tp[0],tp[1]);
        tp[4]=minu(tp[2],tp[1]);
        tp[5]=prod(tp[3],tp[4]);
        tp[5]=unit(tp[5]);
        tp[6]=minu(tp[0],ap[0]);
        tp[8]=minu(tp[2],ap[2]);
        tp[9]=minu(tp[0],ap[2]);
            
        diss[0]=norm(tp[6]);         
        diss[2]=norm(tp[8]); 
        diss[3]=norm(tp[9]); 
            
        diss[5]=phi(tp[0].x,tp[0].y,tp[0].z,tp[1].x,tp[1].y,tp[1].z,tp[2].x,tp[2].y,tp[2].z,
            ap[0].x,ap[0].y,ap[0].z);//phi 
        if(diss[5]>180) diss[5]-=180;
        else diss[5]+=180;
        diss[6]=phi(tp[1].x,tp[1].y,tp[1].z,tp[2].x,tp[2].y,tp[2].z,ap[0].x,ap[0].y,ap[0].z,
            ap[1].x,ap[1].y,ap[1].z);//phi 
        if(diss[6]>180) diss[6]-=180;
        else diss[6]+=180;
        diss[7]=phi(tp[2].x,tp[2].y,tp[2].z,ap[0].x,ap[0].y,ap[0].z,
            ap[1].x,ap[1].y,ap[1].z,ap[2].x,ap[2].y,ap[2].z);//phi 
        if(diss[7]>180) diss[7]-=180;
        else diss[7]+=180;
        
        double tval=1;
        for(k=0;k<8;k++)
        {
            if(k==1 || k==3 || k==5 || k==6 || k==7 || k==0 || k==2)
            {
                tval*=expgaussian(diss[k],tsigma[1][k],tmean[1][k]);
            }
        }
        if(tval>2e-4)
        {
            for(k=0;k<=3;k++)
            {
                bb[i+k].sst=1;
            }
        }
    }
//  FILE *file=fopen("betapair.txt","wt");
    //0         1         2        3       4        5         6       7        8         9       10       11       12        13      14       15    16     17
//  double bmean[ ][18]={
//  5.005258,4.812388,5.021389,8.549857,7.358543,0.313033,0.892573,0.434546,0.225614,1.031090,1.005027,1.557159,1.564989,195.8602,188.2889,2.1134,2.1414,0.2971,
//  4.900784,4.863861,4.999239,7.266375,8.554238,0.348184,2.134800,0.443322,2.842577,2.149802,2.227194,1.587717,1.589364,200.7937,194.2951,2.1542,2.1237,2.9134,
//  8.365075,4.491176,7.735837,5.339629,5.320091,2.691990,0.536371,2.615259,0.327605,0.373315,2.583255,1.674306,1.467345,195.7795,205.6715,2.1604,2.1455,2.7822,
//  6.914873,5.235785,8.768489,4.648967,4.582170,2.336092,2.427004,2.741880,2.699932,2.717415,0.766849,1.690506,1.443870,212.0750,185.1968,2.2088,2.2202,0.4856};
//  double bmean[ ][18]={
//  4.9,4.7,4.9,8.5,7.3,0.313033,0.892573,0.434546,0.12,1.031090,1.005027,1.55,1.58,195.8602,188.2889,2.1134,2.1414,0.12,
//  4.7,4.9,4.9,7.3,8.5,0.348184,2.134800,0.443322,3.00,2.149802,2.227194,1.58,1.58,200.7937,194.2951,2.1542,2.1237,3.05,
//  8.3,4.5,8.1,5.3,5.3,2.691990,0.536371,2.615259,0.18,0.373315,2.583255,1.70,1.43,195.7795,205.6715,2.1604,2.1455,2.97,
//  6.9,5.3,8.9,4.5,4.5,2.336092,2.427004,2.741880,2.86,2.717415,0.766849,1.70,1.41,212.0750,185.1968,2.2088,2.2202,0.28};
//  double bmean[ ][18]={
//  5.005258,4.812388,5.021389,8.549857,7.358543,0.313033,0.892573,0.434546,0.10,1.031090,1.005027,1.56,1.56,195.8602,188.2889,2.1134,2.1414,0.12,
//  4.900784,4.863861,4.999239,7.266375,8.554238,0.348184,2.134800,0.443322,3.03,2.149802,2.227194,1.59,1.59,200.7937,194.2951,2.1542,2.1237,3.06,
//  8.365075,4.491176,7.735837,5.339629,5.320091,2.691990,0.536371,2.615259,0.18,0.373315,2.583255,1.68,1.45,195.7795,205.6715,2.1604,2.1455,2.97,
//  6.914873,5.235785,8.768489,4.648967,4.582170,2.336092,2.427004,2.741880,2.86,2.717415,0.766849,1.70,1.44,212.0750,185.1968,2.2088,2.2202,0.27};
//
//  double bsigma[ ][18]={
//  0.506992,0.277813,0.463843,0.386039,0.769766,0.329624,0.103274,0.229149,0.210566,0.097765,0.157291,0.091175,0.110380,42.3417,62.3685,0.1873,0.2068,0.3450,
//  0.495829,0.281285,0.489831,0.769988,0.366858,0.342107,0.098529,0.242618,0.241717,0.203498,0.160810,0.109177,0.088545,35.8454,45.1044,0.1974,0.1815,0.3165,
//  0.391946,0.316595,1.051334,0.412897,0.388383,0.403889,0.250276,0.278586,0.238441,0.252056,0.364956,0.108726,0.110183,33.6450,45.5233,0.1975,0.1930,0.3910,
//  0.686441,0.297945,0.421034,0.543513,0.527419,0.505260,0.395371,0.216206,0.262468,0.276552,0.503292,0.112229,0.105024,48.3067,51.6512,0.2010,0.1934,0.4407};
//  //12 11 22 21
            double bmean[ ][30]={
            5.005258,4.812388,5.021389,8.549857,7.358543,0.313033,0.892573,0.434546,0.10,1.031090,1.005027,1.56,1.56,195.8602,188.2889,2.1134,2.1414,0.12,0.6143,2.3610,1.6780,0.6083,5.8643,5.8880,6.2619,6.1894,4.8166,5.3186,5.3258,4.8050,
            4.900784,4.863861,4.999239,7.266375,8.554238,0.348184,2.134800,0.443322,3.03,2.149802,2.227194,1.59,1.59,200.7937,194.2951,2.1542,2.1237,3.06,0.5957,1.6530,2.4022,0.6273,6.2469,6.1463,5.8328,5.9151,4.8085,5.2970,5.3377,4.8108,
            8.365075,4.491176,7.735837,5.339629,5.320091,2.691990,0.536371,2.615259,0.18,0.373315,2.583255,1.68,1.45,195.7795,205.6715,2.1604,2.1455,2.97,0.7491,2.4668,2.4778,1.5352,5.7154,6.2214,5.7354,6.2418,5.7361,5.5213,5.5205,6.0509,
            6.914873,5.235785,8.768489,4.648967,4.582170,2.336092,2.427004,2.741880,2.86,2.717415,0.766849,1.70,1.44,212.0750,185.1968,2.2088,2.2202,0.27,1.3055,2.5039,2.4991,0.7340,5.6729,6.3980,5.6743,6.3585,4.4500,4.0924,4.0764,4.9571};
            double bsigma[ ][30]={
            0.506992,0.277813,0.463843,0.386039,0.769766,0.329624,0.103274,0.229149,0.210566,0.097765,0.157291,0.091175,0.110380, 42.3417, 62.3685,0.1873,0.2068,0.3450,0.3303,0.2746,0.3560,0.3367,0.4202,0.4111,0.3630,0.3539,0.2653,0.2306,0.3050,0.2174,
            0.495829,0.281285,0.489831,0.769988,0.366858,0.342107,0.098529,0.242618,0.241717,0.203498,0.160810,0.109177,0.088545, 35.8454, 45.1044,0.1974,0.1815,0.3165,0.3169,0.3472,0.2630,0.3452,0.3145,0.3637,0.4117,0.4182,0.2478,0.3001,0.2290,0.2179,
            0.391946,0.316595,1.051334,0.412897,0.388383,0.403889,0.250276,0.278586,0.238441,0.252056,0.364956,0.108726,0.110183, 33.6450, 45.5233,0.1975,0.1930,0.3910,0.3124,0.3525,0.3375,0.4186,0.3643,0.5255,0.3653,0.5009,0.2996,0.2934,0.2997,0.4543,
            0.686441,0.297945,0.421034,0.543513,0.527419,0.505260,0.395371,0.216206,0.262468,0.276552,0.503292,0.112229,0.105024, 48.3067, 51.6512,0.2010,0.1934,0.4407,0.3769,0.3261,0.3278,0.3086,0.3525,0.4469,0.3832,0.4593,0.3052,0.2829,0.2854,0.2721};
            
    if(btscore)
    {
        delete[]btscore;
        btscore=NULL;
    }
    btscore=new double[numbb*numbb];
    for(j=0;j<numbb*numbb;j++)
        btscore[j]=-10000000;
    for(j=1;j<numbb-4;j++) 
    {
        for(k=j+3;k<numbb-1;k++)
        {   
            if(bb[j-1].indn==-1 || bb[j-1].indca==-1 || bb[j-1].indc==-1
                    || bb[j].indn==-1 || bb[j].indca==-1 || bb[j].indc==-1
                    || bb[j+1].indn==-1 || bb[j+1].indca==-1 || bb[j+1].indc==-1)
                    continue;
            if(bb[k-1].indn==-1 || bb[k-1].indca==-1 || bb[k-1].indc==-1
                    || bb[k].indn==-1 || bb[k].indca==-1 || bb[k].indc==-1
                    || bb[k+1].indn==-1 || bb[k+1].indca==-1 || bb[k+1].indc==-1)
                    continue;
            
            tp[1]=setv(proseq[bb[j].indca].x,proseq[bb[j].indca].y,proseq[bb[j].indca].z);
            ap[1]=setv(proseq[bb[k].indca].x,proseq[bb[k].indca].y,proseq[bb[k].indca].z);
            tp[7]=minu(tp[1],ap[1]);
            if(norm(tp[7])>7.8) continue;
            tp[0]=setv(proseq[bb[j-1].indca].x,proseq[bb[j-1].indca].y,proseq[bb[j-1].indca].z); 
            tp[2]=setv(proseq[bb[j+1].indca].x,proseq[bb[j+1].indca].y,proseq[bb[j+1].indca].z);
            tp[3]=minu(tp[0],tp[1]);
            tp[4]=minu(tp[2],tp[1]);
            tp[5]=prod(tp[3],tp[4]);
            tp[5]=unit(tp[5]);

            ap[0]=setv(proseq[bb[k-1].indca].x,proseq[bb[k-1].indca].y,proseq[bb[k-1].indca].z);         
            ap[2]=setv(proseq[bb[k+1].indca].x,proseq[bb[k+1].indca].y,proseq[bb[k+1].indca].z);
            ap[3]=minu(ap[0],ap[1]);
            ap[4]=minu(ap[2],ap[1]);
            ap[5]=prod(ap[3],ap[4]);
            ap[5]=unit(ap[5]);           
            tp[6]=minu(tp[0],ap[0]);
            
            tp[8]=minu(tp[2],ap[2]);
            tp[9]=minu(tp[0],ap[2]);
            tp[10]=minu(tp[2],ap[0]);
            tp[11]=minu(tp[0],tp[2]);
            tp[12]=minu(ap[0],ap[2]);
            
            np[0]=minu(tp[1],ap[0]);
            np[1]=minu(tp[1],ap[2]);
            np[2]=minu(ap[1],tp[0]);
            np[3]=minu(ap[1],tp[2]);
            np[4]=setv(proseq[bb[j-1].indn].x,proseq[bb[j-1].indn].y,proseq[bb[j-1].indn].z);
            np[5]=setv(proseq[bb[j].indn].x,proseq[bb[j].indn].y,proseq[bb[j].indn].z);
            np[6]=setv(proseq[bb[j+1].indn].x,proseq[bb[j+1].indn].y,proseq[bb[j+1].indn].z);
            np[7]=setv(proseq[bb[j-1].indc].x,proseq[bb[j-1].indc].y,proseq[bb[j-1].indc].z);
            np[8]=setv(proseq[bb[j].indc].x,proseq[bb[j].indc].y,proseq[bb[j].indc].z);
            np[9]=setv(proseq[bb[j+1].indc].x,proseq[bb[j+1].indc].y,proseq[bb[j+1].indc].z);
            np[10]=setv(proseq[bb[k-1].indn].x,proseq[bb[k-1].indn].y,proseq[bb[k-1].indn].z);
            np[11]=setv(proseq[bb[k].indn].x,proseq[bb[k].indn].y,proseq[bb[k].indn].z);
            np[12]=setv(proseq[bb[k+1].indn].x,proseq[bb[k+1].indn].y,proseq[bb[k+1].indn].z);
            np[13]=setv(proseq[bb[k-1].indc].x,proseq[bb[k-1].indc].y,proseq[bb[k-1].indc].z);
            np[14]=setv(proseq[bb[k].indc].x,proseq[bb[k].indc].y,proseq[bb[k].indc].z);
            np[15]=setv(proseq[bb[k+1].indc].x,proseq[bb[k+1].indc].y,proseq[bb[k+1].indc].z);
    
            np[16]=minu(tp[0],np[7]);
            np[17]=minu(np[5],np[7]);
            np[18]=minu(ap[0],np[13]);
            np[19]=minu(np[11],np[13]);
            np[20]=prod(np[16],np[17]);
            np[20]=unit(np[20]);
            np[21]=prod(np[18],np[19]);
            np[21]=unit(np[21]);
            
            np[26]=minu(tp[1],np[8]);
            np[27]=minu(np[6],np[8]);
            np[28]=minu(ap[1],np[14]);
            np[29]=minu(np[12],np[14]);
            np[30]=prod(np[26],np[27]);
            np[30]=unit(np[30]);
            np[31]=prod(np[28],np[29]);
            np[31]=unit(np[31]);
            
            np[32]=minu(np[5],np[11]);//nn
            np[33]=minu(np[5],np[14]);//nc
            np[34]=minu(np[8],np[11]);//cn
            np[35]=minu(np[8],np[14]);//cc



            diss[0]=norm(tp[6]);//dist caca-1
            diss[1]=norm(tp[7]);//dist caca 
            diss[2]=norm(tp[8]);//dist caca+1            
            diss[3]=norm(tp[9]);//dist caca-1 +1
            diss[4]=norm(tp[10]);//dist caca+1 -1
            diss[5]=angv(tp[5],ap[5]);//plane plane good
            diss[6]=angv(tp[5],tp[10]);//plane1 cn good 
            diss[7]=angv(tp[11],tp[12]);//caca caca
            diss[8]=angv(tp[5],tp[7]);//plane1 caca      
            diss[9]=angv(tp[5],tp[9]);//dist nc
            diss[10]=angv(ap[5],tp[9]);//plane2 nc
            diss[11]=angv(tp[11],tp[7]);//ca1 caca  
            diss[12]=angv(tp[12],tp[7]);//ca2 caca  

            ap[11]=setv(proseq[bb[j+2].indca].x,proseq[bb[j+2].indca].y,proseq[bb[j+2].indca].z);
            ap[12]=setv(proseq[bb[k-2].indca].x,proseq[bb[k-2].indca].y,proseq[bb[k-2].indca].z);
            diss[13]=phi(tp[0].x,tp[0].y,tp[0].z,tp[1].x,tp[1].y,tp[1].z,tp[2].x,tp[2].y,tp[2].z,
                ap[11].x,ap[11].y,ap[11].z);//phi 

            diss[14]=phi(ap[12].x,ap[12].y,ap[12].z,ap[0].x,ap[0].y,ap[0].z,ap[1].x,ap[1].y,ap[1].z,
                ap[2].x,ap[2].y,ap[2].z);//phi 

            diss[15]=angv(tp[3],tp[4]);
            diss[16]=angv(ap[3],ap[4]);
            diss[17]=angv(ap[5],tp[7]);//plane2 nc
            
            diss[18]=angv(np[20],np[21]);// 4 pepetide plane
            diss[19]=angv(np[20],np[31]);// 4 pepetide plane
            diss[20]=angv(np[30],np[21]);// 4 pepetide plane
            diss[21]=angv(np[30],np[31]);// 4 pepetide plane
            diss[22]=norm(np[0]);//dist
            diss[23]=norm(np[1]);//dist
            diss[24]=norm(np[2]);//dist
            diss[25]=norm(np[3]);//dist      
            diss[26]=norm(np[32]);//dist
            diss[27]=norm(np[33]);//dist
            diss[28]=norm(np[34]);//dist
            diss[29]=norm(np[35]);//dist

            int delseg=-1;
            if(diss[7]<PI/2.0 && diss[8]<PI/2.0)//parallel1 left
            {
                delseg=0;
            }
            else if(diss[7]<PI/2.0 && diss[8]>=PI/2.0)//parallel2 right
            {
                delseg=1;
            }
            else if(diss[7]>=PI/2.0 && diss[8]<PI/2.0)//anti1 left
            {
                delseg=2;
            }
            else if(diss[7]>=PI/2.0 && diss[8]>=PI/2.0)//anti2 right
            {
                delseg=3;
            }
            
            if(delseg==0 && (diss[15]<1.57 || diss[15]>2.70)) continue;
            if(delseg==1 && (diss[15]<1.57 || diss[15]>2.70)) continue;
            if(delseg==2 && (diss[15]<1.57 || diss[15]>2.70)) continue;
            if(delseg==3 && (diss[15]<1.57 || diss[15]>2.70)) continue;
            
            if(delseg==0 && (diss[16]<1.57 || diss[16]>2.70)) continue;
            if(delseg==1 && (diss[16]<1.57 || diss[16]>2.70)) continue;
            if(delseg==2 && (diss[16]<1.57 || diss[16]>2.70)) continue;
            if(delseg==3 && (diss[16]<1.57 || diss[16]>2.70)) continue;//exclude false

            if(delseg==0 && (diss[1]<3.8 || diss[1]>6.2)) continue;
            if(delseg==1 && (diss[1]<3.8 || diss[1]>6.2)) continue;
            if(delseg==2 && (diss[1]<3.2 || diss[1]>6.0)) continue;
            if(delseg==3 && (diss[1]<3.4 || diss[1]>6.4)) continue;

            if(delseg==0 && (diss[3]<6.6 || diss[3]>10.0)) continue;
            if(delseg==1 && (diss[4]<6.6 || diss[4]>9.8)) continue;
            if(delseg==2 && (diss[0]<6.4 || diss[0]>9.7)) continue;
            if(delseg==3 && (diss[2]<6.4 || diss[2]>9.8)) continue;

            if(delseg==0 && ((diss[13]>30 &&diss[13]<100) || diss[13]>360)) continue;
            if(delseg==1 && ((diss[13]>50 &&diss[13]<100) || diss[13]>310)) continue;
            if(delseg==2 && ((diss[13]>0  &&diss[13]<80)  || diss[13]>360)) continue;
            if(delseg==3 && ((diss[13]>50 &&diss[13]<100) || diss[13]>360)) continue;

            if(delseg==0 && (diss[5]>1.0)) continue;
            if(delseg==1 && (diss[5]>1.2)) continue;
            if(delseg==2 && (diss[5]<1.6)) continue;
            if(delseg==3 && (diss[5]<1.0)) continue;


            double tval=1;
        
            if(delseg==0)
            {
                for(l=0;l<30;l++) //if(l!=5  && l!=7)
                //  if( l==0 || l==1 || l==2 || l==3 || l==4 || l==8 || l==17 || l==26 || l==29 || l==27 || l==28)
                    if(l==11 || l==12 || l==0 || l==1 || l==2 || l==3  || l==26 || l==29 || l==27 || l==28 || l==8 || l==17)
                {
                    tval*=expgaussian2(diss[l],bsigma[delseg][l],bmean[delseg][l]);
                }
            }
            else if(delseg==1)
            {
                for(l=0;l<30;l++)  //if(l!=5  && l!=7)
                //  if( l==0 || l==1 || l==2 || l==3 || l==4 || l==8 || l==17 || l==26 || l==29 || l==27 || l==28)
                    if(l==11 || l==12 || l==0 || l==1 || l==2   || l==4  || l==26 || l==29 || l==27 || l==28 || l==8 || l==17)
                {
                    tval*=expgaussian2(diss[l],bsigma[delseg][l],bmean[delseg][l]);
                }
            }
            else if(delseg==2)
            {
                for(l=0;l<30;l++) //if( l!=5 && l!=7)
                //  if( l==3 || l==1 || l==4 || l==0 || l==2 || l==8 || l==17 || l==27 || l==28 || l==26 || l==29)
                    if(l==11 || l==12 || l==3 || l==1 || l==4 || l==0  || l==26 || l==29 || l==27 || l==28 || l==8 || l==17)
                {
                    tval*=expgaussian2(diss[l],bsigma[delseg][l],bmean[delseg][l]);
                }
            }
            else if(delseg==3)
            {
                for(l=0;l<30;l++) //if(l!=5 && l!=6 && l!=7 && l!=10)
                //  if( l==3 || l==1 || l==4 || l==0 || l==2 || l==8 || l==17 || l==27 || l==28 || l==26 || l==29)
                    if(l==11 || l==12 || l==3 || l==1 || l==4  || l==2  || l==26 || l==29 || l==27 || l==28 || l==8 || l==17)
                {
                    tval*=expgaussian2(diss[l],bsigma[delseg][l],bmean[delseg][l]);
                }
            }
            //24 23 76 75
//          if((delseg==0 && tval>2e-4)|| (delseg==1 && tval>1e-4)||(delseg==2 && tval>2e-4)|| (delseg==3 && tval>3e-4))
//          if((delseg==0 && tval>1e-5)|| (delseg==1 && tval>3e-6)||(delseg==2 && tval>2e-5)|| (delseg==3 && tval>6e-5))
//          if( (delseg==0 && tval>5e-6)|| (delseg==1 && tval>5e-7)||(delseg==2 && tval>1e-5)|| (delseg==3 && tval>3e-5))
            if( (delseg==0 && tval>4e-8)|| (delseg==1 && tval>2e-8)||(delseg==2 && tval>4e-6)|| (delseg==3 && tval>1e-6))
            {
                bb[j].sst=2;
                bb[k].sst=2;
                btscore[j*numbb+k]=log(tval);
                btscore[k*numbb+j]=btscore[j*numbb+k];
//              for(l=0;l<30;l++)
//                  if( l==3 || l==1 || l==4 || l==0 || l==2 || l==8 || l==17 || l==27 || l==28 || l==26 || l==29) 
//                      fprintf(file,"%2d %7.3f ",l,diss[l]);
//              fprintf(file,"\n%2d %2d %2d %f\n",j+1,k+1,delseg,btscore[j*numbb+k]);
            }
        }
    }
//  fclose(file);

    ////
    point3d ptp[2];
    double pr;
//  for(i=0;i<numbb-2;i++)
//  {
//      if( bb[i].indca==-1 || bb[i+1].indca==-1 || bb[i+2].indca==-1 || bb[i+1].sst==1 || bb[i+1].sst==2) continue;
//      ptp[0]=bf.setv(proseq[bb[i+2].indca].x-proseq[bb[i].indca].x,
//          proseq[bb[i+2].indca].y-proseq[bb[i].indca].y,
//          proseq[bb[i+2].indca].z-proseq[bb[i].indca].z);
//      pr=bf.norm(ptp[0]);
//      if(pr<=7.0)
//      {
//          bb[i+1].sst=3;
//      }
//  }
    for(i=0;i<numbb-3;i++)
    {
        if( bb[i].indca==-1 || bb[i+1].indca==-1 || bb[i+2].indca==-1 || bb[i+3].indca==-1
            || bb[i+1].sst==1 || bb[i+1].sst==2 || bb[i+2].sst==1 || bb[i+2].sst==2) continue;
        ptp[0]=setv(proseq[bb[i+3].indca].x-proseq[bb[i].indca].x,
            proseq[bb[i+3].indca].y-proseq[bb[i].indca].y,
            proseq[bb[i+3].indca].z-proseq[bb[i].indca].z);
        pr=norm(ptp[0]);
        if(pr<=7.0)
        {
            bb[i+1].sst=3;
            bb[i+2].sst=3;
        }
    }
    for(i=0;i<numbb-4;i++)
    {
        if( bb[i].indca==-1 || bb[i+1].indca==-1 || bb[i+2].indca==-1 || bb[i+3].indca==-1 || bb[i+4].indca==-1
            || bb[i+1].sst==1 || bb[i+1].sst==2 || bb[i+2].sst==1 || bb[i+2].sst==2
            || bb[i+3].sst==1 || bb[i+3].sst==2) continue;
        ptp[0]=setv(proseq[bb[i+4].indca].x-proseq[bb[i].indca].x,
            proseq[bb[i+4].indca].y-proseq[bb[i].indca].y,
            proseq[bb[i+4].indca].z-proseq[bb[i].indca].z);
        pr=norm(ptp[0]);
        if(pr<=7.0)
        {
            bb[i+1].sst=3;
            bb[i+2].sst=3;
            bb[i+3].sst=3;
        }
    }
    for(i=0;i<numbb-5;i++)
    {
        if( bb[i].indca==-1 || bb[i+1].indca==-1 || bb[i+2].indca==-1 || bb[i+3].indca==-1 || bb[i+4].indca==-1
            || bb[i+5].indca==-1 || bb[i+1].sst==1 || bb[i+1].sst==2 || bb[i+2].sst==1 || bb[i+2].sst==2
            || bb[i+3].sst==1 || bb[i+3].sst==2 || bb[i+4].sst==1 || bb[i+4].sst==2) continue;
        ptp[0]=setv(proseq[bb[i+5].indca].x-proseq[bb[i].indca].x,
            proseq[bb[i+5].indca].y-proseq[bb[i].indca].y,
            proseq[bb[i+5].indca].z-proseq[bb[i].indca].z);
        pr=norm(ptp[0]);
        if(pr<=7.0)
        {
            bb[i+1].sst=3;
            bb[i+2].sst=3;
            bb[i+3].sst=3;
            bb[i+4].sst=3;
        }
    }
}  */

void str2tor(vector<point3f> &decstr,int seqnum,int type)
{
    int i;
    point3d p12,p23;
//    BasicFunc bf;
    if(type==1 || type==6)
    {
        decstr[0].leng=-1.0;
        decstr[0].phi =-1.0;        
        if(decstr[0].leng<0 || decstr[0].phi<0)
        {
            if(type==1)
            {
                decstr[0].leng=lennc;
                decstr[0].angl=angcacn;
                decstr[1].angl=angcnca;
            }
            else
            {
                decstr[0].leng=lencaca;
                decstr[0].angl=106.422f;
                decstr[1].angl=106.422f;
            }
            decstr[0].phi=180.0;
            decstr[1].phi=180.0;
            decstr[2].phi=180.0;
        }

        p12=setv(decstr[0].x-decstr[1].x,decstr[0].y-decstr[1].y,decstr[0].z-decstr[1].z);
        p23=setv(decstr[2].x-decstr[1].x,decstr[2].y-decstr[1].y,decstr[2].z-decstr[1].z);
        decstr[1].leng=norm(p12);
        decstr[2].leng=norm(p23);
        decstr[2].angl=angv(p12,p23);
        for(i=3;i<seqnum;i++)
        {
            p12=setv(decstr[i-2].x-decstr[i-1].x,decstr[i-2].y-decstr[i-1].y,decstr[i-2].z-decstr[i-1].z);
            p23=setv(decstr[i].x-decstr[i-1].x,decstr[i].y-decstr[i-1].y,decstr[i].z-decstr[i-1].z);
            decstr[i].leng=norm(p23);
            decstr[i].angl=angv(p12,p23);
            decstr[i].phi=phi(decstr[i-3].x,decstr[i-3].y,decstr[i-3].z,decstr[i-2].x,decstr[i-2].y,decstr[i-2].z,
                decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,decstr[i].x,decstr[i].y,decstr[i].z);
        }
    }
    else if(type==3)
    {
        decstr[0].tor[0] = -1.0;
        decstr[0].len[0] = -1.0;
        if(decstr[0].tor[0]<0 || decstr[0].len[0]<0)
        {
            decstr[0].tor[0]=180.0;
            decstr[0].tor[1]=180.0;
            decstr[0].tor[2]=180.0;
            decstr[0].len[0]=lennc;
            decstr[0].ang[0]=angcacn;   
            decstr[0].ang[1]=angcnca;
        }
        
        p12=setv(decstr[0].ptn.x-decstr[0].x,decstr[0].ptn.y-decstr[0].y,decstr[0].ptn.z-decstr[0].z);
        p23=setv(decstr[0].ptc.x-decstr[0].x,decstr[0].ptc.y-decstr[0].y,decstr[0].ptc.z-decstr[0].z);
        decstr[0].len[1]=norm(p12);
        decstr[0].len[2]=norm(p23);
        decstr[0].ang[2]=angv(p12,p23)*degrad;
        for(i=1;i<seqnum;i++)
        {
            p12=setv(decstr[i-1].x-decstr[i-1].ptc.x,decstr[i-1].y-decstr[i-1].ptc.y,decstr[i-1].z-decstr[i-1].ptc.z);
            p23=setv(decstr[i].ptn.x-decstr[i-1].ptc.x,decstr[i].ptn.y-decstr[i-1].ptc.y,decstr[i].ptn.z-decstr[i-1].ptc.z);
            decstr[i].len[0]=norm(p23);
            decstr[i].ang[0]=angv(p12,p23)*degrad;
            decstr[i].tor[0]=phi(decstr[i-1].ptn.x,decstr[i-1].ptn.y,decstr[i-1].ptn.z,decstr[i-1].x,decstr[i-1].y,
            decstr[i-1].z,decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z);
            
            p12=setv(decstr[i-1].ptc.x-decstr[i].ptn.x,decstr[i-1].ptc.y-decstr[i].ptn.y,decstr[i-1].ptc.z-decstr[i].ptn.z);
            p23=setv(decstr[i].x-decstr[i].ptn.x,decstr[i].y-decstr[i].ptn.y,decstr[i].z-decstr[i].ptn.z);
            decstr[i].len[1]=norm(p23);
            decstr[i].ang[1]=angv(p12,p23)*degrad;
            decstr[i].tor[1]=phi(decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,decstr[i-1].ptc.x,decstr[i-1].ptc.y,
            decstr[i-1].ptc.z,decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,decstr[i].x,decstr[i].y,decstr[i].z);
        
            p12=setv(decstr[i].ptn.x-decstr[i].x,decstr[i].ptn.y-decstr[i].y,decstr[i].ptn.z-decstr[i].z);
            p23=setv(decstr[i].ptc.x-decstr[i].x,decstr[i].ptc.y-decstr[i].y,decstr[i].ptc.z-decstr[i].z);
            decstr[i].len[2]=norm(p23);
            decstr[i].ang[2]=angv(p12,p23)*degrad;
            decstr[i].tor[2]=phi(decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,decstr[i].ptn.x,decstr[i].ptn.y,
            decstr[i].ptn.z,decstr[i].x,decstr[i].y,decstr[i].z,decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z);
        }
    }
    
}  

void calcssennhoc(vector<boneinfo> &bb, int numbb, vector<point3f> &decstr)//3states
{
    int i,j,k,m;
//    BasicFunc bf;
//    ParseSeq ps;
    point3d tp[20],ap[20],kp[20];
    point3d pd,pd2,pd3;
    double lamda;
    double diss[20];
    for(i=0;i<numbb;i++)
    {
        bb[i].sst=0;
    }
/*    vector<point3f> decstr(numbb);
    for(j=0;j<numbb;j++)
    {
        decstr[j].ss2='C';
        decstr[j].aaa=bb[j].resid;
        decstr[j].ptn.x=proseq[bb[j].indn].x_[0];
        decstr[j].ptn.y=proseq[bb[j].indn].x_[1];
        decstr[j].ptn.z=proseq[bb[j].indn].x_[2];
        decstr[j].x=proseq[bb[j].indca].x_[0];
        decstr[j].y=proseq[bb[j].indca].x_[1];
        decstr[j].z=proseq[bb[j].indca].x_[2];
        decstr[j].ptc.x=proseq[bb[j].indc].x_[0];
        decstr[j].ptc.y=proseq[bb[j].indc].x_[1];
        decstr[j].ptc.z=proseq[bb[j].indc].x_[2];
    }  */
//    str2tor(decstr,numbb,3);
//    tor2stroh(decstr,numbb);

    static  double nhochhmean[22]={
        5.2172,8.7074,230.1397,4.2198,6.5391
        };
    static  double nhochhsigma[22]={
        0.3676,0.4375,10.2234,0.4080,0.3763
        };
    static  double nhocmeanval[][4]={//h-o h-o=c n-h-o n-h-o=c  
        2.85,89.0,110.5,199.5,//i+3
        2.00,147.0,159.0,160.0, //i+4  oldoh 
    //  2.83,89.0,110.0,201.5,//i+3 newoh
    //  2.00,148.0,159.0,155.0, //i+4

        2.00,155.0,164.0,180.0,//0
        2.00,155.0,164.0,180.0,//1
        2.00,151.0,163.0,192.0,//2
        2.00,151.0,163.0,192.0,//3

        };
    static  double nhocstdval[][4]={//h-o h-o=c n-h-o n-h-o=c 
        0.315504,7.697185,8.980366,7.932107,
        0.530859,10.582243,11.249764,25.360054,

        0.299730,11.770196,11.292558,68.955920,
        0.299730,11.770196,11.292558,68.955920,
        0.255088,12.376087,11.020081,69.165282,
        0.255088,12.376087,11.020081,69.165282,
        };    


    double threshval[5]={19.113828,19.113828,20.723266,20.723266,16.118096};
    for(i=0;i<numbb;i++)
    {
        decstr[i].vpos=0;
        decstr[i].tpos=0;
        decstr[i].ssm='C';
        decstr[i].indl=-1;
        decstr[i].indr=-1;
        decstr[i].tpr=-1;
        decstr[i].tpl=-1;
    }
    for(i=1;i<numbb-4;i++)
    {
        j=i+3;  
        lamda=0.2;
        tp[1]=setv(decstr[i].x,decstr[i].y,decstr[i].z);
        ap[1]=setv(decstr[j].x,decstr[j].y,decstr[j].z);
        tp[7]=minu(tp[1],ap[1]);
        diss[1]=norm(tp[7]); 
        if(diss[1]>6.6) continue;

        tp[0]=setv(decstr[i-1].x,decstr[i-1].y,decstr[i-1].z);
        tp[2]=setv(decstr[i+1].x,decstr[i+1].y,decstr[i+1].z);
        ap[0]=setv(decstr[j-1].x,decstr[j-1].y,decstr[j-1].z);
        ap[2]=setv(decstr[j+1].x,decstr[j+1].y,decstr[j+1].z);
        kp[0]=minu(tp[2],ap[0]);
        diss[4]=norm(kp[0]);
        if(diss[4]<3.4 || diss[4]>4.2) continue;
        diss[0]=norm(minu(tp[0],ap[0]));
 
        //nhoc
        pd.x=decstr[j].pth.x-decstr[i].pto.x;
        pd.y=decstr[j].pth.y-decstr[i].pto.y;
        pd.z=decstr[j].pth.z-decstr[i].pto.z;
        diss[5]=norm(pd);
        if(diss[5]>=5.0 || diss[5]<1.6) continue;
 

        pd2.x=decstr[i].ptc.x-decstr[i].pto.x;
        pd2.y=decstr[i].ptc.y-decstr[i].pto.y;
        pd2.z=decstr[i].ptc.z-decstr[i].pto.z;
        diss[6]=angv(pd,pd2)*degrad;
        if(diss[6]<70.0 || diss[6]>140.0) continue;

        pd3.x=decstr[j].pth.x-decstr[j].ptn.x;
        pd3.y=decstr[j].pth.y-decstr[j].ptn.y;
        pd3.z=decstr[j].pth.z-decstr[j].ptn.z;
        diss[7]=angv(pd,pd3)*degrad;
                
        diss[8]=phi(decstr[j].ptn.x,decstr[j].ptn.y,decstr[j].ptn.z,decstr[j].pth.x,decstr[j].pth.y,decstr[j].pth.z,
            decstr[i].pto.x,decstr[i].pto.y,decstr[i].pto.z,decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z);
        if(diss[8]<120.0 || diss[8]>280.0) continue;

        //second
        j=i+4;
        pd.x=decstr[j].pth.x-decstr[i].pto.x;
        pd.y=decstr[j].pth.y-decstr[i].pto.y;
        pd.z=decstr[j].pth.z-decstr[i].pto.z;
        diss[9]=norm(pd);
        if(diss[9]>=5.0 || diss[9]<1.6) continue;
 

        pd2.x=decstr[i].ptc.x-decstr[i].pto.x;
        pd2.y=decstr[i].ptc.y-decstr[i].pto.y;
        pd2.z=decstr[i].ptc.z-decstr[i].pto.z;
        diss[10]=angv(pd,pd2)*degrad;
        if(diss[10]<100.0 || diss[10]>170.0) continue;

        pd3.x=decstr[j].pth.x-decstr[j].ptn.x;
        pd3.y=decstr[j].pth.y-decstr[j].ptn.y;
        pd3.z=decstr[j].pth.z-decstr[j].ptn.z;
        diss[11]=angv(pd,pd3)*degrad;
        if(diss[11]<100.0) continue;
                
        diss[12]=phi(decstr[j].ptn.x,decstr[j].ptn.y,decstr[j].ptn.z,decstr[j].pth.x,decstr[j].pth.y,decstr[j].pth.z,
            decstr[i].pto.x,decstr[i].pto.y,decstr[i].pto.z,decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z);
        if(diss[12]<100.0 || diss[12]>230.0) continue;

        double tval=0;
        for(k=0;k<1;k++)
        {       
            tval+=squgaussian(diss[k],nhochhsigma[k],nhochhmean[k]);
        }
        for(k=0;k<4;k++)
        {
            tval+=squgaussian(diss[5+k],nhocstdval[0][k],nhocmeanval[0][k]);
        }
        for(k=0;k<4;k++)
        {
            tval+=squgaussian(diss[9+k],nhocstdval[1][k],nhocmeanval[1][k]);
        }
        tval=exp(tval);
    //  printf("%3d %.3f\n",i,tval);
        if(tval>1e-7)
        {
            if(tval>1e-3) tval=1e-3;
            tval=lamda*(threshval[4] + log(tval));
            for(k=0;k<=3;k++)
            {
                decstr[i+k].ssm='H';
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////beta
    int inda,indb,indc,indd,inde,indf;
    int delseg;
    
    for(j=1;j<numbb-4;j++) 
    {
        for(k=j+3;k<numbb-1;k++)
        {   
             
            lamda=0.2;
            tp[1]=setv(decstr[j].x,decstr[j].y,decstr[j].z);
            ap[1]=setv(decstr[k].x,decstr[k].y,decstr[k].z);
            tp[7]=minu(tp[1],ap[1]);
            if(norm(tp[7])>8.0) continue;
        
            inda=j;indb=k;
            tp[4]=setv(decstr[inda].ptn.x,decstr[inda].ptn.y,decstr[inda].ptn.z);
            tp[5]=setv(decstr[inda+1].ptn.x,decstr[inda+1].ptn.y,decstr[inda+1].ptn.z);
            tp[9]=setv(decstr[inda-1].pto.x,decstr[inda-1].pto.y,decstr[inda-1].pto.z);
            tp[10]=setv(decstr[inda].pto.x,decstr[inda].pto.y,decstr[inda].pto.z);
            ap[4]=setv(decstr[indb].ptn.x,decstr[indb].ptn.y,decstr[indb].ptn.z);
            ap[5]=setv(decstr[indb+1].ptn.x,decstr[indb+1].ptn.y,decstr[indb+1].ptn.z);
            ap[9]=setv(decstr[indb-1].pto.x,decstr[indb-1].pto.y,decstr[indb-1].pto.z);
            ap[10]=setv(decstr[indb].pto.x,decstr[indb].pto.y,decstr[indb].pto.z);
            kp[0]=minu(tp[9],ap[4]);//oi-1 nj
            kp[1]=minu(tp[5],ap[10]);//ni+1 oj
            kp[2]=minu(tp[4],ap[9]);//ni oj-1
            kp[3]=minu(tp[10],ap[5]);//oi nj+1
            kp[4]=minu(tp[9],ap[5]);//oi-1 nj+1
            kp[5]=minu(tp[5],ap[9]);//ni+1 oj-1
            kp[6]=minu(tp[4],ap[10]);//ni oj
            kp[7]=minu(tp[10],ap[4]);//oi nj
            for(m=0;m<8;m++)
            {
                diss[m]=norm(kp[m]);
            }
            delseg=0;
            diss[9]=diss[0]+diss[1];
            diss[8]=diss[2]+diss[3];
            if(diss[8]<diss[9])
            {
                delseg=1;
                diss[9]=diss[8];
            }
            diss[8]=diss[4]+diss[5];
            if(diss[8]<diss[9])
            {
                delseg=2;
                diss[9]=diss[8];
            }
            diss[8]=diss[6]+diss[7];
            if(diss[8]<diss[9])
            {
                delseg=3;
                diss[9]=diss[8];
            }
            if(diss[2*delseg]>6.8) continue;
            if(diss[2*delseg+1]>6.8) continue;
            if(delseg==0)
            {
                indc=k;indd=j-1;inde=j+1;indf=k;
            }
            else if(delseg==1)
            {
                indc=j;indd=k-1;inde=k+1;indf=j;
            }
            else if(delseg==3)
            {
                indc=k;indd=j;inde=j;indf=k;
            }
            else if(delseg==2)
            {
                indc=k+1;indd=j-1;inde=j+1;indf=k-1;
            }

            pd.x=decstr[indc].pth.x-decstr[indd].pto.x;
            pd.y=decstr[indc].pth.y-decstr[indd].pto.y;
            pd.z=decstr[indc].pth.z-decstr[indd].pto.z;
            diss[0]=norm(pd);
            pd2.x=decstr[indd].ptc.x-decstr[indd].pto.x;
            pd2.y=decstr[indd].ptc.y-decstr[indd].pto.y;
            pd2.z=decstr[indd].ptc.z-decstr[indd].pto.z;
            diss[1]=angv(pd,pd2)*degrad;
            pd3.x=decstr[indc].pth.x-decstr[indc].ptn.x;
            pd3.y=decstr[indc].pth.y-decstr[indc].ptn.y;
            pd3.z=decstr[indc].pth.z-decstr[indc].ptn.z;
            diss[2]=angv(pd,pd3)*degrad;
    //      diss[3]=bf.phi(decstr[indc].ptn.x,decstr[indc].ptn.y,decstr[indc].ptn.z,decstr[indc].pth.x,decstr[indc].pth.y,decstr[indc].pth.z,
    //          decstr[indd].pto.x,decstr[indd].pto.y,decstr[indd].pto.z,decstr[indd].ptc.x,decstr[indd].ptc.y,decstr[indd].ptc.z);
            ////////
            pd.x=decstr[inde].pth.x-decstr[indf].pto.x;
            pd.y=decstr[inde].pth.y-decstr[indf].pto.y;
            pd.z=decstr[inde].pth.z-decstr[indf].pto.z;
            diss[4]=norm(pd);                
            pd2.x=decstr[indf].ptc.x-decstr[indf].pto.x;
            pd2.y=decstr[indf].ptc.y-decstr[indf].pto.y;
            pd2.z=decstr[indf].ptc.z-decstr[indf].pto.z;
            diss[5]=angv(pd,pd2)*degrad;
            pd3.x=decstr[inde].pth.x-decstr[inde].ptn.x;
            pd3.y=decstr[inde].pth.y-decstr[inde].ptn.y;
            pd3.z=decstr[inde].pth.z-decstr[inde].ptn.z;
            diss[6]=angv(pd,pd3)*degrad;
    //      diss[7]=bf.phi(decstr[inde].ptn.x,decstr[inde].ptn.y,decstr[inde].ptn.z,decstr[inde].pth.x,decstr[inde].pth.y,decstr[inde].pth.z,
    //          decstr[indf].pto.x,decstr[indf].pto.y,decstr[indf].pto.z,decstr[indf].ptc.x,decstr[indf].ptc.y,decstr[indf].ptc.z);//n-h-o=c
    
            double tval=0;
            for(m=0;m<3;m++)
            {
                tval+=squgaussian(diss[m],nhocstdval[2+delseg][m],nhocmeanval[2+delseg][m]);
            }
            for(m=0;m<3;m++)
            {
                tval+=squgaussian(diss[4+m],nhocstdval[2+delseg][m],nhocmeanval[2+delseg][m]);
            }
            tval=exp(tval);

            if((delseg==0 && tval>5e-9)|| (delseg==1 && tval>5e-9)||(delseg==2 && tval>1e-9)|| (delseg==3 && tval>1e-9))
            {
                if(tval>1e-3) tval=1e-3;
                tval=lamda*(threshval[delseg]+log(tval))/threshval[delseg]*threshval[4];
                if(delseg==0 || delseg==2)//left of j
                {
                    if(tval>decstr[j].vpos)
                    {
                        decstr[j].vpos=tval;
                        if(decstr[j].indl>0) 
                        {       
                            if(decstr[decstr[j].indl].indr==j)//set to zero
                            {
                                decstr[decstr[j].indl].tpos=0;
                                decstr[decstr[j].indl].indr=-1;
                                if(decstr[decstr[j].indl].indl==-1)
                                    decstr[decstr[j].indl].ssm='C';
                            }
                        }
                        decstr[j].indl=k;
                        decstr[j].tpl=delseg;
                        decstr[k].ssm='E';
                        decstr[j].ssm='E';
                    }
                }
                else if(delseg==1 || delseg==3)//right of j
                {
                    if(tval>decstr[j].tpos)
                    {
                        decstr[j].tpos=tval;
                        if(decstr[j].indr>0) 
                        {
                            if(decstr[decstr[j].indr].indl==j)//set to zero
                            {
                                decstr[decstr[j].indr].vpos=0;
                                decstr[decstr[j].indr].indl=-1;
                                if(decstr[decstr[j].indr].indr==-1)
                                    decstr[decstr[j].indr].ssm='C';
                            }
                        }
                        decstr[j].indr=k;
                        decstr[j].tpr=delseg;
                        decstr[k].ssm='E';
                        decstr[j].ssm='E';
                    }
                }
                if(delseg==1 || delseg==2)//left of k
                {
                    if(tval>decstr[k].vpos)
                    {
                        decstr[k].vpos=tval;
                        if(decstr[k].indl>0) 
                        {
                            if(decstr[decstr[k].indl].indr==k)//set to zero
                            {
                                decstr[decstr[k].indl].tpos=0;
                                decstr[decstr[k].indl].indr=-1;
                                if(decstr[decstr[k].indl].indl==-1)
                                    decstr[decstr[k].indl].ssm='C';
                            }
                        }
                        decstr[k].indl=j;
                        decstr[k].tpl=delseg;
                        decstr[j].ssm='E';
                        decstr[k].ssm='E';
                    }
                }
                else if(delseg==0 || delseg==3)//right of k
                {
                    if(tval>decstr[k].tpos)
                    {
                        decstr[k].tpos=tval;
                        if(decstr[k].indr>0) 
                        {
                            if(decstr[decstr[k].indr].indl==k)//set to zero
                            {
                                decstr[decstr[k].indr].vpos=0;
                                decstr[decstr[k].indr].indl=-1;
                                if(decstr[decstr[k].indr].indr==-1)
                                    decstr[decstr[k].indr].ssm='C';
                            }
                        }
                        decstr[k].indr=j;
                        decstr[k].tpr=delseg;
                        decstr[j].ssm='E';
                        decstr[k].ssm='E';
                    }
                }   
            }
        }
    }
    for(i=0;i<numbb;i++)
    {
        if(decstr[i].ssm=='E') bb[i].sst=2;
        else if(decstr[i].ssm=='H') bb[i].sst=1;
    }
//    delete[]decstr;
}  


char numtoss(int num)
{
    if(num==1) return 'H';
    if(num==2) return 'E';
//    if(num==3) return 'T';
    else
        return 'C';
}

void getLinLen(int Nlen, vector<point3f> xyzCa, vector<int> &a)
{
    int LinLen =0;
    int Ni = Nlen - 1;  //the index of anchor in the N-terminal
    int Ci = Nlen;      //the index of anchor in the C-terminal
    vector<float> tmp_A(3,0.0);
    tmp_A[0] = xyzCa[Ni].x;
    tmp_A[1] = xyzCa[Ni].y;
    tmp_A[2] = xyzCa[Ni].z;
    vector<float> tmp_B(3,0.0);
    tmp_B[0] = xyzCa[Ci].x;
    tmp_B[1] = xyzCa[Ci].y;
    tmp_B[2] = xyzCa[Ci].z;    
    float disanchor = sqrt(dist(tmp_A,tmp_B));
    while(disanchor >=(LinLen +1)*3.8*0.92)
    {
        if(LinLen > 20)
            break;
        if(Ni<=0)
        {
            Ni=0;
            cout<<"Ni = 0"<<endl;
            break;
        }
        if(Ci>=(xyzCa.size()-1)) 
        {
            Ci = xyzCa.size()-1;
            cout<<"Ci = end"<<endl;
            break;
        }
        if(numtoss(xyzCa[Ni-1].ssm)== numtoss(xyzCa[Ci+1].ssm))
        {
            float rand=(randf0and1()*2.0-1.0);
            if(rand>=0) 
                Ni = Ni -1 ;
            else 
                Ci = Ci + 1 ;
            LinLen = LinLen + 1 ;
        } else
        {
            if((numtoss(xyzCa[Ni-1].ssm) == 'H'&& numtoss(xyzCa[Ci+1].ssm)=='E') || (numtoss(xyzCa[Ni-1].ssm) == 'E' && numtoss(xyzCa[Ci+1].ssm)=='H'))
            {
                float rand=(randf0and1()*2.0-1.0);
                if(rand>=0) 
                    Ni = Ni -1 ;
                else 
                    Ci = Ci +1 ;
                LinLen = LinLen + 1 ;                
            } else
            {
                if(numtoss(xyzCa[Ni-1].ssm) == 'C')
                {
                    Ni = Ni -1 ;
                    LinLen = LinLen + 1 ;
                } 
                else if(numtoss(xyzCa[Ci+1].ssm)=='C')
                {
                    Ci = Ci + 1 ;
                    LinLen = LinLen + 1 ;
                }
            }
        }
        tmp_A[0] = xyzCa[Ni].x;
        tmp_A[1] = xyzCa[Ni].y;
        tmp_A[2] = xyzCa[Ni].z;

        tmp_B[0] = xyzCa[Ci].x;
        tmp_B[1] = xyzCa[Ci].y;
        tmp_B[2] = xyzCa[Ci].z;          
        disanchor = sqrt(dist(tmp_A,tmp_B));
    }
    a=vector<int>(2,0);
    a[0] = Ni;
    a[1] = Ci;
}

bool GroupRotationpid(const vector<float> &axisA, const vector<float> &axisB, float angle, vector<point3f> &pointB,int index0,int indexn)
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
   //C
    point_A[0]=pointB[index0].ptc.x-axisA[0];
    point_A[1]=pointB[index0].ptc.y-axisA[1];
    point_A[2]=pointB[index0].ptc.z-axisA[2];

    pointB[index0].ptc.x=axisA[0];
    pointB[index0].ptc.y=axisA[1];
    pointB[index0].ptc.z=axisA[2];
    pointB[index0].ptc.x+=rotmtx[0][0]*point_A[0]+rotmtx[0][1]*point_A[1]+rotmtx[0][2]*point_A[2]+rotmtx[0][3];
    pointB[index0].ptc.y+=rotmtx[1][0]*point_A[0]+rotmtx[1][1]*point_A[1]+rotmtx[1][2]*point_A[2]+rotmtx[1][3];
    pointB[index0].ptc.z+=rotmtx[2][0]*point_A[0]+rotmtx[2][1]*point_A[1]+rotmtx[2][2]*point_A[2]+rotmtx[2][3];  

   for(int i=index0+1; i<indexn; ++i)
   {
    //N
    point_A[0]=pointB[i].ptn.x-axisA[0];
    point_A[1]=pointB[i].ptn.y-axisA[1];
    point_A[2]=pointB[i].ptn.z-axisA[2];

    pointB[i].ptn.x=axisA[0];
    pointB[i].ptn.y=axisA[1];
    pointB[i].ptn.z=axisA[2];
    pointB[i].ptn.x+=rotmtx[0][0]*point_A[0]+rotmtx[0][1]*point_A[1]+rotmtx[0][2]*point_A[2]+rotmtx[0][3];
    pointB[i].ptn.y+=rotmtx[1][0]*point_A[0]+rotmtx[1][1]*point_A[1]+rotmtx[1][2]*point_A[2]+rotmtx[1][3];
    pointB[i].ptn.z+=rotmtx[2][0]*point_A[0]+rotmtx[2][1]*point_A[1]+rotmtx[2][2]*point_A[2]+rotmtx[2][3]; 

    //CA
    point_A[0]=pointB[i].x-axisA[0];
    point_A[1]=pointB[i].y-axisA[1];
    point_A[2]=pointB[i].z-axisA[2];

    pointB[i].x=axisA[0];
    pointB[i].y=axisA[1];
    pointB[i].z=axisA[2];
    pointB[i].x+=rotmtx[0][0]*point_A[0]+rotmtx[0][1]*point_A[1]+rotmtx[0][2]*point_A[2]+rotmtx[0][3];
    pointB[i].y+=rotmtx[1][0]*point_A[0]+rotmtx[1][1]*point_A[1]+rotmtx[1][2]*point_A[2]+rotmtx[1][3];
    pointB[i].z+=rotmtx[2][0]*point_A[0]+rotmtx[2][1]*point_A[1]+rotmtx[2][2]*point_A[2]+rotmtx[2][3]; 

    //C
    point_A[0]=pointB[i].ptc.x-axisA[0];
    point_A[1]=pointB[i].ptc.y-axisA[1];
    point_A[2]=pointB[i].ptc.z-axisA[2];

    pointB[i].ptc.x=axisA[0];
    pointB[i].ptc.y=axisA[1];
    pointB[i].ptc.z=axisA[2];
    pointB[i].ptc.x+=rotmtx[0][0]*point_A[0]+rotmtx[0][1]*point_A[1]+rotmtx[0][2]*point_A[2]+rotmtx[0][3];
    pointB[i].ptc.y+=rotmtx[1][0]*point_A[0]+rotmtx[1][1]*point_A[1]+rotmtx[1][2]*point_A[2]+rotmtx[1][3];
    pointB[i].ptc.z+=rotmtx[2][0]*point_A[0]+rotmtx[2][1]*point_A[1]+rotmtx[2][2]*point_A[2]+rotmtx[2][3]; 
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
    //N
    point_A[0]=pointB[indexn].ptn.x-axisA[0];
    point_A[1]=pointB[indexn].ptn.y-axisA[1];
    point_A[2]=pointB[indexn].ptn.z-axisA[2];

    pointB[indexn].ptn.x=axisA[0];
    pointB[indexn].ptn.y=axisA[1];
    pointB[indexn].ptn.z=axisA[2];
    pointB[indexn].ptn.x+=rotmtx[0][0]*point_A[0]+rotmtx[0][1]*point_A[1]+rotmtx[0][2]*point_A[2]+rotmtx[0][3];
    pointB[indexn].ptn.y+=rotmtx[1][0]*point_A[0]+rotmtx[1][1]*point_A[1]+rotmtx[1][2]*point_A[2]+rotmtx[1][3];
    pointB[indexn].ptn.z+=rotmtx[2][0]*point_A[0]+rotmtx[2][1]*point_A[1]+rotmtx[2][2]*point_A[2]+rotmtx[2][3];    

   return true;
}


void randomwalk(vector<vector<float> > &LinxyzCa, int istart, int iend)
{
//    vector<vector<float> > LinxyzCa;
    int L = LinxyzCa.size(); 
//    int i1=1,i2 = LinxyzCa.size()-2;
    int i1=istart+1,i2 = iend-1;

    //calculate the interval
    float diss = sqrt(dist(LinxyzCa[i2+1],LinxyzCa[i1-1]));
    float adis = diss/(float)((i2+1)-(i1-1));

    // linear connect the gap
    float adis0 = 3.5;
    float bdis0 =3.5;
    float amcheck_dis = 3.1;
    int x_check = 0;
//    double *dexyz = new double[3];
    vector<float> dexyz(3,0.0);
    int flag15_0 = 1;
    int flag13_0 = 1;
    
    bool flag15 =true;
    while (flag15 == true){ 
        x_check ++;
        if(adis >= adis0){
            for(int j=0; j<3; j++)
                dexyz[j] = 3.5*(LinxyzCa[i1-1][j] - LinxyzCa[i2+1][j])/diss;
            for(int i = i2; i>=i1; i--){
                for(int j=0; j<3; j++)
                    LinxyzCa[i][j] = LinxyzCa[i+1][j] + dexyz[j];
            }
        }
        
        //random walk from i1 to i2
        int m_check = 0;
        bool flag13 = true;
        while (flag13 == true){
            if(flag15_0 == 0)
                break;
            m_check ++;
            for(int i= i1; i<=i2; i++){
                if(flag13_0 == 0 || flag15_0 == 0)
                    break;
                int n_check = 0;
                bool flag12 = true;     
                while(flag12 == true){
                    float *axyz = get_bond(3.8);
                    for(int j=0; j<3; j++)
                        LinxyzCa[i][j] = LinxyzCa[i-1][j] + axyz[j];
                    n_check ++;
                    
                    if((int)(n_check/1000)*1000 == n_check){
                        bdis0 = bdis0*1.03;
                        if(bdis0 >= 3.8)
                            bdis0 = 3.8;
                        amcheck_dis = 2.8;
                        if(i1 == i2)
                            amcheck_dis = 2.0;
                    }
                    flag13 = false;
                    if(n_check == 4000){
                        flag13 = true;
                        flag13_0 = 0;
                        break;
                    }
                    flag15 = false;
                    if(m_check >= 2000){
                        if(x_check <= 6){
                            adis0 = adis0*0.995;
                            flag15 = true;
                            flag15_0 = 0;
                            break;
                        }
                    }
                    flag12 = false;
                    if(mcheck(i, LinxyzCa, amcheck_dis) == 3){
                        flag12 = true;
                        continue;
                    }
                    
                    flag12 = false;
                    float bdis = sqrt(dist(LinxyzCa[i], LinxyzCa[i2+1])/(float)(i2+1-i));
                    if(i < i2 && bdis >= bdis0){
                        flag12 = true;
                        continue;
                    }

                    flag12 = false;                 
                    if(i == i2){
                        if(bdis > 4.2 || bdis < 3.4){
                            flag12 = true;
                            continue;
                        }
                    }
                    bdis0 = 3.5;
                    amcheck_dis =3.1;
                }       
            }
        }
    }   
    
//    int tmp=0;
//    for(int i=0;i<L;i++)
//    {
//        if(LinChain[i].elt_ == "CA")
//        {
//            LinChain[i].x_=LinxyzCa[tmp];
//            tmp=tmp+1;
//        }
//    }    
//    vector<float *>().swap(LinxyzCa);  
//    vector<vector<float>>().swap(LinxyzCa); 
    return;

} 

//generate vector with length equa to al
float *get_bond(float al){
    float *xyz = new float[3];
    float athita = acos(1.0-2.0*randf0and1()); //thita angle in random, [0,pi]
    float aphi = 2.0*3.1415926*randf0and1(); //phi angle in random, [0,2pi]
    xyz[0] = al*sin(athita)*cos(aphi);
    xyz[1] = al*sin(athita)*sin(aphi);
    xyz[2] = al*cos(athita);
    return xyz;
}
//generate vector with length equa to al
vector<float> get_bondx(float al){
    vector<float> xyz(3,0.0);
    float athita = acos(1.0-2.0*randf0and1()); //thita angle in random, [0,pi]
    float aphi = 2.0*3.1415926*randf0and1(); //phi angle in random, [0,2pi]
    xyz[0] = al*sin(athita)*cos(aphi);
    xyz[1] = al*sin(athita)*sin(aphi);
    xyz[2] = al*cos(athita);
    return xyz;
}

//check the distance
int mcheck(int i0, vector<vector<float>> LinxyzCa, int amcheck_dis){
    int check = 1;
    for(int i =0; i<LinxyzCa.size(); i++){
        if(i0 != i){
            double dis = sqrt(dist(LinxyzCa[i0], LinxyzCa[i]));
            if (dis <= amcheck_dis)
                check =3;
        }
    }
    return check;
}

void singlemoveLMPf(vector<point3f> &tmstr,int k,int mtype)
{
//  BasicFunc bf;
    point3d rp;
    double rnorm;
    switch(mtype)
    {
        case 0:  
                //1 c-1 n
                rp=setv(tmstr[k].ptn.x-tmstr[k-1].ptc.x,tmstr[k].ptn.y-tmstr[k-1].ptc.y,tmstr[k].ptn.z-tmstr[k-1].ptc.z);
                rnorm=norm(rp);
                if(rnorm<epsilon)//two in the same place
                {
                    tmstr[k].ptn.x+=0.0001f;tmstr[k].ptn.y+=0.0001f;tmstr[k].ptn.z+=0.0001f;
                    rp=setv(tmstr[k].ptn.x-tmstr[k-1].ptc.x,tmstr[k].ptn.y-tmstr[k-1].ptc.y,tmstr[k].ptn.z-tmstr[k-1].ptc.z);
                    rnorm=norm(rp);  
                }
                if(rnorm<lencn-delcn || rnorm>lencn+delcn)
                {
                    rp=scalx(rp,lencn/rnorm);
                    tmstr[k].ptn.x=tmstr[k-1].ptc.x+rp.x;
                    tmstr[k].ptn.y=tmstr[k-1].ptc.y+rp.y;
                    tmstr[k].ptn.z=tmstr[k-1].ptc.z+rp.z;
                } 
            break;
        case 1:  
                //2 ca-1 n
                rp=setv(tmstr[k].ptn.x-tmstr[k-1].x,tmstr[k].ptn.y-tmstr[k-1].y,tmstr[k].ptn.z-tmstr[k-1].z);
                rnorm=norm(rp);
                if(rnorm<epsilon)//two in the same place
                {
                    tmstr[k].ptn.x+=0.0001f;tmstr[k].ptn.y+=0.0001f;tmstr[k].ptn.z+=0.0001f;
                    rp=setv(tmstr[k].ptn.x-tmstr[k-1].x,tmstr[k].ptn.y-tmstr[k-1].y,tmstr[k].ptn.z-tmstr[k-1].z);
                    rnorm=norm(rp);  
                }
                if(rnorm<lencan1-delcan1 || rnorm>lencan1+delcan1)
                {
                    rp=scalx(rp,lencan1/rnorm);
                    tmstr[k].ptn.x=tmstr[k-1].x+rp.x;
                    tmstr[k].ptn.y=tmstr[k-1].y+rp.y;
                    tmstr[k].ptn.z=tmstr[k-1].z+rp.z;
                }
            break;
        case 2:  
                //4 c-1 ca
                rp=setv(tmstr[k].x-tmstr[k-1].ptc.x,tmstr[k].y-tmstr[k-1].ptc.y,tmstr[k].z-tmstr[k-1].ptc.z);
                rnorm=norm(rp);
                if(rnorm<epsilon)//two in the same place
                {
                    tmstr[k].x+=0.0001f;tmstr[k].y+=0.0001f;tmstr[k].z+=0.0001f;
                    rp=setv(tmstr[k].x-tmstr[k-1].ptc.x,tmstr[k].y-tmstr[k-1].ptc.y,tmstr[k].z-tmstr[k-1].ptc.z);
                    rnorm=norm(rp);  
                }
                if(rnorm<lencca1-delcca1 || rnorm>lencca1+delcca1)
                {
                    rp=scalx(rp,lencca1/rnorm);
                    tmstr[k].x=tmstr[k-1].ptc.x+rp.x;
                    tmstr[k].y=tmstr[k-1].ptc.y+rp.y;
                    tmstr[k].z=tmstr[k-1].ptc.z+rp.z;
                }
            break;
        case 3:  
                //5 ca-1 ca
                rp=setv(tmstr[k].x-tmstr[k-1].x,tmstr[k].y-tmstr[k-1].y,tmstr[k].z-tmstr[k-1].z);
                rnorm=norm(rp);
                if(rnorm<epsilon)//two in the same place
                {
                    tmstr[k].x+=0.0001f;tmstr[k].y+=0.0001f;tmstr[k].z+=0.0001f;
                    rp=setv(tmstr[k].x-tmstr[k-1].x,tmstr[k].y-tmstr[k-1].y,tmstr[k].z-tmstr[k-1].z);
                    rnorm=norm(rp);  
                }
                if(rnorm<lencaca-delcaca || rnorm>lencaca+delcaca)
                {
                    rp=scalx(rp,lencaca/rnorm);
                    tmstr[k].x=tmstr[k-1].x+rp.x;
                    tmstr[k].y=tmstr[k-1].y+rp.y;
                    tmstr[k].z=tmstr[k-1].z+rp.z;
                }
            break;
        case 4:  
                //3 n ca
                rp=setv(tmstr[k].x-tmstr[k].ptn.x,tmstr[k].y-tmstr[k].ptn.y,tmstr[k].z-tmstr[k].ptn.z);
                rnorm=norm(rp);
                if(rnorm<epsilon)//two in the same place
                {
                    tmstr[k].x+=0.0001f;tmstr[k].y+=0.0001f;tmstr[k].z+=0.0001f;
                    rp=setv(tmstr[k].x-tmstr[k].ptn.x,tmstr[k].y-tmstr[k].ptn.y,tmstr[k].z-tmstr[k].ptn.z);
                    rnorm=norm(rp);  
                }
                if(rnorm<lennca-delnca || rnorm>lennca+delnca)
                {
                    rp=scalx(rp,lennca/rnorm);
                    tmstr[k].x=tmstr[k].ptn.x+rp.x;
                    tmstr[k].y=tmstr[k].ptn.y+rp.y;
                    tmstr[k].z=tmstr[k].ptn.z+rp.z;
                }
            break;
        case 5:  
                //6 ca c
                rp=setv(tmstr[k].ptc.x-tmstr[k].x,tmstr[k].ptc.y-tmstr[k].y,tmstr[k].ptc.z-tmstr[k].z);
                rnorm=norm(rp);
                if(rnorm<epsilon)//two in the same place
                {
                    tmstr[k].ptc.x+=0.0001f;tmstr[k].ptc.y+=0.0001f;tmstr[k].ptc.z+=0.0001f;
                    rp=setv(tmstr[k].ptc.x-tmstr[k].x,tmstr[k].ptc.y-tmstr[k].y,tmstr[k].ptc.z-tmstr[k].z);
                    rnorm=norm(rp);  
                }
                if(rnorm<lencac-delcac || rnorm>lencac+delcac)
                {
                    rp=scalx(rp,lencac/rnorm);
                    tmstr[k].ptc.x=tmstr[k].x+rp.x;
                    tmstr[k].ptc.y=tmstr[k].y+rp.y;
                    tmstr[k].ptc.z=tmstr[k].z+rp.z;
                }
            break;
        case 6:  
                //7 n c
                rp=setv(tmstr[k].ptc.x-tmstr[k].ptn.x,tmstr[k].ptc.y-tmstr[k].ptn.y,tmstr[k].ptc.z-tmstr[k].ptn.z);
                rnorm=norm(rp);
                if(rnorm<epsilon)//two in the same place
                {
                    tmstr[k].ptc.x+=0.0001f;tmstr[k].ptc.y+=0.0001f;tmstr[k].ptc.z+=0.0001f;
                    rp=setv(tmstr[k].ptc.x-tmstr[k].ptn.x,tmstr[k].ptc.y-tmstr[k].ptn.y,tmstr[k].ptc.z-tmstr[k].ptn.z);
                    rnorm=norm(rp);  
                }
                if(rnorm<lennc-delnc || rnorm>lennc+delnc)
                {
                    rp=scalx(rp,lennc/rnorm);
                    tmstr[k].ptc.x=tmstr[k].ptn.x+rp.x;
                    tmstr[k].ptc.y=tmstr[k].ptn.y+rp.y;
                    tmstr[k].ptc.z=tmstr[k].ptn.z+rp.z;
                }
            break;
        default : 
            break; 
    }
//    cout<<"rp: "<<rp.x<<" "<<rp.y<<" "<<rp.z<<endl;    
}

bool mcfragsweepLMP2(vector<point3f> &tmstr,int numseq,int poss,int pose)//[trandpos, pose)
{
    int trandpos=poss;
    int trand2=pose-poss;
//  point3f *tmstr=new point3f[numseq];
    int numiter;
    int k;
    point3d rp;
//  BasicFunc bf;
//  ParseSeq ps;
    double rnorm;
    bool flagdone;
    numiter=0;
        do{         
            for(k=trandpos;k<trandpos+trand2;k++)
            {   
                singlemoveLMPf(tmstr,k,int(randf0and1()*7.0));
                singlemoveLMPf(tmstr,k,int(randf0and1()*7.0));
                singlemoveLMPf(tmstr,k,int(randf0and1()*7.0));
                singlemoveLMPf(tmstr,k,int(randf0and1()*7.0));
                singlemoveLMPf(tmstr,k,int(randf0and1()*7.0));
                singlemoveLMPf(tmstr,k,int(randf0and1()*7.0));
                singlemoveLMPf(tmstr,k,int(randf0and1()*7.0));
            }
            ////////////////////////////////////////////////////////////////////////////////////////inverse     
            for(k=trandpos+trand2;k>trandpos;k--)
            {
                singlemoveLMPb(tmstr,k,int(randf0and1()*7.0));
                singlemoveLMPb(tmstr,k,int(randf0and1()*7.0));
                singlemoveLMPb(tmstr,k,int(randf0and1()*7.0));
                singlemoveLMPb(tmstr,k,int(randf0and1()*7.0));
                singlemoveLMPb(tmstr,k,int(randf0and1()*7.0));
                singlemoveLMPb(tmstr,k,int(randf0and1()*7.0));
                singlemoveLMPb(tmstr,k,int(randf0and1()*7.0));
                    
            }   

            ///////////////////////////////////////////////////////////////////////////////////////check
            flagdone=true;
            for(k=trandpos;k<=trandpos+trand2;k++)
            {
                //ca-1 ca
                rp=setv(tmstr[k].x-tmstr[k-1].x,tmstr[k].y-tmstr[k-1].y,tmstr[k].z-tmstr[k-1].z);
                rnorm=norm(rp);
                if(rnorm<lencaca-delcaca || rnorm>lencaca+delcaca)
                {
                    //printf("%d %d %d %d ca-1 ca %f lencaca %f\n",k,trandpos,trand2,numiter,rnorm,lencaca);
                    flagdone=false;
                    break;
                }
                //ca-1 n
                rp=setv(tmstr[k].ptn.x-tmstr[k-1].x,tmstr[k].ptn.y-tmstr[k-1].y,tmstr[k].ptn.z-tmstr[k-1].z);
                rnorm=norm(rp);
                if(rnorm<lencan1-delcan1 || rnorm>lencan1+delcan1)
                {
                    //printf("%d %d %d %d ca-1 n %f lencan1 %f\n",k,trandpos,trand2,numiter,rnorm,lencan1);
                    flagdone=false;
                    break;
                }
                //c-1 n
                rp=setv(tmstr[k].ptn.x-tmstr[k-1].ptc.x,tmstr[k].ptn.y-tmstr[k-1].ptc.y,tmstr[k].ptn.z-tmstr[k-1].ptc.z);
                rnorm=norm(rp);
                if(rnorm<lencn-delcn || rnorm>lencn+delcn)
                {
                    //printf("%d %d %d %d c-1 n %f lencn %f\n",k,trandpos,trand2,numiter,rnorm,lencn);
                    flagdone=false;
                    break;
                }
                //c-1 ca
                rp=setv(tmstr[k].x-tmstr[k-1].ptc.x,tmstr[k].y-tmstr[k-1].ptc.y,tmstr[k].z-tmstr[k-1].ptc.z);
                rnorm=norm(rp);
                if(rnorm<lencca1-delcca1 || rnorm>lencca1+delcca1)
                {
                    //printf("%d %d %d %d c-1 ca %f lencca1 %f\n",k,trandpos,trand2,numiter,rnorm,lencca1);
                    flagdone=false;
                    break;
                }
                if(k==trandpos+trand2)
                {
                    break;
                }
                //n ca
                rp=setv(tmstr[k].x-tmstr[k].ptn.x,tmstr[k].y-tmstr[k].ptn.y,tmstr[k].z-tmstr[k].ptn.z);
                rnorm=norm(rp);
                if(rnorm<lennca-delnca || rnorm>lennca+delnca)
                {
                    //printf("%d %d %d %d n ca %f lennca %f\n",k,trandpos,trand2,numiter,rnorm,lennca);
                    flagdone=false;
                    break;
                }
                //n c
                rp=setv(tmstr[k].ptc.x-tmstr[k].ptn.x,tmstr[k].ptc.y-tmstr[k].ptn.y,tmstr[k].ptc.z-tmstr[k].ptn.z);
                rnorm=norm(rp);
                if(rnorm<lennc-delnc || rnorm>lennc+delnc)
                {
                    //printf("%d %d %d %d n c %f lennc %f\n",k,trandpos,trand2,numiter,rnorm,lennc);
                    flagdone=false;
                    break;
                }
                //ca c
                rp=setv(tmstr[k].ptc.x-tmstr[k].x,tmstr[k].ptc.y-tmstr[k].y,tmstr[k].ptc.z-tmstr[k].z);
                rnorm=norm(rp);
                if(rnorm<lencac-delcac || rnorm>lencac+delcac)
                {
                    //printf("%d %d %d %d ca c %f lencac %f\n",k,trandpos,trand2,numiter,rnorm,lencac);
                    flagdone=false;
                    break;
                }
            }
            numiter++;
        //  printf("%d %d %d\n",trandpos,trand2,numiter);
        }while(!flagdone && numiter<200);
        if(!flagdone)
        {
//          printf("failure %d %d %d\n",trandpos,trand2,numiter);
//          if(trandpos+trand2+2>=numseq)
//              str2torp(tmstr,numseq,trandpos,numseq-1);
//          else
//              str2torp(tmstr,numseq,trandpos,trandpos+trand2+2);
//          delete[]tmstr;
            return false;
        }
        else
        {
//          printf("success %d %d %d\n",trandpos,trand2,numiter);
//          if(trandpos+trand2+2>=numseq)
//              str2torp(tmstr,numseq,trandpos,numseq-1);
//          else
//              str2torp(tmstr,numseq,trandpos,trandpos+trand2+2);
//          delete[]tmstr;
            return true;
        }
}

void singlemoveLMPb(vector<point3f> &tmstr,int k,int mtype)
{
    point3d rp;
    double rnorm;
    switch(mtype)
    {
        case 0:  
                //1 c-1 n
                rp=setv(tmstr[k-1].ptc.x-tmstr[k].ptn.x,tmstr[k-1].ptc.y-tmstr[k].ptn.y,tmstr[k-1].ptc.z-tmstr[k].ptn.z);
                rnorm=norm(rp);
                if(rnorm<epsilon)//two in the same place
                {
                    tmstr[k-1].ptc.x-=0.0001f;tmstr[k-1].ptc.y-=0.0001f;tmstr[k-1].ptc.z-=0.0001f;
                    rp=setv(tmstr[k-1].ptc.x-tmstr[k].ptn.x,tmstr[k-1].ptc.y-tmstr[k].ptn.y,tmstr[k-1].ptc.z-tmstr[k].ptn.z);
                    rnorm=norm(rp);  
                }
                if(rnorm<lencn-delcn || rnorm>lencn+delcn)
                {
                    rp=scalx(rp,lencn/rnorm);
                    tmstr[k-1].ptc.x=tmstr[k].ptn.x+rp.x;
                    tmstr[k-1].ptc.y=tmstr[k].ptn.y+rp.y;
                    tmstr[k-1].ptc.z=tmstr[k].ptn.z+rp.z;
                }
            break;
        case 1:  
                //4 c-1 ca
                rp=setv(tmstr[k-1].ptc.x-tmstr[k].x,tmstr[k-1].ptc.y-tmstr[k].y,tmstr[k-1].ptc.z-tmstr[k].z);
                rnorm=norm(rp);
                if(rnorm<epsilon)//two in the same place
                {
                    tmstr[k-1].ptc.x-=0.0001f;tmstr[k-1].ptc.y-=0.0001f;tmstr[k-1].ptc.z-=0.0001f;
                    rp=setv(tmstr[k-1].ptc.x-tmstr[k].x,tmstr[k-1].ptc.y-tmstr[k].y,tmstr[k-1].ptc.z-tmstr[k].z);
                    rnorm=norm(rp);  
                }
                if(rnorm<lencca1-delcca1 || rnorm>lencca1+delcca1)
                {
                    rp=scalx(rp,lencca1/rnorm);
                    tmstr[k-1].ptc.x=tmstr[k].x+rp.x;
                    tmstr[k-1].ptc.y=tmstr[k].y+rp.y;
                    tmstr[k-1].ptc.z=tmstr[k].z+rp.z;
                }
            break;
        case 2:  
                //2 ca-1 n
                rp=setv(tmstr[k-1].x-tmstr[k].ptn.x,tmstr[k-1].y-tmstr[k].ptn.y,tmstr[k-1].z-tmstr[k].ptn.z);
                rnorm=norm(rp);
                if(rnorm<epsilon)//two in the same place
                {
                    tmstr[k-1].x-=0.0001f;tmstr[k-1].y-=0.0001f;tmstr[k-1].z-=0.0001f;
                    rp=setv(tmstr[k-1].x-tmstr[k].ptn.x,tmstr[k-1].y-tmstr[k].ptn.y,tmstr[k-1].z-tmstr[k].ptn.z);
                    rnorm=norm(rp);  
                }
                if(rnorm<lencan1-delcan1 || rnorm>lencan1+delcan1)
                {
                    rp=scalx(rp,lencan1/rnorm);
                    tmstr[k-1].x=tmstr[k].ptn.x+rp.x;
                    tmstr[k-1].y=tmstr[k].ptn.y+rp.y;
                    tmstr[k-1].z=tmstr[k].ptn.z+rp.z;
                }
            break;
        case 3:  
                //5 ca-1 ca
                rp=setv(tmstr[k-1].x-tmstr[k].x,tmstr[k-1].y-tmstr[k].y,tmstr[k-1].z-tmstr[k].z);
                rnorm=norm(rp);
                if(rnorm<epsilon)//two in the same place
                {   
                    tmstr[k-1].x-=0.0001f;tmstr[k-1].y-=0.0001f;tmstr[k-1].z-=0.0001f;
                    rp=setv(tmstr[k-1].x-tmstr[k].x,tmstr[k-1].y-tmstr[k].y,tmstr[k-1].z-tmstr[k].z);
                    rnorm=norm(rp);  
                }
                if(rnorm<lencaca-delcaca || rnorm>lencaca+delcaca)
                {
                    rp=scalx(rp,lencaca/rnorm);
                    tmstr[k-1].x=tmstr[k].x+rp.x;
                    tmstr[k-1].y=tmstr[k].y+rp.y;
                    tmstr[k-1].z=tmstr[k].z+rp.z;
                }
            break;
        case 4:  
                //6 ca c
                rp=setv(tmstr[k-1].x-tmstr[k-1].ptc.x,tmstr[k-1].y-tmstr[k-1].ptc.y,tmstr[k-1].z-tmstr[k-1].ptc.z);
                rnorm=norm(rp);
                if(rnorm<epsilon)//two in the same place
                {
                    tmstr[k-1].x-=0.0001f;tmstr[k-1].y-=0.0001f;tmstr[k-1].z-=0.0001f;
                    rp=setv(tmstr[k-1].x-tmstr[k-1].ptc.x,tmstr[k-1].y-tmstr[k-1].ptc.y,tmstr[k-1].z-tmstr[k-1].ptc.z);
                    rnorm=norm(rp);  
                }
                if(rnorm<lencac-delcac || rnorm>lencac+delcac)
                {
                    rp=scalx(rp,lencac/rnorm);
                    tmstr[k-1].x=tmstr[k-1].ptc.x+rp.x;
                    tmstr[k-1].y=tmstr[k-1].ptc.y+rp.y;
                    tmstr[k-1].z=tmstr[k-1].ptc.z+rp.z;
                }
            break;
        case 5:  
                //3 n ca
                rp=setv(tmstr[k-1].ptn.x-tmstr[k-1].x,tmstr[k-1].ptn.y-tmstr[k-1].y,tmstr[k-1].ptn.z-tmstr[k-1].z);
                rnorm=norm(rp);
                if(rnorm<epsilon)//two in the same place
                {
                    tmstr[k-1].ptn.x-=0.0001f;tmstr[k-1].ptn.y-=0.0001f;tmstr[k-1].ptn.z-=0.0001f;
                    rp=setv(tmstr[k-1].ptn.x-tmstr[k-1].x,tmstr[k-1].ptn.y-tmstr[k-1].y,tmstr[k-1].ptn.z-tmstr[k-1].z);
                    rnorm=norm(rp);  
                }
                if(rnorm<lennca-delnca || rnorm>lennca+delnca)
                {
                    rp=scalx(rp,lennca/rnorm);
                    tmstr[k-1].ptn.x=tmstr[k-1].x+rp.x;
                    tmstr[k-1].ptn.y=tmstr[k-1].y+rp.y;
                    tmstr[k-1].ptn.z=tmstr[k-1].z+rp.z;
                }
            break;
        case 6:  
                //7 n c
                rp=setv(tmstr[k-1].ptn.x-tmstr[k-1].ptc.x,tmstr[k-1].ptn.y-tmstr[k-1].ptc.y,tmstr[k-1].ptn.z-tmstr[k-1].ptc.z);
                rnorm=norm(rp);
                if(rnorm<epsilon)//two in the same place
                {
                    tmstr[k-1].ptn.x-=0.0001f;tmstr[k-1].ptn.y-=0.0001f;tmstr[k-1].ptn.z-=0.0001f;
                    rp=setv(tmstr[k-1].ptn.x-tmstr[k-1].ptc.x,tmstr[k-1].ptn.y-tmstr[k-1].ptc.y,tmstr[k-1].ptn.z-tmstr[k-1].ptc.z);
                    rnorm=norm(rp);  
                }
                if(rnorm<lennc-delnc || rnorm>lennc+delnc)
                {
                    rp=scalx(rp,lennc/rnorm);
                    tmstr[k-1].ptn.x=tmstr[k-1].ptc.x+rp.x;
                    tmstr[k-1].ptn.y=tmstr[k-1].ptc.y+rp.y;
                    tmstr[k-1].ptn.z=tmstr[k-1].ptc.z+rp.z;
                }
            break;
        default : 
            break; 
    }
}

point3d scalx(point3d p1,double f)
{
    point3d temp;
    temp.x=f*p1.x;
    temp.y=f*p1.y;
    temp.z=f*p1.z;
    return temp;
}


int Eenergybondangle(vector<point3f> &decstr,int numseq,float &fene)
{
    int i;
    double totene=0;
    lableproblematic lp1,lp2;
//  double meanstda[][4]={
//  {111.8,2.5*4,111.8-2.5*4,111.8+2.5*4},//ncac
//  {116.9,1.5*4,116.9-1.5*4,116.9+1.5*4},//cacn
//  {122.6,5.0*4,122.6-5.0*4,122.6+5.0*4},//cnca
//  };
    int istart;
    lp1.nn[0]=0;
    lp1.indn[0] = new int [numseq];
    lp1.indn[20] = new int [numseq];    
    memset(lp1.indn[0],0,numseq*sizeof(int));
    memset(lp1.indn[20],0,numseq*sizeof(int));
//    cout<<"rrr"<<endl;
    for(i=1;i<numseq;i++)
    {
        if(decstr[i].iaa==5) istart=3;
        else if(decstr[i].iaa==12) istart=6;
        else istart=0;
        if(decstr[i].ang[0]<stdangle[istart+1][2]) 
        {
            lp1.indn[0][lp1.nn[0]]=i;
            lp1.nn[0]++;
            lp1.indn[20][i]=1;
            totene+=stdangle[istart+1][2]-decstr[i].ang[0];
        }
        else if(decstr[i].ang[0]>stdangle[istart+1][3])
        {
            lp1.indn[0][lp1.nn[0]]=i;
            lp1.nn[0]++;
            lp1.indn[20][i]=1;
            totene+=decstr[i].ang[0]-stdangle[istart+1][3];
        }
    }
//    cout<<"eee"<<endl;
    lp1.nn[1]=0;
    lp1.indn[1] = new int [numseq];
    lp1.indn[21] = new int [numseq];    
    memset(lp1.indn[1],0,numseq*sizeof(int));
    memset(lp1.indn[21],0,numseq*sizeof(int));
    for(i=1;i<numseq;i++)
    {
        if(decstr[i].iaa==5) istart=3;
        else if(decstr[i].iaa==12) istart=6;
        else istart=0;
        if(decstr[i].ang[1]<stdangle[istart+2][2])
        {
            lp1.indn[1][lp1.nn[1]]=i;
            lp1.nn[1]++;
            lp1.indn[21][i]=1;
            totene+=stdangle[istart+2][2]-decstr[i].ang[1];
        }
        else if(decstr[i].ang[1]>stdangle[istart+2][3])
        {
            lp1.indn[1][lp1.nn[1]]=i;
            lp1.nn[1]++;
            lp1.indn[21][i]=1;
            totene+=decstr[i].ang[1]-stdangle[istart+2][3];
        }
    }
//    cout<<"www"<<endl;
    lp1.nn[2]=0;
    lp1.indn[2] = new int [numseq];
    lp1.indn[22] = new int [numseq];    
    memset(lp1.indn[2],0,numseq*sizeof(int));
    memset(lp1.indn[22],0,numseq*sizeof(int));
    for(i=0;i<numseq;i++)
    {
        if(decstr[i].iaa==5) istart=3;
        else if(decstr[i].iaa==12) istart=6;
        else istart=0;
        if(decstr[i].ang[2]<stdangle[istart+0][2])
        {
            lp1.indn[2][lp1.nn[2]]=i;
            lp1.nn[2]++;
            lp1.indn[22][i]=1;
            totene+=stdangle[istart+0][2]-decstr[i].ang[2];
        }
        else if(decstr[i].ang[2]>stdangle[istart+0][3])
        {
            lp1.indn[2][lp1.nn[2]]=i;
            lp1.nn[2]++;
            lp1.indn[22][i]=1;
            totene+=decstr[i].ang[2]-stdangle[istart+0][3];
        }
    }
    int totocn=0;
    double angleocn;
    point3d tp1,tp2;
    for(i=1;i<numseq;i++)
    {
        if(decstr[i].iaa==5) istart=13;
        else if(decstr[i].iaa==12) istart=14;
        else istart=12;
        tp1=setv(decstr[i-1].pto.x-decstr[i-1].ptc.x,decstr[i-1].pto.y-decstr[i-1].ptc.y,decstr[i-1].pto.z-decstr[i-1].ptc.z);
        tp2=setv(decstr[i].ptn.x-decstr[i-1].ptc.x,decstr[i].ptn.y-decstr[i-1].ptc.y,decstr[i].ptn.z-decstr[i-1].ptc.z);
        angleocn=angv(tp1,tp2)*degrad;
        if(angleocn<stdangle[istart+0][2])
        {
            totocn++;
            totene+=stdangle[istart+0][2]-angleocn;
        }
        else if(angleocn>stdangle[istart+0][3])
        {
            totocn++;
            totene+=angleocn-stdangle[istart+0][3];
        }
    }
    fene=totene;
    delete [] lp1.indn[20];
    delete [] lp1.indn[21];
    delete [] lp1.indn[22];
    delete [] lp1.indn[0];
    delete [] lp1.indn[1];
    delete [] lp1.indn[2];
    return lp1.nn[0]+lp1.nn[1]+lp1.nn[2]+totocn;
}

int energybondlength(vector<point3f> &decstr,int numseq,float &fene)
{
    int i;
    double totene=0;
//  double meanstd[][4]={
//  {1.466,0.015*4,1.466-0.015*4,1.466+0.015*4},//nca
//  {1.525,0.021*4,1.525-0.021*4,1.525+0.021*4},//cac
//  {1.341,0.016*4,1.341-0.016*4,1.341+0.016*4},//cn
//  {3.813,0.080*4,3.813-0.080*4,3.813+0.080*4},//caca
//  };
    lableproblematic lp1,lp2;
    int istart;
    lp1.nn[3]=0;
    lp1.indn[23] = new int [numseq];
    lp1.indn[3] = new int [numseq];
    memset(lp1.indn[3],0,numseq*sizeof(int));
    memset(lp1.indn[23],0,numseq*sizeof(int));
    for(i=1;i<numseq;i++)
    {
        if(decstr[i].iaa==5) istart=4;
        else if(decstr[i].iaa==12) istart=7;
        else istart=0;
        if(decstr[i].len[0]<stdlength[istart+2][2])
        {
    //      printf("%3d cn [%.4f %.4f] %.4f\n",i,meanstd[2][2],meanstd[2][3],decstr[i].len[0]);
            lp1.indn[3][lp1.nn[3]]=i;
            lp1.nn[3]++;
            lp1.indn[23][i]=1;
            totene+=stdlength[istart+2][2]-decstr[i].len[0];
        }
        else if(decstr[i].len[0]>stdlength[istart+2][3])
        {
            lp1.indn[3][lp1.nn[3]]=i;
            lp1.nn[3]++;
            lp1.indn[23][i]=1;
            totene+=decstr[i].len[0]-stdlength[istart+2][3];
        }
    }
    lp1.nn[4]=0;
    lp1.indn[24] = new int [numseq];
    lp1.indn[4] = new int [numseq];
    memset(lp1.indn[24],0,numseq*sizeof(int));
    memset(lp1.indn[4],0,numseq*sizeof(int));
//    cout<<"seq: "<<numseq<<endl;
    for(i=0;i<numseq;i++)
    {
        if(decstr[i].iaa==5) istart=4;
        else if(decstr[i].iaa==12) istart=7;
        else istart=0;
        if(decstr[i].len[1]<stdlength[istart+0][2])
        {
    //        cout<<"ll1"<<endl;
    //      printf("%3d nca [%.4f %.4f] %.4f\n",i,meanstd[0][2],meanstd[0][3],decstr[i].len[1]);
            lp1.indn[4][lp1.nn[4]]=i;
            lp1.nn[4]++;
            lp1.indn[24][i]=1;
            totene+=stdlength[istart+0][2]-decstr[i].len[1];
        }
        else if(decstr[i].len[1]>stdlength[istart+0][3])
        {
        //    cout<<"ll2: "<<lp1.nn[4]<<endl;
            lp1.indn[4][lp1.nn[4]]=i;
            lp1.nn[4]++;
            lp1.indn[24][i]=1;
            totene+=decstr[i].len[1]-stdlength[istart+0][3];
        }
    }
    lp1.nn[5]=0;
    lp1.indn[25] = new int [numseq];
    lp1.indn[5] = new int [numseq];
    memset(lp1.indn[25],0,numseq*sizeof(int));
    memset(lp1.indn[5],0,numseq*sizeof(int));
    for(i=0;i<numseq;i++)
    {
        if(decstr[i].iaa==5) istart=4;
        else if(decstr[i].iaa==12) istart=7;
        else istart=0;
        if(decstr[i].len[2]<stdlength[istart+1][2])
        {
    //      printf("%3d cac [%.4f %.4f] %.4f\n",i,meanstd[1][2],meanstd[1][3],decstr[i].len[2]);
            lp1.indn[5][lp1.nn[5]]=i;
            lp1.nn[5]++;
            lp1.indn[25][i]=1;
            totene+=stdlength[istart+1][2]-decstr[i].len[2];
        }
        else if(decstr[i].len[2]>stdlength[istart+1][3])
        {
            lp1.indn[5][lp1.nn[5]]=i;
            lp1.nn[5]++;
            lp1.indn[25][i]=1;
            totene+=decstr[i].len[2]-stdlength[istart+1][3];
        }
    }
    lp1.nn[6]=0;
    lp1.indn[26] = new int [numseq];
    lp1.indn[6] = new int [numseq];
    memset(lp1.indn[26],0,numseq*sizeof(int));
    memset(lp1.indn[6],0,numseq*sizeof(int));
    double tdist;
    point3d tp;
    for(i=1;i<numseq;i++)
    {
        tp.x=decstr[i].x-decstr[i-1].x;
        tp.y=decstr[i].y-decstr[i-1].y;
        tp.z=decstr[i].z-decstr[i-1].z;
        tdist=norm(tp);
        if(tdist<stdlength[3][2])
        {
    //      printf("%3d ca [%.4f %.4f] %.4f\n",i,meanstd[3][2],meanstd[3][3],tdist);
            lp1.indn[6][lp1.nn[6]]=i;
            lp1.nn[6]++;
            lp1.indn[26][i]=1;
            totene+=stdlength[3][2]-tdist;
        }
        else if(tdist>stdlength[3][3])
        {
            lp1.indn[6][lp1.nn[6]]=i;
            lp1.nn[6]++;
            lp1.indn[26][i]=1;
            totene+=tdist-stdlength[3][3];
        }
    }
    fene=totene;

    delete [] lp1.indn[23];
    delete [] lp1.indn[24];
    delete [] lp1.indn[25];
    delete [] lp1.indn[26];
    delete [] lp1.indn[3];
    delete [] lp1.indn[4];
    delete [] lp1.indn[5];
    delete [] lp1.indn[6];
    return lp1.nn[3]+lp1.nn[4]+lp1.nn[5]+lp1.nn[6];
}

int aminoid(vector<char> &aminoname)
{
    int i;
    //empty
    if(aminoname[0]==' ' && aminoname[1]==' ' && aminoname[2]==' ')
    {
    //  printf("empty amino name\n");
        return 23;
    }
    // one in 26
    for(i=0;i<26;i++)
    {
        if(aminoname[0]==aad3[i][0] && aminoname[1]==aad3[i][1] && aminoname[2]==aad3[i][2])
        {
            if(i==24)
            {
                i=23;
            }
            return i;
        }
//      else if(aminoname[0]==' ' && aminoname[1]==' ' && aminoname[2]==aad1[i])
//      {
//          return i;
//      }
    }
    for(i=0;i<8;i++)
    {
        if(aminoname[0]==dnarnares[i][0] && aminoname[1]==dnarnares[i][1] && aminoname[2]==dnarnares[i][2])
        {
            return -i-1;
        }
    }
    // unknown
    //printf("unknown amino name %c%c%c\n",aminoname[0],aminoname[1],aminoname[2]);
    return 23;
}

int aminoid(char aminoname)
{
    int i;
    //empty
    if(aminoname==' ')
    {
    //  printf("empty amino name\n");
        return 23;
    }
    // one in 26
    for(i=0;i<26;i++)
    {
        if(aminoname==aad1[i])
        {
            return i;
        }
    }
    // unknown
    //printf("unknown amino name %c%c%c\n",aminoname[0],aminoname[1],aminoname[2]);
    return 23;
}
int fenergybondangle(vector<residueatoms> &decstr,int numseq,double fene,bool flagacc)
{
    int i;
    double totene=0;
    lableproblematic lp1,lp2;
//  double meanstda[][4]={
//  {111.8,2.5*4,111.8-2.5*4,111.8+2.5*4},//ncac
//  {116.9,1.5*4,116.9-1.5*4,116.9+1.5*4},//cacn
//  {122.6,5.0*4,122.6-5.0*4,122.6+5.0*4},//cnca
//  };
    int istart;
    point3d tp1,tp2;
    double angleocn;
if(!flagacc)
{
    lp1.nn[0]=0;
    lp1.fn[0]=0;
    memset(lp1.indn[0],0,numseq*sizeof(int));
    memset(lp1.indf[0],0,numseq*sizeof(double));
    memset(lp1.flap[0],0,numseq*sizeof(bool));
    for(i=1;i<numseq;i++)
    {
        if(decstr[i].resind==5) istart=3;
        else if(decstr[i].resind==12) istart=6;
        else istart=0;
        if(decstr[i].ang[0]<stdangle[istart+1][2]) 
        {
            lp1.nn[0]++;
            lp1.indn[0][i]=1;
            lp1.indf[0][i]=stdangle[istart+1][2]-decstr[i].ang[0];
            lp1.fn[0]+=lp1.indf[0][i];
        }
        else if(decstr[i].ang[0]>stdangle[istart+1][3])
        {
            lp1.nn[0]++;
            lp1.indn[0][i]=1;
            lp1.indf[0][i]=decstr[i].ang[0]-stdangle[istart+1][3];
            lp1.fn[0]+=lp1.indf[0][i];
        }
    }

    lp1.nn[1]=0;
    lp1.fn[1]=0;
    memset(lp1.indn[1],0,numseq*sizeof(int));
    memset(lp1.indf[1],0,numseq*sizeof(double));
    memset(lp1.flap[1],0,numseq*sizeof(bool));
    for(i=1;i<numseq;i++)
    {
        if(decstr[i].resind==5) istart=3;
        else if(decstr[i].resind==12) istart=6;
        else istart=0;
        if(decstr[i].ang[1]<stdangle[istart+2][2])
        {
            lp1.nn[1]++;
            lp1.indn[1][i]=1;
            lp1.indf[1][i]=stdangle[istart+2][2]-decstr[i].ang[1];
            lp1.fn[1]+=lp1.indf[1][i];
        }
        else if(decstr[i].ang[1]>stdangle[istart+2][3])
        {
            lp1.nn[1]++;
            lp1.indn[1][i]=1;
            lp1.indf[1][i]=decstr[i].ang[1]-stdangle[istart+2][3];
            lp1.fn[1]+=lp1.indf[1][i];
        }
    }
    
    lp1.nn[2]=0;
    lp1.fn[2]=0;
    memset(lp1.indn[2],0,numseq*sizeof(int));
    memset(lp1.indf[2],0,numseq*sizeof(double));
    memset(lp1.flap[2],0,numseq*sizeof(bool));
    for(i=0;i<numseq;i++)
    {
        if(decstr[i].resind==5) istart=3;
        else if(decstr[i].resind==12) istart=6;
        else istart=0;
        if(decstr[i].ang[2]<stdangle[istart+0][2])
        {
            lp1.nn[2]++;
            lp1.indn[2][i]=1;
            lp1.indf[2][i]=stdangle[istart+0][2]-decstr[i].ang[2];
            lp1.fn[2]+=lp1.indf[2][i];
        }
        else if(decstr[i].ang[2]>stdangle[istart+0][3])
        {
            lp1.nn[2]++;
            lp1.indn[2][i]=1;
            lp1.indf[2][i]=decstr[i].ang[2]-stdangle[istart+0][3];
            lp1.fn[2]+=lp1.indf[2][i];
        }
    }
    ////////////////////////////////
    lp1.nn[14]=0;
    lp1.fn[14]=0;
    memset(lp1.indn[14],0,numseq*sizeof(int));
    memset(lp1.indf[14],0,numseq*sizeof(double));
    for(i=1;i<numseq;i++)
    {
        if(decstr[i].resind==5) istart=13;
        else if(decstr[i].resind==12) istart=14;
        else istart=12;
        tp1=setv(decstr[i-1].ptv[3].x-decstr[i-1].ptv[2].x,decstr[i-1].ptv[3].y-decstr[i-1].ptv[2].y,decstr[i-1].ptv[3].z-decstr[i-1].ptv[2].z);
        tp2=setv(decstr[i].ptv[0].x-decstr[i-1].ptv[2].x,decstr[i].ptv[0].y-decstr[i-1].ptv[2].y,decstr[i].ptv[0].z-decstr[i-1].ptv[2].z);
        angleocn=angv(tp1,tp2)*degrad;
        if(angleocn<stdangle[istart][2])
        {
            lp1.nn[14]++;
            lp1.indn[14][i]=1;
            lp1.indf[14][i]=stdangle[istart][2]-angleocn;
            lp1.fn[14]+=lp1.indf[14][i];
        }
        else if(angleocn>stdangle[istart][3])
        {
            lp1.nn[14]++;
            lp1.indn[14][i]=1;
            lp1.indf[14][i]=angleocn-stdangle[istart][3];
            lp1.fn[14]+=lp1.indf[14][i];
        }
    }
}
else
{
    int iold[4];
    int inew[4];
    double fold[4];
    double fnew[4];
    for(i=0;i<4;i++)
    {
        iold[i]=0;
        inew[i]=0;
        fold[i]=0;
        fnew[i]=0;
    }
    ////////////////ocn
    for(i=1;i<numseq;i++) if(lp1.flap[7][i-1] || lp1.flap[0][i])
    {
        if(decstr[i].resind==5) istart=13;
        else if(decstr[i].resind==12) istart=14;
        else istart=12;
        iold[3]+=lp1.indn[14][i];
        fold[3]+=lp1.indf[14][i];
        tp1=setv(decstr[i-1].ptv[3].x-decstr[i-1].ptv[2].x,decstr[i-1].ptv[3].y-decstr[i-1].ptv[2].y,decstr[i-1].ptv[3].z-decstr[i-1].ptv[2].z);
        tp2=setv(decstr[i].ptv[0].x-decstr[i-1].ptv[2].x,decstr[i].ptv[0].y-decstr[i-1].ptv[2].y,decstr[i].ptv[0].z-decstr[i-1].ptv[2].z);
        angleocn=angv(tp1,tp2)*degrad;
        if(angleocn<stdangle[istart][2])
        {
            lp1.indn[14][i]=1;
            lp1.indf[14][i]=stdangle[istart][2]-angleocn;
        }
        else if(angleocn>stdangle[istart][3])
        {
            lp1.indn[14][i]=1;
            lp1.indf[14][i]=angleocn-stdangle[istart][3];
        }
        else
        {
            lp1.indn[14][i]=0;
            lp1.indf[14][i]=0;
        }
        inew[3]+=lp1.indn[14][i];
        fnew[3]+=lp1.indf[14][i];
    }
    lp1.nn[14]+=inew[3]-iold[3];
    lp1.fn[14]+=fnew[3]-fold[3];
    ///////////////////////////////////////////
    for(i=1;i<numseq;i++) if(lp1.flap[0][i])
    {
        if(decstr[i].resind==5) istart=3;
        else if(decstr[i].resind==12) istart=6;
        else istart=0;
        iold[0]+=lp1.indn[0][i];
        fold[0]+=lp1.indf[0][i];
        if(decstr[i].ang[0]<stdangle[istart+1][2]) 
        {
            lp1.indn[0][i]=1;
            lp1.indf[0][i]=stdangle[istart+1][2]-decstr[i].ang[0];
        }
        else if(decstr[i].ang[0]>stdangle[istart+1][3])
        {
            lp1.indn[0][i]=1;
            lp1.indf[0][i]=decstr[i].ang[0]-stdangle[istart+1][3];
        }
        else
        {
            lp1.indn[0][i]=0;
            lp1.indf[0][i]=0;
        }
        inew[0]+=lp1.indn[0][i];
        fnew[0]+=lp1.indf[0][i];
        lp1.flap[0][i]=false;
    }
    lp1.nn[0]+=inew[0]-iold[0];
    lp1.fn[0]+=fnew[0]-fold[0];
    for(i=1;i<numseq;i++) if(lp1.flap[1][i])
    {
        if(decstr[i].resind==5) istart=3;
        else if(decstr[i].resind==12) istart=6;
        else istart=0;
        iold[1]+=lp1.indn[1][i];
        fold[1]+=lp1.indf[1][i];
        if(decstr[i].ang[1]<stdangle[istart+2][2])
        {
            lp1.indn[1][i]=1;
            lp1.indf[1][i]=stdangle[istart+2][2]-decstr[i].ang[1];
        }
        else if(decstr[i].ang[1]>stdangle[istart+2][3])
        {
            lp1.indn[1][i]=1;
            lp1.indf[1][i]=decstr[i].ang[1]-stdangle[istart+2][3];
        }
        else
        {
            lp1.indn[1][i]=0;
            lp1.indf[1][i]=0;
        }
        inew[1]+=lp1.indn[1][i];
        fnew[1]+=lp1.indf[1][i];
        lp1.flap[1][i]=false;
    }
    lp1.nn[1]+=inew[1]-iold[1];
    lp1.fn[1]+=fnew[1]-fold[1];
    for(i=0;i<numseq;i++) if(lp1.flap[2][i])
    {
        if(decstr[i].resind==5) istart=3;
        else if(decstr[i].resind==12) istart=6;
        else istart=0;
        iold[2]+=lp1.indn[2][i];
        fold[2]+=lp1.indf[2][i];
        if(decstr[i].ang[2]<stdangle[istart+0][2])
        {
            lp1.indn[2][i]=1;
            lp1.indf[2][i]=stdangle[istart+0][2]-decstr[i].ang[2];
        }
        else if(decstr[i].ang[2]>stdangle[istart+0][3])
        {
            lp1.indn[2][i]=1;
            lp1.indf[2][i]=decstr[i].ang[2]-stdangle[istart+0][3];
        }
        else
        {
            lp1.indn[2][i]=0;
            lp1.indf[2][i]=0;
        }
        inew[2]+=lp1.indn[2][i];
        fnew[2]+=lp1.indf[2][i];
        lp1.flap[2][i]=false;
    }
    lp1.nn[2]+=inew[2]-iold[2];
    lp1.fn[2]+=fnew[2]-fold[2];
    
}
    fene=lp1.fn[0]+lp1.fn[1]+lp1.fn[2]+lp1.fn[14];
    return lp1.nn[0]+lp1.nn[1]+lp1.nn[2]+lp1.nn[14];
}

double energyhbondnhoc3(vector<point3f> &decstr,int numseq)
{
    int i,j,k,l,m;
    point3d tp[20],ap[20],kp[20];
    point3d pd,pd2,pd3;
    double lamda;
    double diss[20];
    double totenergy=0;
    double totenergy2=0;

    static  double nhochhmean[22]={
        5.2172,8.7074,230.1397,4.2198,6.5391
        };
    static  double nhochhsigma[22]={
        0.3676,0.4375,10.2234,0.4080,0.3763
        };
    static  double nhocmeanval[][4]={//h-o h-o=c n-h-o n-h-o=c  
        2.85,89.0,110.5,199.5,//i+3
        2.00,147.0,159.0,160.0, //i+4  oldoh 
    //  2.83,89.0,110.0,201.5,//i+3 newoh
    //  2.00,148.0,159.0,155.0, //i+4

        2.00,155.0,164.0,180.0,//0
        2.00,155.0,164.0,180.0,//1
        2.00,151.0,163.0,192.0,//2
        2.00,151.0,163.0,192.0,//3

        };
    static  double nhocstdval[][4]={//h-o h-o=c n-h-o n-h-o=c 
        0.315504,7.697185,8.980366,7.932107,
        0.530859,10.582243,11.249764,25.360054,

        0.299730,11.770196,11.292558,68.955920,
        0.299730,11.770196,11.292558,68.955920,
        0.255088,12.376087,11.020081,69.165282,
        0.255088,12.376087,11.020081,69.165282,
        };

    double threshval[5]={19.113828,19.113828,20.723266,20.723266,16.118096};
    for(i=0;i<numseq;i++)
    {
        decstr[i].vpos=0;
        decstr[i].tpos=0;
        decstr[i].ssm='C';
        decstr[i].indl=-1;
        decstr[i].indr=-1;
        decstr[i].tpr=-1;
        decstr[i].tpl=-1;
    }

/*    for(int i=1;i<numseq;i++)
    {
        float pox = 0.0;
        float poy = 0.0;
        float poz = 0.0;
        tor2pos22(decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
                decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
                179.6715f*raddeg,1.229f,2.0961f,&pox,&poy,&poz);
        decstr[i].pto.x = pox ;
        decstr[i].pto.y = poy ;
        decstr[i].pto.z = poz ;
        cout<<"pto: "<<pox<<" "<<poy<<" "<<poz<<endl;
        float phx = 0.0;
        float phy = 0.0;
        float phz = 0.0;
        tor2pos22(decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,decstr[i].x,decstr[i].y,decstr[i].z,
            decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,179.8174f*raddeg,0.987f,2.0814f,&phx,&phy,&phz);//0.9919f,2.0574f
        decstr[i].pth.x = phx ;
        decstr[i].pth.y = phy ;
        decstr[i].pth.z = phz ; 
        cout<<"pth: "<<phx<<" "<<phy<<" "<<phz<<endl;           
    } */
    str2tor(decstr,numseq,3);
    tor2stroh(decstr,numseq);    

    for(i=1;i<numseq-4;i++)
    {
        j=i+3;  
        if((decstr[i].ss2=='E') || (decstr[j].ss2=='E' ))
            continue;
        l=0;
        for(k=0;k<4;k++)
        {
            if(decstr[i+k].ss2=='H') 
                l++;
        } 
    //  m=0;
    //    for(k=0;k<4;k++)
    //    {
    //        if(decstr[i+k].stype=='H') 
    //            m++;
    //    } 
    //    if(l==3 && (decstr[i+1].ss2!='H' || decstr[i+2].ss2!='H'))
    //        continue;
    //    else if(l==2  && 
    //        ((decstr[i+1].ss2=='H' && decstr[i+2].ss2=='H') || (decstr[i+1].ss2!='H' && decstr[i+2].ss2!='H') ))
    //        continue;  
        if(l==0 ) continue;
        else if(l==1 ) continue;
        else if(l==2 ) lamda=0.1;
    //    else if(l==2 ) lamda=0.2;
        else if(l==3) lamda=0.4;
    //    else if(l==3) lamda=0.5;
        else if(l==4) lamda=1.5;
    //  else if(m==3) lamda=0.6;
    //  else if(m==4) lamda=1.0;
        else lamda=0.4;  
        
        tp[1]=setv(decstr[i].x,decstr[i].y,decstr[i].z);
        ap[1]=setv(decstr[j].x,decstr[j].y,decstr[j].z);
        tp[7]=minu(tp[1],ap[1]);
        diss[1]=norm(tp[7]); 
        if(diss[1]>6.6) continue;

        tp[0]=setv(decstr[i-1].x,decstr[i-1].y,decstr[i-1].z);
        tp[2]=setv(decstr[i+1].x,decstr[i+1].y,decstr[i+1].z);
        ap[0]=setv(decstr[j-1].x,decstr[j-1].y,decstr[j-1].z);
        ap[2]=setv(decstr[j+1].x,decstr[j+1].y,decstr[j+1].z);
        kp[0]=minu(tp[2],ap[0]);
        diss[4]=norm(kp[0]);
        if(diss[4]<3.4 || diss[4]>4.2) continue;
//        diss[0]=norm(minu(tp[0],ap[0]));
//      diss[1]=bf.norm(bf.minu(tp[0],ap[2]));          
//      diss[2]=bf.phi(tp[1].x,tp[1].y,tp[1].z,tp[2].x,tp[2].y,tp[2].z,ap[0].x,ap[0].y,ap[0].z,ap[1].x,ap[1].y,ap[1].z);
//      if(diss[2]>180) diss[2]-=180;
//      else diss[2]+=180;


        tp[3]=setv(decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z);
        tp[5]=setv(decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z);
        ap[4]=setv(decstr[j+1].ptn.x,decstr[j+1].ptn.y,decstr[j+1].ptn.z);
        ap[5]=setv(decstr[j].ptc.x,decstr[j].ptc.y,decstr[j].ptc.z); 
        diss[0]=norm(minu(tp[0],ap[0]));
        diss[1]=norm(minu(tp[0],ap[2]));          
        diss[2]=phi(tp[1].x,tp[1].y,tp[1].z,tp[2].x,tp[2].y,tp[2].z,ap[0].x,ap[0].y,ap[0].z,ap[1].x,ap[1].y,ap[1].z); 
        if(diss[2]>180) diss[2]-=180;
        else diss[2]+=180;
        diss[3]=norm(minu(tp[3],ap[4]));
        diss[4]=norm(minu(tp[5],ap[5]));           
//      diss[3]=bf.norm(bf.minu(tp[3],ap[4]));
//      diss[4]=bf.norm(bf.minu(tp[5],ap[5]));

        //nhoc
        pd.x=decstr[j].pth.x-decstr[i].pto.x;
        pd.y=decstr[j].pth.y-decstr[i].pto.y;
        pd.z=decstr[j].pth.z-decstr[i].pto.z;
        diss[5]=norm(pd);
        if(diss[5]>=5.0 || diss[5]<1.6) continue;
 

        pd2.x=decstr[i].ptc.x-decstr[i].pto.x;
        pd2.y=decstr[i].ptc.y-decstr[i].pto.y;
        pd2.z=decstr[i].ptc.z-decstr[i].pto.z;
        diss[6]=angv(pd,pd2)*degrad;
        if(diss[6]<70.0 || diss[6]>140.0) continue;

        pd3.x=decstr[j].pth.x-decstr[j].ptn.x;
        pd3.y=decstr[j].pth.y-decstr[j].ptn.y;
        pd3.z=decstr[j].pth.z-decstr[j].ptn.z;
        diss[7]=angv(pd,pd3)*degrad;
                
        diss[8]=phi(decstr[j].ptn.x,decstr[j].ptn.y,decstr[j].ptn.z,decstr[j].pth.x,decstr[j].pth.y,decstr[j].pth.z,
            decstr[i].pto.x,decstr[i].pto.y,decstr[i].pto.z,decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z);
        if(diss[8]<120.0 || diss[8]>280.0) continue;
/*
        //second
        j=i+4;
        pd.x=decstr[j].pth.x-decstr[i].pto.x;
        pd.y=decstr[j].pth.y-decstr[i].pto.y;
        pd.z=decstr[j].pth.z-decstr[i].pto.z;
        diss[9]=norm(pd);
        if(diss[9]>=5.0 || diss[9]<1.6) continue;
 

        pd2.x=decstr[i].ptc.x-decstr[i].pto.x;
        pd2.y=decstr[i].ptc.y-decstr[i].pto.y;
        pd2.z=decstr[i].ptc.z-decstr[i].pto.z;
        diss[10]=angv(pd,pd2)*degrad;
        if(diss[10]<100.0 || diss[10]>170.0) continue;

        pd3.x=decstr[j].pth.x-decstr[j].ptn.x;
        pd3.y=decstr[j].pth.y-decstr[j].ptn.y;
        pd3.z=decstr[j].pth.z-decstr[j].ptn.z;
        diss[11]=angv(pd,pd3)*degrad;
        if(diss[11]<100.0) continue;
                
        diss[12]= phi(decstr[j].ptn.x,decstr[j].ptn.y,decstr[j].ptn.z,decstr[j].pth.x,decstr[j].pth.y,decstr[j].pth.z,
            decstr[i].pto.x,decstr[i].pto.y,decstr[i].pto.z,decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z);
        if(diss[12]<100.0 || diss[12]>230.0) continue;
 */
        double tval=0;
        for(k=0;k<5;k++)
        {       
            tval+=squgaussian(diss[k],nhochhsigma[k],nhochhmean[k]);
        }
        for(k=0;k<4;k++)
        {
            tval+=squgaussian(diss[5+k],nhocstdval[0][k],nhocmeanval[0][k]);
        }  
//        for(k=0;k<4;k++)
//        {
//            tval+=squgaussian(diss[9+k],nhocstdval[1][k],nhocmeanval[1][k]);
//        } 
        tval=exp(tval);
    //  printf("%3d %.3f\n",i,tval);
        if(tval>2e-4)
        {
            if(tval>1e-2) tval=1e-3;
            tval=lamda*(threshval[4]+log(tval));
            totenergy+=tval;
            for(k=0;k<=3;k++)
            {
                decstr[i+k].ssm='H';
            }
        }
    }
  
 
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////beta
    int inda,indb,indc,indd,inde,indf;
    int delseg;
    
    for(j=1;j<numseq-4;j++) 
    {
        for(k=j+3;k<numseq-1;k++)
        {            
            if( (decstr[j].ss2=='H' ) || (decstr[k].ss2=='H'))//diff
                continue;
    //        lamda=0.2;
            lamda=0.4;
            if((decstr[j].ss2=='E' ) && (decstr[k].ss2=='E' ))
                lamda+=0.1;
            else if((decstr[j].ss2=='E' ) || (decstr[k].ss2=='E' ))
                lamda+=0.2;
            else if(decstr[j].ss2=='E' && decstr[k].ss2=='E') 
                lamda+=1.6;
            else if(decstr[k].ss2=='E' || decstr[j].ss2=='E')
                lamda+=0.6;   
    //        lamda = 1.0; 
            tp[1]=setv(decstr[j].x,decstr[j].y,decstr[j].z);
            ap[1]=setv(decstr[k].x,decstr[k].y,decstr[k].z);
            tp[7]=minu(tp[1],ap[1]);
            if(norm(tp[7])>8.0) continue;
        
            inda=j;indb=k;
            tp[4]=setv(decstr[inda].ptn.x,decstr[inda].ptn.y,decstr[inda].ptn.z);
            tp[5]=setv(decstr[inda+1].ptn.x,decstr[inda+1].ptn.y,decstr[inda+1].ptn.z);
            tp[9]=setv(decstr[inda-1].pto.x,decstr[inda-1].pto.y,decstr[inda-1].pto.z);
            tp[10]=setv(decstr[inda].pto.x,decstr[inda].pto.y,decstr[inda].pto.z);
            ap[4]=setv(decstr[indb].ptn.x,decstr[indb].ptn.y,decstr[indb].ptn.z);
            ap[5]=setv(decstr[indb+1].ptn.x,decstr[indb+1].ptn.y,decstr[indb+1].ptn.z);
            ap[9]=setv(decstr[indb-1].pto.x,decstr[indb-1].pto.y,decstr[indb-1].pto.z);
            ap[10]=setv(decstr[indb].pto.x,decstr[indb].pto.y,decstr[indb].pto.z);
            kp[0]=minu(tp[9],ap[4]);//oi-1 nj
            kp[1]=minu(tp[5],ap[10]);//ni+1 oj
            kp[2]=minu(tp[4],ap[9]);//ni oj-1
            kp[3]=minu(tp[10],ap[5]);//oi nj+1
            kp[4]=minu(tp[9],ap[5]);//oi-1 nj+1
            kp[5]=minu(tp[5],ap[9]);//ni+1 oj-1
            kp[6]=minu(tp[4],ap[10]);//ni oj
            kp[7]=minu(tp[10],ap[4]);//oi nj
            for(m=0;m<8;m++)
            {
                diss[m]=norm(kp[m]);
            }
            delseg=0;
            diss[9]=diss[0]+diss[1];
            diss[8]=diss[2]+diss[3];
            if(diss[8]<diss[9])
            {
                delseg=1;
                diss[9]=diss[8];
            }
            diss[8]=diss[4]+diss[5];
            if(diss[8]<diss[9])
            {
                delseg=2;
                diss[9]=diss[8];
            }
            diss[8]=diss[6]+diss[7];
            if(diss[8]<diss[9])
            {
                delseg=3;
                diss[9]=diss[8];
            }
            if(diss[2*delseg]>6.8) continue;
            if(diss[2*delseg+1]>6.8) continue;
            if(delseg==0)
            {
                indc=k;indd=j-1;inde=j+1;indf=k;
            }
            else if(delseg==1)
            {
                indc=j;indd=k-1;inde=k+1;indf=j;
            }
            else if(delseg==3)
            {
                indc=k;indd=j;inde=j;indf=k;
            }
            else if(delseg==2)
            {
                indc=k+1;indd=j-1;inde=j+1;indf=k-1;
            }

            pd.x=decstr[indc].pth.x-decstr[indd].pto.x;
            pd.y=decstr[indc].pth.y-decstr[indd].pto.y;
            pd.z=decstr[indc].pth.z-decstr[indd].pto.z;
            diss[0]=norm(pd);
            pd2.x=decstr[indd].ptc.x-decstr[indd].pto.x;
            pd2.y=decstr[indd].ptc.y-decstr[indd].pto.y;
            pd2.z=decstr[indd].ptc.z-decstr[indd].pto.z;
            diss[1]=angv(pd,pd2)*degrad;
            pd3.x=decstr[indc].pth.x-decstr[indc].ptn.x;
            pd3.y=decstr[indc].pth.y-decstr[indc].ptn.y;
            pd3.z=decstr[indc].pth.z-decstr[indc].ptn.z;
            diss[2]=angv(pd,pd3)*degrad;
    //      diss[3]=bf.phi(decstr[indc].ptn.x,decstr[indc].ptn.y,decstr[indc].ptn.z,decstr[indc].pth.x,decstr[indc].pth.y,decstr[indc].pth.z,
    //          decstr[indd].pto.x,decstr[indd].pto.y,decstr[indd].pto.z,decstr[indd].ptc.x,decstr[indd].ptc.y,decstr[indd].ptc.z);
            ////////
            pd.x=decstr[inde].pth.x-decstr[indf].pto.x;
            pd.y=decstr[inde].pth.y-decstr[indf].pto.y;
            pd.z=decstr[inde].pth.z-decstr[indf].pto.z;
            diss[4]=norm(pd);               
            pd2.x=decstr[indf].ptc.x-decstr[indf].pto.x;
            pd2.y=decstr[indf].ptc.y-decstr[indf].pto.y;
            pd2.z=decstr[indf].ptc.z-decstr[indf].pto.z;
            diss[5]=angv(pd,pd2)*degrad;
            pd3.x=decstr[inde].pth.x-decstr[inde].ptn.x;
            pd3.y=decstr[inde].pth.y-decstr[inde].ptn.y;
            pd3.z=decstr[inde].pth.z-decstr[inde].ptn.z;
            diss[6]=angv(pd,pd3)*degrad;
    //      diss[7]=bf.phi(decstr[inde].ptn.x,decstr[inde].ptn.y,decstr[inde].ptn.z,decstr[inde].pth.x,decstr[inde].pth.y,decstr[inde].pth.z,
    //          decstr[indf].pto.x,decstr[indf].pto.y,decstr[indf].pto.z,decstr[indf].ptc.x,decstr[indf].ptc.y,decstr[indf].ptc.z);//n-h-o=c
    
            double tval=0;
            for(m=0;m<3;m++)
            {
                tval+=squgaussian(diss[m],nhocstdval[2+delseg][m],nhocmeanval[2+delseg][m]);
            }
            for(m=0;m<3;m++)
            {
                tval+=squgaussian(diss[4+m],nhocstdval[2+delseg][m],nhocmeanval[2+delseg][m]);
            }
            tval=exp(tval);

            if((delseg==0 && tval>5e-9)|| (delseg==1 && tval>5e-9)||(delseg==2 && tval>1e-9)|| (delseg==3 && tval>1e-9))
            {
                if(tval>1e-3) tval=1e-3;
                tval=lamda*(threshval[delseg]+log(tval))/threshval[delseg]*threshval[4];
                if(delseg==0 || delseg==2)//left of j
                {
                    if(tval>decstr[j].vpos)
                    {
                        decstr[j].vpos=tval;
                        if(decstr[j].indl>0) 
                        {       
                            if(decstr[decstr[j].indl].indr==j)//set to zero
                            {
                                decstr[decstr[j].indl].tpos=0;
                                decstr[decstr[j].indl].indr=-1;
                                if(decstr[decstr[j].indl].indl==-1)
                                    decstr[decstr[j].indl].ssm='C';
                            }
                        }
                        decstr[j].indl=k;
                        decstr[j].tpl=delseg;
                        decstr[k].ssm='E';
                        decstr[j].ssm='E';
                    }
                }
                else if(delseg==1 || delseg==3)//right of j
                {
                    if(tval>decstr[j].tpos)
                    {
                        decstr[j].tpos=tval;
                        if(decstr[j].indr>0) 
                        {
                            if(decstr[decstr[j].indr].indl==j)//set to zero
                            {
                                decstr[decstr[j].indr].vpos=0;
                                decstr[decstr[j].indr].indl=-1;
                                if(decstr[decstr[j].indr].indr==-1)
                                    decstr[decstr[j].indr].ssm='C';
                            }
                        }
                        decstr[j].indr=k;
                        decstr[j].tpr=delseg;
                        decstr[k].ssm='E';
                        decstr[j].ssm='E';
                    }
                }
                if(delseg==1 || delseg==2)//left of k
                {
                    if(tval>decstr[k].vpos)
                    {
                        decstr[k].vpos=tval;
                        if(decstr[k].indl>0) 
                        {
                            if(decstr[decstr[k].indl].indr==k)//set to zero
                            {
                                decstr[decstr[k].indl].tpos=0;
                                decstr[decstr[k].indl].indr=-1;
                                if(decstr[decstr[k].indl].indl==-1)
                                    decstr[decstr[k].indl].ssm='C';
                            }
                        }
                        decstr[k].indl=j;
                        decstr[k].tpl=delseg;
                        decstr[j].ssm='E';
                        decstr[k].ssm='E';
                    }
                }
                else if(delseg==0 || delseg==3)//right of k
                {
                    if(tval>decstr[k].tpos)
                    {
                        decstr[k].tpos=tval;
                        if(decstr[k].indr>0) 
                        {
                            if(decstr[decstr[k].indr].indl==k)//set to zero
                            {
                                decstr[decstr[k].indr].vpos=0;
                                decstr[decstr[k].indr].indl=-1;
                                if(decstr[decstr[k].indr].indr==-1)
                                    decstr[decstr[k].indr].ssm='C';
                            }
                        }
                        decstr[k].indr=j;
                        decstr[k].tpr=delseg;
                        decstr[j].ssm='E';
                        decstr[k].ssm='E';
                    }
                }   
            }
        }
    }
    double wthb=2.0;
    double wttot=3.0;
    int protype =0; //attenstion

    for(i=0;i<numseq;i++)
    {
        totenergy2+=decstr[i].vpos;
        totenergy2+=decstr[i].tpos;
        if(decstr[i].indl!=-1 && decstr[i].indr!=-1)
        {
            totenergy2+=wthb;//beta 3.0 alpha 2.0
            if((decstr[i].tpl==0 || decstr[i].tpl==1) && (decstr[i].tpr==0 || decstr[i].tpr==1))
                totenergy2+=1.0;//+decstr[i].vpos+decstr[i].tpos;
            else if((decstr[i].tpl==2 || decstr[i].tpl==3) && (decstr[i].tpr==2 || decstr[i].tpr==3))
                totenergy2+=1.5;//+decstr[i].vpos+decstr[i].tpos;
            else if((decstr[i].tpl==0 || decstr[i].tpl==1) && (decstr[i].tpr==2 || decstr[i].tpr==3))
                totenergy2-=0.0;
            else if((decstr[i].tpl==2 || decstr[i].tpl==3) && (decstr[i].tpr==0 || decstr[i].tpr==1))
                totenergy2-=0.0;
        }
        else if(decstr[i].indl!=-1 || decstr[i].indr!=-1)
        {
            totenergy2+=1.5;
        }
    }  
//    delete [] flagpos;
    return -(1.0*totenergy+wttot*totenergy2);//beta 1,2  alpha 1,3
}

double energyhbondnhoc2(vector<point3f> &decstr,int numseq)
{
    int i,j,k,l,m;
    point3d tp[20],ap[20],kp[20];
    point3d pd,pd2,pd3;
    double lamda;
    double diss[20];
    double totenergy=0;
    double totenergy2=0;

    static  double nhochhmean[22]={
        5.2172,8.7074,230.1397,4.2198,6.5391
        };
    static  double nhochhsigma[22]={
        0.3676,0.4375,10.2234,0.4080,0.3763
        };
    static  double nhocmeanval[][4]={//h-o h-o=c n-h-o n-h-o=c  
        2.85,89.0,110.5,199.5,//i+3
        2.00,147.0,159.0,160.0, //i+4  oldoh 
    //  2.83,89.0,110.0,201.5,//i+3 newoh
    //  2.00,148.0,159.0,155.0, //i+4

        2.00,155.0,164.0,180.0,//0
        2.00,155.0,164.0,180.0,//1
        2.00,151.0,163.0,192.0,//2
        2.00,151.0,163.0,192.0,//3

        };
    static  double nhocstdval[][4]={//h-o h-o=c n-h-o n-h-o=c 
        0.315504,7.697185,8.980366,7.932107,
        0.530859,10.582243,11.249764,25.360054,

        0.299730,11.770196,11.292558,68.955920,
        0.299730,11.770196,11.292558,68.955920,
        0.255088,12.376087,11.020081,69.165282,
        0.255088,12.376087,11.020081,69.165282,
        };

    double threshval[5]={19.113828,19.113828,20.723266,20.723266,16.118096};
    for(i=0;i<numseq;i++)
    {
        decstr[i].vpos=0;
        decstr[i].tpos=0;
        decstr[i].ssm='C';
        decstr[i].indl=-1;
        decstr[i].indr=-1;
        decstr[i].tpr=-1;
        decstr[i].tpl=-1;
    }

/*    for(int i=1;i<numseq;i++)
    {
        float pox = 0.0;
        float poy = 0.0;
        float poz = 0.0;
        tor2pos22(decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
                decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
                179.6715f*raddeg,1.229f,2.0961f,&pox,&poy,&poz);
        decstr[i].pto.x = pox ;
        decstr[i].pto.y = poy ;
        decstr[i].pto.z = poz ;
        cout<<"pto: "<<pox<<" "<<poy<<" "<<poz<<endl;
        float phx = 0.0;
        float phy = 0.0;
        float phz = 0.0;
        tor2pos22(decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,decstr[i].x,decstr[i].y,decstr[i].z,
            decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,179.8174f*raddeg,0.987f,2.0814f,&phx,&phy,&phz);//0.9919f,2.0574f
        decstr[i].pth.x = phx ;
        decstr[i].pth.y = phy ;
        decstr[i].pth.z = phz ; 
        cout<<"pth: "<<phx<<" "<<phy<<" "<<phz<<endl;           
    } */
    str2tor(decstr,numseq,3);
    tor2stroh(decstr,numseq);    

    for(i=1;i<numseq-4;i++)
    {
        j=i+3;  
        if((decstr[i].ss2=='E') || (decstr[j].ss2=='E' ))
            continue;
        l=0;
        for(k=0;k<4;k++)
        {
            if(decstr[i+k].ss2=='H') 
                l++;
        } 
    //  m=0;
    //    for(k=0;k<4;k++)
    //    {
    //        if(decstr[i+k].stype=='H') 
    //            m++;
    //    } 
    /*    if(l==3 && (decstr[i+1].ss2!='H' || decstr[i+2].ss2!='H'))
            continue;
        else if(l==2  && 
            ((decstr[i+1].ss2=='H' && decstr[i+2].ss2=='H') || (decstr[i+1].ss2!='H' && decstr[i+2].ss2!='H') ))
            continue;  */
        if(l==0 ) continue;
        else if(l==1 ) continue;
        else if(l==2 ) lamda=0.05;
    //    else if(l==2 ) lamda=0.2;
        else if(l==3) lamda=0.4;
    //    else if(l==3) lamda=0.5;
        else if(l==4) lamda=1.5;
    //  else if(m==3) lamda=0.6;
    //  else if(m==4) lamda=1.0;
        else lamda=0.4;  
        
        tp[1]=setv(decstr[i].x,decstr[i].y,decstr[i].z);
        ap[1]=setv(decstr[j].x,decstr[j].y,decstr[j].z);
        tp[7]=minu(tp[1],ap[1]);
        diss[1]=norm(tp[7]); 
        if(diss[1]>6.6) continue;

        tp[0]=setv(decstr[i-1].x,decstr[i-1].y,decstr[i-1].z);
        tp[2]=setv(decstr[i+1].x,decstr[i+1].y,decstr[i+1].z);
        ap[0]=setv(decstr[j-1].x,decstr[j-1].y,decstr[j-1].z);
        ap[2]=setv(decstr[j+1].x,decstr[j+1].y,decstr[j+1].z);
        kp[0]=minu(tp[2],ap[0]);
        diss[4]=norm(kp[0]);
        if(diss[4]<3.4 || diss[4]>4.2) continue;
        diss[0]=norm(minu(tp[0],ap[0]));
//      diss[1]=bf.norm(bf.minu(tp[0],ap[2]));          
//      diss[2]=bf.phi(tp[1].x,tp[1].y,tp[1].z,tp[2].x,tp[2].y,tp[2].z,ap[0].x,ap[0].y,ap[0].z,ap[1].x,ap[1].y,ap[1].z);
//      if(diss[2]>180) diss[2]-=180;
//      else diss[2]+=180;


//      tp[3]=bf.setv(decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z);
//      tp[5]=bf.setv(decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z);
//      ap[4]=bf.setv(decstr[j+1].ptn.x,decstr[j+1].ptn.y,decstr[j+1].ptn.z);
//      ap[5]=bf.setv(decstr[j].ptc.x,decstr[j].ptc.y,decstr[j].ptc.z); 
//      diss[3]=bf.norm(bf.minu(tp[3],ap[4]));
//      diss[4]=bf.norm(bf.minu(tp[5],ap[5]));

        //nhoc
        pd.x=decstr[j].pth.x-decstr[i].pto.x;
        pd.y=decstr[j].pth.y-decstr[i].pto.y;
        pd.z=decstr[j].pth.z-decstr[i].pto.z;
        diss[5]=norm(pd);
        if(diss[5]>=5.0 || diss[5]<1.6) continue;
 

        pd2.x=decstr[i].ptc.x-decstr[i].pto.x;
        pd2.y=decstr[i].ptc.y-decstr[i].pto.y;
        pd2.z=decstr[i].ptc.z-decstr[i].pto.z;
        diss[6]=angv(pd,pd2)*degrad;
        if(diss[6]<70.0 || diss[6]>140.0) continue;

        pd3.x=decstr[j].pth.x-decstr[j].ptn.x;
        pd3.y=decstr[j].pth.y-decstr[j].ptn.y;
        pd3.z=decstr[j].pth.z-decstr[j].ptn.z;
        diss[7]=angv(pd,pd3)*degrad;
                
        diss[8]=phi(decstr[j].ptn.x,decstr[j].ptn.y,decstr[j].ptn.z,decstr[j].pth.x,decstr[j].pth.y,decstr[j].pth.z,
            decstr[i].pto.x,decstr[i].pto.y,decstr[i].pto.z,decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z);
        if(diss[8]<120.0 || diss[8]>280.0) continue;

        //second
        j=i+4;
        pd.x=decstr[j].pth.x-decstr[i].pto.x;
        pd.y=decstr[j].pth.y-decstr[i].pto.y;
        pd.z=decstr[j].pth.z-decstr[i].pto.z;
        diss[9]=norm(pd);
        if(diss[9]>=5.0 || diss[9]<1.6) continue;
 

        pd2.x=decstr[i].ptc.x-decstr[i].pto.x;
        pd2.y=decstr[i].ptc.y-decstr[i].pto.y;
        pd2.z=decstr[i].ptc.z-decstr[i].pto.z;
        diss[10]=angv(pd,pd2)*degrad;
        if(diss[10]<100.0 || diss[10]>170.0) continue;

        pd3.x=decstr[j].pth.x-decstr[j].ptn.x;
        pd3.y=decstr[j].pth.y-decstr[j].ptn.y;
        pd3.z=decstr[j].pth.z-decstr[j].ptn.z;
        diss[11]=angv(pd,pd3)*degrad;
        if(diss[11]<100.0) continue;
                
        diss[12]= phi(decstr[j].ptn.x,decstr[j].ptn.y,decstr[j].ptn.z,decstr[j].pth.x,decstr[j].pth.y,decstr[j].pth.z,
            decstr[i].pto.x,decstr[i].pto.y,decstr[i].pto.z,decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z);
        if(diss[12]<100.0 || diss[12]>230.0) continue;

        double tval=0;
        for(k=0;k<1;k++)
        {       
            tval+=squgaussian(diss[k],nhochhsigma[k],nhochhmean[k]);
        }
        for(k=0;k<4;k++)
        {
            tval+=squgaussian(diss[5+k],nhocstdval[0][k],nhocmeanval[0][k]);
        }
        for(k=0;k<4;k++)
        {
            tval+=squgaussian(diss[9+k],nhocstdval[1][k],nhocmeanval[1][k]);
        }
        tval=exp(tval);
    //  printf("%3d %.3f\n",i,tval);
        if(tval>1e-7)
        {
            if(tval>1e-3) tval=1e-3;
            tval=lamda*(threshval[4]+log(tval));
            totenergy+=tval;
            for(k=0;k<=3;k++)
            {
                decstr[i+k].ssm='H';
            }
        }
    }
    
 
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////beta
    int inda,indb,indc,indd,inde,indf;
    int delseg;
    
    for(j=1;j<numseq-4;j++) 
    {
        for(k=j+3;k<numseq-1;k++)
        {   
            if( (decstr[j].ss2=='H' ) || (decstr[k].ss2=='H'))//diff
                continue;
    //        lamda=0.2;
            lamda=0.4;
            if((decstr[j].ss2=='E' ) && (decstr[k].ss2=='E' ))
                lamda+=0.1;
            else if((decstr[j].ss2=='E' ) || (decstr[k].ss2=='E' ))
                lamda+=0.2;
            else if(decstr[j].ss2=='E' && decstr[k].ss2=='E') 
                lamda+=1.6;
            else if(decstr[k].ss2=='E' || decstr[j].ss2=='E')
                lamda+=0.6;   
    //        lamda = 1.0;
            tp[1]=setv(decstr[j].x,decstr[j].y,decstr[j].z);
            ap[1]=setv(decstr[k].x,decstr[k].y,decstr[k].z);
            tp[7]=minu(tp[1],ap[1]);
            if(norm(tp[7])>8.0) continue;
        
            inda=j;indb=k;
            tp[4]=setv(decstr[inda].ptn.x,decstr[inda].ptn.y,decstr[inda].ptn.z);
            tp[5]=setv(decstr[inda+1].ptn.x,decstr[inda+1].ptn.y,decstr[inda+1].ptn.z);
            tp[9]=setv(decstr[inda-1].pto.x,decstr[inda-1].pto.y,decstr[inda-1].pto.z);
            tp[10]=setv(decstr[inda].pto.x,decstr[inda].pto.y,decstr[inda].pto.z);
            ap[4]=setv(decstr[indb].ptn.x,decstr[indb].ptn.y,decstr[indb].ptn.z);
            ap[5]=setv(decstr[indb+1].ptn.x,decstr[indb+1].ptn.y,decstr[indb+1].ptn.z);
            ap[9]=setv(decstr[indb-1].pto.x,decstr[indb-1].pto.y,decstr[indb-1].pto.z);
            ap[10]=setv(decstr[indb].pto.x,decstr[indb].pto.y,decstr[indb].pto.z);
            kp[0]=minu(tp[9],ap[4]);//oi-1 nj
            kp[1]=minu(tp[5],ap[10]);//ni+1 oj
            kp[2]=minu(tp[4],ap[9]);//ni oj-1
            kp[3]=minu(tp[10],ap[5]);//oi nj+1
            kp[4]=minu(tp[9],ap[5]);//oi-1 nj+1
            kp[5]=minu(tp[5],ap[9]);//ni+1 oj-1
            kp[6]=minu(tp[4],ap[10]);//ni oj
            kp[7]=minu(tp[10],ap[4]);//oi nj
            for(m=0;m<8;m++)
            {
                diss[m]=norm(kp[m]);
            }
            delseg=0;
            diss[9]=diss[0]+diss[1];
            diss[8]=diss[2]+diss[3];
            if(diss[8]<diss[9])
            {
                delseg=1;
                diss[9]=diss[8];
            }
            diss[8]=diss[4]+diss[5];
            if(diss[8]<diss[9])
            {
                delseg=2;
                diss[9]=diss[8];
            }
            diss[8]=diss[6]+diss[7];
            if(diss[8]<diss[9])
            {
                delseg=3;
                diss[9]=diss[8];
            }
            if(diss[2*delseg]>6.8) continue;
            if(diss[2*delseg+1]>6.8) continue;
            if(delseg==0)
            {
                indc=k;indd=j-1;inde=j+1;indf=k;
            }
            else if(delseg==1)
            {
                indc=j;indd=k-1;inde=k+1;indf=j;
            }
            else if(delseg==3)
            {
                indc=k;indd=j;inde=j;indf=k;
            }
            else if(delseg==2)
            {
                indc=k+1;indd=j-1;inde=j+1;indf=k-1;
            }

            pd.x=decstr[indc].pth.x-decstr[indd].pto.x;
            pd.y=decstr[indc].pth.y-decstr[indd].pto.y;
            pd.z=decstr[indc].pth.z-decstr[indd].pto.z;
            diss[0]=norm(pd);
            pd2.x=decstr[indd].ptc.x-decstr[indd].pto.x;
            pd2.y=decstr[indd].ptc.y-decstr[indd].pto.y;
            pd2.z=decstr[indd].ptc.z-decstr[indd].pto.z;
            diss[1]=angv(pd,pd2)*degrad;
            pd3.x=decstr[indc].pth.x-decstr[indc].ptn.x;
            pd3.y=decstr[indc].pth.y-decstr[indc].ptn.y;
            pd3.z=decstr[indc].pth.z-decstr[indc].ptn.z;
            diss[2]=angv(pd,pd3)*degrad;
    //      diss[3]=bf.phi(decstr[indc].ptn.x,decstr[indc].ptn.y,decstr[indc].ptn.z,decstr[indc].pth.x,decstr[indc].pth.y,decstr[indc].pth.z,
    //          decstr[indd].pto.x,decstr[indd].pto.y,decstr[indd].pto.z,decstr[indd].ptc.x,decstr[indd].ptc.y,decstr[indd].ptc.z);
            ////////
            pd.x=decstr[inde].pth.x-decstr[indf].pto.x;
            pd.y=decstr[inde].pth.y-decstr[indf].pto.y;
            pd.z=decstr[inde].pth.z-decstr[indf].pto.z;
            diss[4]=norm(pd);               
            pd2.x=decstr[indf].ptc.x-decstr[indf].pto.x;
            pd2.y=decstr[indf].ptc.y-decstr[indf].pto.y;
            pd2.z=decstr[indf].ptc.z-decstr[indf].pto.z;
            diss[5]=angv(pd,pd2)*degrad;
            pd3.x=decstr[inde].pth.x-decstr[inde].ptn.x;
            pd3.y=decstr[inde].pth.y-decstr[inde].ptn.y;
            pd3.z=decstr[inde].pth.z-decstr[inde].ptn.z;
            diss[6]=angv(pd,pd3)*degrad;
    //      diss[7]=bf.phi(decstr[inde].ptn.x,decstr[inde].ptn.y,decstr[inde].ptn.z,decstr[inde].pth.x,decstr[inde].pth.y,decstr[inde].pth.z,
    //          decstr[indf].pto.x,decstr[indf].pto.y,decstr[indf].pto.z,decstr[indf].ptc.x,decstr[indf].ptc.y,decstr[indf].ptc.z);//n-h-o=c
    
            double tval=0;
            for(m=0;m<3;m++)
            {
                tval+=squgaussian(diss[m],nhocstdval[2+delseg][m],nhocmeanval[2+delseg][m]);
            }
            for(m=0;m<3;m++)
            {
                tval+=squgaussian(diss[4+m],nhocstdval[2+delseg][m],nhocmeanval[2+delseg][m]);
            }
            tval=exp(tval);

            if((delseg==0 && tval>5e-9)|| (delseg==1 && tval>5e-9)||(delseg==2 && tval>1e-9)|| (delseg==3 && tval>1e-9))
            {
                if(tval>1e-3) tval=1e-3;
                tval=lamda*(threshval[delseg]+log(tval))/threshval[delseg]*threshval[4];
                if(delseg==0 || delseg==2)//left of j
                {
                    if(tval>decstr[j].vpos)
                    {
                        decstr[j].vpos=tval;
                        if(decstr[j].indl>0) 
                        {       
                            if(decstr[decstr[j].indl].indr==j)//set to zero
                            {
                                decstr[decstr[j].indl].tpos=0;
                                decstr[decstr[j].indl].indr=-1;
                                if(decstr[decstr[j].indl].indl==-1)
                                    decstr[decstr[j].indl].ssm='C';
                            }
                        }
                        decstr[j].indl=k;
                        decstr[j].tpl=delseg;
                        decstr[k].ssm='E';
                        decstr[j].ssm='E';
                    }
                }
                else if(delseg==1 || delseg==3)//right of j
                {
                    if(tval>decstr[j].tpos)
                    {
                        decstr[j].tpos=tval;
                        if(decstr[j].indr>0) 
                        {
                            if(decstr[decstr[j].indr].indl==j)//set to zero
                            {
                                decstr[decstr[j].indr].vpos=0;
                                decstr[decstr[j].indr].indl=-1;
                                if(decstr[decstr[j].indr].indr==-1)
                                    decstr[decstr[j].indr].ssm='C';
                            }
                        }
                        decstr[j].indr=k;
                        decstr[j].tpr=delseg;
                        decstr[k].ssm='E';
                        decstr[j].ssm='E';
                    }
                }
                if(delseg==1 || delseg==2)//left of k
                {
                    if(tval>decstr[k].vpos)
                    {
                        decstr[k].vpos=tval;
                        if(decstr[k].indl>0) 
                        {
                            if(decstr[decstr[k].indl].indr==k)//set to zero
                            {
                                decstr[decstr[k].indl].tpos=0;
                                decstr[decstr[k].indl].indr=-1;
                                if(decstr[decstr[k].indl].indl==-1)
                                    decstr[decstr[k].indl].ssm='C';
                            }
                        }
                        decstr[k].indl=j;
                        decstr[k].tpl=delseg;
                        decstr[j].ssm='E';
                        decstr[k].ssm='E';
                    }
                }
                else if(delseg==0 || delseg==3)//right of k
                {
                    if(tval>decstr[k].tpos)
                    {
                        decstr[k].tpos=tval;
                        if(decstr[k].indr>0) 
                        {
                            if(decstr[decstr[k].indr].indl==k)//set to zero
                            {
                                decstr[decstr[k].indr].vpos=0;
                                decstr[decstr[k].indr].indl=-1;
                                if(decstr[decstr[k].indr].indr==-1)
                                    decstr[decstr[k].indr].ssm='C';
                            }
                        }
                        decstr[k].indr=j;
                        decstr[k].tpr=delseg;
                        decstr[j].ssm='E';
                        decstr[k].ssm='E';
                    }
                }   
            }
        }
    }
    double wthb=2.0;
    double wttot=3.0;
    int protype =0; //attenstion
/*    if(protype==2) 
    {
        wthb=3.0;
        wttot=2.0;
    } */
    for(i=0;i<numseq;i++)
    {
        totenergy2+=decstr[i].vpos;
        totenergy2+=decstr[i].tpos;
        if(decstr[i].indl!=-1 && decstr[i].indr!=-1)
        {
            totenergy2+=wthb;//beta 3.0 alpha 2.0
            if((decstr[i].tpl==0 || decstr[i].tpl==1) && (decstr[i].tpr==0 || decstr[i].tpr==1))
                totenergy2+=1.0;//+decstr[i].vpos+decstr[i].tpos;
            else if((decstr[i].tpl==2 || decstr[i].tpl==3) && (decstr[i].tpr==2 || decstr[i].tpr==3))
                totenergy2+=1.5;//+decstr[i].vpos+decstr[i].tpos;
            else if((decstr[i].tpl==0 || decstr[i].tpl==1) && (decstr[i].tpr==2 || decstr[i].tpr==3))
                totenergy2-=0.0;
            else if((decstr[i].tpl==2 || decstr[i].tpl==3) && (decstr[i].tpr==0 || decstr[i].tpr==1))
                totenergy2-=0.0;
        }
        else if(decstr[i].indl!=-1 || decstr[i].indr!=-1)
        {
            totenergy2+=1.5;
        }
    }
    return -(1.0*totenergy+wttot*totenergy2);//beta 1,2  alpha 1,3
}

double energyhbondnhoc4(vector<point3f> &decstr,int numseq)
{
    int i,j,k,l,m;
    point3d tp[20],ap[20],kp[20];
    point3d pd,pd2,pd3;
    double lamda;
    double diss[20];
    double totenergy=0;
    double totenergy2=0;

    static  double nhochhmean[22]={
        5.2172,8.7074,230.1397,4.2198,6.5391
        };
    static  double nhochhsigma[22]={
        0.3676,0.4375,10.2234,0.4080,0.3763
        };
    static  double nhocmeanval[][4]={//h-o h-o=c n-h-o n-h-o=c  
        2.85,89.0,110.5,199.5,//i+3
        2.00,147.0,159.0,160.0, //i+4  oldoh 
    //  2.83,89.0,110.0,201.5,//i+3 newoh
    //  2.00,148.0,159.0,155.0, //i+4

        2.00,155.0,164.0,180.0,//0
        2.00,155.0,164.0,180.0,//1
        2.00,151.0,163.0,192.0,//2
        2.00,151.0,163.0,192.0,//3

        };
    static  double nhocstdval[][4]={//h-o h-o=c n-h-o n-h-o=c 
        0.315504,7.697185,8.980366,7.932107,
        0.530859,10.582243,11.249764,25.360054,

        0.299730,11.770196,11.292558,68.955920,
        0.299730,11.770196,11.292558,68.955920,
        0.255088,12.376087,11.020081,69.165282,
        0.255088,12.376087,11.020081,69.165282,
        };

    double threshval[5]={19.113828,19.113828,20.723266,20.723266,16.118096};
    for(i=0;i<numseq;i++)
    {
        decstr[i].vpos=0;
        decstr[i].tpos=0;
        decstr[i].ssm='C';
        decstr[i].indl=-1;
        decstr[i].indr=-1;
        decstr[i].tpr=-1;
        decstr[i].tpl=-1;
    }

/*    for(int i=1;i<numseq;i++)
    {
        float pox = 0.0;
        float poy = 0.0;
        float poz = 0.0;
        tor2pos22(decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
                decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
                179.6715f*raddeg,1.229f,2.0961f,&pox,&poy,&poz);
        decstr[i].pto.x = pox ;
        decstr[i].pto.y = poy ;
        decstr[i].pto.z = poz ;
        cout<<"pto: "<<pox<<" "<<poy<<" "<<poz<<endl;
        float phx = 0.0;
        float phy = 0.0;
        float phz = 0.0;
        tor2pos22(decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,decstr[i].x,decstr[i].y,decstr[i].z,
            decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,179.8174f*raddeg,0.987f,2.0814f,&phx,&phy,&phz);//0.9919f,2.0574f
        decstr[i].pth.x = phx ;
        decstr[i].pth.y = phy ;
        decstr[i].pth.z = phz ; 
        cout<<"pth: "<<phx<<" "<<phy<<" "<<phz<<endl;           
    } */
//    str2tor(decstr,numseq,3);
    tor2stroh(decstr,numseq);    

    for(i=1;i<numseq-4;i++)
    {
        j=i+3;  
        if((decstr[i].ss2=='E' &&  decstr[i].stype!='H') || (decstr[j].ss2=='E' &&  decstr[j].stype!='H'))
            continue;
        l=0;
        for(k=0;k<4;k++)
        {
            if(decstr[i+k].ss2=='H') 
                l++;
        }
        m=0;
        for(k=0;k<4;k++)
        {
            if(decstr[i+k].stype=='H') 
                m++;
        }
        if(l==3 && m<2 && (decstr[i+1].ss2!='H' || decstr[i+2].ss2!='H'))
            continue;
        else if(l==2 && m<2 && 
            ((decstr[i+1].ss2=='H' && decstr[i+2].ss2=='H') || (decstr[i+1].ss2!='H' && decstr[i+2].ss2!='H') ))
            continue;
        if(l==0 && m<3) continue;
        else if(l==1 && m<3) lamda=0.05;
        else if(l==2 && m<3) lamda=0.1;
//        else if(l==1 && m<3) continue;
//        else if(l==2 && m<3) lamda=0.05;
        else if(l==3) lamda=0.4;
        else if(l==4) lamda=1.5;
        else if(m==3) lamda=0.6;
        else if(m==4) lamda=1.0;
        else lamda=0.3;

/*        j=i+3;  
        if((decstr[i].ss2=='E') || (decstr[j].ss2=='E' ))
            continue;
        l=0;
        for(k=0;k<4;k++)
        {
            if(decstr[i+k].ss2=='H') 
                l++;
        } 
        if(l==0 ) continue;
        else if(l==1 ) continue;
        else if(l==2 ) lamda=0.05;
    //    else if(l==2 ) lamda=0.2;
        else if(l==3) lamda=0.4;
    //    else if(l==3) lamda=0.5;
        else if(l==4) lamda=1.5;
    //  else if(m==3) lamda=0.6;
    //  else if(m==4) lamda=1.0;
        else lamda=0.4;  */
        
        tp[1]=setv(decstr[i].x,decstr[i].y,decstr[i].z);
        ap[1]=setv(decstr[j].x,decstr[j].y,decstr[j].z);
        tp[7]=minu(tp[1],ap[1]);
        diss[1]=norm(tp[7]); 
        if(diss[1]>6.6) continue;

        tp[0]=setv(decstr[i-1].x,decstr[i-1].y,decstr[i-1].z);
        tp[2]=setv(decstr[i+1].x,decstr[i+1].y,decstr[i+1].z);
        ap[0]=setv(decstr[j-1].x,decstr[j-1].y,decstr[j-1].z);
        ap[2]=setv(decstr[j+1].x,decstr[j+1].y,decstr[j+1].z);
        kp[0]=minu(tp[2],ap[0]);
        diss[4]=norm(kp[0]);
        if(diss[4]<3.4 || diss[4]>4.2) continue;
        diss[0]=norm(minu(tp[0],ap[0]));
//      diss[1]=bf.norm(bf.minu(tp[0],ap[2]));          
//      diss[2]=bf.phi(tp[1].x,tp[1].y,tp[1].z,tp[2].x,tp[2].y,tp[2].z,ap[0].x,ap[0].y,ap[0].z,ap[1].x,ap[1].y,ap[1].z);
//      if(diss[2]>180) diss[2]-=180;
//      else diss[2]+=180;


//      tp[3]=bf.setv(decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z);
//      tp[5]=bf.setv(decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z);
//      ap[4]=bf.setv(decstr[j+1].ptn.x,decstr[j+1].ptn.y,decstr[j+1].ptn.z);
//      ap[5]=bf.setv(decstr[j].ptc.x,decstr[j].ptc.y,decstr[j].ptc.z); 
//      diss[3]=bf.norm(bf.minu(tp[3],ap[4]));
//      diss[4]=bf.norm(bf.minu(tp[5],ap[5]));

        //nhoc
        pd.x=decstr[j].pth.x-decstr[i].pto.x;
        pd.y=decstr[j].pth.y-decstr[i].pto.y;
        pd.z=decstr[j].pth.z-decstr[i].pto.z;
        diss[5]=norm(pd);
        if(diss[5]>=5.0 || diss[5]<1.6) continue;
 

        pd2.x=decstr[i].ptc.x-decstr[i].pto.x;
        pd2.y=decstr[i].ptc.y-decstr[i].pto.y;
        pd2.z=decstr[i].ptc.z-decstr[i].pto.z;
        diss[6]=angv(pd,pd2)*degrad;
        if(diss[6]<70.0 || diss[6]>140.0) continue;

        pd3.x=decstr[j].pth.x-decstr[j].ptn.x;
        pd3.y=decstr[j].pth.y-decstr[j].ptn.y;
        pd3.z=decstr[j].pth.z-decstr[j].ptn.z;
        diss[7]=angv(pd,pd3)*degrad;
                
        diss[8]=phi(decstr[j].ptn.x,decstr[j].ptn.y,decstr[j].ptn.z,decstr[j].pth.x,decstr[j].pth.y,decstr[j].pth.z,
            decstr[i].pto.x,decstr[i].pto.y,decstr[i].pto.z,decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z);
        if(diss[8]<120.0 || diss[8]>280.0) continue;

        //second
        j=i+4;
        pd.x=decstr[j].pth.x-decstr[i].pto.x;
        pd.y=decstr[j].pth.y-decstr[i].pto.y;
        pd.z=decstr[j].pth.z-decstr[i].pto.z;
        diss[9]=norm(pd);
        if(diss[9]>=5.0 || diss[9]<1.6) continue;
 

        pd2.x=decstr[i].ptc.x-decstr[i].pto.x;
        pd2.y=decstr[i].ptc.y-decstr[i].pto.y;
        pd2.z=decstr[i].ptc.z-decstr[i].pto.z;
        diss[10]=angv(pd,pd2)*degrad;
        if(diss[10]<100.0 || diss[10]>170.0) continue;

        pd3.x=decstr[j].pth.x-decstr[j].ptn.x;
        pd3.y=decstr[j].pth.y-decstr[j].ptn.y;
        pd3.z=decstr[j].pth.z-decstr[j].ptn.z;
        diss[11]=angv(pd,pd3)*degrad;
        if(diss[11]<100.0) continue;
                
        diss[12]= phi(decstr[j].ptn.x,decstr[j].ptn.y,decstr[j].ptn.z,decstr[j].pth.x,decstr[j].pth.y,decstr[j].pth.z,
            decstr[i].pto.x,decstr[i].pto.y,decstr[i].pto.z,decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z);
        if(diss[12]<100.0 || diss[12]>230.0) continue;

        double tval=0;
        for(k=0;k<1;k++)
        {       
            tval+=squgaussian(diss[k],nhochhsigma[k],nhochhmean[k]);
        }
        for(k=0;k<4;k++)
        {
            tval+=squgaussian(diss[5+k],nhocstdval[0][k],nhocmeanval[0][k]);
        }
        for(k=0;k<4;k++)
        {
            tval+=squgaussian(diss[9+k],nhocstdval[1][k],nhocmeanval[1][k]);
        }
        tval=exp(tval);
    //  printf("%3d %.3f\n",i,tval);
        if(tval>1e-7)
        {
            if(tval>1e-3) tval=1e-3;
            tval=lamda*(threshval[4]+log(tval));
            totenergy+=tval;
            for(k=0;k<=3;k++)
            {
                decstr[i+k].ssm='H';
            }
        }
    }
    
 
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////beta
    int inda,indb,indc,indd,inde,indf;
    int delseg;
    
    for(j=1;j<numseq-4;j++) 
    {
        for(k=j+3;k<numseq-1;k++)
        {   
            if( (decstr[j].ss2=='H' && decstr[j].stype!='E') || (decstr[k].ss2=='H' && decstr[k].stype!='E'))//diff
                continue;
            lamda=0.2;
            if((decstr[j].ss2=='E' && decstr[j].stype=='H') && (decstr[k].ss2=='E' && decstr[k].stype=='H'))
                lamda+=0.1;
            else if((decstr[j].ss2=='E' && decstr[j].stype=='H') || (decstr[k].ss2=='E' && decstr[k].stype=='H'))
                lamda+=0.2;
            else if(decstr[j].ss2=='E' && decstr[k].ss2=='E') 
                lamda+=1.6;
            else if(decstr[k].ss2=='E' || decstr[j].ss2=='E')
                lamda+=0.6;               
/*            if( (decstr[j].ss2=='H' ) || (decstr[k].ss2=='H'))//diff
                continue;
    //        lamda=0.2;
            lamda=0.4;
            if((decstr[j].ss2=='E' ) && (decstr[k].ss2=='E' ))
                lamda+=0.1;
            else if((decstr[j].ss2=='E' ) || (decstr[k].ss2=='E' ))
                lamda+=0.2;
            else if(decstr[j].ss2=='E' && decstr[k].ss2=='E') 
                lamda+=1.6;
            else if(decstr[k].ss2=='E' || decstr[j].ss2=='E')
                lamda+=0.6;   
    //        lamda = 1.0;  */
            tp[1]=setv(decstr[j].x,decstr[j].y,decstr[j].z);
            ap[1]=setv(decstr[k].x,decstr[k].y,decstr[k].z);
            tp[7]=minu(tp[1],ap[1]);
            if(norm(tp[7])>8.0) continue;
        
            inda=j;indb=k;
            tp[4]=setv(decstr[inda].ptn.x,decstr[inda].ptn.y,decstr[inda].ptn.z);
            tp[5]=setv(decstr[inda+1].ptn.x,decstr[inda+1].ptn.y,decstr[inda+1].ptn.z);
            tp[9]=setv(decstr[inda-1].pto.x,decstr[inda-1].pto.y,decstr[inda-1].pto.z);
            tp[10]=setv(decstr[inda].pto.x,decstr[inda].pto.y,decstr[inda].pto.z);
            ap[4]=setv(decstr[indb].ptn.x,decstr[indb].ptn.y,decstr[indb].ptn.z);
            ap[5]=setv(decstr[indb+1].ptn.x,decstr[indb+1].ptn.y,decstr[indb+1].ptn.z);
            ap[9]=setv(decstr[indb-1].pto.x,decstr[indb-1].pto.y,decstr[indb-1].pto.z);
            ap[10]=setv(decstr[indb].pto.x,decstr[indb].pto.y,decstr[indb].pto.z);
            kp[0]=minu(tp[9],ap[4]);//oi-1 nj
            kp[1]=minu(tp[5],ap[10]);//ni+1 oj
            kp[2]=minu(tp[4],ap[9]);//ni oj-1
            kp[3]=minu(tp[10],ap[5]);//oi nj+1
            kp[4]=minu(tp[9],ap[5]);//oi-1 nj+1
            kp[5]=minu(tp[5],ap[9]);//ni+1 oj-1
            kp[6]=minu(tp[4],ap[10]);//ni oj
            kp[7]=minu(tp[10],ap[4]);//oi nj
            for(m=0;m<8;m++)
            {
                diss[m]=norm(kp[m]);
            }
            delseg=0;
            diss[9]=diss[0]+diss[1];
            diss[8]=diss[2]+diss[3];
            if(diss[8]<diss[9])
            {
                delseg=1;
                diss[9]=diss[8];
            }
            diss[8]=diss[4]+diss[5];
            if(diss[8]<diss[9])
            {
                delseg=2;
                diss[9]=diss[8];
            }
            diss[8]=diss[6]+diss[7];
            if(diss[8]<diss[9])
            {
                delseg=3;
                diss[9]=diss[8];
            }
            if(diss[2*delseg]>6.8) continue;
            if(diss[2*delseg+1]>6.8) continue;
            if(delseg==0)
            {
                indc=k;indd=j-1;inde=j+1;indf=k;
            }
            else if(delseg==1)
            {
                indc=j;indd=k-1;inde=k+1;indf=j;
            }
            else if(delseg==3)
            {
                indc=k;indd=j;inde=j;indf=k;
            }
            else if(delseg==2)
            {
                indc=k+1;indd=j-1;inde=j+1;indf=k-1;
            }

            pd.x=decstr[indc].pth.x-decstr[indd].pto.x;
            pd.y=decstr[indc].pth.y-decstr[indd].pto.y;
            pd.z=decstr[indc].pth.z-decstr[indd].pto.z;
            diss[0]=norm(pd);
            pd2.x=decstr[indd].ptc.x-decstr[indd].pto.x;
            pd2.y=decstr[indd].ptc.y-decstr[indd].pto.y;
            pd2.z=decstr[indd].ptc.z-decstr[indd].pto.z;
            diss[1]=angv(pd,pd2)*degrad;
            pd3.x=decstr[indc].pth.x-decstr[indc].ptn.x;
            pd3.y=decstr[indc].pth.y-decstr[indc].ptn.y;
            pd3.z=decstr[indc].pth.z-decstr[indc].ptn.z;
            diss[2]=angv(pd,pd3)*degrad;
    //      diss[3]=bf.phi(decstr[indc].ptn.x,decstr[indc].ptn.y,decstr[indc].ptn.z,decstr[indc].pth.x,decstr[indc].pth.y,decstr[indc].pth.z,
    //          decstr[indd].pto.x,decstr[indd].pto.y,decstr[indd].pto.z,decstr[indd].ptc.x,decstr[indd].ptc.y,decstr[indd].ptc.z);
            ////////
            pd.x=decstr[inde].pth.x-decstr[indf].pto.x;
            pd.y=decstr[inde].pth.y-decstr[indf].pto.y;
            pd.z=decstr[inde].pth.z-decstr[indf].pto.z;
            diss[4]=norm(pd);               
            pd2.x=decstr[indf].ptc.x-decstr[indf].pto.x;
            pd2.y=decstr[indf].ptc.y-decstr[indf].pto.y;
            pd2.z=decstr[indf].ptc.z-decstr[indf].pto.z;
            diss[5]=angv(pd,pd2)*degrad;
            pd3.x=decstr[inde].pth.x-decstr[inde].ptn.x;
            pd3.y=decstr[inde].pth.y-decstr[inde].ptn.y;
            pd3.z=decstr[inde].pth.z-decstr[inde].ptn.z;
            diss[6]=angv(pd,pd3)*degrad;
    //      diss[7]=bf.phi(decstr[inde].ptn.x,decstr[inde].ptn.y,decstr[inde].ptn.z,decstr[inde].pth.x,decstr[inde].pth.y,decstr[inde].pth.z,
    //          decstr[indf].pto.x,decstr[indf].pto.y,decstr[indf].pto.z,decstr[indf].ptc.x,decstr[indf].ptc.y,decstr[indf].ptc.z);//n-h-o=c
    
            double tval=0;
            for(m=0;m<3;m++)
            {
                tval+=squgaussian(diss[m],nhocstdval[2+delseg][m],nhocmeanval[2+delseg][m]);
            }
            for(m=0;m<3;m++)
            {
                tval+=squgaussian(diss[4+m],nhocstdval[2+delseg][m],nhocmeanval[2+delseg][m]);
            }
            tval=exp(tval);

            if((delseg==0 && tval>5e-9)|| (delseg==1 && tval>5e-9)||(delseg==2 && tval>1e-9)|| (delseg==3 && tval>1e-9))
            {
                if(tval>1e-3) tval=1e-3;
                tval=lamda*(threshval[delseg]+log(tval))/threshval[delseg]*threshval[4];
                if(delseg==0 || delseg==2)//left of j
                {
                    if(tval>decstr[j].vpos)
                    {
                        decstr[j].vpos=tval;
                        if(decstr[j].indl>0) 
                        {       
                            if(decstr[decstr[j].indl].indr==j)//set to zero
                            {
                                decstr[decstr[j].indl].tpos=0;
                                decstr[decstr[j].indl].indr=-1;
                                if(decstr[decstr[j].indl].indl==-1)
                                    decstr[decstr[j].indl].ssm='C';
                            }
                        }
                        decstr[j].indl=k;
                        decstr[j].tpl=delseg;
                        decstr[k].ssm='E';
                        decstr[j].ssm='E';
                    }
                }
                else if(delseg==1 || delseg==3)//right of j
                {
                    if(tval>decstr[j].tpos)
                    {
                        decstr[j].tpos=tval;
                        if(decstr[j].indr>0) 
                        {
                            if(decstr[decstr[j].indr].indl==j)//set to zero
                            {
                                decstr[decstr[j].indr].vpos=0;
                                decstr[decstr[j].indr].indl=-1;
                                if(decstr[decstr[j].indr].indr==-1)
                                    decstr[decstr[j].indr].ssm='C';
                            }
                        }
                        decstr[j].indr=k;
                        decstr[j].tpr=delseg;
                        decstr[k].ssm='E';
                        decstr[j].ssm='E';
                    }
                }
                if(delseg==1 || delseg==2)//left of k
                {
                    if(tval>decstr[k].vpos)
                    {
                        decstr[k].vpos=tval;
                        if(decstr[k].indl>0) 
                        {
                            if(decstr[decstr[k].indl].indr==k)//set to zero
                            {
                                decstr[decstr[k].indl].tpos=0;
                                decstr[decstr[k].indl].indr=-1;
                                if(decstr[decstr[k].indl].indl==-1)
                                    decstr[decstr[k].indl].ssm='C';
                            }
                        }
                        decstr[k].indl=j;
                        decstr[k].tpl=delseg;
                        decstr[j].ssm='E';
                        decstr[k].ssm='E';
                    }
                }
                else if(delseg==0 || delseg==3)//right of k
                {
                    if(tval>decstr[k].tpos)
                    {
                        decstr[k].tpos=tval;
                        if(decstr[k].indr>0) 
                        {
                            if(decstr[decstr[k].indr].indl==k)//set to zero
                            {
                                decstr[decstr[k].indr].vpos=0;
                                decstr[decstr[k].indr].indl=-1;
                                if(decstr[decstr[k].indr].indr==-1)
                                    decstr[decstr[k].indr].ssm='C';
                            }
                        }
                        decstr[k].indr=j;
                        decstr[k].tpr=delseg;
                        decstr[j].ssm='E';
                        decstr[k].ssm='E';
                    }
                }   
            }
        }
    }
    double wthb=2.0;
    double wttot=3.0;
    int protype =0; //attenstion
/*    if(protype==2) 
    {
        wthb=3.0;
        wttot=2.0;
    } */
    for(i=0;i<numseq;i++)
    {
        totenergy2+=decstr[i].vpos;
        totenergy2+=decstr[i].tpos;
        if(decstr[i].indl!=-1 && decstr[i].indr!=-1)
        {
            totenergy2+=wthb;//beta 3.0 alpha 2.0
            if((decstr[i].tpl==0 || decstr[i].tpl==1) && (decstr[i].tpr==0 || decstr[i].tpr==1))
                totenergy2+=1.0;//+decstr[i].vpos+decstr[i].tpos;
            else if((decstr[i].tpl==2 || decstr[i].tpl==3) && (decstr[i].tpr==2 || decstr[i].tpr==3))
                totenergy2+=1.5;//+decstr[i].vpos+decstr[i].tpos;
            else if((decstr[i].tpl==0 || decstr[i].tpl==1) && (decstr[i].tpr==2 || decstr[i].tpr==3))
                totenergy2-=0.0;
            else if((decstr[i].tpl==2 || decstr[i].tpl==3) && (decstr[i].tpr==0 || decstr[i].tpr==1))
                totenergy2-=0.0;
        }
        else if(decstr[i].indl!=-1 || decstr[i].indr!=-1)
        {
            totenergy2+=1.5;
        }
    }
    return -(1.0*totenergy+wttot*totenergy2);//beta 1,2  alpha 1,3
}

double energyhbondnhoc2(vector<point3f> &decstr,int numseq,float &hbalpha,float &hbbeta)
{
    int i,j,k,l,m;
    point3d tp[20],ap[20],kp[20];
    point3d pd,pd2,pd3;
    double lamda;
    double diss[20];
    double totenergy=0;
    double totenergy2=0;
    double totval1=0;
    double totval2=0;    

    static  double nhochhmean[22]={
        5.2172,8.7074,230.1397,4.2198,6.5391
        };
    static  double nhochhsigma[22]={
        0.3676,0.4375,10.2234,0.4080,0.3763
        };
    static  double nhocmeanval[][4]={//h-o h-o=c n-h-o n-h-o=c  
        2.85,89.0,110.5,199.5,//i+3
        2.00,147.0,159.0,160.0, //i+4  oldoh 
    //  2.83,89.0,110.0,201.5,//i+3 newoh
    //  2.00,148.0,159.0,155.0, //i+4

        2.00,155.0,164.0,180.0,//0
        2.00,155.0,164.0,180.0,//1
        2.00,151.0,163.0,192.0,//2
        2.00,151.0,163.0,192.0,//3

        };
    static  double nhocstdval[][4]={//h-o h-o=c n-h-o n-h-o=c 
        0.315504,7.697185,8.980366,7.932107,
        0.530859,10.582243,11.249764,25.360054,

        0.299730,11.770196,11.292558,68.955920,
        0.299730,11.770196,11.292558,68.955920,
        0.255088,12.376087,11.020081,69.165282,
        0.255088,12.376087,11.020081,69.165282,
        };

    double threshval[5]={19.113828,19.113828,20.723266,20.723266,16.118096};
    int *flagpos=new int[numseq];    
    for(i=0;i<numseq;i++)
    {
        decstr[i].vpos=0;
        decstr[i].tpos=0;
        decstr[i].ssm='C';
        decstr[i].indl=-1;
        decstr[i].indr=-1;
        decstr[i].tpr=-1;
        decstr[i].tpl=-1;
        flagpos[i]=0;
    }


    for(i=1;i<numseq-4;i++)
    {
        j=i+3;  
/*        if((decstr[i].ss2=='E') || (decstr[j].ss2=='E' ))
            continue;
        l=0;
        for(k=0;k<4;k++)
        {
            if(decstr[i+k].ss2=='H') 
                l++;
        } 
    //  m=0;
    //    for(k=0;k<4;k++)
    //    {
    //        if(decstr[i+k].stype=='H') 
    //            m++;
    //    } 
    //    if(l==3 && (decstr[i+1].ss2!='H' || decstr[i+2].ss2!='H'))
    //        continue;
    //    else if(l==2  && 
    //        ((decstr[i+1].ss2=='H' && decstr[i+2].ss2=='H') || (decstr[i+1].ss2!='H' && decstr[i+2].ss2!='H') ))
    //        continue;  
        if(l==0 ) continue;
        else if(l==1 ) continue;
        else if(l==2 ) lamda=0.05;
    //    else if(l==2 ) lamda=0.2;
        else if(l==3) lamda=0.4;
    //    else if(l==3) lamda=0.5;
        else if(l==4) lamda=1.5;
    //  else if(m==3) lamda=0.6;
    //  else if(m==4) lamda=1.0;
        else lamda=0.4;  */
        
        tp[1]=setv(decstr[i].x,decstr[i].y,decstr[i].z);
        ap[1]=setv(decstr[j].x,decstr[j].y,decstr[j].z);
        tp[7]=minu(tp[1],ap[1]);
        diss[1]=norm(tp[7]); 
        if(diss[1]>6.6) continue;

        tp[0]=setv(decstr[i-1].x,decstr[i-1].y,decstr[i-1].z);
        tp[2]=setv(decstr[i+1].x,decstr[i+1].y,decstr[i+1].z);
        ap[0]=setv(decstr[j-1].x,decstr[j-1].y,decstr[j-1].z);
        ap[2]=setv(decstr[j+1].x,decstr[j+1].y,decstr[j+1].z);
        kp[0]=minu(tp[2],ap[0]);
        diss[4]=norm(kp[0]);
        if(diss[4]<3.4 || diss[4]>4.2) continue;
        diss[0]=norm(minu(tp[0],ap[0]));
//      diss[1]=bf.norm(bf.minu(tp[0],ap[2]));          
//      diss[2]=bf.phi(tp[1].x,tp[1].y,tp[1].z,tp[2].x,tp[2].y,tp[2].z,ap[0].x,ap[0].y,ap[0].z,ap[1].x,ap[1].y,ap[1].z);
//      if(diss[2]>180) diss[2]-=180;
//      else diss[2]+=180;


//      tp[3]=bf.setv(decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z);
//      tp[5]=bf.setv(decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z);
//      ap[4]=bf.setv(decstr[j+1].ptn.x,decstr[j+1].ptn.y,decstr[j+1].ptn.z);
//      ap[5]=bf.setv(decstr[j].ptc.x,decstr[j].ptc.y,decstr[j].ptc.z); 
//      diss[3]=bf.norm(bf.minu(tp[3],ap[4]));
//      diss[4]=bf.norm(bf.minu(tp[5],ap[5]));

        //nhoc
        pd.x=decstr[j].pth.x-decstr[i].pto.x;
        pd.y=decstr[j].pth.y-decstr[i].pto.y;
        pd.z=decstr[j].pth.z-decstr[i].pto.z;
        diss[5]=norm(pd);
        if(diss[5]>=5.0 || diss[5]<1.6) continue;
 

        pd2.x=decstr[i].ptc.x-decstr[i].pto.x;
        pd2.y=decstr[i].ptc.y-decstr[i].pto.y;
        pd2.z=decstr[i].ptc.z-decstr[i].pto.z;
        diss[6]=angv(pd,pd2)*degrad;
        if(diss[6]<70.0 || diss[6]>140.0) continue;

        pd3.x=decstr[j].pth.x-decstr[j].ptn.x;
        pd3.y=decstr[j].pth.y-decstr[j].ptn.y;
        pd3.z=decstr[j].pth.z-decstr[j].ptn.z;
        diss[7]=angv(pd,pd3)*degrad;
                
        diss[8]=phi(decstr[j].ptn.x,decstr[j].ptn.y,decstr[j].ptn.z,decstr[j].pth.x,decstr[j].pth.y,decstr[j].pth.z,
            decstr[i].pto.x,decstr[i].pto.y,decstr[i].pto.z,decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z);
        if(diss[8]<120.0 || diss[8]>280.0) continue;

        //second
        j=i+4;
        pd.x=decstr[j].pth.x-decstr[i].pto.x;
        pd.y=decstr[j].pth.y-decstr[i].pto.y;
        pd.z=decstr[j].pth.z-decstr[i].pto.z;
        diss[9]=norm(pd);
        if(diss[9]>=5.0 || diss[9]<1.6) continue;
 

        pd2.x=decstr[i].ptc.x-decstr[i].pto.x;
        pd2.y=decstr[i].ptc.y-decstr[i].pto.y;
        pd2.z=decstr[i].ptc.z-decstr[i].pto.z;
        diss[10]=angv(pd,pd2)*degrad;
        if(diss[10]<100.0 || diss[10]>170.0) continue;

        pd3.x=decstr[j].pth.x-decstr[j].ptn.x;
        pd3.y=decstr[j].pth.y-decstr[j].ptn.y;
        pd3.z=decstr[j].pth.z-decstr[j].ptn.z;
        diss[11]=angv(pd,pd3)*degrad;
        if(diss[11]<100.0) continue;
                
        diss[12]= phi(decstr[j].ptn.x,decstr[j].ptn.y,decstr[j].ptn.z,decstr[j].pth.x,decstr[j].pth.y,decstr[j].pth.z,
            decstr[i].pto.x,decstr[i].pto.y,decstr[i].pto.z,decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z);
        if(diss[12]<100.0 || diss[12]>230.0) continue;

        double tval=0;
        for(k=0;k<1;k++)
        {       
            tval+=squgaussian(diss[k],nhochhsigma[k],nhochhmean[k]);
        }
        for(k=0;k<4;k++)
        {
            tval+=squgaussian(diss[5+k],nhocstdval[0][k],nhocmeanval[0][k]);
        }
        for(k=0;k<4;k++)
        {
            tval+=squgaussian(diss[9+k],nhocstdval[1][k],nhocmeanval[1][k]);
        }
        tval=exp(tval);
    //  printf("%3d %.3f\n",i,tval);
        if(tval>1e-7)
        {
        //  if(tval>1e-2) tval=1e-2;
            tval=threshval[4]+log(tval);
            totval1+=tval;
            for(k=0;k<=3;k++)
            {
                decstr[i+k].ssm='H';
                flagpos[i+k]=1;
            }
        }
    }
    for(i=0;i<numseq;i++)
    {
        if(flagpos[i]==1)
            totenergy++;
    }    
 
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////beta
    int inda,indb,indc,indd,inde,indf;
    int delseg;
    
    for(j=1;j<numseq-4;j++) 
    {
        for(k=j+3;k<numseq-1;k++)
        {   
    /*        if( (decstr[j].ss2=='H' ) || (decstr[k].ss2=='H'))//diff
                continue;
    //        lamda=0.2;
            lamda=0.4;
            if((decstr[j].ss2=='E' ) && (decstr[k].ss2=='E' ))
                lamda+=0.1;
            else if((decstr[j].ss2=='E' ) || (decstr[k].ss2=='E' ))
                lamda+=0.2;
            else if(decstr[j].ss2=='E' && decstr[k].ss2=='E') 
                lamda+=1.6;
            else if(decstr[k].ss2=='E' || decstr[j].ss2=='E')
                lamda+=0.6;   */
    //        lamda = 1.0;
            if(flagpos[j]==1 || flagpos[k]==1) continue;
            tp[1]=setv(decstr[j].x,decstr[j].y,decstr[j].z);
            ap[1]=setv(decstr[k].x,decstr[k].y,decstr[k].z);
            tp[7]=minu(tp[1],ap[1]);
            if(norm(tp[7])>8.0) continue;
        
            inda=j;indb=k;
            tp[4]=setv(decstr[inda].ptn.x,decstr[inda].ptn.y,decstr[inda].ptn.z);
            tp[5]=setv(decstr[inda+1].ptn.x,decstr[inda+1].ptn.y,decstr[inda+1].ptn.z);
            tp[9]=setv(decstr[inda-1].pto.x,decstr[inda-1].pto.y,decstr[inda-1].pto.z);
            tp[10]=setv(decstr[inda].pto.x,decstr[inda].pto.y,decstr[inda].pto.z);
            ap[4]=setv(decstr[indb].ptn.x,decstr[indb].ptn.y,decstr[indb].ptn.z);
            ap[5]=setv(decstr[indb+1].ptn.x,decstr[indb+1].ptn.y,decstr[indb+1].ptn.z);
            ap[9]=setv(decstr[indb-1].pto.x,decstr[indb-1].pto.y,decstr[indb-1].pto.z);
            ap[10]=setv(decstr[indb].pto.x,decstr[indb].pto.y,decstr[indb].pto.z);
            kp[0]=minu(tp[9],ap[4]);//oi-1 nj
            kp[1]=minu(tp[5],ap[10]);//ni+1 oj
            kp[2]=minu(tp[4],ap[9]);//ni oj-1
            kp[3]=minu(tp[10],ap[5]);//oi nj+1
            kp[4]=minu(tp[9],ap[5]);//oi-1 nj+1
            kp[5]=minu(tp[5],ap[9]);//ni+1 oj-1
            kp[6]=minu(tp[4],ap[10]);//ni oj
            kp[7]=minu(tp[10],ap[4]);//oi nj
            for(m=0;m<8;m++)
            {
                diss[m]=norm(kp[m]);
            }
            delseg=0;
            diss[9]=diss[0]+diss[1];
            diss[8]=diss[2]+diss[3];
            if(diss[8]<diss[9])
            {
                delseg=1;
                diss[9]=diss[8];
            }
            diss[8]=diss[4]+diss[5];
            if(diss[8]<diss[9])
            {
                delseg=2;
                diss[9]=diss[8];
            }
            diss[8]=diss[6]+diss[7];
            if(diss[8]<diss[9])
            {
                delseg=3;
                diss[9]=diss[8];
            }
            if(diss[2*delseg]>6.8) continue;
            if(diss[2*delseg+1]>6.8) continue;
            if(delseg==0)
            {
                indc=k;indd=j-1;inde=j+1;indf=k;
            }
            else if(delseg==1)
            {
                indc=j;indd=k-1;inde=k+1;indf=j;
            }
            else if(delseg==3)
            {
                indc=k;indd=j;inde=j;indf=k;
            }
            else if(delseg==2)
            {
                indc=k+1;indd=j-1;inde=j+1;indf=k-1;
            }

            pd.x=decstr[indc].pth.x-decstr[indd].pto.x;
            pd.y=decstr[indc].pth.y-decstr[indd].pto.y;
            pd.z=decstr[indc].pth.z-decstr[indd].pto.z;
            diss[0]=norm(pd);
            pd2.x=decstr[indd].ptc.x-decstr[indd].pto.x;
            pd2.y=decstr[indd].ptc.y-decstr[indd].pto.y;
            pd2.z=decstr[indd].ptc.z-decstr[indd].pto.z;
            diss[1]=angv(pd,pd2)*degrad;
            pd3.x=decstr[indc].pth.x-decstr[indc].ptn.x;
            pd3.y=decstr[indc].pth.y-decstr[indc].ptn.y;
            pd3.z=decstr[indc].pth.z-decstr[indc].ptn.z;
            diss[2]=angv(pd,pd3)*degrad;
    //      diss[3]=bf.phi(decstr[indc].ptn.x,decstr[indc].ptn.y,decstr[indc].ptn.z,decstr[indc].pth.x,decstr[indc].pth.y,decstr[indc].pth.z,
    //          decstr[indd].pto.x,decstr[indd].pto.y,decstr[indd].pto.z,decstr[indd].ptc.x,decstr[indd].ptc.y,decstr[indd].ptc.z);
            ////////
            pd.x=decstr[inde].pth.x-decstr[indf].pto.x;
            pd.y=decstr[inde].pth.y-decstr[indf].pto.y;
            pd.z=decstr[inde].pth.z-decstr[indf].pto.z;
            diss[4]=norm(pd);               
            pd2.x=decstr[indf].ptc.x-decstr[indf].pto.x;
            pd2.y=decstr[indf].ptc.y-decstr[indf].pto.y;
            pd2.z=decstr[indf].ptc.z-decstr[indf].pto.z;
            diss[5]=angv(pd,pd2)*degrad;
            pd3.x=decstr[inde].pth.x-decstr[inde].ptn.x;
            pd3.y=decstr[inde].pth.y-decstr[inde].ptn.y;
            pd3.z=decstr[inde].pth.z-decstr[inde].ptn.z;
            diss[6]=angv(pd,pd3)*degrad;
    //      diss[7]=bf.phi(decstr[inde].ptn.x,decstr[inde].ptn.y,decstr[inde].ptn.z,decstr[inde].pth.x,decstr[inde].pth.y,decstr[inde].pth.z,
    //          decstr[indf].pto.x,decstr[indf].pto.y,decstr[indf].pto.z,decstr[indf].ptc.x,decstr[indf].ptc.y,decstr[indf].ptc.z);//n-h-o=c
    
            double tval=0;
            for(m=0;m<3;m++)
            {
                tval+=squgaussian(diss[m],nhocstdval[2+delseg][m],nhocmeanval[2+delseg][m]);
            }
            for(m=0;m<3;m++)
            {
                tval+=squgaussian(diss[4+m],nhocstdval[2+delseg][m],nhocmeanval[2+delseg][m]);
            }
            tval=exp(tval);

            if((delseg==0 && tval>5e-9)|| (delseg==1 && tval>5e-9)||(delseg==2 && tval>1e-9)|| (delseg==3 && tval>1e-9))
            {
                if(tval>1e-3) tval=1e-3;
                tval=lamda*(threshval[delseg]+log(tval))/threshval[delseg]*threshval[4];
                if(delseg==0 || delseg==2)//left of j
                {
                    if(tval>decstr[j].vpos)
                    {
                        decstr[j].vpos=tval;
                        if(decstr[j].indl>0) 
                        {       
                            if(decstr[decstr[j].indl].indr==j)//set to zero
                            {
                                decstr[decstr[j].indl].tpos=0;
                                decstr[decstr[j].indl].indr=-1;
                                if(decstr[decstr[j].indl].indl==-1)
                                    decstr[decstr[j].indl].ssm='C';
                            }
                        }
                        decstr[j].indl=k;
                        decstr[j].tpl=delseg;
                        decstr[k].ssm='E';
                        decstr[j].ssm='E';
                    }
                }
                else if(delseg==1 || delseg==3)//right of j
                {
                    if(tval>decstr[j].tpos)
                    {
                        decstr[j].tpos=tval;
                        if(decstr[j].indr>0) 
                        {
                            if(decstr[decstr[j].indr].indl==j)//set to zero
                            {
                                decstr[decstr[j].indr].vpos=0;
                                decstr[decstr[j].indr].indl=-1;
                                if(decstr[decstr[j].indr].indr==-1)
                                    decstr[decstr[j].indr].ssm='C';
                            }
                        }
                        decstr[j].indr=k;
                        decstr[j].tpr=delseg;
                        decstr[k].ssm='E';
                        decstr[j].ssm='E';
                    }
                }
                if(delseg==1 || delseg==2)//left of k
                {
                    if(tval>decstr[k].vpos)
                    {
                        decstr[k].vpos=tval;
                        if(decstr[k].indl>0) 
                        {
                            if(decstr[decstr[k].indl].indr==k)//set to zero
                            {
                                decstr[decstr[k].indl].tpos=0;
                                decstr[decstr[k].indl].indr=-1;
                                if(decstr[decstr[k].indl].indl==-1)
                                    decstr[decstr[k].indl].ssm='C';
                            }
                        }
                        decstr[k].indl=j;
                        decstr[k].tpl=delseg;
                        decstr[j].ssm='E';
                        decstr[k].ssm='E';
                    }
                }
                else if(delseg==0 || delseg==3)//right of k
                {
                    if(tval>decstr[k].tpos)
                    {
                        decstr[k].tpos=tval;
                        if(decstr[k].indr>0) 
                        {
                            if(decstr[decstr[k].indr].indl==k)//set to zero
                            {
                                decstr[decstr[k].indr].vpos=0;
                                decstr[decstr[k].indr].indl=-1;
                                if(decstr[decstr[k].indr].indr==-1)
                                    decstr[decstr[k].indr].ssm='C';
                            }
                        }
                        decstr[k].indr=j;
                        decstr[k].tpr=delseg;
                        decstr[j].ssm='E';
                        decstr[k].ssm='E';
                    }
                }   
            }
        }
    }
    double wthb=2.0;
    double wttot=3.0;
    int protype =0; //attenstion

/*    for(i=0;i<numseq;i++)
    {
        totenergy2+=decstr[i].vpos;
        totenergy2+=decstr[i].tpos;
        if(decstr[i].indl!=-1 && decstr[i].indr!=-1)
        {
            totenergy2+=wthb;//beta 3.0 alpha 2.0
            if((decstr[i].tpl==0 || decstr[i].tpl==1) && (decstr[i].tpr==0 || decstr[i].tpr==1))
                totenergy2+=1.0;//+decstr[i].vpos+decstr[i].tpos;
            else if((decstr[i].tpl==2 || decstr[i].tpl==3) && (decstr[i].tpr==2 || decstr[i].tpr==3))
                totenergy2+=1.5;//+decstr[i].vpos+decstr[i].tpos;
            else if((decstr[i].tpl==0 || decstr[i].tpl==1) && (decstr[i].tpr==2 || decstr[i].tpr==3))
                totenergy2-=0.0;
            else if((decstr[i].tpl==2 || decstr[i].tpl==3) && (decstr[i].tpr==0 || decstr[i].tpr==1))
                totenergy2-=0.0;
        }
        else if(decstr[i].indl!=-1 || decstr[i].indr!=-1)
        {
            totenergy2+=1.5;
        }
    }  */
    for(i=0;i<numseq;i++)
    {
        totval2+=decstr[i].vpos;
        totval2+=decstr[i].tpos;
        if(decstr[i].indl!=-1)
        {
            totenergy2++;
        }
        if(decstr[i].indr!=-1)
        {
            totenergy2++;
        }
    }
    delete[]flagpos;    
    hbalpha=totval1;
    hbbeta=totval2;
//    return -(1.0*totenergy+wttot*totenergy2);//beta 1,2  alpha 1,3
    return -(totenergy+totenergy2);
}

bool tor2pos22(float xi,float yi,float zi,float xj,float yj,float zj,float xk,
           float yk,float zk,float tang,float tleng,float tinner,float *xl,float *yl,float *zl)
{
    point3d p12,p23,e1,e2,e3;
    point3d p1,p2,p3;
    double tpangle,tcos,tsin,rmat[9];
    double q[4];
    p12.x=xi-xj;p12.y=yi-yj;p12.z=zi-zj;
    p23.x=xk-xj;p23.y=yk-yj;p23.z=zk-zj;
///*
    if(norm(p12)<epsilon*0.00001 || norm(p23)<epsilon*0.00001 || angv(p12,p23)<epsilon*0.001 || (PI-angv(p12,p23))<epsilon*0.001)
    {
        int imax=maxnormal(p23);
        if(imax==0)
        {
            yj+=0.01f;
        }
        else if(imax==1)
        {
            zj+=0.01f;
        }
        else
        {
            xj+=0.01f;
        }
        p12.x=xi-xj;p12.y=yi-yj;p12.z=zi-zj;
        p23.x=xk-xj;p23.y=yk-yj;p23.z=zk-zj;
        printf("make adjustment tor2pos22\n");
    }
//*/
    e2=unit(p23);
    e3=unit(p12);
    e1=prod(e2,e3);
    e1=unit(e1);
    if(norm(e1)<epsilon || norm(e2)<epsilon)
    {
        printf("wrong in tor2pos22 [%f %f %f] [%f %f %f] [%f %f %f]\n",xi,yi,zi,xj,yj,zj,xk,yk,zk);
        *xl=xk;*yl=yk;*zl=zk;
        return false;
    }
    p1=scalx(e2,tleng);
    tpangle=(PI-tinner)/2.0;
    tcos=cos(tpangle);
    tsin=sin(tpangle);
    q[0]=tcos;q[1]=tsin*e1.x;q[2]=tsin*e1.y;q[3]=tsin*e1.z;
    q2rot(q,rmat);
    p2=mmat(rmat,p1);

    tpangle=tang/2.0;
    tcos=cos(tpangle);
    tsin=sin(tpangle);
    q[0]=tcos;q[1]=tsin*e2.x;q[2]=tsin*e2.y;q[3]=tsin*e2.z;
    q2rot(q,rmat);
    p3=mmat(rmat,p2);

    *xl=p3.x+xk;*yl=p3.y+yk;*zl=p3.z+zk;
    return true;
}

int maxnormal(point3d p1)
{
    int i;
    double tem,temp[3]; 
    temp[0]=fabs(p1.x);
    temp[1]=fabs(p1.y);
    temp[2]=fabs(p1.z);
    i=0;
    tem=temp[0];
    if(temp[1]>tem)
    {
        i=1;
        tem=temp[1];
    }
    if(temp[2]>tem)
    {
        i=2;
        tem=temp[2];
    }
    return i;
}

void  q2rot(double *q,double rmax[])
{
    int i,j;
    for(i=0;i<3;i++)
    {
        for(j=0;j<3;j++)
        {
            rmax[i*3+j]=0;
        }
    }
    rmax[0]=q[0]*q[0]+q[1]*q[1]-q[2]*q[2]-q[3]*q[3];
    rmax[1]=2*(q[1]*q[2]-q[0]*q[3]);
    rmax[2]=2*(q[1]*q[3]+q[0]*q[2]);
    rmax[3]=2*(q[1]*q[2]+q[0]*q[3]);
    rmax[4]=q[0]*q[0]-q[1]*q[1]+q[2]*q[2]-q[3]*q[3];
    rmax[5]=2*(q[2]*q[3]-q[0]*q[1]);
    rmax[6]=2*(q[1]*q[3]-q[0]*q[2]);
    rmax[7]=2*(q[2]*q[3]+q[0]*q[1]);
    rmax[8]=q[0]*q[0]-q[1]*q[1]-q[2]*q[2]+q[3]*q[3];
}
point3d mmat(double mat[3][3],point3d p1)
{
    point3d temp;
    temp.x=p1.x*mat[0][0]+p1.y*mat[0][1]+p1.z*mat[0][2];
    temp.y=p1.x*mat[1][0]+p1.y*mat[1][1]+p1.z*mat[1][2];
    temp.z=p1.x*mat[2][0]+p1.y*mat[2][1]+p1.z*mat[2][2];
    return temp;
}
point3d mmat(double mat[9],point3d p1)
{
    point3d temp;
    temp.x=p1.x*mat[0]+p1.y*mat[1]+p1.z*mat[2];
    temp.y=p1.x*mat[3]+p1.y*mat[4]+p1.z*mat[5];
    temp.z=p1.x*mat[6]+p1.y*mat[7]+p1.z*mat[8];
    return temp;
}
bool tor2stroh(vector<point3f> &decstr,int seqnum)
{
    int i;
    int tind;
    point3s pn;
//    BasicFunc bf;
    bool flagts;
    bool flagwhole=true;
    for(i=1;i<seqnum;i++)
    {
    //  /*
        //atom o
//      flagts=bf.tor2pos22(decstr[i-1].ptn.x,decstr[i-1].ptn.y,decstr[i-1].ptn.z,decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
//          decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
//          decstr[i].tor[0]*raddeg-PI,lenco,angcaco*raddeg,&pn.x,&pn.y,&pn.z);
//      flagts=bf.tor2pos22(decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
//          decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
//          179.6715f*raddeg,1.2324f,2.1037f,&pn.x,&pn.y,&pn.z);
        flagts=tor2pos22(decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
            decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
            179.6715f*raddeg,1.229f,2.0961f,&pn.x,&pn.y,&pn.z);
        decstr[i-1].pto.x=pn.x;decstr[i-1].pto.y=pn.y;decstr[i-1].pto.z=pn.z;
        ///*
        //atom h
//      flagts=bf.tor2pos22(decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
//          decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,0,lennh,angcnh*raddeg,&pn.x,&pn.y,&pn.z);
        flagts=tor2pos22(decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,decstr[i].x,decstr[i].y,decstr[i].z,
            decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,179.8174f*raddeg,0.987f,2.0814f,&pn.x,&pn.y,&pn.z);//0.9919f,2.0574f
        if(!flagts)
        {
            flagwhole=false;
            printf("wrong front coordinates d %d\n",i);
        }
        decstr[i].pth.x=pn.x;decstr[i].pth.y=pn.y;decstr[i].pth.z=pn.z;
    //  */
    }
//  /*
    //atom o
    i=seqnum;
//  flagts=bf.tor2pos22(decstr[i-1].ptn.x,decstr[i-1].ptn.y,decstr[i-1].ptn.z,decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
//      decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,0,lenco,angcaco*raddeg,&pn.x,&pn.y,&pn.z);
    flagts=tor2pos22(decstr[i-1].ptn.x,decstr[i-1].ptn.y,decstr[i-1].ptn.z,decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
        decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,0.0f,1.2439f, 2.0855f,&pn.x,&pn.y,&pn.z);
    if(!flagts)
    {
        flagwhole=false;
        printf("wrong front coordinates o %d\n",i);
    }
    decstr[i-1].pto.x=pn.x;decstr[i-1].pto.y=pn.y;decstr[i-1].pto.z=pn.z;
//  */  
//  /*
    //atom h 
    i=0;
//  flagts=bf.tor2pos22(decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,decstr[i].x,decstr[i].y,decstr[i].z,
//      decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,0,lennh,anghnca*raddeg,&pn.x,&pn.y,&pn.z);
    flagts=tor2pos22(decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,decstr[i].x,decstr[i].y,decstr[i].z,
            decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,60*raddeg,0.987f,2.0306f,&pn.x,&pn.y,&pn.z);//0.9972f
    if(!flagts)
    {
        flagwhole=false;
        printf("wrong front coordinates %d\n",i);
    }
    decstr[i].pth.x=pn.x;decstr[i].pth.y=pn.y;decstr[i].pth.z=pn.z;
    //atom cb new
    for(i=0;i<seqnum;i++)
    {
        tind=aminoid(decstr[i].aaa);
        if(tind>19) tind=19;
        flagts=tor2pos22(decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
            decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,decstr[i].x,decstr[i].y,decstr[i].z,
            cbsta[tind][2]*raddeg,cbsta[tind][0],cbsta[tind][1],&pn.x,&pn.y,&pn.z);
        if(!flagts)
        {
            flagwhole=false;
            printf("wrong front coordinates cb2 %d\n",i);
        }
        decstr[i].ptb.x=pn.x;decstr[i].ptb.y=pn.y;decstr[i].ptb.z=pn.z;
        if(decstr[i].aaa=='G')
        {
            decstr[i].ptb.x=decstr[i].x;decstr[i].ptb.y=decstr[i].y;decstr[i].ptb.z=decstr[i].z;
        }
    }
    //ha
    for(i=0;i<seqnum;i++)
    {
        flagts=tor2pos22(decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,
            decstr[i].x,decstr[i].y,decstr[i].z,tornccaha*raddeg,lencaha,angccaha*raddeg,&pn.x,&pn.y,&pn.z);
        if(!flagts)
        {
            flagwhole=false;
            printf("wrong front coordinates f %d\n",i);
        }
        decstr[i].ptha.x=pn.x;decstr[i].ptha.y=pn.y;decstr[i].ptha.z=pn.z;
    }
    return flagwhole;
}

bool movementhel(vector<point3f> &tmstr,int numseq)
{
    vector<sssegment> sse;
    int numsse=0;
    genesse(tmstr,numseq,sse,numsse);
//    cout<<"VV "<<numsse<<" "<<<<endl;
    int i,j;
    int del=3;
    bool flagdir;
    bool flagtor;
    int tothelix=0;
    for(i=0;i<numsse;i++)
    {
        if(sse[i].ss=='H') tothelix++;
    }
//    cout<<"CC: "<<tothelix<<endl;
    if(tothelix==0)
    {
//        nummov[13][1]++;
        return false;
    }
    j=int(tothelix*randf0and1());
    tothelix=0;
    for(i=0;i<numsse;i++)
    {
        if(sse[i].ss=='H' && tothelix==j) break;
        else if(sse[i].ss=='H') tothelix++;
    }
    if(randf0and1()<0.5) flagdir=false;//head
    else flagdir=true;
    if((!flagdir && sse[i].init==0) || (flagdir && sse[i].term==numseq-1)) 
    {
//        nummov[13][1]++;
        return false;
    }
//    vector<point3f> tmstr;
//    tmstr = decstr;
//    memcpy(tmstr,decstr,numseq*sizeof(point3f));
//    memcpy(lp2.nn,lp1.nn,60*sizeof(int));
//    for(j=0;j<60;j++)
//    {
//        memcpy(lp2.indn[j],lp1.indn[j],numseq*sizeof(int));
//        memcpy(lp2.indf[j],lp1.indf[j],numseq*sizeof(double));
//    }
    point3s pt,pn,pc;
    j=-1;
    if(!flagdir) 
    {   
        tothelix=sse[i].init-1;
        if(tothelix<0)
        {
    //        nummov[13][1]++;
            return false;
        }
        if(i>0 && sse[i-1].ss=='C') j=sse[i-1].init;
        tmstr[tothelix+1].tor[2]=2*torideal[tmstr[tothelix+1].iaa][0]+1;
        tmstr[tothelix+1].tor[1]=180.0;
        tmstr[tothelix+1].tor[0]=2*torideal[tmstr[tothelix].iaa][1]+1;
        i=tothelix;
        flagtor=tor2pos22(tmstr[i+1].ptc.x,tmstr[i+1].ptc.y,tmstr[i+1].ptc.z,tmstr[i+1].x,tmstr[i+1].y,tmstr[i+1].z,tmstr[i+1].ptn.x,
        tmstr[i+1].ptn.y,tmstr[i+1].ptn.z,tmstr[i+1].tor[2]*raddeg,tmstr[i+1].len[0],tmstr[i+1].ang[1]*raddeg,&pn.x,&pn.y,&pn.z);
        if(!flagtor)
        {
            printf("wrong %d back c\n",tothelix);
        }
        tmstr[i].ptc.x=pn.x;tmstr[i].ptc.y=pn.y;tmstr[i].ptc.z=pn.z;
        flagtor=tor2pos22(tmstr[i+1].x,tmstr[i+1].y,tmstr[i+1].z,tmstr[i+1].ptn.x,tmstr[i+1].ptn.y,tmstr[i+1].ptn.z,tmstr[i].ptc.x,
        tmstr[i].ptc.y,tmstr[i].ptc.z,tmstr[i+1].tor[1]*raddeg,tmstr[i].len[2],tmstr[i+1].ang[0]*raddeg,&pt.x,&pt.y,&pt.z);
        if(!flagtor)
        {
            printf("wrong %d back ca\n",tothelix);
        }
        tmstr[i].x=pt.x;tmstr[i].y=pt.y;tmstr[i].z=pt.z;
        flagtor=tor2pos22(tmstr[i+1].ptn.x,tmstr[i+1].ptn.y,tmstr[i+1].ptn.z,tmstr[i].ptc.x,tmstr[i].ptc.y,tmstr[i].ptc.z,tmstr[i].x,
        tmstr[i].y,tmstr[i].z,tmstr[i+1].tor[0]*raddeg,tmstr[i].len[1],tmstr[i].ang[2]*raddeg,&pc.x,&pc.y,&pc.z);
        if(!flagtor)
        {
            printf("wrong %d back n\n",tothelix);
        }
        tmstr[i].ptn.x=pc.x;tmstr[i].ptn.y=pc.y;tmstr[i].ptn.z=pc.z;
        //[j,tothelix-1]
        if(tothelix==0) 
        {
            str2torp(tmstr,numseq,0,1);
    //      printf("head0 %d %d\n",0,1);
        }
        else if(tothelix==1)
        {
            str2torp(tmstr,numseq,0,2);
    //      printf("head1 %d %d\n",0,2);
        }
        else 
        {
            if(j==-1)
            {
                j=tothelix-del;
                if(j<1) j=1;
            }
            else
            {
                j=tothelix-j;
                if(j>=del) j=int(j*randf0and1());
                j=tothelix-j;
            }
            mcfragsweepLMP2(tmstr,numseq,j,tothelix);
    //      printf("head2 %d %d\n",j-1,tothelix);
            str2torp(tmstr,numseq,j-1,tothelix);
        }
    }
    else
    {
        tothelix=sse[i].term+1;
        if(tothelix>numseq-1)
        {
    //        nummov[13][1]++;
            return false;
        }
        if(i<numsse-1 && sse[i+1].ss=='C') j=sse[i+1].term;
        if(j==numseq-1) j--;
        tmstr[tothelix].tor[2]=2*torideal[tmstr[tothelix].iaa][0]+1;
        tmstr[tothelix].tor[1]=180.0;
        tmstr[tothelix].tor[0]=2*torideal[tmstr[tothelix-1].iaa][1]+1;
        i=tothelix;
        flagtor=tor2pos22(tmstr[i-1].ptn.x,tmstr[i-1].ptn.y,tmstr[i-1].ptn.z,tmstr[i-1].x,tmstr[i-1].y,tmstr[i-1].z,tmstr[i-1].ptc.x,
        tmstr[i-1].ptc.y,tmstr[i-1].ptc.z,tmstr[i].tor[0]*raddeg,tmstr[i].len[0],tmstr[i].ang[0]*raddeg,&pn.x,&pn.y,&pn.z);
        if(!flagtor)
        {
            printf("wrong %d for n\n",tothelix);
        }
        tmstr[i].ptn.x=pn.x;tmstr[i].ptn.y=pn.y;tmstr[i].ptn.z=pn.z;
        flagtor=tor2pos22(tmstr[i-1].x,tmstr[i-1].y,tmstr[i-1].z,tmstr[i-1].ptc.x,tmstr[i-1].ptc.y,tmstr[i-1].ptc.z,tmstr[i].ptn.x,
        tmstr[i].ptn.y,tmstr[i].ptn.z,tmstr[i].tor[1]*raddeg,tmstr[i].len[1],tmstr[i].ang[1]*raddeg,&pt.x,&pt.y,&pt.z);
        if(!flagtor)
        {
            printf("wrong %d for ca\n",tothelix);
        }
        tmstr[i].x=pt.x;tmstr[i].y=pt.y;tmstr[i].z=pt.z;
        flagtor=tor2pos22(tmstr[i-1].ptc.x,tmstr[i-1].ptc.y,tmstr[i-1].ptc.z,tmstr[i].ptn.x,tmstr[i].ptn.y,tmstr[i].ptn.z,tmstr[i].x,
        tmstr[i].y,tmstr[i].z,tmstr[i].tor[2]*raddeg,tmstr[i].len[2],tmstr[i].ang[2]*raddeg,&pc.x,&pc.y,&pc.z);
        if(!flagtor)
        {
            printf("wrong %d for c\n",tothelix);
        }
        tmstr[i].ptc.x=pc.x;tmstr[i].ptc.y=pc.y;tmstr[i].ptc.z=pc.z;
        //[tothelix+1,j]
        if(tothelix==numseq-1) 
        {
            str2torp(tmstr,numseq,numseq-2,numseq-1);
    //      printf("tail0 %d %d\n",numseq-2,numseq-1);
        }
        else if(tothelix==numseq-2)
        {
            str2torp(tmstr,numseq,numseq-3,numseq-1);
    //      printf("tail1 %d %d\n",numseq-3,numseq-1);
        }
        else 
        {
            if(j==-1)
            {
                j=tothelix+del;
                if(j>numseq-2) j=numseq-2;
            }
            else
            {
                j=j-tothelix;
                if(j>=del) j=int(j*randf0and1());
                j=tothelix+j;
            }
            mcfragsweepLMP2(tmstr,numseq,tothelix+1,j+1);
    //      printf("tail3 %d %d\n",tothelix,j+1);
            str2torp(tmstr,numseq,tothelix,j+1);
        }
    }
    
    return true;
}

bool genesse(vector<point3f> &decstr,int numseq,vector<sssegment> &sse,int &numsse)
{
    int i,j;
//    if(sse)
//    {
//        delete[]sse;
//        sse=NULL;
//    }
//    sse=new sssegment[numseq]; 
    vector<sssegment>().swap(sse);
    sse = vector<sssegment>(numseq);
    numsse=0;
    for(i=0;i<numseq;i++)
    {
        j=i;
        while(j<numseq && decstr[j].ss2==decstr[i].ss2)
        {
            j++;
        }
        sse[numsse].init=i;
        sse[numsse].term=j-1;
        sse[numsse].ss = decstr[i].ss2;
        numsse++;
        i=j-1;
    }
//    if(numsse!=0)
//    sse=(sssegment *)realloc(sse,numsse*sizeof(sssegment));
    return true;
}

bool genessex(vector<point3f> &decstr,int numseq,vector<sssegment> &sse,int &numsse,int &numshelix)
{
    int i,j;
//    if(sse)
//    {
//        delete[]sse;
//        sse=NULL;
//    }
//    sse=new sssegment[numseq]; 
//    vector<int>().swap(alphasind);
 //   alphasind =  vector<int>(numseq,0);
    vector<sssegment>().swap(sse);
    sse = vector<sssegment>(numseq);
    numsse=0;
    int numsh=0;
    for(i=0;i<numseq;i++)
    {
        j=i;
        while(j<numseq && decstr[j].ssm==decstr[i].ssm)
        {
            if(j == numshelix)
            {
                numsh = numsse;
            }
            j++;
        }
        sse[numsse].init=i;
        sse[numsse].term=j-1;
        sse[numsse].ss=decstr[i].ssm;
        numsse++;
        i=j-1;
    }
    numshelix = numsh;
//    if(numsse!=0)
//    sse=(sssegment *)realloc(sse,numsse*sizeof(sssegment));
    return true;
}

void calcabind(vector<int> &alphasind, vector<int> &alphaind, vector<int> &betaind,int numsse,vector<sssegment> &sse)
{
    int i;
/*  if(alphaind)
    {
        delete[]alphaind;
        alphaind=NULL;
    }
    if(betaind)
    {
        delete[]betaind;
        betaind=NULL;
    }
    if(alphasind)
    {
        delete[]alphasind;
        alphasind=NULL;
    } */
    int numsalpha=0;
    int numalpha=0;
    int numbeta=0;
    alphaind = vector<int>(numsse,0);
    betaind = vector<int>(numsse,0);
    alphasind = vector<int>(numsse,0);
//  alphaind=new int[numsse];
//  alphasind=new int[numsse];
//  betaind=new int[numsse];
    for(i=0;i<numsse;i++)
    {
        if(sse[i].ss=='H')
            alphaind[numalpha++]=i;
        else if(sse[i].ss=='E')
            betaind[numbeta++]=i;
    }
    for(i=0;i<numalpha-1;i++)
    {
        if(alphaind[i+1]==alphaind[i]+2)
        {
            alphasind[numsalpha++]=alphaind[i];
        }
    }
}

void str2torp(vector<point3f> &decstr,int seqnum,int istart,int iend)
{
    int i;
    point3d p12,p23;
//    BasicFunc bf;
    int realstart;
    if(istart==0)
        realstart=1;
    else realstart=istart;
    if(realstart==1)
    {
        if(decstr[0].tor[0]<0 || decstr[0].len[0]<0)
        {
            decstr[0].tor[0]=180.0;
            decstr[0].tor[1]=180.0;
            decstr[0].tor[2]=180.0;
            decstr[0].len[0]=lennc;
            decstr[0].ang[0]=angcacn;   
            decstr[0].ang[1]=angcnca;
        }   
        p12=setv(decstr[0].ptn.x-decstr[0].x,decstr[0].ptn.y-decstr[0].y,decstr[0].ptn.z-decstr[0].z);
        p23=setv(decstr[0].ptc.x-decstr[0].x,decstr[0].ptc.y-decstr[0].y,decstr[0].ptc.z-decstr[0].z);
        decstr[0].len[1]=norm(p12);
        decstr[0].len[2]=norm(p23);
        decstr[0].ang[2]=angv(p12,p23)*degrad;
    }
    for(i=realstart;i<=iend;i++)
    {
        p12=setv(decstr[i-1].x-decstr[i-1].ptc.x,decstr[i-1].y-decstr[i-1].ptc.y,decstr[i-1].z-decstr[i-1].ptc.z);
        p23=setv(decstr[i].ptn.x-decstr[i-1].ptc.x,decstr[i].ptn.y-decstr[i-1].ptc.y,decstr[i].ptn.z-decstr[i-1].ptc.z);
        decstr[i].len[0]=norm(p23);
        decstr[i].ang[0]=angv(p12,p23)*degrad;
        decstr[i].tor[0]=phi(decstr[i-1].ptn.x,decstr[i-1].ptn.y,decstr[i-1].ptn.z,decstr[i-1].x,decstr[i-1].y,
        decstr[i-1].z,decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z);
        
        p12=setv(decstr[i-1].ptc.x-decstr[i].ptn.x,decstr[i-1].ptc.y-decstr[i].ptn.y,decstr[i-1].ptc.z-decstr[i].ptn.z);
        p23=setv(decstr[i].x-decstr[i].ptn.x,decstr[i].y-decstr[i].ptn.y,decstr[i].z-decstr[i].ptn.z);
        decstr[i].len[1]=norm(p23);
        decstr[i].ang[1]=angv(p12,p23)*degrad;
        decstr[i].tor[1]=phi(decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,decstr[i-1].ptc.x,decstr[i-1].ptc.y,
        decstr[i-1].ptc.z,decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,decstr[i].x,decstr[i].y,decstr[i].z);
    
        p12=setv(decstr[i].ptn.x-decstr[i].x,decstr[i].ptn.y-decstr[i].y,decstr[i].ptn.z-decstr[i].z);
        p23=setv(decstr[i].ptc.x-decstr[i].x,decstr[i].ptc.y-decstr[i].y,decstr[i].ptc.z-decstr[i].z);
        decstr[i].len[2]=norm(p23);
        decstr[i].ang[2]=angv(p12,p23)*degrad;
        decstr[i].tor[2]=phi(decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,decstr[i].ptn.x,decstr[i].ptn.y,
        decstr[i].ptn.z,decstr[i].x,decstr[i].y,decstr[i].z,decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z);
    }
}

double energyhelixpacking(vector<point3f> &decstr,int numseq,vector<vector<float> > hhda,vector<vector<vector<float>>> haad,float &enehpk)
{
    int j,k,ii,jj,kk,ll;
    int iii,jjj,resval;
    double ainn,ator,dist;
    double hwt=2.879385; 
    point3d ps1,pe1,ps2,pe2,pc1,pc2,tp[3],ap[10];
    int numsse = 0 ;
//    vector<vector<float>> hhda;
//    hhda=loadhelixhelix("helixdistangle3.txt");
    vector<sssegment> sse;
    genesse(decstr,numseq,sse,numsse);
//    cout<<"5"<<endl;
//    vector<alphahelix> ahelix;
//    getahelix(decstr,numseq,ahelix);
    double totene=0;
    for(j=0;j<numsse-1;j++) if(sse[j].ss=='H')
    {
        if(sse[j].term-sse[j].init<3) continue;
            ap[0]=setv(decstr[sse[j].init].x,decstr[sse[j].init].y,decstr[sse[j].init].z);
            ap[1]=setv(decstr[sse[j].init+1].x,decstr[sse[j].init+1].y,decstr[sse[j].init+1].z);
            ap[2]=setv(decstr[sse[j].init+2].x,decstr[sse[j].init+2].y,decstr[sse[j].init+2].z);
            ap[3]=addv(ap[0],ap[2]);
            ap[4]=scalx(ap[3],hwt);
            ap[0]=addv(ap[4],ap[1]);
            ap[2]=scalx(ap[0],1.0/(1.0+2.0*hwt));
            ap[5]=setv(decstr[sse[j].term].x,decstr[sse[j].term].y,decstr[sse[j].term].z);
            ap[6]=setv(decstr[sse[j].term-1].x,decstr[sse[j].term-1].y,decstr[sse[j].term-1].z);
            ap[7]=setv(decstr[sse[j].term-2].x,decstr[sse[j].term-2].y,decstr[sse[j].term-2].z);
            ap[8]=addv(ap[5],ap[7]);
            ap[9]=scalx(ap[8],hwt);
            ap[5]=addv(ap[9],ap[6]);
            ap[7]=scalx(ap[5],1.0/(1.0+2.0*hwt));
            ps1.x=ap[2].x;
            ps1.y=ap[2].y;
            ps1.z=ap[2].z;
            pe1.x=ap[7].x;
            pe1.y=ap[7].y;
            pe1.z=ap[7].z;
//      ps1.x=decstr[sse[j].init].x;
//      ps1.y=decstr[sse[j].init].y;
//      ps1.z=decstr[sse[j].init].z;
//      pe1.x=decstr[sse[j].term].x;
//      pe1.y=decstr[sse[j].term].y;
//      pe1.z=decstr[sse[j].term].z;
        tp[0]=minu(pe1,ps1);
        for(k=j+1;k<numsse;k++) if(sse[k].ss=='H')
        {
            if(sse[k].term-sse[k].init<3) continue;
            ap[0]=setv(decstr[sse[k].init].x,decstr[sse[k].init].y,decstr[sse[k].init].z);
            ap[1]=setv(decstr[sse[k].init+1].x,decstr[sse[k].init+1].y,decstr[sse[k].init+1].z);
            ap[2]=setv(decstr[sse[k].init+2].x,decstr[sse[k].init+2].y,decstr[sse[k].init+2].z);
            ap[3]=addv(ap[0],ap[2]);
            ap[4]=scalx(ap[3],hwt);
            ap[0]=addv(ap[4],ap[1]);
            ap[2]=scalx(ap[0],1.0/(1.0+2.0*hwt));
            ap[5]=setv(decstr[sse[k].term].x,decstr[sse[k].term].y,decstr[sse[k].term].z);
            ap[6]=setv(decstr[sse[k].term-1].x,decstr[sse[k].term-1].y,decstr[sse[k].term-1].z);
            ap[7]=setv(decstr[sse[k].term-2].x,decstr[sse[k].term-2].y,decstr[sse[k].term-2].z);
            ap[8]=addv(ap[5],ap[7]);
            ap[9]=scalx(ap[8],hwt);
            ap[5]=addv(ap[9],ap[6]);
            ap[7]=scalx(ap[5],1.0/(1.0+2.0*hwt));
            ps2.x=ap[2].x;
            ps2.y=ap[2].y;
            ps2.z=ap[2].z;
            pe2.x=ap[7].x;
            pe2.y=ap[7].y;
            pe2.z=ap[7].z;
//          ps2.x=decstr[sse[k].init].x;
//          ps2.y=decstr[sse[k].init].y;
//          ps2.z=decstr[sse[k].init].z;
//          pe2.x=decstr[sse[k].term].x;
//          pe2.y=decstr[sse[k].term].y;
//          pe2.z=decstr[sse[k].term].z;
            tp[1]=minu(pe2,ps2);
            resval=linecross(ps1,pe1,ps2,pe2,&pc1,&pc2,&dist);
            if(resval!=0 && resval!=9 && resval!=14) continue;
//          printf("%.3f %.3f %.3f|%.3f %.3f %.3f|%.3f %.3f %.3f|%.3f %.3f %.3f\n",
//              ps1.x,ps1.y,ps1.z,pe1.x,pe1.y,pe1.z,ps2.x,ps2.y,ps2.z,pe2.x,pe2.y,pe2.z);
//          printf("j %d [%d %d] k %d [%d %d] dist %.3f\n",j,sse[j].init,sse[j].term,k,sse[k].init,sse[k].term,dist);
                //angledist
            kk=int(dist*2.0);
            if(kk<30)  
            {
                ainn=angv(tp[0],tp[1])*degrad;
            //  ator=bf.phi(ps1.x,ps1.y,ps1.z,pe1.x,pe1.y,pe1.z,ps2.x,ps2.y,ps2.z,pe2.x,pe2.y,pe2.z);//old
                ator=phi(ps1.x,ps1.y,ps1.z,pc1.x,pc1.y,pc1.z,pc2.x,pc2.y,pc2.z,ps2.x,ps2.y,ps2.z);
                if(ator>=180.0) ainn=360-ainn;
                ll=int(ainn/10.0);
                if(ll>35) ll=35;
                else if(ll<0) ll=0;
                totene+=hhda[kk][ll];
    //          printf("j %d k %d kk %d ll %d %f\n",j,k,kk,ll,hhda[kk][ll]);
            }
            //pair aa
            for(iii=sse[j].init;iii<=sse[j].term;iii++)
            {       
                if(decstr[iii].ssm!='H') continue;
                tp[0].x=decstr[iii].ptsg.x;
                tp[0].y=decstr[iii].ptsg.y;
                tp[0].z=decstr[iii].ptsg.z;
                ii=aminoid(decstr[iii].aaa);
                if(ii>19) continue;
                for(jjj=sse[k].init;jjj<=sse[k].term;jjj++)
                {
                    if(decstr[jjj].ssm!='H') continue;
                    tp[1].x=decstr[jjj].ptsg.x;
                    tp[1].y=decstr[jjj].ptsg.y;
                    tp[1].z=decstr[jjj].ptsg.z;
                    jj=aminoid(decstr[jjj].aaa);
                    if(jj>19) continue;
                    tp[2]=minu(tp[0],tp[1]);
                    dist=norm(tp[2]);
                    kk=int(dist*2.0);
                    if(kk>29) continue;
                    enehpk+=haad[ii][jj][kk];
    //              printf("ii %d jj %d kk %d %.3f\n",ii,jj,kk,haad[ii][jj][kk]);
                }//jjj
            }//iii
        }
    }
    enehpk*=0.05;
    return totene;
}

int linecross(point3d ps1,point3d pe1,point3d ps2,point3d pe2,point3d *pc1,point3d *pc2,double  *dist)
{
    point3d pse1,pse2,pss12,tfp;
    double k1,k2;
    double ang1,tdist;
    double dp1,dp2,dp3,dp4,minidp;
    int inddp;
    pse1=minu(pe1,ps1);
    pse2=minu(pe2,ps2);
    pss12=minu(ps1,ps2);
    ang1=angv(pse1,pse2);
    ///////////////////////////////////////////////////////////////////////////////////////////parallel
    if(ang1<epsilon || 180-ang1<epsilon)//parallel
    {
        k1=footpoint(ps1,pe1,ps2,&tfp,&tdist);
        if(tdist>=0 && tdist<=norm(pse1))
        {
            *pc1=ps2;
            *pc2=tfp;
            return 9;//good
        }
        k1=footpoint(ps1,pe1,pe2,&tfp,&tdist);
        if(tdist>=0 && tdist<=norm(pse1))
        {
            *pc1=pe2;
            *pc2=tfp;
            return 9;//good
        }
        k1=footpoint(ps2,pe2,ps1,&tfp,&tdist);
        if(tdist>=0 && tdist<=norm(pse2))
        {
            *pc1=ps1;
            *pc2=tfp;
            return 9;//good
        }
        k1=footpoint(ps2,pe2,pe1,&tfp,&tdist);
        if(tdist>=0 && tdist<=norm(pse2))
        {
            *pc1=pe1;
            *pc2=tfp;
            return 9;//good
        }
        dp1=norm(minu(ps1,ps2));
        inddp=1;
        minidp=dp1;
        *pc1=ps1;*pc2=ps2;
        dp2=norm(minu(ps1,pe2));
        dp3=norm(minu(pe1,ps2));
        dp4=norm(minu(pe1,pe2));
        if(dp2<minidp)
        {
            *pc1=ps1;*pc2=pe2;
            inddp=2;
            minidp=dp2;
        }
        if(dp3<minidp)
        {
            *pc1=pe1;*pc2=ps2;
            inddp=3;
            minidp=dp3;
        }
        if(dp4<minidp)
        {
            *pc1=pe1;*pc2=pe2;
            inddp=4;
            minidp=dp4;
        }
        *dist=minidp;
        return 9+inddp;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////coplane
    double tmat[9];
    tmat[0]=pe1.x-ps1.x;
    tmat[1]=pe1.y-ps1.y;
    tmat[2]=pe1.z-ps1.z;
    tmat[3]=ps2.x-ps1.x;
    tmat[4]=ps2.y-ps1.y;
    tmat[5]=ps2.z-ps1.z;
    tmat[6]=pe2.x-ps1.x;
    tmat[7]=pe2.y-ps1.y;
    tmat[8]=pe2.z-ps1.z;
    double tdet=sdet(tmat,3);
    point3d pdi=prod(pse1,pse2);
    int indmax=maxnormal(pdi);
    if(fabs(tdet)<epsilon)//coplane
    {
        if(indmax==2)
            k1=(pss12.x*pse2.y-pss12.y*pse2.x)/(pse1.y*pse2.x-pse1.x*pse2.y);
        else if(indmax==1)
            k1=(pss12.x*pse2.z-pss12.z*pse2.x)/(pse1.z*pse2.x-pse1.x*pse2.z);
        else
            k1=(pss12.z*pse2.y-pss12.y*pse2.z)/(pse1.y*pse2.z-pse1.z*pse2.y);
        *pc1=addv(ps1,scalx(pse1,k1));
        *pc2=*pc1;
        *dist=0;
        k2=norm(minu(*pc1,ps2))/norm(pse2);
        if(k1>=0 && k1<=1 && k2>=0 && k2<=1)
            return 14;//good
        else if(k1<0 && k2>=0 && k2<=1)
            return 15;
        else if(k1>1 && k2>=0 && k2<=1)
            return 16;
        else if(k1>=0 && k1<=1 && k2<0)
            return 17;
        else if(k1<0 && k2<0)
            return 18;
        else if(k1>1 && k2<0)
            return 19;
        else if(k1>=0 && k1<=1 && k2>1)
            return 20;
        else if(k1<0 && k2>1)
            return 21;
        else if(k1>1 && k2>1)
            return 22;
    }
    ////////////////////////////////////////////
    tfp=minu(ps2,ps1);
    *dist=dotv(tfp,pdi)/norm(pdi);
    point3d pfc1=addv(ps1,scalx(pdi,*dist/norm(pdi)));
    point3d pfc2=addv(pe1,scalx(pdi,*dist/norm(pdi)));
    pse1=minu(pfc2,pfc1);
    pss12=minu(pfc1,ps2);
    if(indmax==2)
        k1=(pss12.x*pse2.y-pss12.y*pse2.x)/(pse1.y*pse2.x-pse1.x*pse2.y);
    else if(indmax==1)
        k1=(pss12.x*pse2.z-pss12.z*pse2.x)/(pse1.z*pse2.x-pse1.x*pse2.z);
    else
        k1=(pss12.z*pse2.y-pss12.y*pse2.z)/(pse1.y*pse2.z-pse1.z*pse2.y);
    *pc2=addv(pfc1,scalx(pse1,k1));
    *pc1=minu(*pc2,scalx(pdi,*dist/norm(pdi)));
    *dist=fabs(*dist);
    k2=norm(minu(*pc1,ps2))/norm(pse2);
    if(k1>=0 && k1<=1 && k2>=0 && k2<=1)
        return 0;//good
    else if(k1<0 && k2>=0 && k2<=1)
        return 1;
    else if(k1>1 && k2>=0 && k2<=1)
        return 2;
    else if(k1>=0 && k1<=1 && k2<0)
        return 3;
    else if(k1<0 && k2<0)
        return 4;
    else if(k1>1 && k2<0)
        return 5;
    else if(k1>=0 && k1<=1 && k2>1)
        return 6;
    else if(k1<0 && k2>1)
        return 7;
    else if(k1>1 && k2>1)
        return 8;

    return 100;
}

 double sdet(double a[],int n)//content in a[] will be changed
  { 
    int i,j,k,is,js,l,u,v;
    double f,det,q,d;
    f=1.0; det=1.0;
    for (k=0; k<=n-2; k++)
      { 
        q=0.0;
        for (i=k; i<=n-1; i++)
        for (j=k; j<=n-1; j++)
          {
            l=i*n+j; 
            d=fabs(a[l]);
            if (d>q) 
            { q=d; is=i; js=j;}
          }
        if (q+1.0==1.0)
          { det=0.0; 
        return(det);}
        if (is!=k)
          { f=-f;
            for (j=k; j<=n-1; j++)
              { u=k*n+j;
                v=is*n+j;
                d=a[u]; 
                a[u]=a[v];
                a[v]=d;
              }
          }
        if (js!=k)
          { f=-f;
            for (i=k; i<=n-1; i++)
              {
                u=i*n+js; v=i*n+k;
                d=a[u]; a[u]=a[v]; a[v]=d;
              }
          }
        l=k*n+k;
        det=det*a[l];
        for (i=k+1; i<=n-1; i++)
          { 
            d=a[i*n+k]/a[l];
            for (j=k+1; j<=n-1; j++)
              { u=i*n+j;
                a[u]=a[u]-d*a[k*n+j];
              }
          }
      }
    det=f*det*a[n*n-1];
    return(det);
  }
point3d addv(point3d p1,point3d p2)
{
    point3d temp;
    temp.x=p1.x+p2.x;
    temp.y=p1.y+p2.y;
    temp.z=p1.z+p2.z;
    return temp;
}

//vertical dist tfp foot tdist-to ps1
double footpoint(point3d ps1,point3d pe1, point3d tp,point3d *tfp,double *tdist)
{
    point3d pse,pstp;
    double ang1,k;
    pse=minu(pe1,ps1);
    pstp=minu(tp,ps1);
    ang1=angv(pstp,pse);
    if(ang1<epsilon || 180-ang1<epsilon)
    {
        *tfp=tp;
        if(ang1<epsilon)
            *tdist=1;
        else *tdist=-1;
        *tdist*=norm(pstp);
        return 0.0;
    }
    k=dotv(pse,pstp)/dotv(pse,pse);
    *tfp=addv(ps1,scalx(pse,k));
    *tdist=k*norm(pse);
    pstp=minu(tp,*tfp);
    return norm(pstp);
}

void getahelix(vector<point3f> &decstr,int numseq,vector<alphahelix> &ahelix)
{
//    if(ahelix)
//    {
//        delete[]ahelix;
//        ahelix=NULL;
//    }
    int numahelix=0;
//    ahelix=new alphahelix[numseq];
    ahelix = vector<alphahelix>(numseq);
    int i,j,k;
    int tmptot;
    point3d tp1,td1,tdt;
//  BasicFunc bf;
    for(i=0;i<numseq;i++)
    {
    //    if(decstr[i].stype!='H')
        if(decstr[i].ss2!='H')
        {
            continue;
        }
        j=i;
        while(j<numseq && decstr[j].ss2==decstr[i].ss2)
        {
            j++;
        }
        if(j-i<=1)
        {
            continue;
        }
        ahelix[numahelix].seg.init=i;
        ahelix[numahelix].seg.term=j-1;
        if(j-i==2)
        {
        //  printf("the helix only hava 2\n");
        }
        tp1.x=0;tp1.y=0;tp1.z=0;
        tmptot=0;
        for(k=i;k<j;k++)
        {
            tp1.x+=decstr[k].x;
            tp1.y+=decstr[k].y;
            tp1.z+=decstr[k].z;
            tmptot++;
        }
        tp1.x/=double(tmptot);
        tp1.y/=double(tmptot);
        tp1.z/=double(tmptot);
        ahelix[numahelix].cen=tp1;
        td1.x=0;td1.y=0;td1.z=0;
        for(k=i;k<j;k++)
        {
            tdt.x=decstr[k].x-tp1.x;
            tdt.y=decstr[k].y-tp1.y;
            tdt.z=decstr[k].z-tp1.z;
            double box[3];
            box[0]=fabs(decstr[j-1].x-decstr[i].x);
            box[1]=fabs(decstr[j-1].y-decstr[i].y);
            box[2]=fabs(decstr[j-1].z-decstr[i].z);
            int flagbox=0;
            if(box[1]>box[0]) flagbox=1;
            if(box[2]>box[flagbox]) flagbox=2;

            if( (flagbox==0 && tdt.x*(decstr[j-1].x-decstr[i].x)<0) ||
                (flagbox==1 && tdt.y*(decstr[j-1].y-decstr[i].y)<0) ||
                (flagbox==2 && tdt.z*(decstr[j-1].z-decstr[i].z)<0) )
            {
                tdt.x=-tdt.x;
                tdt.y=-tdt.y;
                tdt.z=-tdt.z;
            }
            td1.x+=tdt.x;td1.y+=tdt.y;td1.z+=tdt.z;         
        }
        ahelix[numahelix].dir=unit(td1);
        tdt.x=decstr[i].x;
        tdt.y=decstr[i].y;
        tdt.z=decstr[i].z;
        ahelix[numahelix].dist1=(dotv(minu(tdt,ahelix[numahelix].cen),ahelix[numahelix].dir));
        tdt.x=decstr[j-1].x;
        tdt.y=decstr[j-1].y;
        tdt.z=decstr[j-1].z;
        ahelix[numahelix].dist2=(dotv(minu(tdt,ahelix[numahelix].cen),ahelix[numahelix].dir));
        if(ahelix[numahelix].dist1>=0 || ahelix[numahelix].dist2<=0) 
            printf("wrong dist1 %f %f\n",ahelix[numahelix].dist1,ahelix[numahelix].dist2);
        i=j-1;
        numahelix++;
    }
//    if(numahelix!=0)
//    ahelix=(alphahelix *)realloc(ahelix,numahelix*sizeof(alphahelix));
}
vector<vector<float>> loadhelixhelix(const char *filename)//distangle
{
    FILE *file;
    vector<vector<float>> hhda;
    int i,j;
    int ii,jj,kk;
    char oneline[300];
    hhda = vector<vector<float> >(30,vector<float>(36,0.0));
    file=fopen(filename,"rt");
    if(!file)
    {
        printf("no helixhelix file %s\n",filename);
        return hhda;
    }
//    if(!hhda)
//    {
//        hhda=new double*[30];
//        for(i=0;i<30;i++)
//            hhda[i]=new double[36];
//    }
    int totnum=0;
    double cutval=1e-5;
    double cutlogval=log(cutval);
    for(i=0;i<30;i++)
    {
        for(j=0;j<36;j++)
        {
            fgets(oneline,300,file);
            sscanf(oneline,"%d %d %d",&ii,&jj,&kk);
            hhda[i][j]=kk;
            totnum+=kk;
            if(i==29 && j==33) printf("helixhelix %d\n",kk);
        }
    }
    fclose(file);
    for(i=0;i<30;i++)
    {
        for(j=0;j<36;j++)
        {
            hhda[i][j]/=double(totnum);
            if(hhda[i][j]<cutval) hhda[i][j]=cutval;
            hhda[i][j]=cutlogval-log(hhda[i][j]);
    //      printf("%d %d %f\n",i,j,hhda[i][j]);
        }
    }
    return hhda;
}

vector<vector<vector<float>>> loadhelixpairsg(const char *filename)
{
    FILE *file;
    vector<vector<vector<float>>> haad;
    int i,j,k;
    int ii,jj,kk,ll;
    char oneline[300];
    haad = vector<vector<vector<float>>> (20, vector<vector<float>>(20, vector<float>(30,0.0)));
    file=fopen(filename,"rt");
    if(!file)
    {
        printf("no helixpairsg file %s\n",filename);
//        return false;
        return haad;
    }
//    if(haad==NULL)
//    {   
//        haad=new double**[20];
//        for(i=0;i<20;i++)
//        {   
//            haad[i]=new double*[20];
//            for(j=0;j<20;j++)
//            {
//                haad[i][j]=new double[30];
//            }   
//        }
//    }
    int totnum=0;
    
    for(i=0;i<20;i++)
    {
        for(j=0;j<20;j++)
        {
            for(k=0;k<30;k++)
            {
                fgets(oneline,300,file);
                sscanf(oneline,"%d %d %d %d",&ii,&jj,&kk,&ll);
                haad[i][j][k]=ll;
        //      if(k==29 && ll==0) printf("zero %d %d %d\n",i,j,k);
        //      if(i==19 && j==19 && k==29) printf("%d\n",ll);
                totnum+=ll;
            }
        }
    }
    fclose(file);
    double cutdist=1e-5;
    double tval1;
    double cutlogdist=log(cutdist);
    double minval=100000;
    for(i=0;i<20;i++)
    {
        for(j=0;j<20;j++)
        {
            for(k=0;k<30;k++)
            {
                haad[i][j][k]/=double(totnum);
                if(haad[i][j][k]<cutdist)
                {
                    haad[i][j][k]=0.0;
                    if(minval>haad[i][j][k]) minval=haad[i][j][k];
            //      printf("%2d %2d %2d %7d %9f\n",i,j,k,totnum,haad[i][j][k]);
                    continue;
                }
                tval1=cutlogdist-log(haad[i][j][k]);
                haad[i][j][k]=tval1;
                if(minval>haad[i][j][k]) minval=haad[i][j][k];
            //  printf("%2d %2d %2d %7d %9f\n",i,j,k,totnum,haad[i][j][k]);
            }
        }
    }
    printf("load helixpairsg minval %f\n",minval);
    /*
    FILE *file2=fopen("helixpaircut.txt","wt");
//  int tottot=0;
        for(i=0;i<20;i++)
        {
            for(j=0;j<20;j++)
            {
                fprintf(file2,"%2d %2d     ",i,j);
                for(k=2;k<28;k++)
                {
                    if(haad[i][j][k]<haad[i][j][k-1] && haad[i][j][k]<haad[i][j][k+1])
                        fprintf(file2,"%2d ",k);
                    else if(haad[i][j][k]==haad[i][j][k-1] && haad[i][j][k]<haad[i][j][k-2] && haad[i][j][k]<haad[i][j][k+1])
                        fprintf(file2,"%2d ",k);
                    else if(haad[i][j][k]<haad[i][j][k-1] && haad[i][j][k]<haad[i][j][k+2] && haad[i][j][k]==haad[i][j][k+1])
                    {
                    //  fprintf(file2,"%2d ",k);
                    }
                }
                fprintf(file2,"\n");
            }
        }
    fclose(file2);
//  */
    return haad;
//    return true;
}

bool mcmovementLMP(vector<point3f> &tmstr,int seqlen,int nums,int nume)
{
    int i,k;
    int trandpos,trand2;
    point3d rp;
    double rnorm;
    double rmax;
    bool flagdone;
    int numiter;

        //[2 7] [trandpos, trandpos+trand2-1] [1 numseq-2]
   //     trand2=2+int(6*Random());
   //     trandpos=1+int((seqlen-trand2-1)*Random());
        trandpos = nums;
        trand2 = nume - nums;

        int threshiter=trand2*200;
        if(threshiter>1000) threshiter=1000;
    /*    double treject=Random();
        double ssrat=0;
        for(i=trandpos;i<trandpos+trand2;i++)
        {
            if(tmstr[i].ssm!='C')
                ssrat+=1;
        }
        ssrat/=double(trand2);
        if(ssrat>epsilon && treject>0.5*(1-ssrat)+0.05)
        {
            summcsep[0][2]++;
            return false;
        } */
        //perturbation
        for(i=trandpos;i<trandpos+trand2;i++)
        {
            rmax=0.25+1.25*randf0and1();
            rp=ranv(rmax);
            tmstr[i].ptn.x+=rp.x;
            tmstr[i].ptn.y+=rp.y;
            tmstr[i].ptn.z+=rp.z;
            rmax=0.25+1.25*randf0and1();
            rp=ranv(rmax);
            tmstr[i].ptc.x+=rp.x;
            tmstr[i].ptc.y+=rp.y;
            tmstr[i].ptc.z+=rp.z;
            rmax=0.25+1.25*randf0and1();
            rp=ranv(rmax);
            tmstr[i].x+=rp.x;
            tmstr[i].y+=rp.y;
            tmstr[i].z+=rp.z;
        }   
        numiter=0;
        do{         
    /*        for(k=trandpos;k<trandpos+trand2;k++)
            {   
                singlemoveLMP(tmstr,k,int(randf0and1()*7.0),0);
                singlemoveLMP(tmstr,k,int(randf0and1()*7.0),0);
                singlemoveLMP(tmstr,k,int(randf0and1()*7.0),0);
                singlemoveLMP(tmstr,k,int(randf0and1()*7.0),0);
                singlemoveLMP(tmstr,k,int(randf0and1()*7.0),0);
                singlemoveLMP(tmstr,k,int(randf0and1()*7.0),0);
                singlemoveLMP(tmstr,k,int(randf0and1()*7.0),0);
            } */
            for(k=trandpos;k<trandpos+trand2;k++)
            {   
                    singlemoveLMPf(tmstr,k,int(randf0and1()*7.0));
                    singlemoveLMPf(tmstr,k,int(randf0and1()*7.0));
                    singlemoveLMPf(tmstr,k,int(randf0and1()*7.0));
                    singlemoveLMPf(tmstr,k,int(randf0and1()*7.0));
                    singlemoveLMPf(tmstr,k,int(randf0and1()*7.0));
                    singlemoveLMPf(tmstr,k,int(randf0and1()*7.0));
                    singlemoveLMPf(tmstr,k,int(randf0and1()*7.0));
            }            
            ////////////////////////////////////////////////////////////////////////////////////////inverse     
        /*    for(k=trandpos+trand2;k>trandpos;k--)
            {
                singlemoveLMP(tmstr,k,int(randf0and1()*7.0),1);
                singlemoveLMP(tmstr,k,int(randf0and1()*7.0),1);
                singlemoveLMP(tmstr,k,int(randf0and1()*7.0),1);
                singlemoveLMP(tmstr,k,int(randf0and1()*7.0),1);
                singlemoveLMP(tmstr,k,int(randf0and1()*7.0),1);
                singlemoveLMP(tmstr,k,int(randf0and1()*7.0),1);
                singlemoveLMP(tmstr,k,int(randf0and1()*7.0),1);             
            }   */
            for(k=trandpos+trand2;k>trandpos;k--)
            {
                    singlemoveLMPb(tmstr,k,int(randf0and1()*7.0));
                    singlemoveLMPb(tmstr,k,int(randf0and1()*7.0));
                    singlemoveLMPb(tmstr,k,int(randf0and1()*7.0));
                    singlemoveLMPb(tmstr,k,int(randf0and1()*7.0));
                    singlemoveLMPb(tmstr,k,int(randf0and1()*7.0));
                    singlemoveLMPb(tmstr,k,int(randf0and1()*7.0));
                    singlemoveLMPb(tmstr,k,int(randf0and1()*7.0));
                    
            }              
            ///////////////////////////////////////////////////////////////////////////////////////check
            flagdone=true;
            for(k=trandpos;k<=trandpos+trand2;k++)
            {
                //ca-1 ca
            //    rp=setv(tmstr[k].ptv[1].x-tmstr[k-1].ptv[1].x,tmstr[k].ptv[1].y-tmstr[k-1].ptv[1].y,tmstr[k].ptv[1].z-tmstr[k-1].ptv[1].z);
                rp=setv(tmstr[k].x-tmstr[k-1].x,tmstr[k].y-tmstr[k-1].y,tmstr[k].z-tmstr[k-1].z);
                rnorm=norm(rp);
                if(rnorm<lencaca-delcaca || rnorm>lencaca+delcaca)
                {
                    //printf("%d %d %d %d ca-1 ca %f lencaca %f\n",k,trandpos,trand2,numiter,rnorm,lencaca);
                    flagdone=false;
                    break;
                }
                //ca-1 n
            //    rp=setv(tmstr[k].ptv[0].x-tmstr[k-1].ptv[1].x,tmstr[k].ptv[0].y-tmstr[k-1].ptv[1].y,tmstr[k].ptv[0].z-tmstr[k-1].ptv[1].z);
                rp=setv(tmstr[k].ptn.x-tmstr[k-1].x,tmstr[k].ptn.y-tmstr[k-1].y,tmstr[k].ptn.z-tmstr[k-1].z);
                rnorm=norm(rp);
                if(rnorm<lencan1-delcan1 || rnorm>lencan1+delcan1)
                {
                    //printf("%d %d %d %d ca-1 n %f lencan1 %f\n",k,trandpos,trand2,numiter,rnorm,lencan1);
                    flagdone=false;
                    break;
                }
                //c-1 n
            //    rp=setv(tmstr[k].ptv[0].x-tmstr[k-1].ptv[2].x,tmstr[k].ptv[0].y-tmstr[k-1].ptv[2].y,tmstr[k].ptv[0].z-tmstr[k-1].ptv[2].z);
                rp=setv(tmstr[k].ptn.x-tmstr[k-1].ptc.x,tmstr[k].ptn.y-tmstr[k-1].ptc.y,tmstr[k].ptn.z-tmstr[k-1].ptc.z);
                rnorm=norm(rp);
                if(rnorm<lencn-delcn || rnorm>lencn+delcn)
                {
                    //printf("%d %d %d %d c-1 n %f lencn %f\n",k,trandpos,trand2,numiter,rnorm,lencn);
                    flagdone=false;
                    break;
                }
                //c-1 ca
            //    rp=setv(tmstr[k].ptv[1].x-tmstr[k-1].ptv[2].x,tmstr[k].ptv[1].y-tmstr[k-1].ptv[2].y,tmstr[k].ptv[1].z-tmstr[k-1].ptv[2].z);
                rp=setv(tmstr[k].x-tmstr[k-1].ptc.x,tmstr[k].y-tmstr[k-1].ptc.y,tmstr[k].z-tmstr[k-1].ptc.z);
                rnorm=norm(rp);
                if(rnorm<lencca1-delcca1 || rnorm>lencca1+delcca1)
                {
                    //printf("%d %d %d %d c-1 ca %f lencca1 %f\n",k,trandpos,trand2,numiter,rnorm,lencca1);
                    flagdone=false;
                    break;
                }
                if(k==trandpos+trand2)
                {
                    break;
                }
                //n ca
            //    rp=setv(tmstr[k].ptv[1].x-tmstr[k].ptv[0].x,tmstr[k].ptv[1].y-tmstr[k].ptv[0].y,tmstr[k].ptv[1].z-tmstr[k].ptv[0].z);
                rp=setv(tmstr[k].x-tmstr[k].ptn.x,tmstr[k].y-tmstr[k].ptn.y,tmstr[k].z-tmstr[k].ptn.z);
                rnorm=norm(rp);
                if(rnorm<lennca-delnca || rnorm>lennca+delnca)
                {
                    //printf("%d %d %d %d n ca %f lennca %f\n",k,trandpos,trand2,numiter,rnorm,lennca);
                    flagdone=false;
                    break;
                }
                //n c
            //    rp=setv(tmstr[k].ptv[2].x-tmstr[k].ptv[0].x,tmstr[k].ptv[2].y-tmstr[k].ptv[0].y,tmstr[k].ptv[2].z-tmstr[k].ptv[0].z);
                rp=setv(tmstr[k].ptc.x-tmstr[k].ptn.x,tmstr[k].ptc.y-tmstr[k].ptn.y,tmstr[k].ptc.z-tmstr[k].ptn.z);
                rnorm=norm(rp);
                if(rnorm<lennc-delnc || rnorm>lennc+delnc)
                {
                    //printf("%d %d %d %d n c %f lennc %f\n",k,trandpos,trand2,numiter,rnorm,lennc);
                    flagdone=false;
                    break;
                }
                //ca c
            //    rp=setv(tmstr[k].ptv[2].x-tmstr[k].ptv[1].x,tmstr[k].ptv[2].y-tmstr[k].ptv[1].y,tmstr[k].ptv[2].z-tmstr[k].ptv[1].z);
                rp=setv(tmstr[k].ptc.x-tmstr[k].x,tmstr[k].ptc.y-tmstr[k].y,tmstr[k].ptc.z-tmstr[k].z);
                rnorm=norm(rp);
                if(rnorm<lencac-delcac || rnorm>lencac+delcac)
                {
                    //printf("%d %d %d %d ca c %f lencac %f\n",k,trandpos,trand2,numiter,rnorm,lencac);
                    flagdone=false;
                    break;
                }
            }
            numiter++;
        //  printf("%d %d %d\n",trandpos,trand2,numiter);
        }while(!flagdone && numiter<threshiter);
/*        if(!flagdone)
        {
            summcsep[0][2]++;
            return false;
            
        } 
        else
        {
            if(trandpos+trand2+2>=seqlen)
            {
                ps.str2torp(tmstr,seqlen,trandpos,seqlen-1);
//              for(i=trandpos-1;i<seqlen;i++)
//              {
//              //  addhvinit(&tmstr[i]);
//                  addheavyinitn5one(tmstr,seqlen,i);
//              }
            }
            else
            {
                ps.str2torp(tmstr,seqlen,trandpos,trandpos+trand2+2);
//              for(i=trandpos-1;i<trandpos+trand2+2;i++)
//              {
//                  addheavyinitn5one(tmstr,seqlen,i);
//              //  addhvinit(&tmstr[i]);
//              }
            }
        }
        for(i=trandpos-1;i<trandpos+trand2+3;i++) if(i>0 && i<seqlen)
        {
            tmstr[i].hd[0]=true;
            tmstr[i].hd[2]=true;
            tmstr[i].hd[3]=true;
        }
        for(i=trandpos-4;i<trandpos+trand2+3;i++) if(i>0 && i<seqlen)
        {
            tmstr[i].hd[4]=true;
        }
        for(i=0;i<seqlen;i++)  
        {
            for(int j=0;j<seqlen;j++) if((j>=trandpos-4 && j<trandpos+trand2+3)|| (i>=trandpos-4 && i<trandpos+trand2+3))
            {
                rppmat[i*seqlen+j].flag[0]=true;
                rppmat[i*seqlen+j].flag[1]=true;
            }
        } */
    /*    memcpy(bkrppmat,rppmat,seqlen*seqlen*sizeof(residuepair));
        memcpy(bkrppene,rppene,40*sizeof(int));
        memcpy(bkrppnum,rppnum,40*sizeof(int));
        memcpy(bkrppval,rppval,40*sizeof(double)); */
    //////////make decision
 /*   *newenergy=calcrmsdenergy(tmstr,seqlen,true);
    if(*newenergy<oldenergy)
    {
        summcsep[0][0]++;
        memcpy(fulstr,tmstr,seqlen*sizeof(residueatoms));
        return true;
    }
    else
    {
        double trand;
        trand=Random();
        if(trand<exp(-tbeta*(*newenergy-oldenergy)))
        {
            summcsep[0][1]++;
            memcpy(fulstr,tmstr,seqlen*sizeof(residueatoms));
            return true;
        }
        else
        {
            summcsep[0][2]++;
            memcpy(rppmat,bkrppmat,seqlen*seqlen*sizeof(residuepair));
            memcpy(rppene,bkrppene,40*sizeof(int));
            memcpy(rppnum,bkrppnum,40*sizeof(int));
            memcpy(rppval,bkrppval,40*sizeof(double));
            return false;
        }
    } */
    return true;
}

point3d ranv(double fac)
{
    point3d tp;
    double dtheta,dphi,dr;
//  dr=0.00001+fac*rand()/double(RAND_MAX+1.0);
//  dtheta=(1.0*rand()/double(RAND_MAX+1.0))*PI;
//  dphi=(2.0*rand()/double(RAND_MAX+1.0)-1.0)*PI;
    dr=0.00001+fac*randf0and1();
    dtheta=(1.0*randf0and1())*PI;
    dphi=(2.0*randf0and1()-1.0)*PI;
    tp.x=dr*sin(dtheta)*cos(dphi);
    tp.y=dr*sin(dtheta)*sin(dphi);
    tp.z=dr*cos(dtheta);
    return(tp);
}

bool mcfragsweepaaa(vector<point3f> &tmstr,int numseq,int numshelix)
{
    int i;
    int trandpos;
    double trand2,tlength;
    int leng1,leng2,leng3,indbin;
    double rmat[9];
    point3d tp[12];
    double tnewenergy;

    vector<sssegment> sse;
    int numsse=0;
    vector<int> alphasind;
    vector<int> alphaind;
    vector<int> betaind;
//    genesse(tmstr,numseq,sse,numsse);
    genessex(tmstr,numseq,sse,numsse,numshelix);
    calcabind(alphasind,alphaind,betaind,numsse,sse);   
//  point3f *tmstr=new point3f[numseq];
//  BasicFunc bf;
//  ParseSeq ps;
    bool flagclash=false;
    bool flagtor;
    double threshclash=0.01;
//  int numclash;
/*    int tothelix=0;
    for(i=0;i<numsse;i++)
    {
        if(sse[i].ss=='H') tothelix++;
    }
//    cout<<"CC: "<<tothelix<<endl;
    if(tothelix==0)
    {
//        nummov[13][1]++;
        return false;
    }   */
    int numsalpha = alphasind.size();
    if(numsalpha==0)
    {
//      summcaaa[2]++;
//      delete[]tmstr;
        return false;
    }
do
{
//  for(i=0;i<numseq;i++)
//  {
//      tmstr[i]=decstr[i];
//  }
//  memcpy(tmstr,decstr,numseq*sizeof(point3f));
//  accdispos[0]=dispos[0];
//  for(i=1;i<numseq;i++)
//  {
//      accdispos[i]=accdispos[i-1]+dispos[i];
//  }   
        //only in coil
//      trandpos=findindpos(accdispos,accdispos[numseq-1]*rand()/double(RAND_MAX+1.0),numseq);
//      if(tmstr[trandpos].stype!='C' || trandpos==numseq-1)
//      {
//          continue;
//      }

//      trandpos=int((numsalpha)*(rand()/double(RAND_MAX+1.0)));
//      trandpos=int((numsalpha)*Random());
        trandpos = numshelix;
        leng1 = sse[alphasind[trandpos]].term-sse[alphasind[trandpos]].init+1;
        leng2 = sse[alphasind[trandpos]+1].term-sse[alphasind[trandpos]+1].init+1;
        leng3 = sse[alphasind[trandpos]+2].term-sse[alphasind[trandpos]+2].init+1;
        if(leng1>4 && leng3>4 && leng3<20)
        {
            indbin=0;
            trand2=randf0and1();
            for(i=0;i<30;i++)
            {
                if(trand2<=angleaa[i])
                {
                    indbin=i;
                    break;
                }
            }
            trand2 = (indbin+randf0and1())/30.0*PI;
        }
        else 
            trand2 = randf0and1()*PI;
        tp[0]=setv(tmstr[sse[alphasind[trandpos]].init].x,tmstr[sse[alphasind[trandpos]].init].y,
            tmstr[sse[alphasind[trandpos]].init].z);
        tp[1]=setv(tmstr[sse[alphasind[trandpos]].term].x,tmstr[sse[alphasind[trandpos]].term].y,
            tmstr[sse[alphasind[trandpos]].term].z);
        tp[2]=setv(tmstr[sse[alphasind[trandpos]+2].init].x,tmstr[sse[alphasind[trandpos]+2].init].y,
            tmstr[sse[alphasind[trandpos]+2].init].z);//for new begin more flexible
        tp[3]=setv(tmstr[sse[alphasind[trandpos]+2].term].x,tmstr[sse[alphasind[trandpos]+2].term].y,
            tmstr[sse[alphasind[trandpos]+2].term].z);
        tp[4]=minu(tp[3],tp[2]);
        tlength=norm(tp[4]);
        tp[5]=rana(trand2);
        tp[6]=setv(0,0,1);
        tp[7]=minu(tp[0],tp[1]);
        tp[8]=unit(tp[7]);
        v2rot(tp[8],rmat);
        tp[9]=rotv(tp[5],rmat);
        tp[10]=scalx(tp[9],tlength);
        tp[11]=addv(tp[2],tp[10]);//for new end
        mcfragsweepCCD4(tmstr,numseq,sse[alphasind[trandpos]+1].init,sse[alphasind[trandpos]+1].term,
            tp[2],sse[alphasind[trandpos]+2].init,tp[11],sse[alphasind[trandpos]+2].term);
        flagtor=tor2str(tmstr,numseq,3);
        if(!flagtor)
        {
            printf("tor2str wrong in aaa\n");
        }
    /*    numclash=calcaclash(tmstr,numseq);
        if(numclash<numseq*threshclash)
        {
          flagclash=true;
        }
        else flagclash=false; */
        flagclash=true;
}while(!flagclash);
    
//  int numclash;
//  numclash=ps.calcaclash(tmstr,numseq);
//  if(numclash>numseq*threshclash)
//  {
//      summcome[2]++;
//      delete[]tmstr;
//      return false;
//  }
    //////////make decision
//  bool flagev=evcheck(tmstr,numseq);
//  if(!flagev) numev3[10][3]++;
//  if(flagev)
//  {
//      summcaaa[2]++;
//      return false;
//  }
    return false;
}

bool mcfragsweepCCD6(vector<point3f> &decstr,int numseq,
        int lps,int lpe,point3d pt1,int ind1,point3d pt2,int ind2)
{
//  int i;
//  BasicFunc bf;
//  ParseSeq ps;
    int trandpos;
//  int indphi,indpsi;
    bool flagphi,flagtor;
    int flagpt,indp[2];
    indp[0]=ind1;indp[1]=ind2;
    point3d tp[2]; 
    double cutdist=5.0;
    tp[0].x=pt1.x;tp[0].y=pt1.y;tp[0].z=pt1.z;
    tp[1].x=pt2.x;tp[1].y=pt2.y;tp[1].z=pt2.z;
//  point3f *tmstr=new point3f[numseq];
    point3d pcur,p12,p13;
    int numiter=0;
    double tdist[2],ttheta,tphi,tpsi,tinner;
//  for(i=0;i<numseq;i++)
//  {
//      tmstr[i]=decstr[i];
//  }
    vector<point3f> tmstr2;
    vector<point3f>().swap(tmstr2);
    tmstr2 = decstr;
//    memcpy(tmstr2,decstr,numseq*sizeof(point3f));
    flagpt=0;
//  trandpos=lps;
    int threshiter=15*(lpe-lps+1);
    if(threshiter>100) threshiter=100;
    do{
        tdist[0]=10000;
        tdist[1]=10000;
//      trandpos=lps+int((lpe-lps+1)*(rand()/double(RAND_MAX+1.0)));
        trandpos=lps+int((lpe-lps+1)*randf0and1());
//      if(trandpos>lpe)
//      {
//          printf("ccd3 pos is large %d %d %d\n",lps,lpe,trandpos);
//          trandpos=lpe;
//      }
//      trandpos=lps+(trandpos+1-lps)%(lpe-lps+1);
    //  if(genrand()%2==0)
        if(randf0and1()<0.5)
        {
            flagphi=true;//change phi n ca
        }
        else
        {
            flagphi=false;//change psi ca c in pos+1
        }
        
        flagpt=(flagpt+1)%2;
        if(flagphi)
        {
            p12=setv(tmstr2[trandpos].x-tmstr2[trandpos].ptn.x,tmstr2[trandpos].y-tmstr2[trandpos].ptn.y,
                tmstr2[trandpos].z-tmstr2[trandpos].ptn.z);
            p13=setv(tmstr2[indp[flagpt]].x-tmstr2[trandpos].ptn.x,tmstr2[indp[flagpt]].y-tmstr2[trandpos].ptn.y,
                tmstr2[indp[flagpt]].z-tmstr2[trandpos].ptn.z);
            tinner=angv(p12,p13);
            if(tinner<epsilon || PI-tinner<epsilon || norm(p12)<epsilon)//in one line no affect when rotate
            {
                numiter++;
            //  printf("oneline 1\n");
                continue;
            }
            pcur.x=tmstr2[indp[flagpt]].x;pcur.y=tmstr2[indp[flagpt]].y;pcur.z=tmstr2[indp[flagpt]].z;
            p13.x=tmstr2[trandpos].x;p13.y=tmstr2[trandpos].y;p13.z=tmstr2[trandpos].z;
            p12.x=tmstr2[trandpos].ptn.x;p12.y=tmstr2[trandpos].ptn.y;p12.z=tmstr2[trandpos].ptn.z;
            ttheta=phi(pcur.x,pcur.y,pcur.z,p12.x,p12.y,p12.z,p13.x,p13.y,p13.z,tp[flagpt].x,tp[flagpt].y,tp[flagpt].z);
            tphi=tmstr2[trandpos].tor[2]+ttheta;
            if(tphi>360.0)
            {
                tphi-=360.0;
            }
//          tpsi=tmstr2[trandpos+1].tor[0];
//          indphi=int(tphi/18.0);
//          indpsi=int(tpsi/18.0);
//          if(!tphipsi[indpsi][indphi])
//          {
//          //  continue;
//          }
            tmstr2[trandpos].tor[2]=tphi;   
        }
        else
        {
            p12=setv(tmstr2[trandpos].ptc.x-tmstr2[trandpos].x,tmstr2[trandpos].ptc.y-tmstr2[trandpos].y,
                tmstr2[trandpos].ptc.z-tmstr2[trandpos].z);
            p13=setv(tmstr2[indp[flagpt]].x-tmstr2[trandpos].x,tmstr2[indp[flagpt]].y-tmstr2[trandpos].y,
                tmstr2[indp[flagpt]].z-tmstr2[trandpos].z);
            tinner=angv(p12,p13);
            if(tinner<epsilon || PI-tinner<epsilon || norm(p12)<epsilon)//in one line no affect when rotate
            {
                numiter++;
            //  printf("oneline 2\n");
                continue;
            }
            pcur.x=tmstr2[indp[flagpt]].x;pcur.y=tmstr2[indp[flagpt]].y;pcur.z=tmstr2[indp[flagpt]].z;
            p12.x=tmstr2[trandpos].x;p12.y=tmstr2[trandpos].y;p12.z=tmstr2[trandpos].z;
            p13.x=tmstr2[trandpos].ptc.x;p13.y=tmstr2[trandpos].ptc.y;p13.z=tmstr2[trandpos].ptc.z;
            ttheta=phi(pcur.x,pcur.y,pcur.z,p12.x,p12.y,p12.z,p13.x,p13.y,p13.z,tp[flagpt].x,tp[flagpt].y,tp[flagpt].z);
            tpsi=tmstr2[trandpos+1].tor[0]+ttheta;
            if(tpsi>360.0)
            {
                tpsi-=360.0;
            }
//          tphi=tmstr2[trandpos].tor[2];
//          indphi=int(tphi/18.0);
//          indpsi=int(tpsi/18.0);
//          if(!tphipsi[indpsi][indphi])
//          {
//          //  continue;//go to the endpos to check tdist and numiter
//          }
            tmstr2[trandpos+1].tor[0]=tpsi;
        }       
        numiter++;
        if(trandpos==0)
            flagtor=tor2str(tmstr2,numseq,3);
        else
            flagtor=tor2strp(tmstr2,numseq,trandpos);
        if(!flagtor)
        {
            printf("tor2str wrong in CCD6 %d\n",trandpos);
        }
    pcur.x=tmstr2[ind1].x-tp[0].x;
    pcur.y=tmstr2[ind1].y-tp[0].y;
    pcur.z=tmstr2[ind1].z-tp[0].z;
    tdist[0]=pcur.x*pcur.x+pcur.y*pcur.y+pcur.z*pcur.z;
    pcur.x=tmstr2[ind2].x-tp[1].x;
    pcur.y=tmstr2[ind2].y-tp[1].y;
    pcur.z=tmstr2[ind2].z-tp[1].z;
    tdist[1]=pcur.x*pcur.x+pcur.y*pcur.y+pcur.z*pcur.z;
//  printf("%d %d %f %f %f\n",numiter,threshiter,tdist[0],tdist[1],tdist[2]);
    }while(numiter<threshiter && (tdist[0]>cutdist || tdist[1]>cutdist));
//  printf("%d [%d %d %d] pos %d %d dist %f theta %f\n",numiter,lps,lpe,lpt,trandpos, flagphi, tdist,ttheta);
//  *newenergy=calcrmsdenergy(tmstr2,numseq);
    if(numiter<threshiter)
    {
//      for(i=0;i<numseq;i++)
//      {
//          decstr[i]=tmstr[i];
//      }
        vector<point3f>().swap(decstr);
        decstr=tmstr2;
//        memcpy(decstr,tmstr2,numseq*sizeof(point3f));
//      delete[]tmstr;
        return true;
    }
    else
    {
//      for(i=0;i<numseq;i++)
//      {
//          decstr[i]=tmstr[i];
//      }
//        vector<point3f>().swap(decstr);
//        decstr=tmstr2;        
//        memcpy(decstr,tmstr2,numseq*sizeof(point3f));
//      delete[]tmstr;
        return false;
    }
}

bool mcfragsweepCCD4(vector<point3f> &tmstr2,int numseq,int lps,int lpe,point3d pt1,int ind1,point3d pt2,int ind2)
{
//  int i;
//  BasicFunc bf;
//  ParseSeq ps;
    int trandpos;
//  int indphi,indpsi;
    bool flagphi,flagtor;
    int flagpt,indp[2];
    indp[0]=ind1;indp[1]=ind2;
    point3d tp[2]; 
    tp[0].x=pt1.x;tp[0].y=pt1.y;tp[0].z=pt1.z;
    tp[1].x=pt2.x;tp[1].y=pt2.y;tp[1].z=pt2.z;
//  point3f *tmstr=new point3f[numseq];
    point3d pcur,p12,p13;
    int numiter=0;
    double tdist[2],ttheta,tphi,tpsi,tinner;
//  for(i=0;i<numseq;i++)
//  {
//      tmstr[i]=decstr[i];
//  }
//  memcpy(tmstr2,decstr,numseq*sizeof(point3f));
    flagpt=0;
//  trandpos=lps;
    int threshiter=25*(lpe-lps+1);
    if(threshiter>600) threshiter=600;
    do{
        tdist[0]=10000;
        tdist[1]=10000;
//      trandpos=lps+int((lpe-lps+1)*(rand()/double(RAND_MAX+1.0)));
        trandpos=lps+int((lpe-lps+1)*randf0and1());
//      if(trandpos>lpe)
//      {
//          printf("ccd3 pos is large %d %d %d\n",lps,lpe,trandpos);
//          trandpos=lpe;
//      }
//      trandpos=lps+(trandpos+1-lps)%(lpe-lps+1);
    //  if(genrand()%2==0)
        if(randf0and1()<0.5)
        {
            flagphi=true;//change phi n ca
        }
        else
        {
            flagphi=false;//change psi ca c in pos+1
        }
        if(tmstr2[trandpos].ss2=='H' && tmstr2[trandpos+1].ss2=='H')
        {
            numiter++;
            continue;
        }
        else if(tmstr2[trandpos].ss2=='H')
        {
            flagphi=false;
        }
        else if(tmstr2[trandpos].ss2=='E' && tmstr2[trandpos+1].ss2=='E')
        {
            numiter++;
            continue;
        }
        else if(tmstr2[trandpos].ss2=='E')
        {
            flagphi=false;
        }
        flagpt=(flagpt+1)%2;
        if(flagphi)
        {
            p12=setv(tmstr2[trandpos].x-tmstr2[trandpos].ptn.x,tmstr2[trandpos].y-tmstr2[trandpos].ptn.y,
                tmstr2[trandpos].z-tmstr2[trandpos].ptn.z);
            p13=setv(tmstr2[indp[flagpt]].x-tmstr2[trandpos].ptn.x,tmstr2[indp[flagpt]].y-tmstr2[trandpos].ptn.y,
                tmstr2[indp[flagpt]].z-tmstr2[trandpos].ptn.z);
            tinner=angv(p12,p13);
            if(tinner<epsilon || PI-tinner<epsilon || norm(p12)<epsilon)//in one line no affect when rotate
            {
                numiter++;
            //  printf("oneline 1\n");
                continue;
            }
            pcur.x=tmstr2[indp[flagpt]].x;pcur.y=tmstr2[indp[flagpt]].y;pcur.z=tmstr2[indp[flagpt]].z;
            p13.x=tmstr2[trandpos].x;p13.y=tmstr2[trandpos].y;p13.z=tmstr2[trandpos].z;
            p12.x=tmstr2[trandpos].ptn.x;p12.y=tmstr2[trandpos].ptn.y;p12.z=tmstr2[trandpos].ptn.z;
            ttheta=phi(pcur.x,pcur.y,pcur.z,p12.x,p12.y,p12.z,p13.x,p13.y,p13.z,tp[flagpt].x,tp[flagpt].y,tp[flagpt].z);
            tphi=tmstr2[trandpos].tor[2]+ttheta;
            if(tphi>360.0)
            {
                tphi-=360.0;
            }
//          tpsi=tmstr2[trandpos+1].tor[0];
//          indphi=int(tphi/18.0);
//          indpsi=int(tpsi/18.0);
//          if(!tphipsi[indpsi][indphi])
//          {
//          //  continue;
//          }
            tmstr2[trandpos].tor[2]=tphi;   
        }
        else
        {
            p12=setv(tmstr2[trandpos].ptc.x-tmstr2[trandpos].x,tmstr2[trandpos].ptc.y-tmstr2[trandpos].y,
                tmstr2[trandpos].ptc.z-tmstr2[trandpos].z);
            p13=setv(tmstr2[indp[flagpt]].x-tmstr2[trandpos].x,tmstr2[indp[flagpt]].y-tmstr2[trandpos].y,
                tmstr2[indp[flagpt]].z-tmstr2[trandpos].z);
            tinner=angv(p12,p13);
            if(tinner<epsilon || PI-tinner<epsilon || norm(p12)<epsilon)//in one line no affect when rotate
            {
                numiter++;
            //  printf("oneline 2\n");
                continue;
            }
            pcur.x=tmstr2[indp[flagpt]].x;pcur.y=tmstr2[indp[flagpt]].y;pcur.z=tmstr2[indp[flagpt]].z;
            p12.x=tmstr2[trandpos].x;p12.y=tmstr2[trandpos].y;p12.z=tmstr2[trandpos].z;
            p13.x=tmstr2[trandpos].ptc.x;p13.y=tmstr2[trandpos].ptc.y;p13.z=tmstr2[trandpos].ptc.z;
            ttheta=phi(pcur.x,pcur.y,pcur.z,p12.x,p12.y,p12.z,p13.x,p13.y,p13.z,tp[flagpt].x,tp[flagpt].y,tp[flagpt].z);
            tpsi=tmstr2[trandpos+1].tor[0]+ttheta;
            if(tpsi>360.0)
            {
                tpsi-=360.0;
            }
//          tphi=tmstr2[trandpos].tor[2];
//          indphi=int(tphi/18.0);
//          indpsi=int(tpsi/18.0);
//          if(!tphipsi[indpsi][indphi])
//          {
//          //  continue;//go to the endpos to check tdist and numiter
//          }
            tmstr2[trandpos+1].tor[0]=tpsi;
        }       
        numiter++;
        if(trandpos==0)
            flagtor=tor2str(tmstr2,numseq,3);
        else
            flagtor=tor2strp(tmstr2,numseq,trandpos);
        if(!flagtor)
        {
            printf("tor2str wrong in CCD4 %d\n",trandpos);
        }
    pcur.x=tmstr2[ind1].x-tp[0].x;
    pcur.y=tmstr2[ind1].y-tp[0].y;
    pcur.z=tmstr2[ind1].z-tp[0].z;
    tdist[0]=pcur.x*pcur.x+pcur.y*pcur.y+pcur.z*pcur.z;
    pcur.x=tmstr2[ind2].x-tp[1].x;
    pcur.y=tmstr2[ind2].y-tp[1].y;
    pcur.z=tmstr2[ind2].z-tp[1].z;
    tdist[1]=pcur.x*pcur.x+pcur.y*pcur.y+pcur.z*pcur.z;
//  printf("%d %d %f %f %f\n",numiter,threshiter,tdist[0],tdist[1],tdist[2]);
    }while(numiter<threshiter && (tdist[0]>6.50 || tdist[1]>2.50));
//  printf("%d [%d %d %d] pos %d %d dist %f theta %f\n",numiter,lps,lpe,lpt,trandpos, flagphi, tdist,ttheta);
//  *newenergy=calcrmsdenergy(tmstr2,numseq);
/*  if(numiter<threshiter)
    {

//
//      for(i=0;i<numseq;i++)
//      {
//          decstr[i]=tmstr[i];
//      }
        memcpy(decstr,tmstr2,numseq*sizeof(point3f));
//      delete[]tmstr;
        return true;
    }
    else
    {
//      for(i=0;i<numseq;i++)
//      {
//          decstr[i]=tmstr[i];
//      }
        memcpy(decstr,tmstr2,numseq*sizeof(point3f));
//      delete[]tmstr;
        return false;
    }  */
    return true;
}

bool tor2str(vector<point3f> &decstr,int seqnum,int type)
{
    int i;
    point3s pt,pn,pc;
//    BasicFunc bf;
    bool flagts;
    bool flagwhole=true;
if(type==1)//all 0 n 1 ca 2 c
{
    if(decstr[1].leng<0)
    {
        decstr[1].leng=float(lennca);
    }
    if(decstr[2].leng<0)
    {
        decstr[2].leng=float(lencac);
    }
    if(decstr[2].angl<0)
    {
        decstr[2].angl=float(angncac);
    }
    decstr[0].x=0;decstr[0].y=0;decstr[0].z=0;
    decstr[1].x=decstr[1].leng;decstr[1].y=0;decstr[1].z=0;
    decstr[2].x=decstr[1].leng-decstr[2].leng*cos(decstr[2].angl*raddeg);
    decstr[2].y=decstr[2].leng*sin(decstr[2].angl*raddeg);decstr[2].z=0;
    for(i=3;i<seqnum;i++)
    {
        flagts=tor2pos22(decstr[i-3].x,decstr[i-3].y,decstr[i-3].z,decstr[i-2].x,decstr[i-2].y,decstr[i-2].z,
            decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,decstr[i].phi*raddeg,decstr[i].leng,decstr[i].angl*raddeg,&pt.x,&pt.y,&pt.z);
        if(!flagts)
        {
            flagwhole=false;
            printf("wrong front coordinates %d\n",i);
        }
        decstr[i].x=pt.x;decstr[i].y=pt.y;decstr[i].z=pt.z;
    }
}
else if(type==6)//ca 
{
    if(decstr[1].leng<0)
    {
        decstr[1].leng=float(lencaca);
    }
    if(decstr[2].leng<0)
    {
        decstr[2].leng=float(lencaca);
    }
    if(decstr[2].angl<0)
    {
        decstr[2].angl=106.422f;
    }
    decstr[0].x=0;decstr[0].y=0;decstr[0].z=0;
    decstr[1].x=decstr[1].leng;decstr[1].y=0;decstr[1].z=0;
    decstr[2].x=decstr[1].leng-decstr[2].leng*cos(decstr[2].angl*raddeg);
    decstr[2].y=decstr[2].leng*sin(decstr[2].angl*raddeg);decstr[2].z=0;
    for(i=3;i<seqnum;i++)
    {
        flagts=tor2pos22(decstr[i-3].x,decstr[i-3].y,decstr[i-3].z,decstr[i-2].x,decstr[i-2].y,decstr[i-2].z,
            decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,decstr[i].phi*raddeg,decstr[i].leng,decstr[i].angl*raddeg,&pt.x,&pt.y,&pt.z);
        if(!flagts)
        {
            flagwhole=false;
            printf("wrong front coordinates %d\n",i);
        }
        decstr[i].x=pt.x;decstr[i].y=pt.y;decstr[i].z=pt.z;
    }
}
else if(type==3)
{
    if(decstr[0].len[1]<0)
    {
        decstr[0].len[1]=float(lennca);
    }
    if(decstr[0].len[2]<0)
    {
        decstr[0].len[2]=float(lencac);
    }
    if(decstr[0].ang[2]<0)
    {
        decstr[0].ang[2]=float(angncac);
    }
    decstr[0].ptn.x=0;decstr[0].ptn.y=0;decstr[0].ptn.z=0;
    decstr[0].x=decstr[0].len[1];decstr[0].y=0;decstr[0].z=0;
    decstr[0].ptc.x=decstr[0].len[1]-decstr[0].len[2]*cos(decstr[0].ang[2]*raddeg);
    decstr[0].ptc.y=decstr[0].len[2]*sin(decstr[0].ang[2]*raddeg);decstr[0].ptc.z=0;
    if(decstr[0].tor[0]>=0 && decstr[0].tor[2]>=0)
    {
        flagts=tor2pos22(0,0,0,lennca,0,0,lennca+lencac*sin(angncac*raddeg),-lencac*cos(angncac*raddeg),0,
            decstr[0].tor[0]*raddeg,float(lencn),float(angcacn)*raddeg,&pn.x,&pn.y,&pn.z);
        if(!flagts)
        {
            flagwhole=false;
            printf("wrong front coordinates n %d\n",0);
        }
        decstr[0].ptn.x=pn.x;decstr[0].ptn.y=pn.y;decstr[0].ptn.z=pn.z;
        if(decstr[0].tor[1]<0) 
            decstr[0].tor[1]=180;
        flagts=tor2pos22(lennca,0,0,lennca+lencac*sin(angncac*raddeg),-lencac*cos(angncac*raddeg),0,
            decstr[0].ptn.x,decstr[0].ptn.y,decstr[0].ptn.z,
            decstr[0].tor[1]*raddeg,float(lennca),angcnca*raddeg,&pt.x,&pt.y,&pt.z);
        if(!flagts)
        {
            flagwhole=false;
            printf("wrong front coordinates ca %d\n",0);
        }
        decstr[0].x=pt.x;decstr[0].y=pt.y;decstr[0].z=pt.z;
        flagts=tor2pos22(lennca+lencac*sin(angncac*raddeg),-lencac*cos(angncac*raddeg),0,
            decstr[0].ptn.x,decstr[0].ptn.y,decstr[0].ptn.z,decstr[0].x,decstr[0].y,decstr[0].z,
            decstr[0].tor[2]*raddeg,float(lencac),angncac*raddeg,&pc.x,&pc.y,&pc.z);
        if(!flagts)
        {
            flagwhole=false;
            printf("wrong front coordinates c %d\n",0);
        }
        decstr[0].ptc.x=pc.x;decstr[0].ptc.y=pc.y;decstr[0].ptc.z=pc.z;
    }
    for(i=1;i<seqnum;i++)
    {
        if(decstr[i].tor[0]<0) 
        {
            if(decstr[i].stype=='E')
                decstr[i].tor[0]=120;
            else if(decstr[i].stype=='H')
                decstr[i].tor[0]=300;
            else
                decstr[i].tor[0]=120;
        }
        if(decstr[i].tor[1]<0) 
            decstr[i].tor[1]=180;
        if(decstr[i].tor[2]<0) 
            decstr[i].tor[2]=290;
        if(decstr[i].len[0]<0) 
            decstr[i].len[0]=float(lencn);
        if(decstr[i].len[1]<0) 
            decstr[i].len[1]=float(lennca);
        if(decstr[i].len[2]<0) 
            decstr[i].len[2]=float(lencac);
        if(decstr[i].ang[0]<0) 
            decstr[i].ang[0]=float(angcacn);
        if(decstr[i].ang[1]<0) 
            decstr[i].ang[1]=float(angcnca);
        if(decstr[i].ang[2]<0) 
            decstr[i].ang[2]=float(angncac);
        //0 1 2 n ca c
        //original phi i-1 psi i-1 omega i
        flagts=tor2pos22(decstr[i-1].ptn.x,decstr[i-1].ptn.y,decstr[i-1].ptn.z,decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
            decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
            decstr[i].tor[0]*raddeg,decstr[i].len[0],decstr[i].ang[0]*raddeg,&pn.x,&pn.y,&pn.z);
        if(!flagts)
        {
            flagwhole=false;
            printf("wrong front coordinates n %d %f %f %f\n",i,decstr[i].tor[0],decstr[i].len[0],decstr[i].ang[0]);
        }
        decstr[i].ptn.x=pn.x;decstr[i].ptn.y=pn.y;decstr[i].ptn.z=pn.z;
        flagts=tor2pos22(decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
            decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
            decstr[i].tor[1]*raddeg,decstr[i].len[1],decstr[i].ang[1]*raddeg,&pt.x,&pt.y,&pt.z);
        if(!flagts)
        {
            flagwhole=false;
            printf("wrong front coordinates ca %d %f %f %f\n",i,decstr[i].tor[1],decstr[i].len[1],decstr[i].ang[1]);
        }
        decstr[i].x=pt.x;decstr[i].y=pt.y;decstr[i].z=pt.z;
        flagts=tor2pos22(decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
            decstr[i].x,decstr[i].y,decstr[i].z,
            decstr[i].tor[2]*raddeg,decstr[i].len[2],decstr[i].ang[2]*raddeg,&pc.x,&pc.y,&pc.z);
        if(!flagts)
        {
            flagwhole=false;
            printf("wrong front coordinates c %d %f %f %f\n",i,decstr[i].tor[2],decstr[i].len[2],decstr[i].ang[2]);
        }
        decstr[i].ptc.x=pc.x;decstr[i].ptc.y=pc.y;decstr[i].ptc.z=pc.z;
        /*
        //atom o
        flagts=bf.tor2pos22(decstr[i-1].ptn.x,decstr[i-1].ptn.y,decstr[i-1].ptn.z,decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
            decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
            decstr[i].tor[0]*raddeg-PI,lenco,angcaco*raddeg,&pn.x,&pn.y,&pn.z);
        if(!flagts)
        {
            flagwhole=false;
            printf("wrong front coordinates o %d %f %f %f\n",i,decstr[i].tor[0],decstr[i].len[0],decstr[i].ang[0]);
        }
        decstr[i-1].pto.x=pn.x;decstr[i-1].pto.y=pn.y;decstr[i-1].pto.z=pn.z;
        */
        /*
        //atom h
        flagts=bf.tor2pos22(decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
            decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,0,lennh,angcnh*raddeg,&pn.x,&pn.y,&pn.z);
        if(!flagts)
        {
            flagwhole=false;
            printf("wrong front coordinates d %d\n",i);
        }
        decstr[i].pth.x=pn.x;decstr[i].pth.y=pn.y;decstr[i].pth.z=pn.z;
        */
    }
    
    /*
    //atom o
    i=seqnum;
    flagts=bf.tor2pos22(decstr[i-1].ptn.x,decstr[i-1].ptn.y,decstr[i-1].ptn.z,decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
        decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,0,lenco,angcaco*raddeg,&pn.x,&pn.y,&pn.z);
    if(!flagts)
    {
        flagwhole=false;
        printf("wrong front coordinates o %d\n",i);
    }
    decstr[i-1].pto.x=pn.x;decstr[i-1].pto.y=pn.y;decstr[i-1].pto.z=pn.z;
    */
    /*
    //atom h 
    i=0;
    flagts=bf.tor2pos22(decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,decstr[i].x,decstr[i].y,decstr[i].z,
        decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,0,lennh,anghnca*raddeg,&pn.x,&pn.y,&pn.z);
    if(!flagts)
    {
        flagwhole=false;
        printf("wrong front coordinates %d\n",i);
    }
    decstr[i].pth.x=pn.x;decstr[i].pth.y=pn.y;decstr[i].pth.z=pn.z;
    //atom cb
    for(i=0;i<seqnum;i++)
    {
        flagts=bf.tor2pos22(decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,decstr[i].x,decstr[i].y,decstr[i].z,
            decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,
            torncaccb*raddeg,lenccb,angcaccb*raddeg,&pn.x,&pn.y,&pn.z);
        if(!flagts)
        {
            flagwhole=false;
            printf("wrong front coordinates e %d\n",i);
        }
        decstr[i].ptb.x=pn.x;decstr[i].ptb.y=pn.y;decstr[i].ptb.z=pn.z;
        if(decstr[i].aaa=='G')
        {
            decstr[i].ptb.x=decstr[i].x;decstr[i].ptb.y=decstr[i].y;decstr[i].ptb.z=decstr[i].z;
        }
    }
    //atom ha
    for(i=0;i<seqnum;i++)
    {
        flagts=bf.tor2pos22(decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,
            decstr[i].x,decstr[i].y,decstr[i].z,tornccaha*raddeg,lencaha,angccaha*raddeg,&pn.x,&pn.y,&pn.z);
        if(!flagts)
        {
            flagwhole=false;
            printf("wrong front coordinates f %d\n",i);
        }
        decstr[i].ptha.x=pn.x;decstr[i].ptha.y=pn.y;decstr[i].ptha.z=pn.z;
    }   
    */
}
    return flagwhole;
}

bool itor2strp(vector<point3f> &decstr,int seqnum,int iss,int istart)
{
    if(istart>=seqnum-1) return itor2str(decstr,seqnum);
    int i;
    point3s pt,pn,pc;
//    BasicFunc bf;
    bool flagts;
    bool flagwhole=true;
    for(i=istart;i>=iss;i--)
    {
        if(decstr[i+1].tor[0]<0) 
        {
            if(decstr[i+1].stype=='E')
                decstr[i+1].tor[0]=120;
            else if(decstr[i+1].stype=='H')
                decstr[i+1].tor[0]=300;
            else
                decstr[i+1].tor[0]=120;
        }
        if(decstr[i+1].tor[1]<0) 
            decstr[i+1].tor[1]=180;
        if(decstr[i+1].tor[2]<0) 
            decstr[i+1].tor[2]=290;
        if(decstr[i+1].len[0]<0) 
            decstr[i+1].len[0]=float(lencn);
        if(decstr[i].len[1]<0) 
            decstr[i].len[1]=float(lennca);
        if(decstr[i].len[2]<0) 
            decstr[i].len[2]=float(lencac);
        if(decstr[i+1].ang[0]<0) 
            decstr[i+1].ang[0]=float(angcacn);
        if(decstr[i+1].ang[1]<0) 
            decstr[i+1].ang[1]=float(angcnca);
        if(decstr[i].ang[2]<0) 
            decstr[i].ang[2]=float(angncac);
        //0 1 2 n ca c
        //original phi i-1 psi i-1 omega i
        flagts=tor2pos22(decstr[i+1].ptc.x,decstr[i+1].ptc.y,decstr[i+1].ptc.z,decstr[i+1].x,decstr[i+1].y,decstr[i+1].z,
            decstr[i+1].ptn.x,decstr[i+1].ptn.y,decstr[i+1].ptn.z,
            decstr[i+1].tor[2]*raddeg,decstr[i+1].len[0],decstr[i+1].ang[1]*raddeg,&pn.x,&pn.y,&pn.z);
        if(!flagts)
        {
            flagwhole=false;
            printf("wrong front coordinates c %d %f %f %f\n",i,decstr[i+1].tor[2],decstr[i+1].len[0],decstr[i+1].ang[1]);
        }
        decstr[i].ptc.x=pn.x;decstr[i].ptc.y=pn.y;decstr[i].ptc.z=pn.z;
        flagts=tor2pos22(decstr[i+1].x,decstr[i+1].y,decstr[i+1].z,decstr[i+1].ptn.x,decstr[i+1].ptn.y,decstr[i+1].ptn.z,
            decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,
            decstr[i+1].tor[1]*raddeg,decstr[i].len[2],decstr[i+1].ang[0]*raddeg,&pt.x,&pt.y,&pt.z);
        if(!flagts)
        {
            flagwhole=false;
            printf("wrong front coordinates ca %d %f %f %f\n",i,decstr[i+1].tor[1],decstr[i].len[2],decstr[i+1].ang[0]);
        }
        decstr[i].x=pt.x;decstr[i].y=pt.y;decstr[i].z=pt.z;
        flagts=tor2pos22(decstr[i+1].ptn.x,decstr[i+1].ptn.y,decstr[i+1].ptn.z,decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,
            decstr[i].x,decstr[i].y,decstr[i].z,
            decstr[i+1].tor[0]*raddeg,decstr[i].len[1],decstr[i].ang[2]*raddeg,&pc.x,&pc.y,&pc.z);
        if(!flagts)
        {
            flagwhole=false;
            printf("wrong front coordinates n %d %f %f %f\n",i,decstr[i+1].tor[0],decstr[i].len[1],decstr[i].ang[2]);
        }
        decstr[i].ptn.x=pc.x;decstr[i].ptn.y=pc.y;decstr[i].ptn.z=pc.z;
    }
    return true;
}

bool itor2str(vector<point3f> &decstr,int seqnum)
{
    int i;
    point3s pt,pn,pc;
//    BasicFunc bf;
    bool flagts;
    bool flagwhole=true;
    if(decstr[seqnum-1].len[2]<0)
    {
        decstr[seqnum-1].len[2]=float(lencac);
    }
    if(decstr[seqnum-1].len[1]<0)
    {
        decstr[seqnum-1].len[1]=float(lennca);
    }
    if(decstr[seqnum-1].ang[2]<0)
    {
        decstr[seqnum-1].ang[2]=float(angncac);
    }
    decstr[seqnum-1].ptc.x=0;decstr[seqnum-1].ptc.y=0;decstr[seqnum-1].ptc.z=0;
    decstr[seqnum-1].x=decstr[seqnum-1].len[2];decstr[seqnum-1].y=0;decstr[seqnum-1].z=0;
    decstr[seqnum-1].ptn.x=decstr[seqnum-1].len[2]-decstr[seqnum-1].len[1]*cos(decstr[seqnum-1].ang[2]*raddeg);
    decstr[seqnum-1].ptn.y=decstr[seqnum-1].len[1]*sin(decstr[seqnum-1].ang[2]*raddeg);decstr[seqnum-1].ptn.z=0;
    for(i=seqnum-2;i>=0;i--)
    {
        if(decstr[i+1].tor[0]<0) 
        {
            if(decstr[i+1].stype=='E')
                decstr[i+1].tor[0]=120;
            else if(decstr[i+1].stype=='H')
                decstr[i+1].tor[0]=300;
            else
                decstr[i+1].tor[0]=120;
        }
        if(decstr[i+1].tor[1]<0) 
            decstr[i+1].tor[1]=180;
        if(decstr[i+1].tor[2]<0) 
            decstr[i+1].tor[2]=290;
        if(decstr[i+1].len[0]<0) 
            decstr[i+1].len[0]=float(lencn);
        if(decstr[i].len[1]<0) 
            decstr[i].len[1]=float(lennca);
        if(decstr[i].len[2]<0) 
            decstr[i].len[2]=float(lencac);
        if(decstr[i+1].ang[0]<0) 
            decstr[i+1].ang[0]=float(angcacn);
        if(decstr[i+1].ang[1]<0) 
            decstr[i+1].ang[1]=float(angcnca);
        if(decstr[i].ang[2]<0) 
            decstr[i].ang[2]=float(angncac);
        //0 1 2 n ca c
        //original phi i-1 psi i-1 omega i
        flagts=tor2pos22(decstr[i+1].ptc.x,decstr[i+1].ptc.y,decstr[i+1].ptc.z,decstr[i+1].x,decstr[i+1].y,decstr[i+1].z,
            decstr[i+1].ptn.x,decstr[i+1].ptn.y,decstr[i+1].ptn.z,
            decstr[i+1].tor[2]*raddeg,decstr[i+1].len[0],decstr[i+1].ang[1]*raddeg,&pn.x,&pn.y,&pn.z);
        if(!flagts)
        {
            flagwhole=false;
            printf("wrong front coordinates c %d %f %f %f\n",i,decstr[i+1].tor[2],decstr[i+1].len[0],decstr[i+1].ang[1]);
        }
        decstr[i].ptc.x=pn.x;decstr[i].ptc.y=pn.y;decstr[i].ptc.z=pn.z;
        flagts=tor2pos22(decstr[i+1].x,decstr[i+1].y,decstr[i+1].z,decstr[i+1].ptn.x,decstr[i+1].ptn.y,decstr[i+1].ptn.z,
            decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,
            decstr[i+1].tor[1]*raddeg,decstr[i].len[2],decstr[i+1].ang[0]*raddeg,&pt.x,&pt.y,&pt.z);
        if(!flagts)
        {
            flagwhole=false;
            printf("wrong front coordinates ca %d %f %f %f\n",i,decstr[i+1].tor[1],decstr[i].len[2],decstr[i+1].ang[0]);
        }
        decstr[i].x=pt.x;decstr[i].y=pt.y;decstr[i].z=pt.z;
        flagts=tor2pos22(decstr[i+1].ptn.x,decstr[i+1].ptn.y,decstr[i+1].ptn.z,decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,
            decstr[i].x,decstr[i].y,decstr[i].z,
            decstr[i+1].tor[0]*raddeg,decstr[i].len[1],decstr[i].ang[2]*raddeg,&pc.x,&pc.y,&pc.z);
        if(!flagts)
        {
            flagwhole=false;
            printf("wrong front coordinates n %d %f %f %f\n",i,decstr[i+1].tor[0],decstr[i].len[1],decstr[i].ang[2]);
        }
        decstr[i].ptn.x=pc.x;decstr[i].ptn.y=pc.y;decstr[i].ptn.z=pc.z;
    }
    return flagwhole;
}

bool tor2strp(vector<point3f> &decstr,int seqnum,int istart)
{
    int i;
    point3s pt,pn,pc;
//  BasicFunc bf;
    bool flagts;
    bool flagwhole=true;
    for(i=istart;i<seqnum;i++)
    {
        if(decstr[i].tor[0]<0) 
        {
            if(decstr[i].stype=='E')
                decstr[i].tor[0]=120;
            else if(decstr[i].stype=='H')
                decstr[i].tor[0]=300;
            else
                decstr[i].tor[0]=120;
        }
        if(decstr[i].tor[1]<0) 
            decstr[i].tor[1]=180;
        if(decstr[i].tor[2]<0) 
            decstr[i].tor[2]=290;
        if(decstr[i].len[0]<0) 
            decstr[i].len[0]=float(lencn);
        if(decstr[i].len[1]<0) 
            decstr[i].len[1]=float(lennca);
        if(decstr[i].len[2]<0) 
            decstr[i].len[2]=float(lencac);
        if(decstr[i].ang[0]<0) 
            decstr[i].ang[0]=float(angcacn);
        if(decstr[i].ang[1]<0) 
            decstr[i].ang[1]=float(angcnca);
        if(decstr[i].ang[2]<0) 
            decstr[i].ang[2]=float(angncac);
        //0 1 2 n ca c
        //original phi i-1 psi i-1 omega i
        flagts=tor2pos22(decstr[i-1].ptn.x,decstr[i-1].ptn.y,decstr[i-1].ptn.z,decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
            decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
            decstr[i].tor[0]*raddeg,decstr[i].len[0],decstr[i].ang[0]*raddeg,&pn.x,&pn.y,&pn.z);
        if(!flagts)
        {
            flagwhole=false;
            printf("wrong front coordinatesp n %d %f %f %f\n",i,decstr[i].tor[0],decstr[i].len[0],decstr[i].ang[0]);
        }
        decstr[i].ptn.x=pn.x;decstr[i].ptn.y=pn.y;decstr[i].ptn.z=pn.z;
        //ca
        flagts=tor2pos22(decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
            decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
            decstr[i].tor[1]*raddeg,decstr[i].len[1],decstr[i].ang[1]*raddeg,&pt.x,&pt.y,&pt.z);
        if(!flagts)
        {
            flagwhole=false;
            printf("wrong front coordinatesp ca %d %f %f %f\n",i,decstr[i].tor[1],decstr[i].len[1],decstr[i].ang[1]);
        }
        decstr[i].x=pt.x;decstr[i].y=pt.y;decstr[i].z=pt.z;
        //c
        flagts=tor2pos22(decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
            decstr[i].x,decstr[i].y,decstr[i].z,
            decstr[i].tor[2]*raddeg,decstr[i].len[2],decstr[i].ang[2]*raddeg,&pc.x,&pc.y,&pc.z);
        if(!flagts)
        {
            flagwhole=false;
            printf("wrong front coordinatesp c %d %f %f %f\n",i,decstr[i].tor[2],decstr[i].len[2],decstr[i].ang[2]);
        }
        decstr[i].ptc.x=pc.x;decstr[i].ptc.y=pc.y;decstr[i].ptc.z=pc.z;
    }
    return flagwhole;
}

bool tor2strp(vector<point3f> &decstr,int seqnum,int istart,int iend)
{
    int i;
    point3s pt,pn,pc;
//  BasicFunc bf;
    bool flagts;
    bool flagwhole=true;
    for(i=istart;i<=iend;i++)
    {
        if(decstr[i].tor[0]<0) 
        {
            if(decstr[i].stype=='E')
                decstr[i].tor[0]=120;
            else if(decstr[i].stype=='H')
                decstr[i].tor[0]=300;
            else
                decstr[i].tor[0]=120;
        }
        if(decstr[i].tor[1]<0) 
            decstr[i].tor[1]=180;
        if(decstr[i].tor[2]<0) 
            decstr[i].tor[2]=290;
        if(decstr[i].len[0]<0) 
            decstr[i].len[0]=float(lencn);
        if(decstr[i].len[1]<0) 
            decstr[i].len[1]=float(lennca);
        if(decstr[i].len[2]<0) 
            decstr[i].len[2]=float(lencac);
        if(decstr[i].ang[0]<0) 
            decstr[i].ang[0]=float(angcacn);
        if(decstr[i].ang[1]<0) 
            decstr[i].ang[1]=float(angcnca);
        if(decstr[i].ang[2]<0) 
            decstr[i].ang[2]=float(angncac);
        //0 1 2 n ca c
        //original phi i-1 psi i-1 omega i
        flagts=tor2pos22(decstr[i-1].ptn.x,decstr[i-1].ptn.y,decstr[i-1].ptn.z,decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
            decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
            decstr[i].tor[0]*raddeg,decstr[i].len[0],decstr[i].ang[0]*raddeg,&pn.x,&pn.y,&pn.z);
        if(!flagts)
        {
            flagwhole=false;
            printf("wrong front coordinatesp n %d %f %f %f\n",i,decstr[i].tor[0],decstr[i].len[0],decstr[i].ang[0]);
        }
        decstr[i].ptn.x=pn.x;decstr[i].ptn.y=pn.y;decstr[i].ptn.z=pn.z;
        //ca
        flagts=tor2pos22(decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
            decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
            decstr[i].tor[1]*raddeg,decstr[i].len[1],decstr[i].ang[1]*raddeg,&pt.x,&pt.y,&pt.z);
        if(!flagts)
        {
            flagwhole=false;
            printf("wrong front coordinatesp ca %d %f %f %f\n",i,decstr[i].tor[1],decstr[i].len[1],decstr[i].ang[1]);
        }
        decstr[i].x=pt.x;decstr[i].y=pt.y;decstr[i].z=pt.z;
        //c
        flagts=tor2pos22(decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
            decstr[i].x,decstr[i].y,decstr[i].z,
            decstr[i].tor[2]*raddeg,decstr[i].len[2],decstr[i].ang[2]*raddeg,&pc.x,&pc.y,&pc.z);
        if(!flagts)
        {
            flagwhole=false;
            printf("wrong front coordinatesp c %d %f %f %f\n",i,decstr[i].tor[2],decstr[i].len[2],decstr[i].ang[2]);
        }
        decstr[i].ptc.x=pc.x;decstr[i].ptc.y=pc.y;decstr[i].ptc.z=pc.z;
    }
    return flagwhole;
}


bool tor2strsg2(vector<point3f> &decstr,int seqnum,vector<vector<vector<vector<float>>>> &sgposdat)
{
    int i;
    int tind;
    point3s pn;
//    BasicFunc bf;
    bool flagts;
    bool flagwhole=true;
    for(i=1;i<seqnum;i++)
    {
    //  /*
        //atom o
//      flagts=bf.tor2pos22(decstr[i-1].ptn.x,decstr[i-1].ptn.y,decstr[i-1].ptn.z,decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
//          decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
//          decstr[i].tor[0]*raddeg-PI,lenco,angcaco*raddeg,&pn.x,&pn.y,&pn.z);
//      flagts=bf.tor2pos22(decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
//          decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
//          179.6715f*raddeg,1.2324f,2.1037f,&pn.x,&pn.y,&pn.z);
        flagts=tor2pos22(decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
            decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
            179.6715f*raddeg,1.229f,2.0961f,&pn.x,&pn.y,&pn.z);
        decstr[i-1].pto.x=pn.x;decstr[i-1].pto.y=pn.y;decstr[i-1].pto.z=pn.z;
        ///*
        //atom h
//      flagts=bf.tor2pos22(decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
//          decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,0,lennh,angcnh*raddeg,&pn.x,&pn.y,&pn.z);
        flagts=tor2pos22(decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,decstr[i].x,decstr[i].y,decstr[i].z,
            decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,179.8174f*raddeg,0.987f,2.0814f,&pn.x,&pn.y,&pn.z);//0.9919f,2.0574f
        if(!flagts)
        {
            flagwhole=false;
            printf("wrong front coordinates d %d\n",i);
        }
        decstr[i].pth.x=pn.x;decstr[i].pth.y=pn.y;decstr[i].pth.z=pn.z; 
    //  */
    }
//  /*
    //atom o
    i=seqnum;
//  flagts=bf.tor2pos22(decstr[i-1].ptn.x,decstr[i-1].ptn.y,decstr[i-1].ptn.z,decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
//      decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,0,lenco,angcaco*raddeg,&pn.x,&pn.y,&pn.z);
    flagts=tor2pos22(decstr[i-1].ptn.x,decstr[i-1].ptn.y,decstr[i-1].ptn.z,decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
        decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,0.0f,1.2439f, 2.0855f,&pn.x,&pn.y,&pn.z);
    if(!flagts)
    {
        flagwhole=false;
        printf("wrong front coordinates o %d\n",i);
    }
    decstr[i-1].pto.x=pn.x;decstr[i-1].pto.y=pn.y;decstr[i-1].pto.z=pn.z;
//  */  
//  /*
    //atom h 
    i=0;
//  flagts=bf.tor2pos22(decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,decstr[i].x,decstr[i].y,decstr[i].z,
//      decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,0,lennh,anghnca*raddeg,&pn.x,&pn.y,&pn.z);
    flagts=tor2pos22(decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,decstr[i].x,decstr[i].y,decstr[i].z,
            decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,60*raddeg,0.987f,2.0306f,&pn.x,&pn.y,&pn.z);//0.9972f
    if(!flagts)
    {
        flagwhole=false;
        printf("wrong front coordinates %d\n",i);
    }
    decstr[i].pth.x=pn.x;decstr[i].pth.y=pn.y;decstr[i].pth.z=pn.z;
//  */
//  /*
    //atom cb new
/*    for(i=0;i<seqnum;i++)
    {
        tind=aminoid(decstr[i].aaa);
        if(tind>19) tind=19;
        flagts=tor2pos22(decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
            decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,decstr[i].x,decstr[i].y,decstr[i].z,
            cbsta[tind][2]*raddeg,cbsta[tind][0],cbsta[tind][1],&pn.x,&pn.y,&pn.z);
        if(!flagts)
        {
            flagwhole=false;
            printf("wrong front coordinates cb2 %d\n",i);
        }
        decstr[i].ptb.x=pn.x;decstr[i].ptb.y=pn.y;decstr[i].ptb.z=pn.z;
        if(decstr[i].aaa=='G')
        {
            decstr[i].ptb.x=decstr[i].x;decstr[i].ptb.y=decstr[i].y;decstr[i].ptb.z=decstr[i].z;
        }
    }  */
//  */
    /*
    //atom cb
    for(i=0;i<seqnum;i++)
    {
        flagts=bf.tor2pos22(decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,decstr[i].x,decstr[i].y,decstr[i].z,
            decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,
            torncaccb*raddeg,lenccb,angcaccb*raddeg,&pn.x,&pn.y,&pn.z);
        if(!flagts)
        {
            flagwhole=false;
            printf("wrong front coordinates cb %d\n",i);
        }
        decstr[i].ptb.x=pn.x;decstr[i].ptb.y=pn.y;decstr[i].ptb.z=pn.z;
        if(decstr[i].aaa=='G')
        {
            decstr[i].ptb.x=decstr[i].x;decstr[i].ptb.y=decstr[i].y;decstr[i].ptb.z=decstr[i].z;
        }
    }
    */
    //atom sg
    int ti,tj;
    int cutnum=72;
    double delta=5.0;
    for(i=0;i<seqnum;i++)
    {
        if(decstr[i].aaa=='G' || decstr[i].aaa=='A')
        {
            decstr[i].ptsg.x=decstr[i].x;decstr[i].ptsg.y=decstr[i].y;decstr[i].ptsg.z=decstr[i].z;
            continue;
        }   
        tind=aminoid(decstr[i].aaa);
        if(tind>19) tind=19;    
        if(i<seqnum-1)
        {
            ti=int(decstr[i].tor[2]/delta);
            tj=int(decstr[i+1].tor[0]/delta);
            if(ti>=0 && ti<cutnum && tj>=0 && tj<cutnum)
            {
                flagts=tor2pos22(decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
                decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,
                decstr[i].x,decstr[i].y,decstr[i].z,
                sgposdat[2][tind][ti][tj],sgposdat[0][tind][ti][tj],sgposdat[1][tind][ti][tj],&pn.x,&pn.y,&pn.z);
                if(!flagts)
                {
                    flagwhole=false;
                    printf("wrong front coordinates full sg %d\n",i);
                }
                decstr[i].ptsg.x=pn.x;decstr[i].ptsg.y=pn.y;decstr[i].ptsg.z=pn.z;      
                continue;
            }
        }
        flagts=tor2pos22(decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
            decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,
            decstr[i].x,decstr[i].y,decstr[i].z,
            sglatavg2[tind][2]*raddeg,sglatavg2[tind][0],sglatavg2[tind][1],&pn.x,&pn.y,&pn.z);
        if(!flagts)
        {
            flagwhole=false;
            printf("wrong front coordinates sg %d\n",i);
        }
        decstr[i].ptsg.x=pn.x;decstr[i].ptsg.y=pn.y;decstr[i].ptsg.z=pn.z;      
    }
    //atom ct
/*    for(i=0;i<seqnum;i++)
    {
        tind=bf.aminoid(decstr[i].aaa);
        pn.x=0;pn.y=0;pn.z=0;
        pn.x+=decstr[i].ptn.x;
        pn.y+=decstr[i].ptn.y;
        pn.z+=decstr[i].ptn.z;
        pn.x+=decstr[i].x;
        pn.y+=decstr[i].y;
        pn.z+=decstr[i].z;
        pn.x+=decstr[i].ptc.x;
        pn.y+=decstr[i].ptc.y;
        pn.z+=decstr[i].ptc.z;
        pn.x+=decstr[i].pto.x;
        pn.y+=decstr[i].pto.y;
        pn.z+=decstr[i].pto.z;
        if(tind!=5)
        {
            pn.x+=decstr[i].ptb.x;
            pn.y+=decstr[i].ptb.y;
            pn.z+=decstr[i].ptb.z;
            pn.x+=decstr[i].ptsg.x*(sgatomnum[tind]-1);
            pn.y+=decstr[i].ptsg.y*(sgatomnum[tind]-1);
            pn.z+=decstr[i].ptsg.z*(sgatomnum[tind]-1);
        }
        else
        {
        }
        decstr[i].ptg.x=pn.x/double(sgatomnum[tind]+4.0);
        decstr[i].ptg.y=pn.y/double(sgatomnum[tind]+4.0);
        decstr[i].ptg.z=pn.z/double(sgatomnum[tind]+4.0);
    }  */

    return flagwhole;
}

bool loadsgpos2(char *filename,int ndim,vector<vector<vector<vector<float>>>> &sgpos2x)
{
    FILE *file;
    file=fopen(filename,"rt");
    if(!file)
    {
        printf("no sgpos2 file %s\n",filename);
        return false;
    }
    int i,j,k;
    char oneline[300];
    double ****sgpos2;
    if(!sgpos2)
    {
        sgpos2=new double***[6];
    
        for(i=0;i<6;i++)
        {
            sgpos2[i]=new double**[20];
            for(j=0;j<20;j++)
            {
                sgpos2[i][j]=new double*[180];
                for(k=0;k<180;k++)
                {
                    sgpos2[i][j][k]=new double[180];
                }       
            }       
        }   
    }  
    for(i=0;i<20;i++)
    {
        for(j=0;j<ndim;j++)
        {
            for(k=0;k<ndim;k++)
            {
                fgets(oneline,300,file);
                sscanf(oneline,"%lf %lf %lf %lf %lf %lf",&sgpos2[0][i][j][k],&sgpos2[1][i][j][k],&sgpos2[2][i][j][k]
                    ,&sgpos2[3][i][j][k],&sgpos2[4][i][j][k],&sgpos2[5][i][j][k]);
            }
        }
        if(i==19)
        printf("sgposition %d %f %f %f %f %f %f\n",i,sgpos2[0][i][ndim-1][ndim-1],sgpos2[1][i][ndim-1][ndim-1],sgpos2[2][i][ndim-1][ndim-1]
                    ,sgpos2[3][i][ndim-1][ndim-1],sgpos2[4][i][ndim-1][ndim-1],sgpos2[5][i][ndim-1][ndim-1]);
    }
    fclose(file);
    sgpos2x = vector<vector<vector<vector<float>>>> (6,vector<vector<vector<float>>>(20,vector<vector<float>>(ndim,vector<float>(ndim,0.0))));
    for(int i=0;i<6;i++)
    {
        for(int j=0;j<20;j++)
        {
            for(int k=0;k<ndim;k++)
            {
                for(int m=0;m<ndim;m++)
                {
                    sgpos2x[i][j][k][m] = (float) sgpos2[i][j][k][m];
                }
            }
        }
    }  
    for(int i=0;i<6;i++)
    {
        for(int j=0;j<20;j++)
        {
            for(int k=0;k<ndim;k++)
            {

                delete [] sgpos2[i][j][k];
            }
            delete [] sgpos2[i][j];
        }
        delete [] sgpos2[i];
    }       
    return true;
}

double energytorsion2(vector<point3f> &decstr,int numseq,vector<vector<vector<double>>> phipsidisss)
{
    int i,ti,tj;
    int tind,tss;
    double ttot=0;
//  double sf[3]={6.455290/7.702060,1,1};
    for(i=0;i<numseq-1;i++) if(decstr[i].ss2!='C')
    {
            if(decstr[i].ss2=='H' && decstr[i].stype!='E') tss=0;
            else if(decstr[i].ss2=='E' && decstr[i].stype!='H') tss=1;
            else tss=2;
            tind=decstr[i].iaa;
            if(tind>19) continue;
            ti=int(decstr[i].tor[2]/2.0);//phi
            ti=ti%180;
            if(ti<0) ti+=180;
  
            tj=int(decstr[i+1].tor[0]/2.0);//psi
            tj=tj%180;
            if(tj<0) tj+=180;
 
        //  ttot+=phipsidis[tind][ti][tj];
            ttot+=phipsidisss[tss*20+tind][ti][tj];
    }
    return ttot;
}

bool loadphipsidisss(char *filename,vector<vector<vector<double>>> &phipsidisss)
{
    FILE *file;
    if((file=fopen(filename,"rt"))==NULL)
    {
      printf("Unable to open phipsidisss %s\n",filename);
      return false;
    }
    phipsidisss = vector<vector<vector<double>>>(60,vector<vector<double>>(180,vector<double>(180,0.0)));
    int i,j,k;
    char oneline[200];
/*    for(i=0;i<60;i++)
    {
        if(!phipsidisss[i])
        {
            phipsidisss[i]=new double*[180];
            for(j=0;j<180;j++)
            {
                phipsidisss[i][j]=new double[180];
            }
        }
    } */

    for(i=0;i<60;i++)
    {
        for(j=0;j<180;j++)
        {
            for(k=0;k<180;k++)
            {
                fgets(oneline,200,file);
                sscanf(oneline,"%lf",&phipsidisss[i][j][k]);
                phipsidisss[i][j][k]-=11.51293;
        //      if(phipsidisss[i][j][k]<-5.5) phipsidisss[i][j][k]=-5.5;//useless
            }
        }
//      printf("%d %f\n",i,phipsidisss[i][0][0]);
    }
    fclose(file);
    return true;
}
int energyrama(vector<point3f> &decstr,int numseq,float &fene,vector<vector<vector<float>>> &rlogduke,vector<vector<vector<int>>> &ramaduke)
{
    int i;
    int ti,tj;
    int aatype;
    double totene=0;
    lableproblematic lp1,lp2;
    lp1.nn[7]=0;
    lp1.indn[7] = new int [numseq];
    lp1.indn[27] = new int [numseq];
    memset(lp1.indn[7],0,numseq*sizeof(int));    
    memset(lp1.indn[27],0,numseq*sizeof(int));
    for(i=1;i<numseq-1;i++)
    {
        if(decstr[i].iaa==5) aatype=1;
        else if(decstr[i].iaa==12) aatype=2;
        else if(decstr[i+1].iaa==12) aatype=3;
        else aatype=0;
        ti=int(decstr[i].tor[2]);
        if(ti<0) ti=0;
        else if(ti>359) ti=359;
        tj=int(decstr[i+1].tor[0]);
        if(tj<0) tj=0;
        else if(tj>359) tj=359;
//      if(!ramaout[ti][tj])
        if(ramaduke[aatype][ti][tj]>0)
        {
            lp1.indn[7][lp1.nn[7]]=i;
            lp1.nn[7]++;
            lp1.indn[27][i]=1;
        }
        totene+=rlogduke[aatype][ti][tj];
    }
    fene=totene;
    delete [] lp1.indn[7];
    delete [] lp1.indn[27];
    return lp1.nn[7];
}
bool loadca2ncbins(char *filename,vector<vector<vector<vector<float> > > > &cancbinsx)
{
    cancbinsx = vector<vector<vector<vector<float> > > >(13,vector<vector<vector<float> > >(105,vector<vector<float> >(105,vector<float>(160,0.0))));
    float ***cancbins[13];
    FILE *file;
    file=fopen(filename,"rt");
    if(!file)
    {
        printf("no ca2nc file %s\n",filename);
        return false;
    }
    printf("loading ca2nc bins\n");
    int i,j,k,l;
    int ii,jj,kk;
    char oneline[300];
    for(i=0;i<13;i++)
    {
        if(!cancbins[i])
        {
            cancbins[i]=new float**[105];
            for(j=0;j<105;j++)
            {
                cancbins[i][j]=new float *[105];
                for(k=0;k<105;k++)
                    cancbins[i][j][k]=new float[160];
            }
        }
    }
    for(j=0;j<105;j++)
    {
            for(k=0;k<105;k++)
            {
                for(l=0;l<160;l++)
                {
                    fgets(oneline,300,file);
                    sscanf(oneline,"%d %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f",
                        &ii,&jj,&kk,&cancbins[0][j][k][l],&cancbins[1][j][k][l],&cancbins[2][j][k][l],&cancbins[3][j][k][l]
                        ,&cancbins[4][j][k][l],&cancbins[5][j][k][l],&cancbins[6][j][k][l],&cancbins[7][j][k][l],&cancbins[8][j][k][l]
                        ,&cancbins[9][j][k][l],&cancbins[10][j][k][l],&cancbins[11][j][k][l],&cancbins[12][j][k][l]);
                }
            }
    }
    fclose(file);

    for(int m1=0;m1<13;m1++)
    {
        for(int m2=0;m2<105;m2++)
        {
            for(int m3=0;m3<105;m3++)
            {
                for(int m4=0;m4<160;m4++)
                {
                    cancbinsx[m1][m2][m3][m4] = cancbins[m1][m2][m3][m4];
                }
            }
        }
    }

    cout<<"yyy"<<endl;
    for(i=0;i<13;i++)
    {
         for(j=0;j<105;j++)
         {
            for(k=0;k<105;k++)
            {
                delete [] cancbins[i][j][k];
            }
            delete [] cancbins[i][j];
         }
         delete [] cancbins[i];    
    }    
    return true;
}

bool loadca2ncbins2(char *filename,vector<vector<vector<vector<float> > > > &cancbinsx)
{
//    if(cancbins[0]) return true;
    cancbinsx = vector<vector<vector<vector<float> > > >(13,vector<vector<vector<float> > >(105,vector<vector<float> >(105,vector<float>(160,0.0))));
    float ***cancbins[13];
    FILE *file;
    file=fopen(filename,"rb");
    if(!file)
    {
        printf("no ca2nc file %s\n",filename);
        return false;
    }
    printf("loading ca2nc bins\n");
    int i,j,k;
    for(i=0;i<13;i++)
    {
//        if(!cancbins[i])
//        {
            cancbins[i]=new float**[105];
            for(j=0;j<105;j++)
            {
                cancbins[i][j]=new float *[105];
                for(k=0;k<105;k++)
                    cancbins[i][j][k]=new float[160];
            }
//        }
    }

    for(i=0;i<13;i++)
    {
        if(i%2==0) continue;
        for(j=0;j<105;j++)
        {
            for(k=0;k<105;k++)
            {
                fread(cancbins[i][j][k],sizeof(float),160,file);
            }
        }
    }

    for(int m1=0;m1<13;m1++)
    {
        for(int m2=0;m2<105;m2++)
        {
            for(int m3=0;m3<105;m3++)
            {
                for(int m4=0;m4<160;m4++)
                {
                    cancbinsx[m1][m2][m3][m4] = cancbins[m1][m2][m3][m4];
                }
            }
        }
    }

    for(i=0;i<13;i++)
    {
         for(j=0;j<105;j++)
         {
            for(k=0;k<105;k++)
            {
                delete [] cancbins[i][j][k];
            }
            delete [] cancbins[i][j];
         }
         delete [] cancbins[i];    
    }
    fclose(file);
    return true;
}
void ca2nc(vector<point3f> &decstr,int numseq,vector<vector<vector<vector<float> > > > &cancbins)
{
    int j;
    point3d tp[5];
    point3s pout;
    int ii,jj,kk;
//    BasicFunc bf;
    double tdist;
    ////////////////////////////head
    ii=26;jj=26;kk=21;
    tp[1]=setv(decstr[0].x,decstr[0].y,decstr[0].z);
    tp[2]=setv(decstr[1].x,decstr[1].y,decstr[1].z);
    tp[3]=setv(decstr[2].x,decstr[2].y,decstr[2].z);
    tor2pos22(tp[3].x,tp[3].y,tp[3].z,tp[2].x,tp[2].y,tp[2].z,tp[1].x,tp[1].y,tp[1].z,
        float(kk)/25.0f,3.813f,float(ii+53.0)/50.0f,&pout.x,&pout.y,&pout.z);
    tp[0]=setv(pout.x,pout.y,pout.z);
    tor2pos22(tp[0].x,tp[0].y,tp[0].z,tp[2].x,tp[2].y,tp[2].z,tp[1].x,tp[1].y,tp[1].z,
        cancbins[5][ii][jj][kk]*raddeg,cancbins[1][ii][jj][kk],cancbins[3][ii][jj][kk],&pout.x,&pout.y,&pout.z);
    decstr[0].ptc.x=pout.x; 
    decstr[0].ptc.y=pout.y;
    decstr[0].ptc.z=pout.z;

    tor2pos22(tp[3].x,tp[3].y,tp[3].z,tp[1].x,tp[1].y,tp[1].z,tp[2].x,tp[2].y,tp[2].z,
        cancbins[11][ii][jj][kk]*raddeg,cancbins[7][ii][jj][kk],cancbins[9][ii][jj][kk],&pout.x,&pout.y,&pout.z);
    decstr[1].ptn.x=pout.x;
    decstr[1].ptn.y=pout.y;
    decstr[1].ptn.z=pout.z; 
            
    tor2pos22(decstr[1].ptn.x,decstr[1].ptn.y,decstr[1].ptn.z,decstr[0].ptc.x,decstr[0].ptc.y,decstr[0].ptc.z,
        tp[1].x,tp[1].y,tp[1].z,180.0f*raddeg,1.460f,111.008f*raddeg,&pout.x,&pout.y,&pout.z);
    decstr[0].ptn.x=pout.x;
    decstr[0].ptn.y=pout.y;
    decstr[0].ptn.z=pout.z;

    ////////////////////////////////////tail
    tp[0]=setv(decstr[numseq-3].x,decstr[numseq-3].y,decstr[numseq-3].z);
    tp[1]=setv(decstr[numseq-2].x,decstr[numseq-2].y,decstr[numseq-2].z);
    tp[2]=setv(decstr[numseq-1].x,decstr[numseq-1].y,decstr[numseq-1].z);
    tor2pos22(tp[0].x,tp[0].y,tp[0].z,tp[1].x,tp[1].y,tp[1].z,tp[2].x,tp[2].y,tp[2].z,
        float(kk)/25.0f,3.813f,float(ii+53.0)/50.0f,&pout.x,&pout.y,&pout.z);
    tp[3]=setv(pout.x,pout.y,pout.z);
    tor2pos22(tp[0].x,tp[0].y,tp[0].z,tp[2].x,tp[2].y,tp[2].z,tp[1].x,tp[1].y,tp[1].z,
        cancbins[5][ii][jj][kk]*raddeg,cancbins[1][ii][jj][kk],cancbins[3][ii][jj][kk],&pout.x,&pout.y,&pout.z);
    decstr[numseq-2].ptc.x=pout.x;
    decstr[numseq-2].ptc.y=pout.y;
    decstr[numseq-2].ptc.z=pout.z;

    tor2pos22(tp[3].x,tp[3].y,tp[3].z,tp[1].x,tp[1].y,tp[1].z,tp[2].x,tp[2].y,tp[2].z,
        cancbins[11][ii][jj][kk]*raddeg,cancbins[7][ii][jj][kk],cancbins[9][ii][jj][kk],&pout.x,&pout.y,&pout.z);
    decstr[numseq-1].ptn.x=pout.x;
    decstr[numseq-1].ptn.y=pout.y;
    decstr[numseq-1].ptn.z=pout.z;
            
    tor2pos22(decstr[numseq-2].ptc.x,decstr[numseq-2].ptc.y,decstr[numseq-2].ptc.z,
        decstr[numseq-1].ptn.x,decstr[numseq-1].ptn.y,decstr[numseq-1].ptn.z,
        tp[2].x,tp[2].y,tp[2].z,180.0f*raddeg,1.525f,111.008f*raddeg,&pout.x,&pout.y,&pout.z);
    decstr[numseq-1].ptc.x=pout.x;
    decstr[numseq-1].ptc.y=pout.y;
    decstr[numseq-1].ptc.z=pout.z;
    for(j=0;j<numseq-3;j++)
    {
        tp[0]=setv(decstr[j].x,decstr[j].y,decstr[j].z);
        tp[1]=setv(decstr[j+1].x,decstr[j+1].y,decstr[j+1].z);
        tp[2]=setv(decstr[j+2].x,decstr[j+2].y,decstr[j+2].z);
        tp[3]=setv(decstr[j+3].x,decstr[j+3].y,decstr[j+3].z);
        tdist=angv(minu(tp[0],tp[1]),minu(tp[2],tp[1]));
        ii=int(tdist*50.0)-53;
        if(ii<0) ii=0;
        else if(ii>102) ii=102;
            
        tdist=angv(minu(tp[3],tp[2]),minu(tp[1],tp[2]));
        jj=int(tdist*50.0)-53;
        if(jj<0) jj=0;
        else if(jj>102) jj=102;

        tdist=phi(tp[0].x,tp[0].y,tp[0].z,tp[1].x,tp[1].y,tp[1].z,tp[2].x,tp[2].y,tp[2].z,tp[3].x,tp[3].y,tp[3].z)*raddeg;
        kk=int(tdist*25.0);
            
        tor2pos22(tp[0].x,tp[0].y,tp[0].z,tp[2].x,tp[2].y,tp[2].z,tp[1].x,tp[1].y,tp[1].z,
            cancbins[5][ii][jj][kk]*raddeg,cancbins[1][ii][jj][kk],cancbins[3][ii][jj][kk],&pout.x,&pout.y,&pout.z);
//      bf.tor2pos22(tp[0].x,tp[0].y,tp[0].z,tp[2].x,tp[2].y,tp[2].z,tp[1].x,tp[1].y,tp[1].z,
//          totbin[5][aa+ii][bb+jj][cc+kk]*raddeg,1.525,0.359,&pout.x,&pout.y,&pout.z);
        decstr[j+1].ptc.x=pout.x;
        decstr[j+1].ptc.y=pout.y;
        decstr[j+1].ptc.z=pout.z;

        tor2pos22(tp[3].x,tp[3].y,tp[3].z,tp[1].x,tp[1].y,tp[1].z,tp[2].x,tp[2].y,tp[2].z,
            cancbins[11][ii][jj][kk]*raddeg,cancbins[7][ii][jj][kk],cancbins[9][ii][jj][kk],&pout.x,&pout.y,&pout.z);
        tor2pos22(decstr[j+1].ptc.x,decstr[j+1].ptc.y,decstr[j+1].ptc.z,tp[1].x,tp[1].y,tp[1].z,tp[2].x,tp[2].y,tp[2].z,
            180.0f*raddeg,1.460f,0.275f,&pout.x,&pout.y,&pout.z);
        decstr[j+2].ptn.x=pout.x;
        decstr[j+2].ptn.y=pout.y;
        decstr[j+2].ptn.z=pout.z;   
    }
    //get tor[0]
    int i;
    for(i=1;i<numseq;i++)
    {
        decstr[i].tor[0]=phi(decstr[i-1].ptn.x,decstr[i-1].ptn.y,decstr[i-1].ptn.z,decstr[i-1].x,decstr[i-1].y,
        decstr[i-1].z,decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z);
    }
    decstr[0].tor[0]=180;

    //add o
    bool flagts;
    point3s pn;
    for(i=1;i<numseq;i++)
    {
//      flagts=bf.tor2pos22(decstr[i-1].ptn.x,decstr[i-1].ptn.y,decstr[i-1].ptn.z,decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
//          decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
//          decstr[i].tor[0]*raddeg-PI,lenco,angcaco*raddeg,&pn.x,&pn.y,&pn.z);
//      flagts=bf.tor2pos22(decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
//          decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
//          179.6715f*raddeg,1.2324f,2.1037f,&pn.x,&pn.y,&pn.z);
        flagts=tor2pos22(decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
            decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
            179.6715f*raddeg,1.229f,2.0961f,&pn.x,&pn.y,&pn.z);
        decstr[i-1].pto.x=pn.x;decstr[i-1].pto.y=pn.y;decstr[i-1].pto.z=pn.z;
        
    }
    i=numseq;
    flagts=tor2pos22(decstr[i-1].ptn.x,decstr[i-1].ptn.y,decstr[i-1].ptn.z,decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
        decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,0.0f,1.2439f,2.0855f,&pn.x,&pn.y,&pn.z);
    decstr[i-1].pto.x=pn.x;decstr[i-1].pto.y=pn.y;decstr[i-1].pto.z=pn.z;

}
bool loadramadukebn(char *datafile,vector<vector<vector<float>>> &rlogdukex,vector<vector<vector<int>>> &ramadukex)
{
    int i,j;
    FILE *file=fopen(datafile,"rb");
    if(!file)
    {
        printf(" no rama binary file %s\n",datafile);
        return false;
    }
    rlogdukex = vector<vector<vector<float>>>(4,vector<vector<float>>(360,vector<float>(360,0.0)));
    ramadukex = vector<vector<vector<int>>>(4,vector<vector<int>>(360,vector<int>(360,0)));
    double **rlogduke[4];
    int **ramaduke[4];
    for(i=0;i<4;i++)
    {
//        if(rlogduke[i]==NULL)
//        {
    //      printf("malloc log %d\n",i);
            rlogduke[i]=new double*[360];
            for(j=0;j<360;j++)
            {
                rlogduke[i][j]=new double[360];
            }
//        }
    }  
    for(i=0;i<4;i++)
    {
//        if(ramaduke[i]==NULL)
//        {
    //      printf("malloc rama %d\n",i);
            ramaduke[i]=new int*[360];
            for(j=0;j<360;j++)
            {
                ramaduke[i][j]=new int[360];
            }
//        }   
    }  
//    int *ram;
//    double *rlog;
    for(i=0;i<4;i++)
    {
        for(j=0;j<360;j++)
        {
            fread(ramaduke[i][j],sizeof(int),360,file);
//            fread(ram,sizeof(int),360,file);
        }
//      printf("%d %d %d\n",ramaduke[i][60][220],ramaduke[i][72][200],ramaduke[i][50][420]);
    }
    for(i=0;i<4;i++)
    {
        for(j=0;j<360;j++)
        {
            fread(rlogduke[i][j],sizeof(double),360,file);
        }
//      printf("%f %f %f\n",rlogduke[i][60][220],rlogduke[i][72][200],rlogduke[i][50][420]);
    }
    for(i=0;i<4;i++)
    {
        for(j=0;j<360;j++)
        {
            for(int k=0;k<360;k++)
            {
                ramadukex[i][j][k] = ramaduke[i][j][k];
                rlogdukex[i][j][k] = float (rlogduke[i][j][k]);
            }
        } 
    }
    for(i=0;i<4;i++)
    {
         for(j=0;j<360;j++)
         {
            delete [] ramaduke[i][j];
            delete [] rlogduke[i][j];
         }
         delete [] ramaduke[i];
         delete [] rlogduke[i];    
    }
//    delete [] ramaduke;
//    delete [] rlogduke;
    fclose(file);
    printf("loading four rama\n");
    return true;
}

int getmovetype(double tmov[],int totmov,double trandnum)
{
    int i;
    for(i=0;i<totmov;i++)
    {
        if(trandnum>=tmov[i] && trandnum<tmov[i+1])
        {
            return i;
        }
    }
    return 0;
}

bool mcfragsweepome(vector<point3f> &decstr,int numseq,float *phipsiprob[8][20])//omega 
{
//  int i;
    int trandpos;
    double trand2;
//  point3f *tmstr=new point3f[numseq];
//  ParseSeq ps;
    bool flagclash=false;
    bool flagtor;
    double threshclash=0.01;
    vector<point3f> tmstr;
//  int numclash;
do
{
//  for(i=0;i<numseq;i++)
//  {
    vector<point3f>().swap(tmstr);
    tmstr = decstr;
//      tmstr[i]=decstr[i];
//  }
//    memcpy(tmstr,decstr,numseq*sizeof(point3f));
//  accdispos[0]=dispos[0];
//  for(i=1;i<numseq;i++)
//  {
//      accdispos[i]=accdispos[i-1]+dispos[i];
//  }   
        //only in coil
//      trandpos=findindpos(accdispos,accdispos[numseq-1]*rand()/double(RAND_MAX+1.0),numseq);
//      if(tmstr[trandpos].stype!='C' || trandpos==numseq-1)
//      {
//          continue;
//      }
//      trandpos=int((numseq)*(rand()/double(RAND_MAX+1.0)));
        trandpos=int((numseq)*randf0and1());
        if(tmstr[trandpos].ssm!='C' && randf0and1()>0.1)
        {
//            summcome[2]++;
//          delete[]tmstr;
            return false;
        }
//      trandpos=bf.posinarray(probaccseq,numseq,Random());
//      if(trandpos<0) trandpos=0;
//      else if(trandpos>=numseq) trandpos=numseq-1;
//      trand2=double(8*(rand()/double(RAND_MAX+1.0)))-4;
        trand2=double(16*randf0and1())-8;
//      if(trand2==0)
//      {
//          continue;
//      }
        tmstr[trandpos].tor[1]+=trand2;
        if(tmstr[trandpos].tor[1]<0)
        {
            tmstr[trandpos].tor[1]+=360;
        }
        else if(tmstr[trandpos].tor[1]>=360)
        {
            tmstr[trandpos].tor[1]-=360;
        }
        if(tmstr[trandpos].tor[1]>190 || tmstr[trandpos].tor[1]<170)
        {
//            summcome[2]++;
//          delete[]tmstr;
            return false;
        }
        if(trandpos==0)
            flagtor=tor2str(tmstr,numseq,3);
        else
            flagtor=tor2strp(tmstr,numseq,trandpos);

//      if(trandpos<numseq-trandpos)
//      {
//          flagtor=ps.itor2strp(tmstr,numseq,0,trandpos+2);
//      }
//      else
//      {
//          flagtor=ps.tor2strp(tmstr,numseq,trandpos);
//      }
        if(!flagtor)
        {
            printf("tor2str wrong in ome %d\n",trandpos);
        }
//  numclash=ps.calcaclash(tmstr,numseq);
//  if(numclash<numseq*threshclash)
//  {
//      flagclash=true;
//  }
//  else flagclash=false;
        flagclash=true;
}while(!flagclash);
    
//  int numclash;
//  numclash=ps.calcaclash(tmstr,numseq);
//  if(numclash>numseq*threshclash)
//  {
//      summcome[2]++;
//      delete[]tmstr;
//      return false;
//  }
    //////////make decision
//  bool flagev=evcheck(tmstr,numseq);
//  if(!flagev) numev3[7][3]++;
//  if(flagev)
//  {
//      summcome[2]++;
//      return false;
//  }
    return true;

}

bool mcfragsweepphi(vector<point3f> &decstr,int numseq,float *phipsiprob[8][20])//phi 
{
//  int i;
    int trandpos;
    double trand2;
//  point3f *tmstr=new point3f[numseq];
//  ParseSeq ps;
    bool flagclash=false;
    bool flagtor;
    double threshclash=0.01;
    int inda,indp;
    vector<point3f> tmstr;
//  int numclash;
do
{
//  for(i=0;i<numseq;i++)
//  {
//      tmstr[i]=decstr[i];
//  }
    vector<point3f>().swap(tmstr);
    tmstr = decstr;
//    memcpy(tmstr,decstr,numseq*sizeof(point3f));
//  accdispos[0]=dispos[0];
//  for(i=1;i<numseq;i++)
//  {
//      accdispos[i]=accdispos[i-1]+dispos[i];
//  }   
        //only in coil
//      trandpos=findindpos(accdispos,accdispos[numseq-1]*rand()/double(RAND_MAX+1.0),numseq);
//      if(tmstr[trandpos].stype!='C' || trandpos==numseq-1)
//      {
//          continue;
//      }
//      trandpos=int((numseq)*(rand()/double(RAND_MAX+1.0)));
        trandpos=int((numseq)*randf0and1());
        if(tmstr[trandpos].ssm!='C' && randf0and1()>0.2)
        {
//            summcphi[2]++;
//          delete[]tmstr;
            return false;
        }
        //newposition
//      trandpos=bf.posinarray(probaccseq,numseq,Random());
//      if(trandpos<0) trandpos=0;
//      else if(trandpos>=numseq) trandpos=numseq-1;

//      trand2=double(20*(rand()/double(RAND_MAX+1.0)))-10;
//      trand2=double(20*Random())-10;
//      if(trand2==0)
//      {
//          continue;
//      }
//      tmstr[trandpos].tor[2]+=trand2;
        trand2=randf0and1();
        inda=tmstr[trandpos].iaa; 
        indp=findpos2(phipsiprob[3][inda],0,359,trand2);
        tmstr[trandpos].tor[2]=indp+0.5;
        if(tmstr[trandpos].tor[2]<0)
        {
            tmstr[trandpos].tor[2]+=360;
        }
        else if(tmstr[trandpos].tor[2]>=360)
        {
            tmstr[trandpos].tor[2]-=360;
        }
        if(trandpos==0)
            flagtor=tor2str(tmstr,numseq,3);
        else
            flagtor=tor2strp(tmstr,numseq,trandpos);
//      if(trandpos<numseq-trandpos)
//      {
//          flagtor=ps.itor2strp(tmstr,numseq,0,trandpos+2);
//      }
//      else
//      {
//          flagtor=ps.tor2strp(tmstr,numseq,trandpos);
//      }
        if(!flagtor)
        {
            printf("tor2str wrong in phi %d\n",trandpos);
        }
//  numclash=ps.calcaclash(tmstr,numseq);
//  if(numclash<numseq*threshclash)
//  {
//      flagclash=true;
//  }
//  else flagclash=false;
        flagclash=true;
}while(!flagclash);

//  int numclash;
//  numclash=ps.calcaclash(tmstr,numseq);
//  if(numclash>numseq*threshclash)
//  {
//      summcphi[2]++;
//      delete[]tmstr;
//      return false;
//  }
    //////////make decision
//  bool flagev=evcheck(tmstr,numseq);
//  if(!flagev) numev3[1][3]++;
//  if(flagev)
//  {
//      summcphi[2]++;
//      return false;
//  }
    return true;
}
bool mcfragsweeppsi(vector<point3f> &decstr,int numseq,float *phipsiprob[8][20])//psi 
{
//  int i;
    int trandpos;
    double trand2;
//  point3f *tmstr=new point3f[numseq];
//  ParseSeq ps;
    bool flagclash=false;
    bool flagtor;
    double threshclash=0.01;
    int inda,indp;
    vector<point3f> tmstr;
//  int numclash;
do
{
//  for(i=0;i<numseq;i++)
//  {
//      tmstr[i]=decstr[i];
//  }
    vector<point3f>().swap(tmstr);
    tmstr = decstr;
//  accdispos[0]=dispos[0];
//  for(i=1;i<numseq;i++)
//  {
//      accdispos[i]=accdispos[i-1]+dispos[i];
//  }   
        //only in coil
//      trandpos=findindpos(accdispos,accdispos[numseq-1]*rand()/double(RAND_MAX+1.0),numseq);
//      if(tmstr[trandpos].stype!='C')
//      {
//          continue;
//      }

//      trandpos=int((numseq)*(rand()/double(RAND_MAX+1.0)));
        trandpos=int((numseq)*randf0and1());
        if(tmstr[trandpos].ssm!='C' && randf0and1()>0.2)
        {
//            summcpsi[2]++;
//          delete[]tmstr;
            return false;
        }

//      trandpos=bf.posinarray(probaccseq,numseq,Random());
//      if(trandpos<0) trandpos=0;
//      else if(trandpos>=numseq) trandpos=numseq-1;

//      trand2=double(20*(rand()/double(RAND_MAX+1.0)))-10;
//      trand2=double(20*Random())-10;
//      if(trand2==0)
//      {
//          continue;
//      }
//      tmstr[trandpos].tor[0]+=trand2;
        if(trandpos==0)
        {
//            summcpsi[2]++;
//          delete[]tmstr;
            return false;
        }
        trand2=randf0and1();
        inda=tmstr[trandpos-1].iaa; 
        indp=findpos2(phipsiprob[7][inda],0,359,trand2);
        tmstr[trandpos].tor[0]=indp+0.5;
        
        if(tmstr[trandpos].tor[0]<0)
        {
            tmstr[trandpos].tor[0]+=360;
        }
        else if(tmstr[trandpos].tor[0]>=360)
        {
            tmstr[trandpos].tor[0]-=360;
        }
        if(trandpos==0)
            flagtor=tor2str(tmstr,numseq,3);
        else
            flagtor=tor2strp(tmstr,numseq,trandpos);
//      if(trandpos<numseq-trandpos)
//      {
//          flagtor=ps.itor2strp(tmstr,numseq,0,trandpos+2);
//      }
//      else
//      {
//          flagtor=ps.tor2strp(tmstr,numseq,trandpos);
//      }
        if(!flagtor)
        {
            printf("tor2str wrong in psi %d\n",trandpos);
        }
//  numclash=ps.calcaclash(tmstr,numseq);
//  if(numclash<numseq*threshclash)
//  {
//      flagclash=true;
//  }
//  else flagclash=false;
    flagclash=true;
}while(!flagclash);
    
//  int numclash;
//  numclash=ps.calcaclash(tmstr,numseq);
//  if(numclash>numseq*threshclash)
//  {
//      summcpsi[2]++;
//      delete[]tmstr;
//      return false;
//  }
    //////////make decision
//  bool flagev=evcheck(tmstr,numseq);
//  if(!flagev) numev3[2][3]++;
//  if(flagev)
//  {
//      summcpsi[2]++;
//      return false;
//  }
    return true;
}

int findpos2(double *p,int is,int ie,double pt)
{
    if(is==ie)
    {
        return is;
    }
    else if(is==ie-1)
    {
        if(pt<=p[is])
        {
            return is;
        }
        else
        {
            return ie;
        }
    }
    else//is!=im!=ie
    {
        int im=(is+ie)/2;
        if(pt==p[im])
        {
            return(im);
        }
        else if(pt<p[im])
        {
            return(findpos2(p,is,im,pt));
        }
        else
        {
            return(findpos2(p,im+1,ie,pt));
        }
    }
}
int findpos2(float *p,int is,int ie,float pt)
{
    if(is==ie)
    {
        return is;
    }
    else if(is==ie-1)
    {
        if(pt<=p[is])
        {
            return is;
        }
        else
        {
            return ie;
        }
    }
    else//is!=im!=ie
    {
        int im=(is+ie)/2;
        if(pt==p[im])
        {
            return(im);
        }
        else if(pt<p[im])
        {
            return(findpos2(p,is,im,pt));
        }
        else
        {
            return(findpos2(p,im+1,ie,pt));
        }
    }
}
/*
bool loadphipsiprob(char *filename,double *phipsiprob[8][20])
{
    FILE *file; 
    if((file=fopen(filename,"rt"))==NULL)
    {
      printf("Unable to open phipsidisprob %s\n",filename);
      return false;
    }  
    int i,j,k;
    char oneline[200];  
    for(i=0;i<8;i++)
    {
        for(j=0;j<20;j++)
            phipsiprob[i][j]=NULL;
    }       
    for(i=0;i<8;i++)
    {
        for(j=0;j<20;j++)
        {
            if(!phipsiprob[i][j])
                phipsiprob[i][j]=new [360];
        }
    } 
    for(i=0;i<8;i++)
    {
        for(j=0;j<20;j++)
        {
            for(k=0;k<360;k++)
            {
                fgets(oneline,200,file);
                sscanf(oneline,"%lf",&phipsiprob[i][j][k]);
            }
        }

    }
    fclose(file); 
    for(i=0;i<8;i++)
    {
        for(j=0;j<20;j++)
        {
            for(k=1;k<360;k++)
            {
                phipsiprob[i][j][k]+=float(phipsiprob[i][j][k-1]);
            }
        }
    }
    return true;
}
*/
int findpos2(vector<float> p,int is,int ie,float pt)
{
    if(is==ie)
    {
        return is;
    }
    else if(is==ie-1)
    {
        if(pt<=p[is])
        {
            return is;
        }
        else
        {
            return ie;
        }
    }
    else//is!=im!=ie
    {
        int im=float(is+ie)/2.0;
        if(pt==p[im])
        {
            return(im);
        }
        else if(pt<p[im])
        {
            return(findpos2(p,is,im,pt));
        }
        else
        {
            return(findpos2(p,im+1,ie,pt));
        }
    }
}
bool loadphipsiprob(char *filename,vector<vector<vector<float>>> &phipsiprob)
{
    FILE *file; 
    if((file=fopen(filename,"rt"))==NULL)
    {
      printf("Unable to open phipsidisprob %s\n",filename);
      return false;
    }  
    int i,j,k;
    char oneline[200]; 
    
    phipsiprob=vector<vector<vector<float>>> (8,vector<vector<float>>(20,vector<float>(360,0.0))); 
    double *phipsiprobx[8][20];
//    cout<<"TTT"<<endl;
    for(i=0;i<8;i++)
    {
        for(j=0;j<20;j++)
            phipsiprobx[i][j]=new double[360];
    }     
//    cout<<"FFF"<<endl;  
/*    for(i=0;i<8;i++)
    {
        for(j=0;j<20;j++)
        {
//            if(!phipsiprob[i][j])
                phipsiprobx[i][j]=new double[360];
        }
    } */

    int aw1,aw2,aw3;
    for(i=0;i<8;i++)
    {
        for(j=0;j<20;j++)
        {
            for(k=0;k<360;k++)
            {
                fgets(oneline,200,file);
            //    sscanf(oneline,"%lf",&phipsiprobx[i][j][k]);
                sscanf(oneline,"%lf %d %d %d",&phipsiprobx[i][j][k],&aw1,&aw2,&aw3);
            }
        }

    }
//    cout<<"ZZZ"<<endl;
    for(i=0;i<8;i++)
    {
        for(j=0;j<20;j++)
        {
            for(k=0;k<360;k++)
            {
                phipsiprob[i][j][k]=float(phipsiprobx[i][j][k]);
            }
        }

    }
    fclose(file); 
//    cout<<"GG"<<endl;
    for(i=0;i<8;i++)
    {
         for(j=0;j<20;j++)
         {
            delete [] phipsiprobx[i][j];
         }
//         delete [] phipsiprobx[i];   
    }    
    for(i=0;i<8;i++)
    {
        for(j=0;j<20;j++)
        {
            for(k=1;k<360;k++)
            {
                phipsiprob[i][j][k]+=float(phipsiprob[i][j][k-1]);
            }
        }
    }
//    cout<<"HH"<<endl;
    return true;
} 
bool mcfragsweeplen(vector<point3f> &decstr,int numseq)  //length 
{
//  int i;
    int trandpos;
    double trand2;
//  point3f *tmstr=new point3f[numseq];
//  ParseSeq ps;
    bool flagclash=false;
    bool flagtor;
    double threshclash=0.01;
    vector<point3f> tmstr;
//  int numclash;
do
{
//  for(i=0;i<numseq;i++)
//  {
//      tmstr[i]=decstr[i];
//  }
    vector<point3f>().swap(tmstr);
    tmstr = decstr;
//    memcpy(tmstr,decstr,numseq*sizeof(point3f));
//  accdispos[0]=dispos[0];
//  for(i=1;i<numseq;i++)
//  {
//      accdispos[i]=accdispos[i-1]+dispos[i];
//  }   
        //only in coil
//      trandpos=findindpos(accdispos,accdispos[numseq-1]*rand()/double(RAND_MAX+1.0),numseq);
//      if(tmstr[trandpos].stype!='C' || trandpos==numseq-1)
//      {
//          continue;
//      }
//      trandpos=int((numseq)*(rand()/double(RAND_MAX+1.0)));
        trandpos=int((numseq)*randf0and1());
        if(tmstr[trandpos].ssm!='C' && randf0and1()>0.1)
        {
//            summclen[2]++;
//          delete[]tmstr;
            return false;
        }
//      trandpos=bf.posinarray(probaccseq,numseq,Random());
//      if(trandpos<0) trandpos=0;
//      else if(trandpos>=numseq) trandpos=numseq-1;
//      int lentype=rand()%3;
    //  int lentype=genrand()%3;
        int lentype=int(randf0and1()*3.0);
//      trand2=double(0.056*(rand()/double(RAND_MAX+1.0)))-0.028;
        trand2=double(0.48*randf0and1())-0.24;
//      if(trand2==0)
//      {
//          continue;
//      }
        tmstr[trandpos].len[lentype]+=trand2;
        if((lentype==0 && fabs(tmstr[trandpos].len[lentype]-lencn)>3*delcn) ||
            (lentype==1 && fabs(tmstr[trandpos].len[lentype]-lennca)>3*delnca) ||
            (lentype==2 && fabs(tmstr[trandpos].len[lentype]-lencac)>3*delcac) )
        {
//            summclen[2]++;
//          delete[]tmstr;
            return false;
        }
        if(trandpos==0)
            flagtor=tor2str(tmstr,numseq,3);
        else
            flagtor=tor2strp(tmstr,numseq,trandpos);
//      if(trandpos<numseq-trandpos)
//      {
//          flagtor=ps.itor2strp(tmstr,numseq,0,trandpos+2);
//      }
//      else
//      {
//          flagtor=ps.tor2strp(tmstr,numseq,trandpos);
//      }
        if(!flagtor)
        {
            printf("tor2str wrong in len %d\n",trandpos);
        }
//  numclash=ps.calcaclash(tmstr,numseq);
//  if(numclash<numseq*threshclash)
//  {
//      flagclash=true;
//  }
//  else flagclash=false;
        flagclash=true;
}while(!flagclash);
    
//  int numclash;
//  numclash=ps.calcaclash(tmstr,numseq);
//  if(numclash>numseq*threshclash)
//  {
//      summclen[2]++;
//      delete[]tmstr;
//      return false;
//  }
    //////////make decision
//  bool flagev=evcheck(tmstr,numseq);
//  if(!flagev) numev3[8][3]++;
//  if(flagev)
//  {
//      summclen[2]++;
//      return false;
//  }
    return true;
}
bool mcfragsweepang(vector<point3f> &decstr,int numseq)//angle 
{
//  int i;
    int trandpos;
    double trand2;
//  point3f *tmstr=new point3f[numseq];
//  ParseSeq ps;
    bool flagclash=false;
    bool flagtor;
    double threshclash=0.01;
    vector<point3f> tmstr;
//  int numclash;
do
{
//  for(i=0;i<numseq;i++)
//  {
//      tmstr[i]=decstr[i];
//  }
    vector<point3f>().swap(tmstr);
    tmstr = decstr;
//    memcpy(tmstr,decstr,numseq*sizeof(point3f));
//  accdispos[0]=dispos[0];
//  for(i=1;i<numseq;i++)
//  {
//      accdispos[i]=accdispos[i-1]+dispos[i];
//  }   
        //only in coil
//      trandpos=findindpos(accdispos,accdispos[numseq-1]*rand()/double(RAND_MAX+1.0),numseq);
//      if(tmstr[trandpos].stype!='C' || trandpos==numseq-1)
//      {
//          continue;
//      }
//      trandpos=int((numseq)*(rand()/double(RAND_MAX+1.0)));
        trandpos=int((numseq)*randf0and1());
        if(tmstr[trandpos].ssm!='C' && randf0and1()>0.1)
        {
//            summcang[2]++;
//          delete[]tmstr;
            return false;
        }
//      trandpos=bf.posinarray(probaccseq,numseq,Random());
//      if(trandpos<0) trandpos=0;
//      else if(trandpos>=numseq) trandpos=numseq-1;
//      int angtype=rand()%3;
    //  int angtype=genrand()%3;
        int angtype=int(randf0and1()*3.0);
//      trand2=double(10*(rand()/double(RAND_MAX+1.0)))-5;
        trand2=double(20*randf0and1())-10;
//      if(trand2==0)
//      {
//          continue;
//      }
        tmstr[trandpos].ang[angtype]+=trand2;
        if( (angtype==0 && fabs(tmstr[trandpos].ang[angtype]-angcacn)>2.5*2.009)||
            (angtype==1 && fabs(tmstr[trandpos].ang[angtype]-angcnca)>2.5*2.227)||
            (angtype==2 && fabs(tmstr[trandpos].ang[angtype]-angncac)>2*2.818) )
        {
//            summcang[2]++;
//          delete[]tmstr;
            return false;
        }
        if(trandpos==0)
            flagtor=tor2str(tmstr,numseq,3);
        else
            flagtor=tor2strp(tmstr,numseq,trandpos);
//      if(trandpos<numseq-trandpos)
//      {
//          flagtor=ps.itor2strp(tmstr,numseq,0,trandpos+2);
//      }
//      else
//      {
//          flagtor=ps.tor2strp(tmstr,numseq,trandpos);
//      }
        if(!flagtor)
        {
            printf("tor2str wrong in ang %d\n",trandpos);
        }
//  numclash=ps.calcaclash(tmstr,numseq);
//  if(numclash<numseq*threshclash)
//  {
//      flagclash=true;
//  }
//  else flagclash=false;
        flagclash=true;
}while(!flagclash);
    
//  int numclash;
//  numclash=ps.calcaclash(tmstr,numseq);
//  if(numclash>numseq*threshclash)
//  {
//      summcang[2]++;
//      delete[]tmstr;
//      return false;
//  }
    //////////make decision
//  bool flagev=evcheck(tmstr,numseq);
//  if(!flagev) numev3[9][3]++;
//  if(flagev)
//  {
//      summcang[2]++;
//      return false;
//  }
    return true;
}
/*bool movementang(vector<point3f> &decstr,int numseq,int resst,int resend)
{
    float meanstda[][4]={
        //  {116.9,1.5*4,116.9-1.5*4,116.9+1.5*4},//cacn
        //  {122.6,5.0*4,122.6-5.0*4,122.6+5.0*4},//cnca
        //  {111.8,2.5*4,111.8-2.5*4,111.8+2.5*4},//ncac    
        {117.2,2.2*4.0,117.2-2.2*4.0,117.2+2.2*4.0},//cacn
        {121.7,2.5*4.0,121.7-2.5*4.0,121.7+2.5*4.0},//cnca
        {111.0,2.7*4.0,111.0-2.7*4.0,111.0+2.7*4.0},//ncac
    };    
    int angtype=int(randf0and1()*3.0);
    decstr[resst].ang[angtype] = meanstda[angtype][0] + meanstda[angtype][1]*(2*randf0and1()-1);
    tor2strp(decstr,numseq,resst,resend);
}  */
/*
bool genesse(vector<point3f> &decstr,int numseq,vector<sssegment> &sse,int &numsse)
{
    int i,j;
//    if(sse)
//    {
//        delete[]sse;
//        sse=NULL;
//    }
//    sse=new sssegment[numseq];
    sse = vector<sssegment>(numseq);
    numsse=0;
    for(i=0;i<numseq;i++)
    {
        j=i;
        while(j<numseq && decstr[j].ss2==decstr[i].ss2)
        {
            j++;
        }
        sse[numsse].init=i;
        sse[numsse].term=j-1;
        sse[numsse].ss=decstr[i].ss2;
        numsse++;
        i=j-1;
    }
//    if(numsse!=0)
//    sse=(sssegment *)realloc(sse,numsse*sizeof(sssegment));
    return true;
} */
/*
double calcassenergy(point3f *decstr,int numseq)
{
    double trmsd=0; 
    int i;
    for(i=0;i<20;i++)
    {
        enelist[i]=0;
    }
//    ps.tor2strsg(decstr,numseq,sgpos);
    // 1   2    3    4   5    6    7   8    9         10  11  12   13     14
    //ev   hb   tor  df  sg   sol  rg  ct   bep  bbab hp  lj   cl   hlaa  caca
    //0.03 0.03 0.40 1.0 0.10 4.00 1.0 0.60 0.10 1.0 0.05 0.30 0.50 0.20  
    //pr        si   pr  pr   si   si  pr   0.074 
    enelist[0]=1.0*energypartialrmsd(decstr,numseq);
    trmsd+=enelist[0];
    enelist[2]=0.40*energyhbondnhoc2ass(decstr,numseq);//put in the font 0.03
    trmsd+=enelist[2];//hb
    calcallenergyass(decstr,numseq,gsoldat,enelist);
    enelist[0]=1.0*energydisres(decstr,numseq);//put after calcallenergyass
    trmsd+=enelist[0];
    enelist[1]*=10.0;
    trmsd+=enelist[1];//ev
    trmsd+=enelist[4];//df
    trmsd+=enelist[5];//sg
//  trmsd+=enelist[6];//sol
    trmsd+=enelist[7];//rg
    trmsd+=enelist[14];
    enelist[3]=0.40*energytorsion2ass(decstr,numseq);//coil2helix and extended if large  0.40
    trmsd+=enelist[3];//tor
    return trmsd;
}

float energypartialrmsd(vector<point3f> &decstr,int numseq,vector<boneinfo> &bb)
{
    float eneval=0;
    int i;
    int j=0;
    int lentmp = numseq;
    vector<double> ntstr(3*numseq);
    vector<double> tpstr(3*numseq);
    if(lentmp==0) return eneval;
    for(i=0;i<numseq;i++)
    {
        tpstr[j]=decstr[i].x;
        tpstr[lentmp+j]=decstr[i].y;
        tpstr[2*lentmp+j]=decstr[i].z;
        ntstr[j]= decstr[i].x;
        ntstr[lentmp+j]= decstr[i].y;
        ntstr[2*lentmp+j]=decstr[i].z;
        j++;
    }
    double pmat[9],ptrans[3];
    eneval=lsfrmsd(ntstr,tpstr,lentmp,pmat,ptrans);
//  eneval=pp.lsfrmsd(ntstr,tpstr,lentmp,5.0,pmat,ptrans);//more clash
    eneval=5.0*eneval*eneval;
    return eneval;
} 
 
double lsfrmsd(vector<double> pin, vector<double> pout, int np, double pmat[], double ptrans[])
{
    int i;
    transpara(pin, pout, np, pmat, ptrans);
    double *ptmp=new double[3*np];
    double rmsd=0;
//    BasicFunc bf;
    trmul(pmat,pout,3,3,np,ptmp);
    for(i=0;i<np;i++)
    {
        ptmp[i]+=ptrans[0];
        ptmp[np+i]+=ptrans[1];
        ptmp[2*np+i]+=ptrans[2];
    }
    for(i=0;i<np;i++)
    {
        rmsd+=(ptmp[i]-pin[i])*(ptmp[i]-pin[i])+(ptmp[np+i]-pin[np+i])*(ptmp[np+i]-pin[np+i])+
            (ptmp[2*np+i]-pin[2*np+i])*(ptmp[2*np+i]-pin[2*np+i]);
    }
    rmsd=sqrt(rmsd/double(np));
    delete[]ptmp;
    return rmsd;
}
//3*np 3*np 3*3 3*1
int transpara(double pin[], double pout[], int np, double pmat[], double ptrans[])
{
    //Rpout+t=pin
    //cpin cpout
    //h=qout*qinT=uav
    //R=vT*uT
    //t=cpin-R*cpout
    int i,j;
    int rval;
//    BasicFunc bf;
    //center
    double cpin[3],cpout[3];
    for(i=0;i<3;i++)
    {
        cpin[i]=0;
        cpout[i]=0;
    }   
    for(i=0;i<np;i++)
    {
        for(j=0;j<3;j++)
        {
            cpin[j]+=pin[j*np+i];
            cpout[j]+=pout[j*np+i];
        }
    }
    for(i=0;i<3;i++)
    {
        cpin[i]/=double(np);
        cpout[i]/=double(np);
    }
    for(i=0;i<np;i++)
    {
        for(j=0;j<3;j++)
        {
            pin[j*np+i]-=cpin[j];
            pout[j*np+i]-=cpout[j];
        }
    }
    ///start
    double *tpin=new double[3*np];
    
    tranmat(pin,3,np,tpin);
    double hmat[9],umat[9],vmat[9],vt[9],ut[9];
    double eps=0.000001;

    trmul(pout,tpin,3,np,3,hmat);
    rval=muav(hmat,3,3,umat,vmat,eps,4);
    if(rval==-1) 
    {
        delete[]tpin;
        return rval;
    }
    tranmat(umat,3,3,ut);
    tranmat(vmat,3,3,vt);
    //begin svd
    double sign=sdet(hmat,3);// Kabsch algorithm
    int bb;
    if(sign>=0) 
    {
        bb=1;
    }
    else
    {
        bb=-1;
//      printf("mirror happens\n");
    }
    double bmat[9]={1,0,0,0,1,0,0,0,bb};
    double tmpmat[9];
    trmul(vt,bmat,3,3,3,tmpmat);
    trmul(tmpmat,ut,3,3,3,pmat);
    //end svd

    trmul(pmat,cpout,3,3,1,ptrans);
    for(i=0;i<3;i++)
    {
        ptrans[i]=cpin[i]-ptrans[i];
    }
    
    //recover 
    for(i=0;i<np;i++)
    {
        for(j=0;j<3;j++)
        {
            pin[j*np+i]+=cpin[j];
            pout[j*np+i]+=cpout[j];
        }
    }
    delete[]tpin;
    return rval;
} 
  void tranmat(double a[],int m,int n,double b[])
  {
      int i,j;
      for(i=0;i<m;i++)
      {
          for(j=0;j<n;j++)
          {
              b[j*m+i]=a[i*n+j];
          }
      }
  }
void trmul(double **a,double *b,int m,int n,double *c)
{
    int i,l;
    for (i=0; i<=m-1; i++)
    {
            c[i]=0.0;
            for (l=0; l<=n-1; l++)
            {
                c[i]=c[i]+a[i][l]*b[l];
            }
    }
    return;
}
void trmul(double **a,double **b,int m,int n,int k,double **c)
{
    int i,j,l;
    for (i=0; i<=m-1; i++)
    {
        for (j=0; j<=k-1; j++)
        { 
        //  u=i*k+j; 
            c[i][j]=0.0;
            for (l=0; l<=n-1; l++)
            {
                c[i][j]=c[i][j]+a[i][l]*b[l][j];
            }
        }
    }
    return;
}
//a[m*n] b[n*k] c[m*k]
void trmul(double a[],double b[],int m,int n,int k,double c[])
{
    int i,j,l,u;
    for (i=0; i<=m-1; i++)
    {
        for (j=0; j<=k-1; j++)
        { 
            u=i*k+j; c[u]=0.0;
            for (l=0; l<=n-1; l++)
            {
                c[u]=c[u]+a[i*n+l]*b[l*k+j];
            }
        }
    }
    return;

} 

int muav(double a[],int m,int n,double u[],double v[],double eps,int ka)
{
    int i,j,k,l,it,ll,kk,ix,iy,mm,nn,iz,m1,ks;
    double d,dd,t,sm,sm1,em1,sk,ek,b,c,shh,fg[2],cs[2];
    double *s,*e,*w;
    s=(double *)malloc(ka*sizeof(double));
    e=(double *)malloc(ka*sizeof(double));
    w=(double *)malloc(ka*sizeof(double));
    it=60; k=n;
    if (m-1<n) k=m-1;
    l=m;
    if (n-2<m) l=n-2;
    if (l<0) l=0;
    ll=k;
    if (l>k) ll=l;
    if (ll>=1)
    { 
        for (kk=1; kk<=ll; kk++)
        { 
            if (kk<=k)
            {
                d=0.0;
                for (i=kk; i<=m; i++)
                { 
                    ix=(i-1)*n+kk-1; d=d+a[ix]*a[ix];
                }
                s[kk-1]=sqrt(d);
                if (s[kk-1]!=0.0)
                {
                    ix=(kk-1)*n+kk-1;
                    if (a[ix]!=0.0)
                    { 
                        s[kk-1]=fabs(s[kk-1]);
                        if (a[ix]<0.0) s[kk-1]=-s[kk-1];
                    }
                    for (i=kk; i<=m; i++)
                    {
                        iy=(i-1)*n+kk-1;
                        a[iy]=a[iy]/s[kk-1];
                    }
                    a[ix]=1.0+a[ix];
                }
                s[kk-1]=-s[kk-1];
            }
            if (n>=kk+1)
            { 
                for (j=kk+1; j<=n; j++)
                { 
                    if ((kk<=k)&&(s[kk-1]!=0.0))
                    { 
                        d=0.0;
                        for (i=kk; i<=m; i++)
                        { 
                            ix=(i-1)*n+kk-1;
                            iy=(i-1)*n+j-1;
                            d=d+a[ix]*a[iy];
                        }
                        d=-d/a[(kk-1)*n+kk-1];
                        for (i=kk; i<=m; i++)
                        {
                            ix=(i-1)*n+j-1;
                            iy=(i-1)*n+kk-1;
                            a[ix]=a[ix]+d*a[iy];
                        }
                    }
                    e[j-1]=a[(kk-1)*n+j-1];
                }
            }
            if (kk<=k)
            {
                for (i=kk; i<=m; i++)
                { 
                    ix=(i-1)*m+kk-1; iy=(i-1)*n+kk-1;
                    u[ix]=a[iy];
                }
            }
            if (kk<=l)
            { 
                d=0.0;
                for (i=kk+1; i<=n; i++)
                    d=d+e[i-1]*e[i-1];
                e[kk-1]=sqrt(d);
                if (e[kk-1]!=0.0)
                { 
                    if (e[kk]!=0.0)
                    { 
                        e[kk-1]=fabs(e[kk-1]);
                        if (e[kk]<0.0) e[kk-1]=-e[kk-1];
                    }
                    for (i=kk+1; i<=n; i++)
                      e[i-1]=e[i-1]/e[kk-1];
                    e[kk]=1.0+e[kk];
                }
                e[kk-1]=-e[kk-1];
                if ((kk+1<=m)&&(e[kk-1]!=0.0))
                { 
                    for (i=kk+1; i<=m; i++) w[i-1]=0.0;
                     for (j=kk+1; j<=n; j++)
                         for (i=kk+1; i<=m; i++)
                              w[i-1]=w[i-1]+e[j-1]*a[(i-1)*n+j-1];
                    for (j=kk+1; j<=n; j++)
                      for (i=kk+1; i<=m; i++)
                      {
                          ix=(i-1)*n+j-1;
                          a[ix]=a[ix]-w[i-1]*e[j-1]/e[kk];
                      }
                }
                for (i=kk+1; i<=n; i++)
                  v[(i-1)*n+kk-1]=e[i-1];
              }
          }
      }
    mm=n;
    if (m+1<n) mm=m+1;
    if (k<n) s[k]=a[k*n+k];
    if (m<mm) s[mm-1]=0.0;
    if (l+1<mm) e[l]=a[l*n+mm-1];
    e[mm-1]=0.0;
    nn=m;
    if (m>n) nn=n;
    if (nn>=k+1)
      { for (j=k+1; j<=nn; j++)
          { for (i=1; i<=m; i++)
              u[(i-1)*m+j-1]=0.0;
            u[(j-1)*m+j-1]=1.0;
          }
      }
    if (k>=1)
      { for (ll=1; ll<=k; ll++)
          { kk=k-ll+1; iz=(kk-1)*m+kk-1;
            if (s[kk-1]!=0.0)
              { if (nn>=kk+1)
                  for (j=kk+1; j<=nn; j++)
                    { d=0.0;
                      for (i=kk; i<=m; i++)
                        { ix=(i-1)*m+kk-1;
                          iy=(i-1)*m+j-1;
                          d=d+u[ix]*u[iy]/u[iz];
                        }
                      d=-d;
                      for (i=kk; i<=m; i++)
                        { ix=(i-1)*m+j-1;
                          iy=(i-1)*m+kk-1;
                          u[ix]=u[ix]+d*u[iy];
                        }
                    }
                  for (i=kk; i<=m; i++)
                    { ix=(i-1)*m+kk-1; u[ix]=-u[ix];}
                  u[iz]=1.0+u[iz];
                  if (kk-1>=1)
                    for (i=1; i<=kk-1; i++)
                      u[(i-1)*m+kk-1]=0.0;
              }
            else
              { for (i=1; i<=m; i++)
                  u[(i-1)*m+kk-1]=0.0;
                u[(kk-1)*m+kk-1]=1.0;
              }
          }
      }
    for (ll=1; ll<=n; ll++)
      { kk=n-ll+1; iz=kk*n+kk-1;
        if ((kk<=l)&&(e[kk-1]!=0.0))
          { for (j=kk+1; j<=n; j++)
              { d=0.0;
                for (i=kk+1; i<=n; i++)
                  { ix=(i-1)*n+kk-1; iy=(i-1)*n+j-1;
                    d=d+v[ix]*v[iy]/v[iz];
                  }
                d=-d;
                for (i=kk+1; i<=n; i++)
                  { ix=(i-1)*n+j-1; iy=(i-1)*n+kk-1;
                    v[ix]=v[ix]+d*v[iy];
                  }
              }
          }
        for (i=1; i<=n; i++)
          v[(i-1)*n+kk-1]=0.0;
        v[iz-n]=1.0;
      }
    for (i=1; i<=m; i++)
    for (j=1; j<=n; j++)
      a[(i-1)*n+j-1]=0.0;
    m1=mm; it=60;
    while (1==1)
    { 
        if (mm==0)
        {
            ppp(a,e,s,v,m,n);
            free(s); free(e); free(w); return(1);
        }
        if (it==0)
        {
            ppp(a,e,s,v,m,n);
            free(s); free(e); free(w); return(-1);
        }
        kk=mm-1;
        while ((kk!=0)&&(fabs(e[kk-1])!=0.0))
        {
            d=fabs(s[kk-1])+fabs(s[kk]);
            dd=fabs(e[kk-1]);
            if (dd>eps*d) kk=kk-1;
            else e[kk-1]=0.0;
        }
        if (kk==mm-1)
        { 
            kk=kk+1;
            if (s[kk-1]<0.0)
            { 
                s[kk-1]=-s[kk-1];
                for (i=1; i<=n; i++)
                  { ix=(i-1)*n+kk-1; v[ix]=-v[ix];}
            }
            while ((kk!=m1)&&(s[kk-1]<s[kk]))
            { 
                d=s[kk-1]; s[kk-1]=s[kk]; s[kk]=d;
                if (kk<n)
                  for (i=1; i<=n; i++)
                    { ix=(i-1)*n+kk-1; iy=(i-1)*n+kk;
                      d=v[ix]; v[ix]=v[iy]; v[iy]=d;
                    }
                if (kk<m)
                  for (i=1; i<=m; i++)
                    { ix=(i-1)*m+kk-1; iy=(i-1)*m+kk;
                      d=u[ix]; u[ix]=u[iy]; u[iy]=d;
                    }
                kk=kk+1;
            }
            it=60;
            mm=mm-1;
          }
        else
          { ks=mm;
            while ((ks>kk)&&(fabs(s[ks-1])!=0.0))
              { d=0.0;
                if (ks!=mm) d=d+fabs(e[ks-1]);
                if (ks!=kk+1) d=d+fabs(e[ks-2]);
                dd=fabs(s[ks-1]);
                if (dd>eps*d) ks=ks-1;
                else s[ks-1]=0.0;
              }
            if (ks==kk)
              { kk=kk+1;
                d=fabs(s[mm-1]);
                t=fabs(s[mm-2]);
                if (t>d) d=t;
                t=fabs(e[mm-2]);
                if (t>d) d=t;
                t=fabs(s[kk-1]);
                if (t>d) d=t;
                t=fabs(e[kk-1]);
                if (t>d) d=t;
                sm=s[mm-1]/d; sm1=s[mm-2]/d;
                em1=e[mm-2]/d;
                sk=s[kk-1]/d; ek=e[kk-1]/d;
                b=((sm1+sm)*(sm1-sm)+em1*em1)/2.0;
                c=sm*em1; c=c*c; shh=0.0;
                if ((b!=0.0)||(c!=0.0))
                  { shh=sqrt(b*b+c);
                    if (b<0.0) shh=-shh;
                    shh=c/(b+shh);
                  }
                fg[0]=(sk+sm)*(sk-sm)-shh;
                fg[1]=sk*ek;
                for (i=kk; i<=mm-1; i++)
                  { sss(fg,cs);
                    if (i!=kk) e[i-2]=fg[0];
                    fg[0]=cs[0]*s[i-1]+cs[1]*e[i-1];
                    e[i-1]=cs[0]*e[i-1]-cs[1]*s[i-1];
                    fg[1]=cs[1]*s[i];
                    s[i]=cs[0]*s[i];
                    if ((cs[0]!=1.0)||(cs[1]!=0.0))
                      for (j=1; j<=n; j++)
                        { ix=(j-1)*n+i-1;
                          iy=(j-1)*n+i;
                          d=cs[0]*v[ix]+cs[1]*v[iy];
                          v[iy]=-cs[1]*v[ix]+cs[0]*v[iy];
                          v[ix]=d;
                        }
                    sss(fg,cs);
                    s[i-1]=fg[0];
                    fg[0]=cs[0]*e[i-1]+cs[1]*s[i];
                    s[i]=-cs[1]*e[i-1]+cs[0]*s[i];
                    fg[1]=cs[1]*e[i];
                    e[i]=cs[0]*e[i];
                    if (i<m)
                      if ((cs[0]!=1.0)||(cs[1]!=0.0))
                        for (j=1; j<=m; j++)
                          { ix=(j-1)*m+i-1;
                            iy=(j-1)*m+i;
                            d=cs[0]*u[ix]+cs[1]*u[iy];
                            u[iy]=-cs[1]*u[ix]+cs[0]*u[iy];
                            u[ix]=d;
                          }
                  }
                e[mm-2]=fg[0];
                it=it-1;
              }
            else
              { if (ks==mm)
                  { kk=kk+1;
                    fg[1]=e[mm-2]; e[mm-2]=0.0;
                    for (ll=kk; ll<=mm-1; ll++)
                      { i=mm+kk-ll-1;
                        fg[0]=s[i-1];
                        sss(fg,cs);
                        s[i-1]=fg[0];
                        if (i!=kk)
                          { fg[1]=-cs[1]*e[i-2];
                            e[i-2]=cs[0]*e[i-2];
                          }
                        if ((cs[0]!=1.0)||(cs[1]!=0.0))
                          for (j=1; j<=n; j++)
                            { ix=(j-1)*n+i-1;
                              iy=(j-1)*n+mm-1;
                              d=cs[0]*v[ix]+cs[1]*v[iy];
                              v[iy]=-cs[1]*v[ix]+cs[0]*v[iy];
                              v[ix]=d;
                            }
                      }
                  }
                else
                  { kk=ks+1;
                    fg[1]=e[kk-2];
                    e[kk-2]=0.0;
                    for (i=kk; i<=mm; i++)
                      { fg[0]=s[i-1];
                        sss(fg,cs);
                        s[i-1]=fg[0];
                        fg[1]=-cs[1]*e[i-1];
                        e[i-1]=cs[0]*e[i-1];
                        if ((cs[0]!=1.0)||(cs[1]!=0.0))
                          for (j=1; j<=m; j++)
                            { ix=(j-1)*m+i-1;
                              iy=(j-1)*m+kk-2;
                              d=cs[0]*u[ix]+cs[1]*u[iy];
                              u[iy]=-cs[1]*u[ix]+cs[0]*u[iy];
                              u[ix]=d;
                            }
                      }
                  }
              }
          }
      }
    return(1);
  }
 void ppp(double a[],double e[],double s[],double v[],int m,int n)
{ 
    int i,j,p,q;
    double d;
    if (m>=n) 
        i=n;
    else i=m;
    for (j=1; j<=i-1; j++)
    {
        a[(j-1)*n+j-1]=s[j-1];
        a[(j-1)*n+j]=e[j-1];
    }
    a[(i-1)*n+i-1]=s[i-1];
    if (m<n)
        a[(i-1)*n+i]=e[i-1];
    for (i=1; i<=n-1; i++)
    {
        for (j=i+1; j<=n; j++)
        {
            p=(i-1)*n+j-1; 
            q=(j-1)*n+i-1;
            d=v[p]; 
            v[p]=v[q];
            v[q]=d;
        }
    }
    return;
}
void sss(double fg[2],double cs[2])
{
    double r,d;
    if ((fabs(fg[0])+fabs(fg[1]))==0.0)
    { 
        cs[0]=1.0; 
        cs[1]=0.0; 
        d=0.0;
    }
    else 
    { 
        d=sqrt(fg[0]*fg[0]+fg[1]*fg[1]);
        if (fabs(fg[0])>fabs(fg[1]))
        { 
            d=fabs(d);
            if (fg[0]<0.0)
            {
                d=-d;
            }
        }
        if (fabs(fg[1])>=fabs(fg[0]))
        {
            d=fabs(d);
            if (fg[1]<0.0)
            {
                d=-d;
            }
        }
        cs[0]=fg[0]/d; 
        cs[1]=fg[1]/d;
    }
    r=1.0;
    if (fabs(fg[0])>fabs(fg[1]))
    {
        r=cs[1];
    }
    else
    {
        if (cs[0]!=0.0) 
        {
            r=1.0/cs[0];
        }
    }
    fg[0]=d; fg[1]=r;
    return;
}

double sdet(double a[],int n)//content in a[] will be changed
  { 
    int i,j,k,is,js,l,u,v;
    double f,det,q,d;
    f=1.0; det=1.0;
    for (k=0; k<=n-2; k++)
      { 
        q=0.0;
        for (i=k; i<=n-1; i++)
        for (j=k; j<=n-1; j++)
          {
            l=i*n+j; 
            d=fabs(a[l]);
            if (d>q) 
            { q=d; is=i; js=j;}
          }
        if (q+1.0==1.0)
          { det=0.0; 
        return(det);}
        if (is!=k)
          { f=-f;
            for (j=k; j<=n-1; j++)
              { u=k*n+j;
                v=is*n+j;
                d=a[u]; 
                a[u]=a[v];
                a[v]=d;
              }
          }
        if (js!=k)
          { f=-f;
            for (i=k; i<=n-1; i++)
              {
                u=i*n+js; v=i*n+k;
                d=a[u]; a[u]=a[v]; a[v]=d;
              }
          }
        l=k*n+k;
        det=det*a[l];
        for (i=k+1; i<=n-1; i++)
          { 
            d=a[i*n+k]/a[l];
            for (j=k+1; j<=n-1; j++)
              { u=i*n+j;
                a[u]=a[u]-d*a[k*n+j];
              }
          }
      }
    det=f*det*a[n*n-1];
    return(det);
  }
  double sdet2(double a[],int n)//content in a[] will be same
  { 
    int i,j,k,is,js,l,u,v;
    double f,det,q,d;
    double *abk=(double *)malloc(n*n*sizeof(double));
    memcpy(abk,a,n*n*sizeof(double));
    f=1.0; det=1.0;
    for (k=0; k<=n-2; k++)
      { 
        q=0.0;
        for (i=k; i<=n-1; i++)
        for (j=k; j<=n-1; j++)
          {
            l=i*n+j; 
            d=fabs(a[l]);
            if (d>q) 
            { q=d; is=i; js=j;}
          }
        if (q+1.0==1.0)
          { det=0.0; 
        return(det);}
        if (is!=k)
          { f=-f;
            for (j=k; j<=n-1; j++)
              { u=k*n+j;
                v=is*n+j;
                d=a[u]; 
                a[u]=a[v];
                a[v]=d;
              }
          }
        if (js!=k)
          { f=-f;
            for (i=k; i<=n-1; i++)
              {
                u=i*n+js; v=i*n+k;
                d=a[u]; a[u]=a[v]; a[v]=d;
              }
          }
        l=k*n+k;
        det=det*a[l];
        for (i=k+1; i<=n-1; i++)
          { 
            d=a[i*n+k]/a[l];
            for (j=k+1; j<=n-1; j++)
              { u=i*n+j;
                a[u]=a[u]-d*a[k*n+j];
              }
          }
      }
    det=f*det*a[n*n-1];
    memcpy(a,abk,n*n*sizeof(double));
    free(abk);
    return(det);
  }  */

