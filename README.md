# DockEM
How to compile this program?

    Users must compile the program before implementation by following step a) and b) sequentially.

    a)  > cd $pkgdir/DockEM-main/CaDensity/ObjexxFCL
        > g++ -c *.cc -I/$pkgdir/DockEM-main/CaDensity -std=c++11 -lstdc++
        > cd $pkgdir/DockEM-main/CaDensity/fourier
        > g++ -c *.cc -I/$pkgdir/DockEM-main/CaDensity -std=c++11 -lstdc++
        > cp $pkgdir/DockEM-main/CaDensity/ObjexxFCL/*.o  $pkgdir/DockEM-main/CaDensity/fourier
        > cd $pkgdir/DockEM-main/CaDensity/fourier
        > ar cru libfourierx.a  *.o

    b) use the relative path to compile in $pkgdir/DockEM-main
        > cd $pkgdir/DockEM-main
        > g++ edock_flexible_new1.cpp MathTools.cpp mol.cpp energy.cpp base_score.cpp match.cpp output.cpp translate_rotate.cpp REMC1.cpp  MC_simulation2.cpp initialmol.cpp flex_dock.cpp precon.cpp generaterandnum.cpp CaDensityx.cpp GeometryTools.cpp GetVdwRadius.cpp randomx.cpp SplineInterp.cpp xray_scattering.cpp -o DockEM -I ./CaDensity -L ./CaDensity/fourier -lfourierx -std=c++11 -lstdc++ -O3

How to run this program?
  
        An example of implementation:
        > cd $pkgdir/DockEM-main/1afkA_BS01_PAP/af2_edock
        > ../../DockEM af2_protein.mol2 super_lig_af2.mol2 bindingsite_pre.dat test1225 1afkA_BS01_PAP_map.mrc res 8.21 sam 0
       
	5 input paramerters:
	af2_protein.mol2 : protein mol2 file
	super_lig_af2.mol2 : ligand mol2 file, add the H
	bindingsite_pre.dat : binding sites file, amino acid residue index according to the protein index
	test1225 : energy output file name, save the docking results
	1afkA_BS01_PAP_map.mrc: density map mrc file
	res 8.21 : Resolution of the density map
                sam 0：sampling rate, set to 0 by default.
	 
	2 output file：
	test_protein.mol2：protein structure aligned with the density map
	af2_protein.mol2.mol2：docked ligand result



