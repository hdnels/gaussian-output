#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// This function can:
// 	Calculate all of the cation-anion bonds.
// 	Choose whether to use an atom position or the center as the origin.
// 	Optionally, specify a type of cation.

// Wave names must be basename_num, basename_atom, basename_x, basename_y, basename_z.

// Example:
// 	BondLengths("CIS_gs") outputs CIS_gs_bondlength and CIS_gs_bonddist
//		which includes all bond lengths in the NC and their positions relative to the NC center (average position of all atoms)
//	BondLengths("CIS_gs", type="Cu") outputs CIS_gs_bondlength_Cu and CIS_gs_bonddist_Cu
//		which includes all Cu-S bond lengths in the NC and their positions relative to the NC center
//	BondLengths("CIS_gs", center=5) outputs CIS_gs_bondlength and CIS_gs_bonddist
//		which includes all bond lengths in the NC and their positions relative to cation site 5
//		(no change in the output wave names, but it does print the center definition and coordinates)

Function BondLengths(basename, [type, center])
	string basename
	string type // type of cation to calculate
	variable center // atom number of the atom in the center
	
	wave numw = $(basename + "_num")
	wave/t atomw = $(basename + "_atom")
	wave xw = $(basename + "_x")
	wave yw = $(basename + "_y")
	wave zw = $(basename + "_z")
	
	// find the number of non-H atoms in the data set
	variable num_atoms = 0
	variable h = 0
	for (h = 0; h < numpnts(atomw); h += 1)
		if (StringMatch(atomw[h], "!H"))	
			num_atoms += 1
		endif
	endfor
	//print num_atoms
	
	variable maxbonds = num_atoms*2 // maximum bonds in the NC is 4 * number of cations
	make /n=(maxbonds)/O out_cation, out_anion, out_distance, out_cent_dist
	
	// define the center
	variable cent_x, cent_y, cent_z
	if (ParamIsDefault(center)) // if the user did not specify a central atom
		cent_x = mean(xw, 0, (num_atoms - 1))
		cent_y = mean(yw, 0, (num_atoms - 1))
		cent_z = mean(zw, 0, (num_atoms - 1))
		print "Centered at the NC center"
	else // if the user did specify a central atom
		cent_x = xw[center-1]
		cent_y = yw[center-1]
		cent_z = zw[center-1]
		print "Centered at atom ", atomw[center-1], center
	endif
	print cent_x, cent_y, cent_z
	
	
	variable i, j, k
	string atom_i, atom_j
	variable distance
	variable bond_x, bond_y, bond_z
	
	// To determine which atoms are bonded, find the minimum distance from any atom to atom 1.
	// Take 1.3 x that distance as the threshold for a bond.
	// This may need modification if the NC contains significantly different atom sizes.
	variable min_bondlength1 = 10
	for (i = 1; i < num_atoms; i += 1)
		distance = sqrt( (xw[0] - xw[i])^2 + (yw[0] - yw[i])^2 + (zw[0] - zw[i])^2)
		if (distance < min_bondlength1)
			min_bondlength1 = distance
		endif
	endfor
	variable bondcutoff = 1.55 * min_bondlength1
	print "Using cutoff ", bondcutoff
	
	
	// i: current atom 1 (cation)
	// j: current atom 2 (anion)
	// k: current bond
	if (ParamIsDefault(type)) // calculate for all atoms if no type is specified
		for (i=0; i<num_atoms; i+=1) // for all atoms in the NC
			atom_i = atomw[i]
			if (!isAnion(atom_i)) // if the atom is a cation
				for (j=0; j<num_atoms; j+=1) // go through all atoms in the NC
					if (isAnion(atomw[j])) // if the second atom is an anion
						distance = sqrt( (xw[i] - xw[j])^2 + (yw[i] - yw[j])^2 + (zw[i] - zw[j])^2) // calculate the distance
						if (distance < bondcutoff) // if the distance is < cutoff, it's a bond
							out_cation[k] = numw[i] // save the cation, anion, and bond length
							out_anion[k] = numw[j]
							out_distance[k] = distance
					
							bond_x = (xw[i] + xw[j])/2 // the xyz coordinates of the bond center
							bond_y = (yw[i] + yw[j])/2
							bond_z = (zw[i] + zw[j])/2
							
							// save the bond distance
							out_cent_dist[k] = sqrt( (bond_x - cent_x)^2 + (bond_y - cent_y)^2 + (bond_z - cent_z)^2)				
							k+=1 // move to the next bond
						endif
					endif
				endfor
			endif
		endfor
	
	else // calculate for a specific type of atom
		for (i=0; i<num_atoms; i+=1) // for all atoms in the NC
			atom_i = atomw[i]
			if (StringMatch(atom_i, type)) // if the atom matches the specified type
				for (j=0; j<num_atoms; j+=1) // go through all atoms in the NC
					if (isAnion(atomw[j])) // if the second atom is an anion
						distance = sqrt( (xw[i] - xw[j])^2 + (yw[i] - yw[j])^2 + (zw[i] - zw[j])^2) // calculate the distance
						if (distance < bondcutoff) // if the distance is < cutoff, it's a bond
							out_cation[k] = numw[i] // save the cation, anion, and bond length
							out_anion[k] = numw[j]
							out_distance[k] = distance
					
							bond_x = (xw[i] + xw[j])/2 // the xyz coordinates of the bond center
							bond_y = (yw[i] + yw[j])/2
							bond_z = (zw[i] + zw[j])/2
							
							// save the bond distance
							out_cent_dist[k] = sqrt( (bond_x - cent_x)^2 + (bond_y - cent_y)^2 + (bond_z - cent_z)^2)				
							k+=1 // move to the next bond
						endif
					endif
				endfor
			endif
		endfor
		
	endif
	
	
	variable c = 0
	// delete extra points in the bond-related waves
	for (c=0; c<maxbonds; c+=1) // go through the bond length wave
		if (out_distance[c] == 0) // when you find a zero, delete the following points from all bond waves
			DeletePoints c, (maxbonds - c), out_distance, out_anion, out_cation, out_cent_dist
			Break
		endif
	endfor
	
	
	// rename the output waves (overwriting them if they already exist)
	if (ParamIsDefault(type)) // if no cation type was specified
		Duplicate /O out_cation, $(basename + "_cation")
		Duplicate /O out_anion, $(basename + "_anion")
		Duplicate /O out_distance, $(basename + "_bondlength")
		Duplicate /O out_cent_dist, $(basename + "_bonddist")
		print "Calculated ", (numpnts(out_cation)), " bonds"
	else // if a cation type was specified
		Duplicate /O out_cation, $(basename + "_cation_" + type)
		Duplicate /O out_anion, $(basename + "_anion_" + type)
		Duplicate /O out_distance, $(basename + "_bondlength_" + type)
		Duplicate /O out_cent_dist, $(basename + "_bonddist_" + type)
		print "Calculated ", (numpnts(out_cation)), " bonds to ", type
	endif
	
	killwaves out_cation, out_anion, out_distance, out_cent_dist
	
end





function isAnion(atomtype) // currently covers S, Se, and O
	string atomtype
	if (StringMatch(atomtype, "S"))
		return 1
	elseif (StringMatch(atomtype, "Se"))
		return 1
	elseif (StringMatch(atomtype, "O"))
		return 1
	else
		return 0
	endif
end





// Load and rename waves from file
Function LoadXYZ(filename, basename)
	string filename, basename
	NewPath /M="Select data folder."/Q filepath
	LoadWave /A/D/J/K=0/N/O/W/P=filepath filename
	// A auto-name and go
	// D double-precision
	// N overwrite with auto-name and go
	// J delimited text
	// W wave names in file
	// O overwrite
	// K=0

	wave Num, Atom, xW, yW, zW
	
	duplicate /O Num $(basename + "_num")
	duplicate /O/T Atom $(basename + "_atom")
	duplicate /O xW $(basename + "_x")
	duplicate /O yW $(basename + "_y")
	duplicate /O zW $(basename + "_z")
	
	killwaves Num, Atom, xW, yW, zW
end



Function bonddiffs(basename1, basename2, outname)
	string basename1, basename2, outname
	
	if (waveexists($(basename1 + "_bondlength_Cu")))
		wave Cu_bonds1 = $(basename1 + "_bondlength_Cu")
		wave Cu_bonds2 = $(basename2 + "_bondlength_Cu")
		make /O/n=(numpnts(Cu_bonds1)) out_Cu = Cu_bonds2 - Cu_bonds1
		duplicate /O out_Cu $(outname + "_Cu_diff")
		killwaves out_Cu
	endif
		
	if (waveexists($(basename1 + "_bondlength_Zn")))
		wave Zn_bonds1 = $(basename1 + "_bondlength_Zn")
		wave Zn_bonds2 = $(basename2 + "_bondlength_Zn")
		make /O/n=(numpnts(Zn_bonds1)) out_Zn = Zn_bonds2 - Zn_bonds1
		duplicate /O out_Zn $(outname + "_Zn_diff")
		killwaves out_Zn
	endif
	
	if (waveexists($(basename1 + "_bondlength_In")))
		wave In_bonds1 = $(basename1 + "_bondlength_In")
		wave In_bonds2 = $(basename2 + "_bondlength_In")
		make /O/n=(numpnts(In_bonds1)) out_In = In_bonds2 - In_bonds1
		duplicate /O out_In $(outname + "_In_diff")
		killwaves out_In
	endif
	

end