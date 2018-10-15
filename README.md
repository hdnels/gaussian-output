# gaussian-output

Scripts for working with Gaussian input and output files, specifically for DFT calculations on semiconductor nanocrystals.

**getDOS.py**: Extract density of states and atomic-orbital contributions to molecular orbitals from a Gaussian output file. Loosely based on mkDOS.py from the Li research group in the Department of Chemistry at the University of Washington, but entirely rewritten and restructured for speed with added options and cleaner formatting. Part of these calculations may be redundant with the Gaussian keyword Population=Orbitals.

**opt_to_dos.py**: Take the output geometry from a Gaussian optimization output file and use it to set up an input file for a density of states calculation. Input file settings are based on the stf queue on the Hyak supercomputing cluster at the University of Washington.

**BondLengths_general.ipf**: Igor Pro 6 file for calculating bond lengths and geometric distortions in a nanocrystal, given the Cartesian coordinates of its atoms. Used in references 1 and 2.

**CuInZnS_random_generator.ipynb**: Scripts used for generating a series of alloyed 34-cation nanocrystals ranging from ZnS to CuInS<sub>2</sub>, as in reference 3.

**cube_remove_atoms.py**: Removes all except the specified atoms from the structure in a Gaussian cube file, for visualization (in GaussView) of structure fragments with localized molecular orbitals as in reference 3.




### References and examples:
[1. Nelson, Li, Gamelin. "Computational Studies of the Electronic Structures of Copper-Doped CdSe Nanocrystals: Oxidation States, Jahn-Teller Distortions, Vibronic Bandshapes, and Singlet-Triplet Splittings." J. Phys. Chem. C 2016, 120, 5714-5723.](https://pubs.acs.org/doi/10.1021/acs.jpcc.5b11319)

[2. Nelson, Hinterding, Fainblat, Creutz, Li, Gamelin. "Mid-Gap States and Normal vs Inverted Bonding in Luminescent Cu+- and Ag+-Doped CdSe Nanocrystals." J. Am. Chem. Soc. 2017, 139, 6411-6421.](https://pubs.acs.org/doi/10.1021/jacs.7b01924)

[3. Nelson, Gamelin. "Valence-Band Electronic Structures of Cu<sup>+</sup>-Doped ZnS, Alloyed Cu−In−Zn−S, and Ternary CuInS<sub>2</sub> Nanocrystals: A Unified Description of Photoluminescence Across Compositions." J. Phys. Chem. C 2018, 122, 18124-18133.](https://pubs.acs.org/doi/10.1021/acs.jpcc.8b05286)
