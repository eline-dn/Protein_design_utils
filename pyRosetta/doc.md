# A bunch of useful links to tutorials and documentation:

## ligands with pyrosetta:
https://www.pyrosetta.org/documentation/obtaining-and-preparing-ligand-pdb-files
https://docs.rosettacommons.org/docs/latest/rosetta_basics/preparation/preparing-ligands
https://docs.rosettacommons.org/demos/latest/tutorials/prepare_ligand/prepare_ligand_tutorial

Quick overview: 

. If starting from an RCSB crystal structure, you can use PyMOL's "Save Molecule" feature to produce an .mdl file of a ligand (the file extension appears as ".mol"). 
Converting to Params Files
An additional script (and other necessary scripts) are provided for converting an .mdl file to a .params file (required for the PyRosetta database) and .pdb files. Execute this script from the commandline providing the .mdl file as the first argument and the ResidueType name as option "-n". For the .mdl file provided with this script, the example commandline call would be:

>python molfile_to_params.py ATP.mdl -n ATP

Preparing Ligand PDB Files
Now that the ResidueType is defined, the PDB file for ligand interface prediction can be made. If the PDB file already has the ligand present, ensure that its ResisueType column (PDB file format) is set to the ligand ResidueType ("ATP" for the example case). It is common practice to rename the chain to "X". 

Inside the relevant script or interpreter, create a non-standard ResidueTypeSet using the method generate_nonstandard_residue_set and use use pose_from_pdb to load data into to a pose object. The method pose_from_pdb is overloaded such that it can accept a Pose (poses), a ResidueTypeSet (residue_set), and a string (filename) and load into the poses the data in the PDB file filename using residue_set to define any unknown residues.

Also possible to convert smiles into param file with this https://pypi.org/project/rdkit-to-params/ 
