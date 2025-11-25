# TO DO: sort out the usued packages
from Bio.PDB import PDBParser
import os
import pandas as pd
import colabdesign
import sys
import numpy as np
import jax
from colabdesign.shared.utils import copy_dict
from Bio.PDB import PDBParser, PDBIO

from collections import defaultdict
from scipy.spatial import cKDTree
from Bio import BiopythonWarning
from Bio.PDB import DSSP, Selection, Polypeptide, Select, Chain, Superimposer
from Bio.PDB.SASA import ShrakeRupley
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.PDB.Selection import unfold_entities


from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB import PDBParser, PDBIO, Model, Chain, Structure
from Bio.PDB import StructureBuilder
from Bio.PDB.Polypeptide import is_aa # Assuming is_aa is needed and available

def chain_length(pdb_path, chain_id=str):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("x", pdb_path)
    model = structure[0]             # first model
    chain = model[chain_id]               # chain of interest
    # count only standard residues
    residues = [res for res in chain.get_residues() if res.id[0] == " "]
    return len(residues)


def chain_seq(pdb_path, chain_id=str,make_str=False): 
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("x", pdb_path)
    model = structure[0]             # first model
    chain = model[chain_id]               # chain
    # count only standard residues
    residues = [res for res in chain.get_residues() if res.id[0] == " "] #res.id[0] == " " filters out HETATM residues.
    if make_str==False:
	return residues
    else:
	return "".join(residues)


def trim_pdb(input_pdb_path, output_pdb_path, trim_length=26):
    """
    Trims the first N amino acids from a PDB file. 

    Args:
        input_pdb_path (str): Path to the input PDB file.
        output_pdb_path (str): Path to save the trimmed PDB file.
        trim_length (int): The number of amino acids to trim from the beginning.
    """
    parser = PDBParser()
    structure = parser.get_structure("protein", input_pdb_path)

    for model in structure:
        for chain in model:
            # Get all residues in the chain
            residues = list(chain.get_residues())
            # Keep residues from index trim_length onwards
            for i, residue in enumerate(residues):
                if i < trim_length:
                    chain.detach_child(residue.get_id())

    io = PDBIO(use_model_flag=1)
    io.set_structure(structure)
    io.save(output_pdb_path)

# Example usage:
# Assuming you have a PDB file named 'input.pdb' in the current directory
# trim_pdb('input.pdb', 'trimmed_output.pdb')
# print("Trimmed PDB file saved as 'trimmed_output.pdb'")



def _copy_structure_with_only_chain(structure, chain_id):
    """
	From BindCraft's Biopyhton_utils (https://github.com/martinpacesa/BindCraft) 
	Return a new Structure containing only model 0 and a deep copy of chain `chain_id`."""
    # Build a tiny structure hierarchy: Structure -> Model(0) -> Chain(chain_id) -> Residues/Atoms

    sb = StructureBuilder.StructureBuilder()
    sb.init_structure("single")
    sb.init_model(1)
    sb.init_chain(chain_id)
    # Set segment ID, padded to 4 characters
    sb.init_seg(chain_id.ljust(4))    
    model0 = structure[0]
    if chain_id not in [c.id for c in model0.get_chains()]:
        raise ValueError(f"Chain '{chain_id}' not found.")
    chain = model0[chain_id]
    for res in chain:
        # Keep only amino-acid residues
        # Assuming is_aa is defined elsewhere and available
        if not is_aa(res, standard=False):
            continue
        hetflag, resseq, icode = res.id
        sb.init_residue(res.resname, hetflag, resseq, icode)

        for atom in res:
            sb.init_atom(atom.name, atom.coord, atom.bfactor, atom.occupancy,
                         atom.altloc, atom.fullname, element=atom.element)
    return sb.get_structure()


def extract_chain(input_pdb_path: str, output_pdb_path: str, chain_id: str):
    """
	
    Extracts a specific chain from a PDB file using _copy_structure_with_only_chain
    and saves it to a new PDB file with explicit MODEL/ENDMDL records.

    Args:
        input_pdb_path (str): Path to the input PDB file (complex).
        output_pdb_path (str): Path to save the extracted chain PDB file.
        chain_id (str): The identifier of the chain to extract (e.g., "A", "B").
    """
    parser = PDBParser()
    structure = parser.get_structure("protein", input_pdb_path)
    io = PDBIO(use_model_flag=1)

    # Use the helper function to get a new structure with only the desired chain
    new_structure = _copy_structure_with_only_chain(structure, chain_id)

    # --- Debug Print Statements ---
    print(f"--- Debug: Saving structure for {output_pdb_path} ---")
    print(f"Number of models in structure to save: {len(new_structure)}")
    for i, model in enumerate(new_structure):
        print(f"  Model {i}: Number of chains = {len(model)}")
        for j, chain in enumerate(model):
            print(f"    Chain {chain.id}: Number of residues = {len(chain)}")
            # Optional: print a few residue IDs and segids to confirm content
            print(f"      First few residues (ID, SegID): {[(r.id, r.segid) for r in list(chain.get_residues())[:5]]}")
    print("----------------------------------------------------")
    # --- End Debug Print Statements ---

    # Save the new structure, explicitly writing model records
    io.set_structure(new_structure)
    io.save(output_pdb_path)




# renumlbering a pdb and its chains:

import sys, os
import math
from Bio.PDB import *
from Bio import PDB

### Define Functions ###

def renumber_pdb(pdb):

    '''
    Switches chains A and B and renumbers the pdb before saving it in a given directory.
    '''
    
    # Create a PDB parser object
    parser = PDB.PDBParser()
    pdb_io = PDB.PDBIO()

    # Parse PDB file
    pdb_id = pdb.split('.')[0].split('/')[-1]
    structure = parser.get_structure(pdb_id, pdb)
    
    # Get all chains
    chains = []
    for model in structure:
        for chain in model:
            chains.append(chain)
    
    # Switch chain ID
    chains[0].id = 'A'

    for j, residue in enumerate(chains[0].get_residues()):
        id = list(residue.id)
        id[1] = j+1
        residue.id = tuple(id)

    # Add chains to the 'new' structure
    for model in structure:
        # Detach old chains
        model.detach_child('A')
        # Add new chain
        model.add(chains[0])

    # Save the merged chains
    pdb_io.set_structure(structure)
    pdb_io.save(pdb.split('.')[0]+'_renumbered.pdb')

    return


import os

## not functionnal:
def addB_pdb(pdb):

    '''
    add chain B on chain A's N-terminus as chain A, and renumber chain A (fuse two chains and renumber, starting from chan B ) 
    '''
    
    # Create a PDB parser object
    parser = PDB.PDBParser()
    pdb_io = PDB.PDBIO()

    # Parse PDB file
    pdb_id = os.path.basename(pdb)
    structure = parser.get_structure(pdb_id, pdb)
    
    # Get all chains
    chains = []
    for model in structure:
        for chain in model:
            chains.append(chain)

	
    # Switch chain ID
    #chains[0].id = 'A'
	# renumber chain B from 1
	len_B=0
	for j, residue in enumerate(chains[1].get_residues()):
		len_B+=1
		residue.id=j+1
		

    for j, residue in enumerate(chains[0].get_residues()):
        id = list(residue.id)
        id[1] = j+1
        residue.id = tuple(id)

    # Add chains to the 'new' structure
    for model in structure:
		for chain in model
        # Detach old chains
        model.detach_child('A')
		model.detach_child('B')
        # Add new chain
        model.add(chains[1])
		model.add(chains[0])

    # Save the merged chains
    pdb_io.set_structure(structure)
    pdb_io.save(pdb.split('.')[0]+'_renumbered.pdb')

    return




def merge_chains_with_structurebuilder(input_pdb):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("x", input_pdb)
    model = structure[0]

    if "A" not in model or "B" not in model:
        raise ValueError("The PDB must contain both chain A and chain B")

    chainA = model["A"]
    chainB = model["B"]

    # Initialize structure builder
    sb = StructureBuilder()
    sb.init_structure("merged")
    sb.init_model(0)

    # create one merged chain 'A'
    sb.init_chain("A")

    # --- IMPORTANT: init_seg sets self.segid used by init_residue internally ---
    sb.init_seg("    ")  # 4-space segid, could be any 1-4 char string

    residue_counter = 1

    def add_chain(chain, residue_counter):
        for res in chain:
			if not is_aa(res, standard=False):
				continue
            # res.id is (hetflag, resseq, icode)
            hetflag, _, icode = res.id

            # Use documented signature: init_residue(resname, field, resseq, icode)
            sb.init_residue(res.get_resname(), hetflag, residue_counter, icode or " ")
            residue_counter += 1

            for atom in res:
                # keep only primary altlocs
                alt = atom.get_altloc()
                if alt not in (" ", "A"):
                    continue

                # StructureBuilder.init_atom(name, coord, bfactor, occupancy,
                #                            altloc, fullname, serial_number, element)
                sb.init_atom(
                    atom.get_name(),
                    atom.get_coord(),
                    atom.get_bfactor(),
                    atom.get_occupancy(),
                    atom.get_altloc() if atom.get_altloc() != " " else " ",
                    atom.get_fullname(),
                    None,              # let builder assign serial number
                    atom.element
                )
        return residue_counter

    # Append chain B then chain A
    residue_counter = add_chain(chainB, residue_counter)
    residue_counter = add_chain(chainA, residue_counter)

    # Retrieve built structure and write to an in-memory buffer
    new_structure = sb.get_structure()
    io = PDBIO()
    io.set_structure(new_structure)
    buffer = StringIO()
    io.save(buffer)
    return buffer.getvalue()
