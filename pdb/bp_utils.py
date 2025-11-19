from Bio.PDB import PDBParser

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
