import argparse
from Bio import PDB

from consts import *

def filter_atoms(structure):
    cg_atoms = set()

    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname()
                if resname in PURINES:
                    for atom in residue:
                        if atom.get_name() in BACKBONE_ATOMS:
                            cg_atoms.add(atom)
                        elif atom.get_name() in PURINE_ATOMS:
                            cg_atoms.add(atom)
                elif resname in PYRIMIDINES:
                    for atom in residue:
                        if atom.get_name() in BACKBONE_ATOMS:
                            cg_atoms.add(atom)
                        elif atom.get_name() in PYRIMIDINE_ATOMS:
                            cg_atoms.add(atom)

    return cg_atoms


def save_cg_structure(structure, atoms, output_fname):
    io = PDB.PDBIO()
    io.set_structure(structure)

    class Select(PDB.Select):
        def accept_atom(self, atom):
            return atom in atoms

    io.save(output_fname, select=Select())
    

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, help="PDB input file")
    parser.add_argument('--output', required=True, help="Coarse-grained PDB output file")

    args = parser.parse_args()

    pdb_parser = PDB.PDBParser(QUIET=True)
    structure = pdb_parser.get_structure('RNA', args.input)

    cg_atoms = filter_atoms(structure)

    save_cg_structure(structure, cg_atoms, args.output)

if __name__ == "__main__":
    main()
