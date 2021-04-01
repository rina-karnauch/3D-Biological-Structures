import sys
import Bio.PDB
from Bio.PDB import MMCIFIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBList import PDBList


def alignment_method():
    # given arguments readr
    arguments = sys.argv
    pdb_id1, chain1, pdb_id2, chain2 = arguments[1], arguments[2], arguments[
        3], arguments[4]
    # parser of pdb files

    downloader = PDBList()
    downloader.retrieve_pdb_file(pdb_id1, file_format="bundle", pdir="./")
    downloader.retrieve_pdb_file(pdb_id2, file_format="bundle", pdir="./")
    parser = PDBParser(PERMISSIVE=1)

    # given structures
    structure_1 = parser.get_structure("pdb1", "pdb" + pdb_id1 + ".ent")[0]
    structure_2 = parser.get_structure("pdb2", "pdb" + pdb_id2 + ".ent")[0]

    # CA chains arrays
    CA_in_1 = get_CAs(structure_1, chain1)
    CA_in_2 = get_CAs(structure_2, chain2)

    alignment_object = Bio.PDB.Superimposer()

    # align structure_1 on top of structure_2
    alignment_object.set_atoms(CA_in_1, CA_in_2)
    alignment_object.apply(structure_2.get_atoms())
    # RMSD of first alignment
    print(alignment_object.rms)

    # align structure_2 on top of structure_1
    alignment_object.set_atoms(CA_in_2, CA_in_1)
    alignment_object.apply(structure_1.get_atoms())
    # RMSD of second alignment
    print(alignment_object.rms)

    mmCIF_output = MMCIFIO()

    mmCIF_output.set_structure(structure_1)
    mmCIF_output.save(pdb_id1 + ".cif")

    mmCIF_output.set_structure(structure_2)
    mmCIF_output.save(pdb_id2 + ".cif")


def get_CAs(structure, chain):
    # CA atoms list of structure
    CA_in_struct = []
    # go through structure chains
    for ch in structure:
        # if chain id fits the chain we go through
        if ch.id == chain:
            for residue in ch:
                # and if it includes an CA atom we append it
                if 'CA' in residue:
                    CA_in_struct.append(residue['CA'])
    # return CA atoms
    return CA_in_struct


if __name__ == '__main__':
    alignment_method()
