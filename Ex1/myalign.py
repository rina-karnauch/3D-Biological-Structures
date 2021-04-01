import sys
import Bio.PDB
from Bio.PDB import MMCIFIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBList import PDBList


def alignment_method():
    """
    method to align 2 structures from pdb together
    one time 1 on top of 2, and then 2 on top of 1
    :return: none
    """
    # given arguments readr
    arguments = sys.argv
    pdb_id1, chain1, pdb_id2, chain2 = arguments[1], arguments[2], arguments[
        3], arguments[4]
    # parser of pdb files

    # download and get structure methods
    download_structures(pdb_id1, pdb_id2)
    structure_1, structure_2 = get_structures(pdb_id1, pdb_id1)

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

    # save alignment on top of structure1 as a sif file
    mmCIF_output.set_structure(structure_1)
    mmCIF_output.save(pdb_id1 + ".cif")

    # save alignment on top of structure2 as a sif file
    mmCIF_output.set_structure(structure_2)
    mmCIF_output.save(pdb_id2 + ".cif")


def download_structures(code1: str, code2: str) -> None:
    """
    method to download structures from pdb sitr
    :param code1: structure 1 to download
    :param code2: structure 2 to download
    :return: none
    """
    # download object to download from pdb net and saving in current folder
    downloader = PDBList()
    downloader.retrieve_pdb_file(code1, file_format="bundle", pdir="./")
    downloader.retrieve_pdb_file(code2, file_format="bundle", pdir="./")


def get_structures(pdb_id1: str, pdb_id2: str):
    """
    method to get structures of pdb files which we downloaded
    :param pdb_id1: code of structure1
    :param pdb_id2: code of structure2
    :return: structure object of 1 and 2
    """
    # parsing object
    parser = PDBParser(PERMISSIVE=1)
    # given structures
    structure_1 = parser.get_structure("pdb1", "pdb" + pdb_id1 + ".ent")[0]
    structure_2 = parser.get_structure("pdb2", "pdb" + pdb_id2 + ".ent")[0]
    return structure_1, structure_2


def get_CAs(structure, chain):
    """
    method to get CA chains of current given structure
    :param structure: pdb structure
    :param chain: the chain we want to read the CA's from
    :return: array of CAs in structure
    """
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
