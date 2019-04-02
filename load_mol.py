from sys import stderr
from io import StringIO
import os
import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem


class MoleculeLoadException(Exception):
    def __init__(self, *args, **kwargs):
        Exception.__init__(*args, **kwargs)


def get_xyz_from_mol(mol):
    xyz = np.zeros((mol.GetNumAtoms(), 3))
    conf = mol.GetConformer()
    for i in range(conf.GetNumAtoms()):
        position = conf.GetAtomPosition(i)
        xyz[i, 0] = position.x
        xyz[i, 1] = position.y
        xyz[i, 2] = position.z
    return (xyz)


def add_hydrogens_to_mol(mol):
    """
    Add hydrogens to a molecule object
    TODO (LESWING) see if there are more flags to add here for default
    :param mol: Rdkit Mol
    :return: Rdkit Mol
    """
    molecule_file = None
    try:
        pdbblock = Chem.MolToPDBBlock(mol)
        pdb_stringio = StringIO()
        pdb_stringio.write(pdbblock)
        pdb_stringio.seek(0)
        import pdbfixer
        fixer = pdbfixer.PDBFixer(pdbfile=pdb_stringio)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.4)

        hydrogenated_io = StringIO()
        import simtk
        simtk.openmm.app.PDBFile.writeFile(fixer.topology, fixer.positions, hydrogenated_io)
        hydrogenated_io.seek(0)
        return Chem.MolFromPDBBlock(hydrogenated_io.read(), sanitize=False, removeHs=False)
    except ValueError as e:
        print('Unable to add hydrogens %s' % e, file=stderr)
        raise MoleculeLoadException(e)
    finally:
        try:
            os.remove(molecule_file)
        except (OSError, TypeError):
            pass


def compute_charges(mol):
    try:
        AllChem.ComputeGasteigerCharges(mol)
    except Exception as e:
        print('Unable to compute charges for mol', file=stderr)
        raise MoleculeLoadException(e)
    return mol


def load_molecule(molecule_file, add_hydrogens=True, calc_charges=True, sanitize=False):
    if '.mol2' in molecule_file:
        my_mol = Chem.MolFromMol2File(molecule_file, sanitize=False, removeHs=False)
    elif '.sdf' in molecule_file:
        suppl = Chem.SDMolSupplier(str(molecule_file), sanitize=False)
        my_mol = suppl[0]
    elif '.pdb' in molecule_file:
        my_mol = Chem.MolFromPDBFile(str(molecule_file), sanitize=False, removeHs=False)
    else:
        raise ValueError('Unrecognized file type')

    if my_mol is None:
        raise ValueError('Unable to read non None Molecule Object')

    if add_hydrogens or calc_charges:
        my_mol = add_hydrogens_to_mol(my_mol)
    if sanitize:
        Chem.SanitizeMol(my_mol)
    if calc_charges:
        compute_charges(my_mol)

    xyz = get_xyz_from_mol(my_mol)

    return xyz, my_mol
