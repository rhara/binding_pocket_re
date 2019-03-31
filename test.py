from rdkit import Chem
import numpy as np
import mdtraj
from deepchem.feat import CircularFingerprint
from deepchem.utils.save import log
from deepchem.feat import Featurizer
import binding_pocket as bp

class BindingPocketFeaturizer(Featurizer):
    """
    Featurizes binding pockets with information about chemical environments.
    """
    residues = [
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'PYL', 'SER', 'SEC', 'THR', 'TRP',
        'TYR', 'VAL', 'ASX', 'GLX'
    ]

    n_features = len(residues)

    def featurize(self, protein_fname, pockets, pocket_atoms_map, pocket_coords, verbose=False):
        protein = mdtraj.load(protein_fname)
        n_pockets = len(pockets)
        n_residues = len(self.residues)
        res_map = {r: i for i, r in enumerate(self.residues)}
        all_features = np.zeros((n_pockets, n_residues))
        for pocket_num, (pocket, coords) in enumerate(zip(pockets, pocket_coords)):
            pocket_atoms = pocket_atoms_map[pocket]
            for ind, atom in enumerate(pocket_atoms):
                atom_name = str(protein.top.atom(atom))
                residue = atom_name[:3]
                if residue not in res_map:
                    log('Warning: Non-stardard residue in PDB file "%s"' % residue, verbose)
                    continue
                atomtype = atom_name.split('-')[1]
                all_features[pocket_num, res_map[residue]] += 1
        return all_features



def compute_binding_pocket_features(protein_fname, ligand_fname, threshold=.3):
    ligand_mol2 = ligand_fname.replace('.sdf', '.mol2')

    active_site_box, active_site_atoms, active_site_coords = bp.extract_active_site(protein_fname, ligand_fname)
    mol = Chem.MolFromMol2File(ligand_mol2, removeHs=False)

    n_ligand_features = 1024
    ligand_features = ligand_featurizer.featurize([mol])

    finder = bp.ConvexHullPocketFinder()
    pockets, pocket_atoms, pocket_coords = finder.find_pockets(protein_fname, ligand_fname)
    # pockets = [((x, x), (x, x), (x,x)), ((x, x), (x, x), (x,x)), ((x, x), (x, x), (x,x))] size 3
    # len(pocket_atoms) = 173
    # pocket_coords shape (774, 3) (312, 3) (1709, 3)
    n_pockets = len(pockets)
    n_pocket_features = BindingPocketFeaturizer.n_features #24

    features = np.zeros((n_pockets, n_pocket_features+n_ligand_features)) # shape (3, 1048)
    pocket_features = pocket_featurizer.featurize(protein_fname, pockets, pocket_atoms, pocket_coords, verbose=True)
    features[:, :n_pocket_features] = pocket_features
    features[:, n_pocket_features:] = ligand_features

    labels = np.zeros(n_pockets)
    pocket_atoms[active_site_box] = active_site_atoms
    for ind, pocket in enumerate(pockets):
        overlap = bp.compute_overlap(pocket_atoms, active_site_box, pocket)
        labels[ind] = 0 if overlap <= threshold else 1

    return features, labels

  
pocket_featurizer = BindingPocketFeaturizer()
ligand_featurizer = CircularFingerprint(size=1024)

compute_binding_pocket_features('v2015/3zzf/3zzf_protein.pdb', 'v2015/3zzf/3zzf_ligand.sdf')
