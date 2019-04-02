from sys import stderr
import numpy as np
from scipy.spatial import ConvexHull
from rdkit import Chem
import mdtraj
from fingerprints import CircularFingerprint
from load_mol import load_molecule

import sys

np.set_printoptions(threshold=np.inf)

def mdtraj_count_atoms(protein):
    n_hvys = 0
    for atom in protein.top.atoms:
        if 1 < atom.element.atomic_number:
            n_hvys += 1
    return protein.n_atoms, n_hvys

def mdtraj_write_pocket_atoms(fname, protein, pocket_atoms):
    xyzs = protein.xyz[0]*10
    syms = [atom.element.symbol for atom in protein.top.atoms]

    with open(fname, 'wt') as out:
        out.write('%d\n\n' % len(pocket_atoms))

        for idx in pocket_atoms:
            sym = syms[idx]
            x, y, z = xyzs[idx]
            out.write('%3s%15.5f%15.5f%15.5f\n' % (sym, x, y, z))

def check_residues(fname, residue_list):
    protein = mdtraj.load(fname)
    res_map = {r: i for i, r in enumerate(residue_list)}
    atomic_hist = [0]*(len(residue_list)+1)
    coords = protein.xyz
    for atom in protein.top.atoms:
        el = atom.element.symbol
        resname = atom.residue.name
        if resname in residue_list:
            atomic_hist[res_map[resname]] += 1
        else:
            atomic_hist[-1] += 1
    return atomic_hist

class BindingPocketFeaturizer:
    """
    Featurizes binding pockets with information about chemical environments.
    """
    residues = [
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'PYL', 'SER', 'SEC', 'THR', 'TRP',
        'TYR', 'VAL', 'ASX', 'GLX'
    ]
    n_features = len(residues)

    def featurize(self, protein_fname, pockets, pocket_atoms_map, pocket_coords):
        protein = mdtraj.load(protein_fname)

        n_atoms, n_hvys = mdtraj_count_atoms(protein)
        print('protein: %s n_hvys/n_atoms=%d/%d' % (protein_fname, n_hvys, protein.n_atoms))

        n_pockets = len(pockets)
        n_residues = len(self.residues)
        res_map = {r: i for i, r in enumerate(self.residues)}
        all_features = np.zeros((n_pockets, n_residues))
        protein_coords = protein.xyz[0]
        count = 1
        for pocket_num, (pocket, coords) in enumerate(zip(pockets, pocket_coords)):
            pocket_atoms = pocket_atoms_map[pocket]

            mdtraj_write_pocket_atoms('out%d.xyz' % count, protein, pocket_atoms)
            count += 1

            print('%d:  pocket_atoms: %d' % (pocket_num, len(pocket_atoms)))

            for atom in pocket_atoms:
                atom_name = str(protein.top.atom(atom))
                residue = atom_name[:3]
                if residue not in res_map:
                    # print('Warning: Non-stardard residue in PDB file "%s"' % residue, file=stderr)
                    continue
                all_features[pocket_num, res_map[residue]] += 1

        return all_features

    def multi_featurize(self, protein_fname, ligand_fname, threshold=.3):
        active_site_box, active_site_atoms, active_site_coords = extract_active_site(protein_fname, ligand_fname)
        print('active_site_box:', active_site_box)

        finder = ConvexHullPocketFinder()
        pockets, pocket_atoms, pocket_coords = finder.find_pockets(protein_fname, ligand_fname)
        n_pockets = len(pockets)
        n_pocket_features = BindingPocketFeaturizer.n_features
        pocket_features = self.featurize(protein_fname, pockets, pocket_atoms, pocket_coords)

        ligand_mol2 = ligand_fname.replace('.sdf', '.mol2')
        mol = Chem.MolFromMol2File(ligand_mol2, removeHs=False)
        # for mol in Chem.rdmolfiles.SDMolSupplier(ligand_fname, removeHs=False):
        #     break
        n_ligand_features = 1024
        ligand_featurizer = CircularFingerprint(size=n_ligand_features)
        ligand_features = ligand_featurizer.featurize(mol)

        labels = np.zeros(n_pockets)
        pocket_atoms[active_site_box] = active_site_atoms
        for i, pocket in enumerate(pockets):
            overlap = compute_overlap(pocket_atoms, active_site_box, pocket)
            labels[i] = 0 if overlap <= threshold else 1

        return pocket_features, ligand_features, labels


def extract_active_site(protein_file, ligand_file, cutoff=4):
    """Extracts a box for the active site."""
    pro_coords = load_molecule(protein_file, add_hydrogens=False)[0]
    lig_coords = load_molecule(ligand_file, add_hydrogens=True, calc_charges=True)[0]
    n_lig_atoms = len(lig_coords)
    n_pro_atoms = len(pro_coords)
    pocket_atoms = set()
    for i in range(n_lig_atoms):
        lig_atom = np.tile(lig_coords[i], n_pro_atoms).reshape((n_pro_atoms, 3))
        D = np.linalg.norm(lig_atom - pro_coords, axis=1) < cutoff
        pocket_atoms = pocket_atoms.union(np.where(D == True)[0].flat)
    # Should be an array of size (n_pocket_atoms, 3)
    pocket_atoms = sorted(pocket_atoms)
    n_pocket_atoms = len(pocket_atoms)
    pocket_coords = pro_coords[pocket_atoms]

    x_min, y_min, z_min = np.floor(np.amin(pocket_coords, axis=0)).astype(int)
    x_max, y_max, z_max = np.ceil(np.amax(pocket_coords, axis=0)).astype(int)
    return (((x_min, x_max), (y_min, y_max), (z_min, z_max)), pocket_atoms, pocket_coords)


def compute_overlap(mapping, box1, box2):
    """Computes overlap between the two boxes.

    Overlap is defined as % atoms of box1 in box2. Note that
    overlap is not a symmetric measurement.
    """
    atom1 = set(mapping[box1])
    atom2 = set(mapping[box2])
    return len(atom1.intersection(atom2)) / float(len(atom1))


def get_all_boxes(coords, pad=5):
    """Get all pocket boxes for protein coords.

    We pad all boxes the prescribed number of angstroms.
    """
    hull = ConvexHull(coords)
    boxes = []
    for triangle in hull.simplices:
        # coords[triangle, 0] gives the x-dimension of all triangle points
        # Take transpose to make sure rows correspond to atoms.
        points = coords[triangle,:]
        # We voxelize so all grids have integral coordinates (convenience)
        x_min, y_min, z_min = np.floor(np.amin(points, axis=0)) - pad
        x_max, y_max, z_max = np.ceil(np.amax(points, axis=0)) + pad
        boxes.append(((x_min, x_max), (y_min, y_max), (z_min, z_max)))
    return boxes


def boxes_to_atoms(atom_coords, boxes):
    """Maps each box to a list of atoms in that box.

    TODO(rbharath): This does a num_atoms x num_boxes computations. Is
    there a reasonable heuristic we can use to speed this up?
    """
    mapping = {}
    for i, box in enumerate(boxes):
        box_atoms = []
        (x_min, x_max), (y_min, y_max), (z_min, z_max) = box
        # print("Handing box %d/%d" % (i, len(boxes)), file=stderr)
        for j in range(atom_coords.shape[0]):
            atom = atom_coords[j]
            if not (x_min <= atom[0] <= x_max):
                continue
            if not (y_min <= atom[1] <= y_max):
                continue
            if not (z_min <= atom[2] <= z_max):
                continue
            box_atoms.append(j)
        mapping[box] = box_atoms
    return mapping


def merge_boxes(box1, box2):
    """Merges two boxes."""
    (x_min1, x_max1), (y_min1, y_max1), (z_min1, z_max1) = box1
    (x_min2, x_max2), (y_min2, y_max2), (z_min2, z_max2) = box2
    x_min = min(x_min1, x_min2)
    y_min = min(y_min1, y_min2)
    z_min = min(z_min1, z_min2)
    x_max = max(x_max1, x_max2)
    y_max = max(y_max1, y_max2)
    z_max = max(z_max1, z_max2)
    return ((x_min, x_max), (y_min, y_max), (z_min, z_max))


def merge_overlapping_boxes(mapping, boxes, threshold=.8):
    """Merge boxes which have an overlap greater than threshold.

    TODO(rbharath): This merge code is terribly inelegant. It's also quadratic
    in number of boxes. It feels like there ought to be an elegant divide and
    conquer approach here. Figure out later...
    """
    num_boxes = len(boxes)
    outputs = []
    for i in range(num_boxes):
        box = boxes[0]
        new_boxes = []
        new_mapping = {}
        # If overlap of box with previously generated output boxes, return
        contained = False
        for output_box in outputs:
            # Carry forward mappings
            new_mapping[output_box] = mapping[output_box]
            if compute_overlap(mapping, box, output_box) == 1:
                contained = True
        if contained:
            continue
        # We know that box has at least one atom not in outputs
        unique_box = True
        for merge_box in boxes[1:]:
            overlap = compute_overlap(mapping, box, merge_box)
            if overlap < threshold:
                new_boxes.append(merge_box)
                new_mapping[merge_box] = mapping[merge_box]
            else:
                # Current box has been merged into box further down list.
                # No need to output current box
                unique_box = False
                merged = merge_boxes(box, merge_box)
                new_boxes.append(merged)
                new_mapping[merged] = list(set(mapping[box]).union(set(mapping[merge_box])))
        if unique_box:
            outputs.append(box)
            new_mapping[box] = mapping[box]
        boxes = new_boxes
        mapping = new_mapping
    return outputs, mapping


class BindingPocketFinder(object):
    """Abstract superclass for binding pocket detectors"""
    def find_pockets(self, protein_file, ligand_file):
      """Finds potential binding pockets in proteins."""
      raise NotImplementedError


class ConvexHullPocketFinder(BindingPocketFinder):
    """Implementation that uses convex hull of protein to find pockets.

    Based on https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4112621/pdf/1472-6807-14-18.pdf
    """

    def __init__(self, pad=5):
        self.pad = pad

    def find_all_pockets(self, protein_file):
        """Find list of binding pockets on protein."""
        # protein_coords is (N, 3) tensor
        coords = load_molecule(protein_file)[0]
        return get_all_boxes(coords, self.pad)

    def find_pockets(self, protein_file, ligand_file):
        """Find list of suitable binding pockets on protein."""
        protein_coords = load_molecule(protein_file, add_hydrogens=False, calc_charges=False)[0]
        ligand_coords = load_molecule(ligand_file, add_hydrogens=False, calc_charges=False)[0]
        boxes = get_all_boxes(protein_coords, self.pad)
        mapping = boxes_to_atoms(protein_coords, boxes)
        pockets, pocket_atoms_map = merge_overlapping_boxes(mapping, boxes)
        pocket_coords = []
        for pocket in pockets:
            atoms = pocket_atoms_map[pocket]
            coords = np.zeros((len(atoms), 3))
            for ind, atom in enumerate(atoms):
                coords[ind] = protein_coords[atom]
            pocket_coords.append(coords)
        return pockets, pocket_atoms_map, pocket_coords
