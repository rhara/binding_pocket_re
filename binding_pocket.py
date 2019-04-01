from sys import stderr
import numpy as np
from scipy.spatial import ConvexHull
from rdkit import Chem
import mdtraj
from deepchem.feat.fingerprints import CircularFingerprint
from deepchem.utils import rdkit_util


def mdtraj_count_atoms(trj_protein):
    n_hvys = 0
    for atom in trj_protein.top.atoms:
        if 1 < atom.element.atomic_number:
            n_hvys += 1
    return trj_protein.n_atoms, n_hvys

def write_pocket_atoms(fname, rdk_protein, protein_coords, pocket_atoms):
    out = open(fname, 'wt')
    out.write('%d\n\n' % len(pocket_atoms))

    for idx in pocket_atoms:
        atom = rdk_protein.GetAtomWithIdx(idx)
        sym = atom.GetSymbol()
        x, y, z = protein_coords[idx]
        out.write('%3s%15.5f%15.5f%15.5f\n' % (sym, x, y, z))
    out.close()

def check_residues(protein_fname, residue_list):
    rdk_protein = Chem.rdmolfiles.MolFromPDBFile(protein_fname, removeHs=False)
    res_map = {r: i for i, r in enumerate(residue_list)}
    atomic_hist = [0]*(len(residue_list)+1)
    for atom in rdk_protein.GetAtoms():
        idx = atom.GetIdx()
        sym = atom.GetSymbol()
        resinfo = atom.GetPDBResidueInfo()
        chainID = resinfo.GetChainId()
        resno = resinfo.GetResidueNumber()
        resname = resinfo.GetResidueName()
        # print('%d %s %s%d%s' % (idx, sym, resname, resno, chainID))
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
        rdk_protein = Chem.rdmolfiles.MolFromPDBFile(protein_fname, removeHs=False)
        protein_coords = rdk_protein.GetConformer(0).GetPositions()
        count = 1
        for pocket_num, (pocket, coords) in enumerate(zip(pockets, pocket_coords)):
            pocket_atoms = pocket_atoms_map[pocket]

            write_pocket_atoms('out%d.xyz' % count, rdk_protein, protein_coords, pocket_atoms)
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
        n_ligand_features = 1024
        ligand_featurizer = CircularFingerprint(size=n_ligand_features)
        ligand_features = ligand_featurizer.featurize([mol])

        labels = np.zeros(n_pockets)
        pocket_atoms[active_site_box] = active_site_atoms
        for i, pocket in enumerate(pockets):
            overlap = compute_overlap(pocket_atoms, active_site_box, pocket)
            labels[i] = 0 if overlap <= threshold else 1

        return pocket_features, ligand_features, labels


def extract_active_site(protein_file, ligand_file, cutoff=4):
    """Extracts a box for the active site."""
    protein_coords = rdkit_util.load_molecule(protein_file, add_hydrogens=False)[0]
    ligand_coords = rdkit_util.load_molecule(ligand_file, add_hydrogens=True, calc_charges=True)[0]
    num_ligand_atoms = len(ligand_coords)
    num_protein_atoms = len(protein_coords)
    pocket_inds = []
    pocket_atoms = set([])
    for lig_atom_ind in range(num_ligand_atoms):
        lig_atom = ligand_coords[lig_atom_ind]
        for protein_atom_ind in range(num_protein_atoms):
            protein_atom = protein_coords[protein_atom_ind]
            if np.linalg.norm(lig_atom - protein_atom) < cutoff:
                if protein_atom_ind not in pocket_atoms:
                    pocket_atoms = pocket_atoms.union(set([protein_atom_ind]))
    # Should be an array of size (n_pocket_atoms, 3)
    pocket_atoms = list(pocket_atoms)
    n_pocket_atoms = len(pocket_atoms)
    pocket_coords = np.zeros((n_pocket_atoms, 3))
    for ind, pocket_ind in enumerate(pocket_atoms):
        pocket_coords[ind] = protein_coords[pocket_ind]

    x_min = int(np.floor(np.amin(pocket_coords[:, 0])))
    x_max = int(np.ceil(np.amax(pocket_coords[:, 0])))
    y_min = int(np.floor(np.amin(pocket_coords[:, 1])))
    y_max = int(np.ceil(np.amax(pocket_coords[:, 1])))
    z_min = int(np.floor(np.amin(pocket_coords[:, 2])))
    z_max = int(np.ceil(np.amax(pocket_coords[:, 2])))
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

    TODO(rbharath): It looks like this may perhaps be non-deterministic?
    """
    hull = ConvexHull(coords)
    boxes = []
    for triangle in hull.simplices:
        # coords[triangle, 0] gives the x-dimension of all triangle points
        # Take transpose to make sure rows correspond to atoms.
        points = np.array([coords[triangle, 0], coords[triangle, 1], coords[triangle, 2]]).T
        # We voxelize so all grids have integral coordinates (convenience)
        x_min, x_max = np.amin(points[:, 0]), np.amax(points[:, 0])
        x_min, x_max = int(np.floor(x_min)) - pad, int(np.ceil(x_max)) + pad
        y_min, y_max = np.amin(points[:, 1]), np.amax(points[:, 1])
        y_min, y_max = int(np.floor(y_min)) - pad, int(np.ceil(y_max)) + pad
        z_min, z_max = np.amin(points[:, 2]), np.amax(points[:, 2])
        z_min, z_max = int(np.floor(z_min)) - pad, int(np.ceil(z_max)) + pad
        boxes.append(((x_min, x_max), (y_min, y_max), (z_min, z_max)))
    return boxes


def boxes_to_atoms(atom_coords, boxes):
    """Maps each box to a list of atoms in that box.

    TODO(rbharath): This does a num_atoms x num_boxes computations. Is
    there a reasonable heuristic we can use to speed this up?
    """
    mapping = {}
    for box_ind, box in enumerate(boxes):
        box_atoms = []
        (x_min, x_max), (y_min, y_max), (z_min, z_max) = box
        # print("Handing box %d/%d" % (box_ind, len(boxes)), file=stderr)
        for atom_ind in range(len(atom_coords)):
            atom = atom_coords[atom_ind]
            x_cont = x_min <= atom[0] and atom[0] <= x_max
            y_cont = y_min <= atom[1] and atom[1] <= y_max
            z_cont = z_min <= atom[2] and atom[2] <= z_max
            if x_cont and y_cont and z_cont:
                box_atoms.append(atom_ind)
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
        coords = rdkit_util.load_molecule(protein_file)[0]
        return get_all_boxes(coords, self.pad)

    def find_pockets(self, protein_file, ligand_file):
        """Find list of suitable binding pockets on protein."""
        protein_coords = rdkit_util.load_molecule(protein_file, add_hydrogens=False, calc_charges=False)[0]
        ligand_coords = rdkit_util.load_molecule(ligand_file, add_hydrogens=False, calc_charges=False)[0]
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
