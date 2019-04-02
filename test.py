from sys import stderr
import re
import numpy as np
import pandas as pd

from binding_pocket import BindingPocketFeaturizer

def load_pdbbind_labels(labels_file):
    """Loads pdbbind labels as dataframe"""
    # Some complexes have labels but no PDB files. Filter these manually
    missing_pdbs = ["1d2v", "1jou", "1s8j", "1cam", "4mlt", "4o7d"]
    contents = []
    with open(labels_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            else:
                # Some of the ligand-names are of form (FMN ox). Use regex
                # to merge into form (FMN-ox)
                p = re.compile('\(([^\)\s]*) ([^\)\s]*)\)')
                line = p.sub('(\\1-\\2)', line)
                elts = line.split()
                # Filter if missing PDB files
                if elts[0] in missing_pdbs:
                    continue
                contents.append(elts)
    contents_df = pd.DataFrame(contents,
                               columns=('PDB code', 'resolution', 'release year', '-logKd/Ki', 'Kd/Ki',
                                        'ignore-this-field', 'reference', 'ligand name'))
    return contents_df


np.random.seed(123)

split = 'random'
subset = 'core'
tasks = ['active-site']
LABEL_FILES = {
    'full': 'v2015/INDEX_general_PL_data.2015',
    'refined': 'v2015/INDEX_refined_data.2015',
    'core': 'v2015/INDEX_core_data.2013',
}
labels_file = LABEL_FILES[subset]
contents_df = load_pdbbind_labels(labels_file)
ids = contents_df['PDB code'].values
y = np.array([float(val) for val in contents_df['-logKd/Ki'].values])

all_features = []
all_labels = []
all_ids = []
featurizer = BindingPocketFeaturizer()
for i, pdb_code in enumerate(ids):
    try:
        print('Processing complex %d, %s' % (i, pdb_code))
        protein_fname = 'v2015/%s/%s_protein.pdb' % (pdb_code, pdb_code)
        ligand_fname = 'v2015/%s/%s_ligand.sdf' % (pdb_code, pdb_code)
        pocket_features, ligand_features, labels = featurizer.multi_featurize(protein_fname, ligand_fname)
        n_pockets = pocket_features.shape[0]
        n_pocket_features = len(featurizer.residues)
        n_ligand_features = 1024
        features = np.zeros((n_pockets, n_pocket_features+n_ligand_features))
        features[:, :n_pocket_features] = pocket_features
        features[:, n_pocket_features:] = ligand_features
    except:
        print('Warning: failed %s' % pdb_code, file=stderr)
        continue
    print(labels)
    all_features.append(features)
    all_labels.append(labels)
    ids = np.array(['%s%d' % (pdb_code, j) for j in range(len(labels))])
    all_ids.append(ids)
    if 10 <= len(all_ids):
        break

X = np.vstack(all_features)
y = np.concatenate(all_labels)
w = np.ones_like(y)
ids = np.concatenate(all_ids)

np.savez('out.npz', X=X, y=y, w=w, ids=ids)
