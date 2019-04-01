# binding_pocket_re

This is to do reverse engineering the original code.

Code is quite much improved by using numpy arrays.

The goal is close to get rid of deepchem, but for now, deepchem.feat.fingerprints.CircularFingerprint and deepchem.utils.rdkit_util.load_molecule remain in the code.

```
python pocket_descriptor.py -l <ligand.sdf> -p <protein.pdb>
```
