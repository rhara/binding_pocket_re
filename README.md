# binding_pocket_re

This is to do reverse engineering the original code.

This code does not depend on deepchem package. Code is quite much improved by using numpy arrays.

Part of the publication: "Massively Multitask Networks for Drug Discovery," B Ramsundar, S Kearnes, P Riley, D Webster, D Konerding, V Pande arXiv preprint arXiv:1502.02072 (2015)

https://arxiv.org/abs/1502.02072

### Requisites:

- numpy (no specific version)
- rdkit (no specific version)
- mdtraj (no specific version)
- pdbfixer (no specific version)

```
python pocket_descriptor.py -l <ligand.sdf> -p <protein.pdb>
```
