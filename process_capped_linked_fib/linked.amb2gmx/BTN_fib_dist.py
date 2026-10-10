import MDAnalysis as mda
from MDAnalysis.analysis import distances
import numpy as np

u = mda.Universe("fib_only_prod.gro", "linked_nojump.xtc")

# select CYL residues, then take the LAST N atoms of each (biotin is at the
# far end of the chain from the M-F junction, so it's consistently the
# terminal segment of the residue regardless of reordering — as long as
# residue membership is intact, which .tpr/.xtc preserves even if the .gro
# print-order looks scrambled)
cyl = u.select_atoms("resname CYL")
protein = u.select_atoms("protein and not resname CYL")

# split each CYL residue and take its last ~15 heavy atoms as "biotin"
biotin_atoms = []
for res in cyl.residues:
    heavy = res.atoms.select_atoms("not name H*")
    biotin_atoms.append(heavy[-13:])   # adjust count to match biotin's heavy-atom size
biotin = sum(biotin_atoms)

dists = []
for ts in u.trajectory[::10]:   # every 10th frame, adjust as needed
    d = mda.analysis.distances.distance_array(biotin.positions, protein.positions)
    dists.append(d.min())

dists = np.array(dists)
print("mean:", dists.mean(), "min:", dists.min(), "max:", dists.max())
