import MDAnalysis as mda
from MDAnalysis.analysis import distances
import numpy as np

u = mda.Universe("fib_only_prod.gro", "linked_nojump.xtc")

BTN = "C19 C20 C21 C22 C23 C24 C25 C26 C27 C28 C29 N3 N4 N5 O7 O8 S1"
PEG = "C11 C12 C13 C14 C15 C16 C17 C18 O3 O4 O5 O6"
hets = u.select_atoms("protein and not resname CYL and not name H*")

for res in u.select_atoms("resname CYL").residues:
    btn = res.atoms.select_atoms(f"name {BTN}")
    peg = res.atoms.select_atoms(f"name {PEG}")
    sg  = res.atoms.select_atoms("name SG")
    tip = res.atoms.select_atoms("name C28")      # far end of biotin
    print(f"\nlinker at resid {res.resid}: BTN {btn.n_atoms} atoms, PEG {peg.n_atoms} atoms")

    dmin_b, dmin_p, ext, cont = [], [], [], []
    for ts in u.trajectory[::10]:
        db = distances.distance_array(btn.positions, hets.positions)
        dp = distances.distance_array(peg.positions, hets.positions)
        dmin_b.append(db.min() / 10)                       # A -> nm
        dmin_p.append(dp.min() / 10)
        ext.append(np.linalg.norm(tip.positions[0] - sg.positions[0]) / 10)
        cont.append(int((db < 4.5).any(axis=1).sum()))     # biotin atoms in contact

    for name, v in (("BTN-HETs min dist (nm)", dmin_b),
                    ("PEG-HETs min dist (nm)", dmin_p),
                    ("SG-biotin extension (nm)", ext),
                    ("biotin atoms in contact", cont)):
        v = np.array(v)
        print(f"  {name:<26} mean {v.mean():.2f}  min {v.min():.2f}  max {v.max():.2f}")
