import MDAnalysis as mda
from MDAnalysis.analysis import distances
import numpy as np

u = mda.Universe("fib_only_prod.gro", "linked_nojump.xtc")
hets = u.select_atoms("protein and not resname CYL and not name H*")

groups = {
    "SUC":    "C4 C5 C6 O1 N1 C30 O9",
    "EO1":    "C10 C11 O3",
    "amide1": "N2 C9 O2",
    "EO2":    "C12 C13 O4",
    "EO3":    "C14 C15 O5",
    "EO4":    "C16 C17 O6 C18 C19",
    "amide2": "N3 C20 O7",
    "BTN":    "C21 C22 C23 C24 C25 S1 C26 C27 N4 C28 O8 N5 C29",
}

CONTACT_CUTOFF = 0.5   # nm — typical vdW contact threshold

for res in u.select_atoms("resname CYL").residues:
    print(f"\n=== linker at resid {res.resid} ===")
    for label, names in groups.items():
        sel = res.atoms.select_atoms(f"name {names}")
        if sel.n_atoms == 0:
            print(f"  {label}: no atoms matched — check names")
            continue

        d = []
        for ts in u.trajectory[::10]:
            dmat = distances.distance_array(sel.positions, hets.positions)
            d.append(dmat.min() / 10)   # A -> nm
        d = np.array(d)

        n = len(d)
        frac_contact = (d < CONTACT_CUTOFF).mean() * 100
        # split into thirds to see if behaviour changes over the run
        thirds = np.array_split(d, 3)
        third_means = [t.mean() for t in thirds]

        print(f"  {label:<9} mean {d.mean():.2f}  std {d.std():.2f}  "
              f"min {d.min():.2f}  max {d.max():.2f}  "
              f"%<{CONTACT_CUTOFF}nm {frac_contact:5.1f}%  "
              f"thirds {[round(x,2) for x in third_means]}")
