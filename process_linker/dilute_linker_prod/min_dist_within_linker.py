import MDAnalysis as mda
from MDAnalysis.analysis import distances
import numpy as np

u = mda.Universe("prod_dil.tpr", "dilute_linker_whole.xtc")
lig = u.select_atoms("resname CYL")

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

sels = {label: lig.select_atoms(f"name {names}") for label, names in groups.items()}
for label, sel in sels.items():
    if sel.n_atoms == 0:
        print(f"  {label}: no atoms matched — check names")

CONTACT_CUTOFF = 0.5   # nm

labels = list(groups.keys())
print(f"\n{'pair':<20}{'mean':>8}{'std':>8}{'min':>8}{'max':>8}{'%<0.5nm':>10}")

for i in range(len(labels)):
    for j in range(i + 1, len(labels)):
        a, b = labels[i], labels[j]
        if sels[a].n_atoms == 0 or sels[b].n_atoms == 0:
            continue

        d = []
        for ts in u.trajectory[::10]:
            dmat = distances.distance_array(sels[a].positions, sels[b].positions)
            d.append(dmat.min() / 10)
        d = np.array(d)

        frac = (d < CONTACT_CUTOFF).mean() * 100
        pair = f"{a}-{b}"
        print(f"{pair:<20}{d.mean():>8.2f}{d.std():>8.2f}{d.min():>8.2f}"
              f"{d.max():>8.2f}{frac:>9.1f}%")
