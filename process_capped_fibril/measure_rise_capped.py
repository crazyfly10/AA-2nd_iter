#!/usr/bin/env python3
"""
measure_rise_capped.py — inter-monomer rise for a fibril with a CAP on one end.

TWO CHANGES FROM measure_rise.py

  1. ALL consecutive pairs are measured, not just interior ones.
     The original hardcoded INTERIOR = range(3,9) to exclude fraying termini,
     because end effects propagate ~4 monomers deep and would bias the bulk
     value used for k_theta. That is still the right choice for deriving CG
     constants — but it is the wrong choice when the question is "how does
     the cap sit relative to chain K", which requires the terminal pairs.

     Interior pairs are still flagged separately in the output so you can
     read the bulk value without the ends contaminating it.

  2. Chains may have DIFFERENT residue numbering.
     The original indexed arithmetically (chain_index * RES_PER_CHAIN + offset),
     which assumes every chain is identical. A cap with 130 residues numbered
     1-130 breaks that completely. This version selects by chain ID with a
     per-chain numbering offset instead.

Requires: pip install MDAnalysis

Usage:
    python3 measure_rise_capped.py hets11_capped.pdb
    python3 measure_rise_capped.py prod.tpr nojump.xtc
"""

import sys
import numpy as np

try:
    import MDAnalysis as mda
except ImportError:
    sys.exit("pip install MDAnalysis")

# --------------------------------------------------------------- config
# The four beta-strands of the solenoid core, in HET-s AUTHOR numbering.
# Fitting on the core only matters: ~30 of 73 residues are disordered and a
# whole-chain centroid is dominated by them.
CORE_HETS = [(226, 234), (237, 245), (262, 270), (273, 282)]

# Per-chain residue-numbering offset relative to HET-s author numbering.
# Fibril chains already use 217-289, so offset 0.
# The cap is numbered 1-130 where cap residue n == HET-s residue n + 216,
# so to express the core in CAP numbering we subtract 216.
CHAIN_OFFSET = {
    "A": 0, "B": 0, "C": 0, "D": 0, "E": 0, "F": 0,
    "G": 0, "H": 0, "I": 0, "J": 0, "K": 0,
    "L": -216,          # <-- the cap
}

# Chains in axial order, first to last along the fibril.
CHAIN_ORDER = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"]

# Which pairs are bulk-representative. End effects reach ~4 monomers deep,
# so in a 12-unit stack only the middle pairs give the value that should
# feed b = 4 x rise and therefore k_theta.
INTERIOR_PAIRS = [("F", "G"), ("G", "H"), ("H", "I"), ("I", "J")]

STRIDE = 1


def core_selection(u, chain, offset):
    """core CA atoms of one chain, with its own numbering offset applied"""
    parts = []
    for lo, hi in CORE_HETS:
        parts.append(f"resid {lo + offset}:{hi + offset}")
    sel = f"(segid {chain} or chainID {chain}) and ({' or '.join(parts)}) and name CA"
    try:
        g = u.select_atoms(sel)
    except Exception:
        # some topologies expose chains only as segid
        sel = f"segid {chain} and ({' or '.join(parts)}) and name CA"
        g = u.select_atoms(sel)
    return g


# ------------------------------------------------------------------ load
args = sys.argv[1:]
if not args:
    sys.exit("usage: measure_rise_capped.py structure [trajectory]")

u = mda.Universe(*args)
n_frames = len(u.trajectory)
print(f"{n_frames} frame(s), {len(u.atoms)} atoms\n")

# --------------------------------------------------- build core selections
cores = {}
print("core atom counts per chain  (all should match — 37 for this core)")
for ch in CHAIN_ORDER:
    off = CHAIN_OFFSET.get(ch, 0)
    g = core_selection(u, ch, off)
    cores[ch] = g
    flag = "" if g.n_atoms == 37 else "   <-- MISMATCH"
    print(f"  {ch}  offset {off:>5}  {g.n_atoms:3d} atoms{flag}")

missing = [c for c in CHAIN_ORDER if cores[c].n_atoms == 0]
if missing:
    sys.exit(f"\nno atoms selected for chains {missing} — "
             f"check CHAIN_OFFSET and whether chains are segid or chainID")

# ---------------------------------------------------------------- measure
pairs = list(zip(CHAIN_ORDER, CHAIN_ORDER[1:]))
series = {p: [] for p in pairs}

for ts in u.trajectory[::STRIDE]:
    for a, b in pairs:
        d = cores[b].center_of_mass() - cores[a].center_of_mass()
        series[(a, b)].append(np.linalg.norm(d))

n_used = len(series[pairs[0]])
print(f"\nanalysed {n_used} frame(s)\n")

# ----------------------------------------------------------------- report
print(f"{'pair':<8}{'mean (A)':>11}{'std':>9}{'min':>9}{'max':>9}   region")
print("-" * 60)

interior_vals = []
for p in pairs:
    v = np.array(series[p])
    is_int = p in INTERIOR_PAIRS
    is_cap = "L" in p
    region = "INTERIOR" if is_int else ("cap junction" if is_cap else "end effect")
    if is_int:
        interior_vals.append(v)
    print(f"{p[0]}-{p[1]:<6}{v.mean():>11.3f}{v.std():>9.3f}"
          f"{v.min():>9.3f}{v.max():>9.3f}   {region}")

if interior_vals:
    allv = np.concatenate(interior_vals)
    print("-" * 60)
    print(f"{'INTERIOR':<8}{allv.mean():>11.3f}{allv.std():>9.3f}"
          f"{allv.min():>9.3f}{allv.max():>9.3f}")
    b = 4 * allv.mean() / 10.0
    print(f"\n  bulk rise      {allv.mean():.3f} A")
    print(f"  b = 4 x rise   {b:.3f} nm")
    print(f"  k_theta        {7500 * 2.478 / b:.0f} kJ/mol/rad^2   "
          f"(l_p = 7.5 um)")

cap_pair = ("K", "L")
if cap_pair in series:
    v = np.array(series[cap_pair])
    print(f"""
  CAP JUNCTION  K-L = {v.mean():.3f} A

  Compare against the interior value. The cap's first two rungs are ordinary
  HET-s, so if it docked in register this spacing should be close to the bulk
  rise. A large deviation means the cap is sitting too far out or too close —
  worth knowing before the CG bead for that end is defined, since the capped
  bead's axial length depends on it.
""")

print("""
Reading this
  INTERIOR pairs are the ones that should set b and k_theta. End-effect pairs
  near A and near the cap are real physics but not bulk-representative.

  A std of ~0.1 A on a trajectory is normal thermal breathing. On a single
  static frame std will be 0 and only the mean is meaningful.
""")
