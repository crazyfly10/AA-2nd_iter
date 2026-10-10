TYPES = """  OW  OW  8  16.00000  0.00000  A  3.15061e-01  6.36386e-01
  HW  HW  1   1.00800  0.00000  A  0.00000e+00  0.00000e+00
  Na  Na 11  22.98977  0.00000  A  3.32840e-01  1.15897e-02
  Cl  Cl 17  35.45300  0.00000  A  4.40104e-01  4.18400e-01
"""
MOLS = """
[ moleculetype ]
SOL     2
[ atoms ]
1   OW  1  SOL  OW   1  -0.834  16.00000
2   HW  1  SOL  HW1  1   0.417   1.00800
3   HW  1  SOL  HW2  1   0.417   1.00800
[ settles ]
1   1   0.09572   0.15139
[ exclusions ]
1 2 3
2 1 3
3 1 2

[ moleculetype ]
NA      1
[ atoms ]
1   Na  1  NA  NA  1   1.00000  22.98977

[ moleculetype ]
CL      1
[ atoms ]
1   Cl  1  CL  CL  1  -1.00000  35.45300

"""
lines, out, in_at = open("linked_GMX.top").readlines(), [], False
for l in lines:
    if l.strip().startswith("[ atomtypes ]"): in_at = True; out.append(l); continue
    if in_at and l.strip().startswith("["):   out.append(TYPES); in_at = False
    if l.strip().startswith("[ system ]"):    out.append(MOLS)
    out.append(l)
open("linked_GMX.top", "w").writelines(out)
