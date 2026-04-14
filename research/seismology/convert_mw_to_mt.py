"""
Convert Mw -> Mrr, Mtt, Mpp for an explosion

Simple conversion script to go from a declared magnitude value to an isotropic
moment tensor (explosion), assuiming Mrr = Mtt = Mpp and zero off-diagnoals.
This is very simplifying because magnitudes are not always given in Mw or M0,
but this gives a first order approximation, to be refined later if necessary.
"""
import sys
from math import sqrt

magnitude = float(sys.argv[1])  # Mw

# Convert moment magnitude to seismic moment (Mw -> M0; Hanks & Kanamori 1979)
c = 10.7  # for units of N*m
m0 = 10 ** ((3/2) * (magnitude + c))

# Distribute seismic moment, M0, across 3 moment tensor components Mrr, Mtt, Mpp
m_ii = sqrt(2/3 * m0 ** 2)

print(f"Mw={magnitude} -> m0={m0:.3e} -> m_ii={m_ii:.3e}")

m0_check = (1 / sqrt(2)) * (sqrt(3 * m_ii ** 2))

print(f"(CHECK: {m0_check:.3e} == {m0:.3e})")


