# Scan of the mass of planet B (supplementary note, not part of the Zenodo record V_1.0)

Question: is planet B (0.95 Earth mass in the monograph) required to be that massive for the retrograde flyby to keep the disk bound?

Method: the phase B tracer code of the monograph (700 tracers, disk from 2.90 to 5 R_E, 2 km/s at infinity, inclination 25 deg) is run with the mass of planet B and the periapsis distance q as free parameters; the radius of B follows R_B = 0.95 R_E (M_B / 0.95 M_E)^0.27, and passes inside the contact distance (R_eq + R_B) are excluded. Three seeds per cell, rotation 5 h (R_eq = 1.52 R_E) and 6 h (R_eq = 1.70 R_E, as in the original phase B script).

Result (5 h, q = 3.0 R_E): fraction of the disk kept bound, retrograde versus prograde flyby: M_B = 0.95 -> 0.47 versus 0.14 (ratio 3.5); 0.7 -> 0.50 / 0.23 (2.2); 0.5 -> 0.57 / 0.36 (1.6); 0.3 -> 0.74 / 0.65 (1.1); 0.1 -> 0.98 / 0.91 (1.1). The retrograde advantage disappears for a light planet B, and the inclination of the bound orbits falls from 31 deg to 15-20 deg.

Limits: test particles only (no braking of the Earth's rotation, no assembly of the body); mass-radius law chosen by the author; three seeds; 5 h and 6 h cases only.

Files: scanMB.py (one run), grid.py and grid2.py (grids), stripping.py and phaseB_band.py (geometry and integrator, unchanged from the monograph code), balayage_masse_B.csv (all results).
Run: python3 grid.py (about 5 minutes on two cores).
