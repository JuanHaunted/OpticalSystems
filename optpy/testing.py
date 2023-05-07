import OpticalInstruments as opi
from test import ABCD_3

primary_mirror = opi.SphericalMirrror(-178)
secondary_mirror = opi.SphericalMirrror(57.15)
eyepiece = opi.ThinLens(40)
eye = opi.ThinLens(18.78, 1.337)
thick_eyepiece = opi.ThickLens(1.5, 26, 80, 10)
deg_mirror = opi.PlaneMirror(45)

system = [primary_mirror, secondary_mirror, deg_mirror, eyepiece]
system_2 = [primary_mirror, secondary_mirror, deg_mirror, thick_eyepiece, eye]
dist = [-59.7745, 2990, 960]

cassegrain_telescope = opi.OpticalSystem(system, dist)
print(cassegrain_telescope.ABCD_matrix)
