import OpticalInstruments as opi
from PIL import Image

primary_mirror = opi.SphericalMirrror(-178)
secondary_mirror = opi.SphericalMirrror(57.15)
eyepiece = opi.ThinLens(40)
eye = opi.ThinLens(18.78, 1.337)
thick_eyepiece = opi.ThickLens(1.5, 26, 80, 10)
deg_mirror = opi.PlaneMirror(45)

#two lens eyepice 
thin_eyep_one = opi.ThinLens(-12.65)
thin_eyep_two = opi.ThinLens(12.65)

system = [primary_mirror, secondary_mirror, deg_mirror, thin_eyep_one, thin_eyep_two]
dist = [-59.7745, 2990, 960, 4]

cassegrain_telescope = opi.OpticalSystem(system, dist)