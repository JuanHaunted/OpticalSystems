import OpticalInstruments as opi
from PIL import Image


primary = opi.SphericalMirrror(2400)
secondary = opi.PlaneMirror(angle=45)
eyepiece = opi.ThinLens(42)
eye = opi.ThinLens(10)

elements = [primary, secondary, eyepiece]

distances = [-1067, 175]

newtonian_telescope = opi.OpticalSystem(elements, distances)

obj_path = 'images/moon.png'
so = 34400000
si = 0.010
res =  49.66 * 1e+6

data = newtonian_telescope.observe_spherical_aberration(obj_path, so, si, eye, res, 1)

