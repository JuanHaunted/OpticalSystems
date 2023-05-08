import OpticalInstruments as opi
from test import ABCD_3
from PIL import Image

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

obj_path = 'images/moon.png'
resolution = 49.66 * 1e+6
so = 384400 * 1e+6
si = 18.78 * 1e-6

data = cassegrain_telescope.observe(obj_path, so, si, eye, resolution)

pixels = opi.interpolation(data[0], data[1], data[2])

#image.save('output/moon_out.png', format='PNG')







print(cassegrain_telescope.ABCD_matrix)
