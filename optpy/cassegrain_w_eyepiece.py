import OpticalInstruments as opi
from PIL import Image

# We define the components of the Optical System using OPI
# Notice that not all elements go into the system and they are included purely
# for testing the eyepiece
primary_mirror = opi.SphericalMirrror(-178)
secondary_mirror = opi.SphericalMirrror(57.15)
eyepiece = opi.ThinLens(40)
eye = opi.ThinLens(18.78, 1.337)
thick_eyepiece = opi.ThickLens(1.5, 26, 80, 10)
deg_mirror = opi.PlaneMirror(45)


# Testing for generating spherical aberration
aberrating_plate = opi.AberratingPlaneRefraction(200, 1.5, 1)

#two lens eyepice simplified to thin lens
#It's designed to mitigate the chromatic aberration
thin_eyep_one = opi.ThinLens(-12.65)
thin_eyep_two = opi.ThinLens(12.65)

#We define 
system = [primary_mirror, secondary_mirror, deg_mirror, thin_eyep_one, thin_eyep_two]
dist = [-59.7745, 2990, 960, 4]

obj_path = 'images/moon.png' #Path of the object
resolution = 49.66 * 1e+6 #Pixel size in mm
so = 384400 * 1e+6 #Distance to object
si = 18.78 * 1e-6 #Distance to image
atten = 13
max_angle = 0.04548328 # Max acceptation angle in radians


cassegrain_telescope = opi.OpticalSystem(system, dist)


data = cassegrain_telescope.observe(obj_path, so, si, eye, resolution, atten, max_angle)