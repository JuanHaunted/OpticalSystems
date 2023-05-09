import OpticalInstruments as opi
from PIL import Image

primary_mirror = opi.SphericalMirrror(-178)
secondary_mirror = opi.SphericalMirrror(57.15)
eyepiece = opi.ThinLens(40)#40
eye = opi.ThinLens(18.78, 1.337)
thick_eyepiece = opi.ThickLens(1.5, 26, 80, 10)
deg_mirror = opi.PlaneMirror(45)

system = [primary_mirror, secondary_mirror, deg_mirror, eyepiece]

system_2 = [primary_mirror, secondary_mirror, deg_mirror, thick_eyepiece, eye]
dist = [-59.7745, 2990, 960] #960

cassegrain_telescope = opi.OpticalSystem(system, dist)

obj_path = 'images/moon.png' #Path of the object
resolution = 49.66 * 1e+6 #Pixel size in mm
so = 384400 * 1e+6 #Distance to object
si = 18.78 * 1e-6 #Distance to image
max_angle = 0.04548328 # Max acceptation angle in radians

atten = 14
#data = cassegrain_telescope.observe(obj_path, so, si, eye, resolution, atten)

#Aberrated DATA
aberrated_data = cassegrain_telescope.observe(obj_path, so, si, eye, resolution, atten, max_angle)



#We interpolate and visualize results
#cassegrain_telescope.visualize_results('ouput/moon_out.png')

#image.save('output/moon_out.png', format='PNG')







print(cassegrain_telescope.ABCD_matrix)
