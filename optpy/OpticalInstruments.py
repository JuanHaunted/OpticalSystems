import numpy as np
from PIL import Image
import math
from scipy.interpolate import griddata

class ThinLens:
    def __init__(self, focal_length, n = 1):
        self.focal_lenght = focal_length
        self.n = n
        self.transference_matrix = np.array([[1, 0],[-n/focal_length, 1]])

    def get_transference_matrix(self):
        return self.transference_matrix
    
    def get_focal_length(self):
        return self.focal_lenght
    
    def set_focal_length(self, new_focal_length):
        self.focal_lenght = new_focal_length
        self.transference_matrix = np.array([[1, 0],[-self.n/self.focal_length, 0]])

class SphericalMirrror:
    def __init__(self, radius):
        self.radius = radius
        self.transference_matrix = np.array([[1, 0],[2/radius, 1]])

    def get_transference_matrix(self):
        return self.transference_matrix
    
    def get_radius(self):
        return self.radius
    
    def set_focal_length(self, new_radius):
        self.focal_lenght = new_radius
        self.transference_matrix = np.array([[1, 0],[2/self.radius, 1]])

class PlaneMirror:
    def __init__(self, angle=0):
        if angle == 0:
            self.transference_matrix = np.identity()
        else:
            self.transference_matrix = np.array([[1, 0],[0, -1]])

    def get_transference_matrix(self):
        return self.transference_matrix

class PlaneRefraction:
    def __init__(self, n1, n2):
        self.n1 = n1
        self.n2 = n2
        self.transferene_matrix = np.array([[1, 0],[0, n1/n2]])

    def get_transference_matrix(self):
        return self.transference_matrix
    
    def get_n1(self):
        return self.n1

    def get_n2(self):
        return self.n2

    def set_n1(self, new_n1):
        self.n1 = new_n1
        self.transferene_matrix = np.array([[1, 0],[0, new_n1/self.n2]])

    def set_n2(self, new_n2):
        self.n2 = new_n2
        self.transferene_matrix = np.array([[1, 0],[0, self.n1/new_n2]])

class SphericalRefraction:
    def __init__(self, n1, n2, radius):
        self.n1 = n1
        self.n2 = n2
        self.radius = radius
        self.transferene_matrix = np.array([[1, 0],[(1/radius)*((n1/n2)-1), n1/n2]])

    def get_n1(self):
        return self.n2

    def get_n2(self):
        return self.n2
    
    def get_radius(self):
        return self.radius
    
    def set_n1(self, new_n1):
        self.n1 = new_n1
        self.transferene_matrix = np.array([[1, 0],[(1/self.radius)*((new_n1/self.n2)-1), new_n1/self.n2]])

    def set_n2(self, new_n2):
        self.n1 = new_n2
        self.transferene_matrix = np.array([[1, 0],[(1/self.radius)*((self.n1/new_n2)-1), self.n1/new_n2]])

    
    def set_radius(self, new_radius):
        self.radius = new_radius
        self.transferene_matrix = np.array([[1, 0],[(1/new_radius)*((self.n1/self.n2)-1), self.n1/self.n2]])

class ThickLens:
    def __init__(self, nl, R1, R2, dl):
        D1 = (nl - 1) / R1
        D2 = (nl - 1) / (-R2)

    # Lens matrix
        a1 = (1 - (D2 * dl) / nl)
        a2 = -D1 - D2 + (D1 * D2 * dl / nl)
        a3 = dl / nl
        a4 = (1 - (D1 * dl) / nl)
        
        self.r1 = R1
        self.r2 = -R2
        self.d = dl
        self.n = nl
        self.transference_matrix = np.array([[a1, a3], [a2, a4]])
    
    def get_transference_matrix(self):
        return self.transferene_matrix
    
    def get_radius(self):
        return (self.r1, self.r2)
    
    def get_specs(self):
        print(f'r1 = {self.r1}\n r2 = {self.r2} \n d = {self.d} \n n = {self.n}')

    

class OpticalSystem:

    def __init__(self, optical_elements = [] , distances = []):
        self.optical_elements = optical_elements
        self.distances = distances
        self.rev_optical_elements = list(reversed(optical_elements))
        self.rev_distances = list(reversed(distances))
        self.ABCD_matrix = np.identity(2)

        
        for i in range(len(self.optical_elements)):
            self.ABCD_matrix = self.ABCD_matrix @ self.rev_optical_elements[i].transference_matrix
            if i != len(self.optical_elements) - 1:
                self.ABCD_matrix = self.ABCD_matrix @ np.array([[1, self.rev_distances[i]], [0, 1]])
            


    def observe(self, object_path, so, si, sensor, res, attenuator, n_so = 1 ,type = 'eye', dist_eyepice_sensor = 20):
        object = Image.open(object_path, "r")
        width, height = object.size

        mt = self.ABCD_matrix[0, 0] 

        # It is needed because most images would get extremely big if normal mt was used
        usable_mt = mt/attenuator #20 is an arbitrary atenuator of transverse magnification factor
        
        # Output weight 
        width_output = int(width*(abs(usable_mt)))
        height_output = int(height*(abs(usable_mt)))

        # Create empty pixel map in order to build the new image
        image = Image.new("RGB", (width_output, height_output), "white")
        pixels = image.load()

        for i in range(width):
            for j in range(height):

                pos_x = i
                pos_y = j

                pixel = object.getpixel((pos_x, pos_y))
               

                #Generate approximate center coordinates for the image
                x = pos_x - width/2
                y = pos_y - height/2


                #Calculate distance from center to pixel
                r = math.sqrt((x*x) + (y*y)) + 1

                y_object = r * res #Resolution equals size of pixel in mm

                #For objects in infinite all rays enter parallel
                #alpha = 0
                #ray_e = np.array([y_object, alpha]).reshape(2, 1) #Enter ray

                #Principal_ray
                entry_angles = np.linspace(0, 0.04548328, 12) #2.606 grados #0.04548328

                for alpha in entry_angles:
                    ray_ie = np.array([y_object, alpha]).reshape(2, 1)
                    d = dist_eyepice_sensor
                    ray_is = np.array([[1, si],[0, 1]]) @ sensor.transference_matrix @ np.array([[1, d],[0, 1]]) @ self.ABCD_matrix @ np.array([[1, n_so/so],[0, 1]]) @ ray_ie

                    y_image = ray_is[0]
                    Mt = y_image / y_object


                    x_prime = Mt*x
                    y_prime = Mt*y

                    pos_x_prime = int(x_prime + width_output/2)
                    pos_y_prime = int(y_prime + height_output/2)

                    if pos_y_prime <= 0 or pos_y_prime >= height_output:   
                        continue 
                    elif pos_x_prime <= 0 or pos_x_prime >= width_output:
            	        continue
                    
                    rays_lost = 0
                    try:
                        new_gray = (int(pixel) + pixels[pos_x_prime, pos_y_prime][0])/2
                    except:
                        rays_lost += 1

                    pix_fin = ( int(new_gray), int(new_gray), int(new_gray) )        
                    pixels[pos_x_prime, pos_y_prime] = pix_fin


        print(rays_lost)

        image.save('output/moon_out.png', format='PNG')

        return pixels, width_output, height_output
    

    def observe_spherical_aberration(self, object_path, so, si, sensor, res, attenuator, n_so = 1 ,type = 'eye', dist_eyepice_sensor = 20):
        object = Image.open(object_path, "r")
        width, height = object.size

        mt = self.ABCD_matrix[0, 0] 

        # It is needed because most images would get extremely big if normal mt was used
        usable_mt = mt/attenuator #20 is an arbitrary atenuator of transverse magnification factor
        
        # Output weight 
        width_output = int(width*(abs(usable_mt)))
        height_output = int(height*(abs(usable_mt)))

        # Create empty pixel map in order to build the new image
        image = Image.new("RGB", (width_output, height_output), "white")
        pixels = image.load()

        for i in range(width):
            for j in range(height):

                pos_x = i
                pos_y = j

                pixel = object.getpixel((pos_x, pos_y))
               

                #Generate approximate center coordinates for the image
                x = pos_x - width/2
                y = pos_y - height/2


                #Calculate distance from center to pixel
                r = math.sqrt((x*x) + (y*y)) + 1

                y_object = r * res #Resolution equals size of pixel in mm

                #We map the real height for the object to a height in the telescope aperture
                apperture = 356
                y_percent = y_object/((height/2)*res)
                h = (apperture/2)*y_percent

                temporal_ABCD = self.gen_aberrated_ABCD_mat(h, so, si)                

                #For objects in infinite all rays enter parallel
                #alpha = 0
                #ray_e = np.array([y_object, alpha]).reshape(2, 1) #Enter ray

                #Principal_ray
                entry_angles = np.linspace(0, 0.04548328, 20) #2.606 grados #0.04548328

                for alpha in entry_angles:
                    ray_ie = np.array([y_object, alpha]).reshape(2, 1)
                    d = dist_eyepice_sensor
                    ray_is = np.array([[1, si],[0, 1]]) @ sensor.transference_matrix @ np.array([[1, d],[0, 1]]) @ temporal_ABCD @ np.array([[1, n_so/so],[0, 1]]) @ ray_ie

                    y_image = ray_is[0]
                    Mt = y_image / y_object


                    x_prime = Mt*x
                    y_prime = Mt*y

                    pos_x_prime = int(x_prime + width_output/2)
                    pos_y_prime = int(y_prime + height_output/2)

                    if pos_y_prime <= 0 or pos_y_prime >= height_output:   
                        continue 
                    elif pos_x_prime <= 0 or pos_x_prime >= width_output:
            	        continue
                    
                    rays_lost = 0
                    try:
                        new_gray = (int(pixel) + pixels[pos_x_prime, pos_y_prime][0])/2
                    except:
                        rays_lost += 1

                    pix_fin = ( int(new_gray), int(new_gray), int(new_gray) )        
                    pixels[pos_x_prime, pos_y_prime] = pix_fin


        print(rays_lost)

        image.save('output/moon_out.png', format='PNG')

        return pixels, width_output, height_output

    def gen_aberrated_ABCD_mat(self, h, so, si):
        primary_radius = abs(self.optical_elements[0].get_radius())
        f = primary_radius/2
        f_prime = f + (h*h)*((1/2*so)*((1/so)+(1/primary_radius))**2 + (1.337/2*si)*((1/primary_radius)+(1/si))**2)

        new_radius = -2*f_prime
        temp_mirror = SphericalMirrror(new_radius)
        temp_elements = list(self.rev_optical_elements.copy())
        temp_elements[-1] = temp_mirror
        

        temp_ABCD = np.identity(2)

        for i in range(len(temp_elements)):
            temp_ABCD = temp_ABCD @ temp_elements[i].transference_matrix
            if i != len(temp_elements) - 1:
                temp_ABCD = temp_ABCD @ np.array([[1, self.rev_distances[i]], [0, 1]])

        return temp_ABCD

    

    def visualize_results(self, image_path):
        obj = Image.open(image_path, "r")
        width, height = obj.size
        pixels = obj.load()
        pixels = interpolation(pixels, width, height)
        for col in range(width):
            for row in range(height):
                if pixels[col,row] == (255, 0, 0):
                    pixels[col,row] = (0, 0 ,0)

        obj.show()


def interpolation (pixels, width_output, height_output):

  arry = np.zeros((width_output, height_output))
  for i in range(width_output):
    for j in range(height_output):
      arry[i,j] = pixels[i,j][0]

  # Get the coordinates of the non-white pixels
  nonwhite_coords = np.argwhere(arry != 255)
  # Get the coordinates of the white pixels
  white_coords = np.argwhere(arry == 255)

  # Get the pixel values of the non-white pixels
  nonwhite_pixels = arry[nonwhite_coords[:,0], nonwhite_coords[:,1]]

  # Interpolate the pixel values of the white pixels
  interpolated_pixels = griddata(nonwhite_coords, nonwhite_pixels, white_coords, method='linear', rescale=True)

  #Change resulting NaN values with some value (zero, for instance)
  interpolated_pixels = np.nan_to_num(interpolated_pixels, nan=127.0)

  #Round interpolated values to integers
  int_out = np.round(interpolated_pixels).astype(int)

  #Fill white pixels locations with interpolated values
  for i in range(int_out.shape[0]):
    pixels[white_coords[i,0], white_coords[i,1] ] = ( int_out[i], int_out[i], int_out[i] )

  #Returne pixels array with interpolated values
  return pixels
                


            

            
                

            
        
    