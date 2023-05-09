#General imports
import numpy as np
from PIL import Image
import math
from scipy.interpolate import griddata

# We define a class for each possible optical elements
# All classes include a self.transference_matrix generated from the specifications of the elements
# Most classes include setters and getters that also redefine the transference matrix  
# Elements alone do not posses any method that processes images

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

class AberratingPlaneRefraction:
    def __init__(self, f, n, t):
        self.f = f
        self.n = n
        self.t = t
        self.transference_matrix = np.array([[1 + (t/f)*(math.sqrt(n*n + 1)/n*n), -f],[1/f, 0]])

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

    
# Main class, used to define the optical system.
# This class has the ability of ray tracing and post processing of images
class OpticalSystem:

    def __init__(self, optical_elements = [] , distances = []):
        # We store the elements and distances between them
        self.optical_elements = optical_elements
        self.distances = distances

        # Reversing the elements and distances list to facilitate multiplication
        self.rev_optical_elements = list(reversed(optical_elements))
        self.rev_distances = list(reversed(distances))

        # Define an identity matrix in order to initialize the multiplication
        self.ABCD_matrix = np.identity(2)

        # Get the ABCD matrix of the system by multiplying all elements and translation matrix
        # Notice there must always be one less propagation than elements in the telescope
        # This means first and last matrix will always correspond to optical elements. 
        for i in range(len(self.optical_elements)):
            self.ABCD_matrix = self.ABCD_matrix @ self.rev_optical_elements[i].transference_matrix
            if i != len(self.optical_elements) - 1:
                self.ABCD_matrix = self.ABCD_matrix @ np.array([[1, self.rev_distances[i]], [0, 1]])
        
        # Verify the matrix just in case
        print(self.ABCD_matrix)


    # Ray tracing
    # Notice this method lets you select a sensor (ej. Camera sensor, eye, etc...)
    # Notice also that there's an attenuator factor. 
    # For low magnification telescopes set attenuation to 1
    # For high magnification telescopes it's recommended to use between 13 and 20 depending on the original image size
    def observe(self, object_path, so, si, sensor, res, attenuator, max_angle ,n_so = 1 ,type = 'eye', dist_eyepice_sensor = 20):

        # Load image
        object = Image.open(object_path, "r")
        width, height = object.size

        # Get transversal magnification of the system (Depending on the telescope it may be mt > 0 or mt < 0)
        # If mt < 0 it is recommended to set it to one, otherwise, image resolution may be lost due to lack of pixels
        mt = self.ABCD_matrix[0, 0] 

        usable_mt = mt/attenuator #arbitrary atenuator of transverse magnification factor
        
        # Gen output pixel dimensions based on the transversal magnification
        width_output = int(width*(abs(usable_mt)))
        height_output = int(height*(abs(usable_mt)))

        # Create empty pixel map in order to build the new image with white background
        image = Image.new("RGB", (width_output, height_output), "white")

        # Load the image as a pixel map
        pixels = image.load()


        # Go through the pixel map, going pixel by pixel
        # Notice this algorith has complexity n^3, therefore it may take up to 10 minutes for large images
        for i in range(width):
            for j in range(height):

                pos_x = i
                pos_y = j

                # Get pixel position
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

                # Generate m rays from 0 to the max acceptation angle of the telescope in radians
                entry_angles = np.linspace(0, max_angle, 16) #0.04548328
                
                # For each pixel propagate the m rays 
                for alpha in entry_angles:

                    # Generate entry bector
                    ray_ie = np.array([y_object, alpha]).reshape(2, 1)
                    # Just notation
                    d = dist_eyepice_sensor

                    # Compute exit ray using ABCD matrix plus the senor and the other elements
                    ray_is = np.array([[1, si],[0, 1]]) @ sensor.transference_matrix @ np.array([[1, d],[0, 1]]) @ self.ABCD_matrix @ np.array([[1, n_so/so],[0, 1]]) @ ray_ie

                    # Get y_image for the exit ray and calculate the transversal magnification for the ray
                    y_image = ray_is[0]
                    Mt = y_image / y_object

                    # Calculate the exit ray postion in pixels (absolute)
                    x_prime = -Mt*x
                    y_prime = Mt*y

                    # Calculate the exit ray postion in pixels relative to the center
                    pos_x_prime = int(x_prime + width_output/2)
                    pos_y_prime = int(y_prime + height_output/2)

                    # Verify if the rays position is inside the exit frame
                    # If not go to the next ray
                    if pos_y_prime <= 0 or pos_y_prime >= height_output:   
                        continue 
                    elif pos_x_prime <= 0 or pos_x_prime >= width_output:
            	        continue
                        

                    # Catch exception for bug that occurs when a ray outside the range goes through the verification method
                    rays_lost = 0
                    try:
                        new_gray = (int(pixel) + pixels[pos_x_prime, pos_y_prime][0])/2
                    except:
                        rays_lost += 1
                    
                    # Position RGB values from the incoming pixel in the resulting pixel
                    pix_fin = ( int(new_gray), int(new_gray), int(new_gray) )        
                    pixels[pos_x_prime, pos_y_prime] = pix_fin


        # Verify how many rays were lost
        print(rays_lost)

        # Save the image
        image.save('output/moon_out.png', format='PNG')

        return pixels, width_output, height_output
    

    # This method is identical to observe, the only difference is that the focal lenght for the primary mirror varies 
    # in function of the mirror height
    def observe_spherical_aberration(self, object_path, so, si, sensor, res, attenuator, max_angle, n_so = 1 ,type = 'eye', dist_eyepice_sensor = 20):
        object = Image.open(object_path, "r")
        width, height = object.size

        mt = self.ABCD_matrix[0, 0] 

        # It is needed because most images would get extremely big if normal mt was used
        usable_mt = mt/attenuator #20 is an arbitrary signifying the mag already in the picture atenuator of transverse magnification factor
        
        # Output dimensions
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

                #Generate an ABCD matrix with spherical aberration in the primary mirror
                temporal_ABCD = self.gen_aberrated_ABCD_mat(h, so, si)                

                #For objects in infinite all rays enter parallel
                #alpha = 0
                #ray_e = np.array([y_object, alpha]).reshape(2, 1) #Enter ray

                #Generate a vector with m angles between 0 and the max acceptance angle for the telescope
                entry_angles = np.linspace(0, max_angle, 20) #2.606 grados #0.04548328

                for alpha in entry_angles:
                    ray_ie = np.array([y_object, alpha]).reshape(2, 1)
                    d = dist_eyepice_sensor

                    # Propagate the entrance ray through the system with modified ABCD_matrix
                    ray_is = np.array([[1, si],[0, 1]]) @ sensor.transference_matrix @ np.array([[1, d],[0, 1]]) @ temporal_ABCD @ np.array([[1, n_so/so],[0, 1]]) @ ray_ie

                    y_image = ray_is[0]
                    Mt = y_image / y_object


                    # Same as observe for the rest of the method
                    x_prime = -Mt*x
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

        # The function that models shperical aberration is a prototype derived from the formula used 
        # for a refraction in a spherical surface
        # h is taken proportional to the real height of the incoming ray
        f_prime = f + (h*h)*((1/2*so)*((1/so)+(1/primary_radius))**2 + (1.337/2*si)*((1/primary_radius)+(1/si))**2)

        # Redefine the spherical mirror
        new_radius = -2*f_prime
        temp_mirror = SphericalMirrror(new_radius)
        temp_elements = list(self.rev_optical_elements.copy())
        temp_elements[-1] = temp_mirror
        

        temp_ABCD = np.identity(2)

        # Generate a New system matrix with the modified mirror
        for i in range(len(temp_elements)):
            temp_ABCD = temp_ABCD @ temp_elements[i].transference_matrix
            if i != len(temp_elements) - 1:
                temp_ABCD = temp_ABCD @ np.array([[1, self.rev_distances[i]], [0, 1]])

        return temp_ABCD

    

    def visualize_results(self, image_path):
        # Open image with pilloe
        obj = Image.open(image_path, "r")
        width, height = obj.size
        pixels = obj.load()
        # Bicuadratic interpolation
        pixels = interpolation(pixels, width, height)
        for col in range(width):
            for row in range(height):
                if pixels[col,row] == (255, 0, 0):
                    pixels[col,row] = (0, 0 ,0)
        #Visualize the results 
        obj.show()

#Interpolation function to recover the image
def interpolation (pixels, width_output, height_output):

  arry = np.zeros((width_output, height_output))

  # Go through the first channel of the Pixel Map
  # Only the R channel is needed sine RGB values are equal for all the images
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
                


            

            
                

            
        
    