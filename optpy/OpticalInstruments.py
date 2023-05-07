import numpy as np
from PIL import Image

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
            



    def observe(self, object_path, so, si, sensor, res, n_sensor = 1,type = 'eye', dist_eyepice_sensor = 20):
        object = Image.open(object_path, "r")
        width, height = object.size

        mt = self.ABCD_matrix[0, 0] 

        # It is needed because most images would get extremely big if normal mt was used
        usable_mt = abs(mt/20) #20 is an arbitrary atenuator of transverse magnification factor
        
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
                



            
                

            
        
    