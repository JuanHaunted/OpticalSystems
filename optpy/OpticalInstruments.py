import numpy as np

class ThinLens:
    def __init__(self, focal_length, n = 1):
        self.focal_lenght = focal_length
        self.n = n
        self.transference_matrix = np.array([[1, 0],[-n/focal_length, 0]])

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
        self.transferene_matrix = np.array([[a1, a2], [a3, a4]])
    
    def get_transference_matrix(self):
        return self.transferene_matrix
    
    def get_radius(self):
        return (self.r1, self.r2)
    
    def get_specs(self):
        print(f'r1 = {self.r1}\n r2 = {self.r2} \n d = {self.d} \n n = {self.n}')

    
class OpticalSystem:
    def __init__(self, si, so, optical_elements = [], distances = []):
        self.optical_elements = optical_elements
        self.distances = distances .
        distances = distances.reverse()
        optical_elements = optical_elements.reverse()
        initial_propagation = np.array([[1, so],[0, 1]])
        final_propagation = np.array([[1, si],[0, 1]])
        self.telescope_matrix = np.identity()
        self.ABCD_matrix = np.identity()
        

        for i in range(len(optical_elements)):
            self.ABCD_matrix = self.ABCD_matrix @ optical_elements[i].transference_matrix
            if i != len(optical_elements) - 1:
                self.ABCD_matrix @ np.array([[1, distances[i]], [0, 1]])


        self.propagation_matrix = final_propagation @ self.ABCD_matrix @ initial_propagation
    