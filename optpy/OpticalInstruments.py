import numpy as np

class ThinLens:
    def __init__(self, focal_length):
        self.focal_lenght = focal_length
        self.transference_matrix = np.array([[1, 0],[2/focal_length, 0]])

    def get_transference_matrix(self):
        return self.transference_matrix
    
    def get_focal_length(self):
        return self.focal_lenght
    
    def set_focal_length(self, new_focal_length):
        self.focal_lenght = new_focal_length
        self.transference_matrix = np.array([[1, 0],[2/self.focal_length, 0]])

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
        self.transferene_matrix = np.array([[1, 0],[(1/self.radius)*((n1/n2)-1), n1/n2]])

    def set_n2(self, new_n2):
        self.n1 = new_n2
        self.transferene_matrix = np.array([[1, 0],[(1/self.radius)*((n1/n2)-1), n1/n2]])

    
    def set_radius(self, radius):
        pass
    
