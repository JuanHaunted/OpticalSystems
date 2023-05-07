import numpy as np

eye_dist = np.array([[1, 18.78/1.337],[0, 1]]) #1.337
eye = np.array([[1, 0],[-1/17, 1]])
eye_relief = np.array([[1, 20],[0, 1]])
eyepice = np.array([[1, 0],[-1/40, 1]])
sec_eye_dist = np.array([[1, 3950],[0, 1]])
sec_mirror = np.array([[1, 0], [2/(57.15), 1]])
prim_sec_trans = np.array([[1, -59.7745],[0, 1]]) #116.82
primary_mirror = np.array([[1, 0],[-2/(178), 1]]) #1.398e-4
aperture_mirror = np.array([[1, 500],[0, 1]])

imaging = np.array([[1, 3900],[0, 1]])
fifty_tran = np.array([[1, 50],[0, 1]])

plane_mirror = np.array([[1, 0],[0, -1]])

primary_eyepie_abcd = eyepice @ sec_eye_dist @ sec_mirror @ prim_sec_trans @ primary_mirror
ABCD = eye_dist @ eye @ eye_relief @ eyepice @ sec_eye_dist @ sec_mirror @ prim_sec_trans @ primary_mirror @ aperture_mirror
no_eyepiece = sec_mirror @ prim_sec_trans @ primary_mirror @ aperture_mirror
ABCD_2 = eye_dist @ eye @ eye_relief @ eyepice @ fifty_tran @ plane_mirror @ imaging @ sec_mirror @ prim_sec_trans @ primary_mirror @ aperture_mirror
ABCD_3 = eyepice @ np.array([[1, 450],[0, 1]]) @ plane_mirror @ np.array([[1, 3500],[0, 1]]) @ sec_mirror @ prim_sec_trans @ primary_mirror



#print(primary_eyepie_abcd)
#print(np.linalg.det(ABCD))
#print(np.linalg.det(no_eyepiece))
#print(np.linalg.det(primary_eyepie_abcd))
#print(ABCD_3)