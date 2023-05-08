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

def observe(self, object_path, so, si, sensor, res, n_sensor = 1 ,type = 'eye', dist_eyepice_sensor = 20):
        object = Image.open(object_path, "r")
        width, height = object.size

        mt = self.ABCD_matrix[0, 0] 

        # It is needed because most images would get extremely big if normal mt was used
        usable_mt = mt/20 #20 is an arbitrary atenuator of transverse magnification factor
        
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
                alpha = 0
                ray_e = np.array([y_object, alpha]).reshape(2, 1) #Enter ray

                #Principal_ray
                alpha_p = math.atan(y_object/so)
                ray_ep = np.array([y_object, alpha_p]).reshape(2, 1)
                
                d = dist_eyepice_sensor
                ray_s = np.array([[1, si],[0, 1]]) @ sensor.transference_matrix @ np.array([[1, d],[0, 1]]) @ self.ABCD_matrix @ np.array([[1, n_sensor/so],[0, 1]]) @ ray_e
                ray_sp = np.array([[1, si],[0, 1]]) @ sensor.transference_matrix @ np.array([[1, d],[0, 1]]) @ self.ABCD_matrix @ np.array([[1, n_sensor/so],[0, 1]]) @ ray_ep

                y_image = ray_s[0]
                y_image_p = ray_sp[0]

                Mt = y_image / y_object
                Mtp = (y_image_p / y_object)


                #Conversion from image coordinates to mirror coordinates        
                x_prime = Mt*x
                y_prime = Mt*y

                pos_x_prime = int(x_prime + width_output/2)
                pos_y_prime = int(y_prime + height_output/2)


                converged = 0

                if pos_y_prime < 0 or pos_y_prime >= height_output:   
                    continue 
                elif pos_x_prime < 0 or pos_x_prime >= width_output:
            	    continue

                converged += 1
            

                new_gray = (int(pixel) + pixels[pos_x_prime, pos_y_prime][0])/2
                if i == 1:
                    print(new_gray)

                pix_fin = ( int(new_gray), int(new_gray), int(new_gray) )        
                pixels[pos_x_prime, pos_y_prime] = pix_fin


        image.save('output/moon_out.png', format='PNG')

        return pixels