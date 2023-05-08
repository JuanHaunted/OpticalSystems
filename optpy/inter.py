import OpticalInstruments as opi
from PIL import Image

obj = Image.open("output/moon_out.png", "r")
width, height = obj.size
pixels = obj.load()
pixels = opi.interpolation(pixels, width, height)
for col in range(width):
    for row in range(height):
        if pixels[col,row] == (255, 0, 0):
            pixels[col,row] = (0, 0 ,0)

obj.show()