
import Image
import ImageDraw
import ImageFont
from math import *
import matplotlib.pyplot as plt
import numpy as np

def plotMesh(mesh, bit_size = 500):

    # create image
    bit_x = float(bit_size) / mesh.width
    bit_y = float(bit_size) / mesh.height
    img = Image.new('RGB', (bit_size,bit_size), 'white')
    draw = ImageDraw.Draw(img)
    
    cw = mesh.cells_x
    y_val = 0.0
    
    # draw cells
    for y in range(mesh.cells_y):
        x_val = 0.0
        for x in range(mesh.cells_x):
            
            cell = mesh.cells[y*cw+x]
        
            # fuel red; moderator blue
            if cell.material.id == 'fuel':
                draw.rectangle([x_val, y_val, x_val + mesh.widths[x] * bit_x, y_val + mesh.heights[y] * bit_y], (255,0,0))
            else:
                draw.rectangle([x_val, y_val, x_val + mesh.widths[x] * bit_x, y_val + mesh.heights[y] * bit_y], (0,0,255))
                
            x_val += mesh.widths[x] * bit_x
        
        y_val += mesh.heights[y] * bit_y
                

    # draw vertical grid lines
    x_val = 0.0
    for x in range(mesh.cells_x-1):
        x_val += mesh.widths[x]*bit_x
        draw.line((x_val, 0, x_val, bit_size), fill=(0,0,0), width=2)

    # draw horizontal grid lines
    y_val = 0.0
    for y in range(mesh.cells_y-1):
        y_val += mesh.heights[y]*bit_y
        draw.line((0, y_val, bit_size, y_val), fill=(0,0,0), width=2)
        

    # save image
    img.save('material.png')
    
    