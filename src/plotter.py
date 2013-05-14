
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
    
    
def plotFlux(solver):
    
    xf = []
    xm = []
    phi_fuel = []
    phi_mod = []
    phi_fuel0 = []
    phi_mod0 = []    
    for i in range(11):
        xf.append(i)
        xm.append(10 + i)
        phi_fuel.append(solver.phi[1] + solver.coeffs[4] * solver.P1(i/10.0) + solver.coeffs[5] * solver.P2(i/10.0) + solver.coeffs[6] * solver.P3(i/10.0) + solver.coeffs[7] * solver.P3(i/10.0))
        phi_mod.append(solver.phi[3] + solver.coeffs[12] * solver.P1(i/10.0) + solver.coeffs[13] * solver.P2(i/10.0) + solver.coeffs[14] * solver.P3(i/10.0) + solver.coeffs[15] * solver.P3(i/10.0))         
        phi_fuel0.append(solver.phi[0] + solver.coeffs[0] * solver.P1(i/10.0) + solver.coeffs[1] * solver.P2(i/10.0) + solver.coeffs[2] * solver.P3(i/10.0) + solver.coeffs[3] * solver.P3(i/10.0))
        phi_mod0.append(solver.phi[2] + solver.coeffs[8] * solver.P1(i/10.0) + solver.coeffs[9] * solver.P2(i/10.0) + solver.coeffs[10] * solver.P3(i/10.0) + solver.coeffs[11] * solver.P3(i/10.0))        
    
    plt.figure()
    plt.plot(xf,phi_fuel, 'r')
    plt.plot(xf,phi_fuel0, 'b')
    plt.plot(xm,phi_mod, 'r')
    plt.plot(xm,phi_mod0, 'b')
    plt.savefig('NEM4_flux.png')
    
def plotCurrent(solver):
    
    xf = []
    xm = []
    J_fuel = []
    J_mod = []
    J_fuel0 = []
    J_mod0 = []    
    
    for i in range(11):
        xf.append(i)
        xm.append(10 + i)
        J_fuel.append(solver.phi[1] + solver.coeffs[4] * solver.P1(i/10.0) + solver.coeffs[5] * solver.P2(i/10.0) + solver.coeffs[6] * solver.P3(i/10.0) + solver.coeffs[7] * solver.P3(i/10.0))
        J_mod.append(solver.phi[3] + solver.coeffs[12] * solver.P1(i/10.0) + solver.coeffs[13] * solver.P2(i/10.0) + solver.coeffs[14] * solver.P3(i/10.0) + solver.coeffs[15] * solver.P3(i/10.0))         
        J_fuel0.append(solver.phi[0] + solver.coeffs[0] * solver.P1(i/10.0) + solver.coeffs[1] * solver.P2(i/10.0) + solver.coeffs[2] * solver.P3(i/10.0) + solver.coeffs[3] * solver.P3(i/10.0))
        J_mod0.append(solver.phi[2] + solver.coeffs[8] * solver.P1(i/10.0) + solver.coeffs[9] * solver.P2(i/10.0) + solver.coeffs[10] * solver.P3(i/10.0) + solver.coeffs[11] * solver.P3(i/10.0))        
    
    plt.figure()
    plt.plot(xf,J_fuel, 'r')
    plt.plot(xf,J_fuel0, 'b')
    plt.plot(xm,J_mod, 'r')
    plt.plot(xm,J_mod0, 'b')
    plt.savefig('NEM4_flux.png')
    
    
    
    
    