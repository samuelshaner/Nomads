
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
    
    x = np.zeros(solver.mesh.cells_x*100)
    ng = solver.mesh.cells[0].material.num_groups
    flux = np.zeros((ng, solver.mesh.cells_x*100))

    for xc in range(solver.mesh.cells_x):
        cell = solver.mesh.cells[xc]
        x_val = sum(solver.mesh.widths[0:xc])
        
        for i in range (100):
            
            x[xc*100+i] = x[xc*100+i-1] + cell.width/100.0
            
            for e in range(ng):
                if solver.method == 'NEM4':
                    flux[e,xc*100+i] = solver.phi[xc*ng+e] + solver.coeffs[4*xc*ng+4*e] * solver.P1(i/100.0) + solver.coeffs[4*xc*ng+4*e+1] * solver.P2(i/100.0) + solver.coeffs[4*xc*ng+4*e+2] * solver.P3(i/100.0) + solver.coeffs[4*xc*ng+4*e+3] * solver.P4(i/100.0)
                else:
                    flux[e,xc*100+i] = solver.phi[xc*ng+e] + solver.coeffs[2*xc*ng+2*e] * solver.P1(i/100.0) + solver.coeffs[2*xc*ng+2*e+1] * solver.P2(i/100.0) 
                    
    
    plt.figure()
    for e in range(ng):
        plt.plot(x,flux[e,:])
    
    plt.savefig(solver.method + '_flux.png')
    
    
def plotCurrent(solver):
    
    x = np.zeros(solver.mesh.cells_x*100)
    ng = solver.mesh.cells[0].material.num_groups
    current = np.zeros((ng, solver.mesh.cells_x*100))

    for xc in range(solver.mesh.cells_x):
        cell = solver.mesh.cells[xc]
        x_val = sum(solver.mesh.widths[0:xc])
        
        for i in range (100):
            
            x[xc*100+i] = x[xc*100+i-1] + cell.width/100.0
            
            for e in range(ng):
                if solver.method == 'NEM4':
                    current[e,xc*100+i] = - cell.material.D[e] / cell.width * (solver.coeffs[4*xc*ng+4*e] * solver.dP1(i/100.0) + solver.coeffs[4*xc*ng+4*e+1] * solver.dP2(i/100.0) + solver.coeffs[4*xc*ng+4*e+2] * solver.dP3(i/100.0) + solver.coeffs[4*xc*ng+4*e+3] * solver.dP4(i/100.0))    
                else:
                    current[e,xc*100+i] = - cell.material.D[e] / cell.width * (solver.coeffs[2*xc*ng+2*e] * solver.dP1(i/100.0) + solver.coeffs[2*xc*ng+2*e+1] * solver.dP2(i/100.0))                    
    
    plt.figure()
    for e in range(ng):
        plt.plot(x,current[e,:])
    
    plt.savefig(solver.method + '_current.png')    
    
    
    