
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

    # find number of cell types
    num_cell_types = 0
    cell_types = np.array([])
    for y in range(mesh.cells_y):
        for x in range(mesh.cells_x):
            cell = mesh.cells[y*cw+x]
            
            if cell.material.id not in cell_types:
                cell_types = np.append(cell_types, cell.material.id)
                num_cell_types += 1
    
    # draw cells
    val_max = len(cell_types)
    for y in range(mesh.cells_y):
        x_val = 0.0
        for x in range(mesh.cells_x):
            
            red = 0.0
            blue = 0.0
            green = 0.0
            
            cell = mesh.cells[y*cw+x]
            
            val = np.where(cell_types == cell.material.id)[0]
            
            # get color
            if (float(val) / val_max <= 1.0/3.0):
                red = 0.0
                green = 3.0 * val / val_max
                blue = 1.0
            elif (float(val) / val_max <= 2.0/3.0):
                red = 3.0 * val / val_max - 1.0
                green = 1.0
                blue = -3.0 * val / val_max + 2.0
            else:
                red = 1.0
                green = -3.0 * val / val_max + 3.0
                blue = 0.0
            
            # convert color to RGB triplet
            red = int(255*red)
            green = int(255*green)
            blue = int(255*blue)

            draw.rectangle([x_val, y_val, x_val + mesh.widths[x] * bit_x, y_val + mesh.heights[y] * bit_y], (red,green,blue))
                
            x_val += mesh.widths[x] * bit_x
        
        y_val += mesh.heights[y] * bit_y
                

    if mesh.cells_x <= 100 and mesh.cells_y <= 100:
        # draw vertical grid lines
        x_val = 0.0
        for x in range(mesh.cells_x-1):
            x_val += mesh.widths[x]*bit_x
            draw.line((x_val, 0, x_val, bit_size), fill=(0,0,0), width=1)

        # draw horizontal grid lines
        y_val = 0.0
        for y in range(mesh.cells_y-1):
            y_val += mesh.heights[y]*bit_y
            draw.line((0, y_val, bit_size, y_val), fill=(0,0,0), width=1)
        
    # save image
    img.save('material.png')
    
    
def plotFluxProb(solvers):
    
    plots = [plt.figure(), plt.figure()]
    
    for j in range(len(solvers)):
    
        solver = solvers[j]
        x = np.zeros(solver.mesh.cells_x * 100)
        ng = solver.mesh.cells[0].material.num_groups
        flux = np.zeros((ng, len(solvers), solver.mesh.cells_x * 100))
    
    
        if j == 3:
            solver.method = 'benchmark'
            flux_sum = np.zeros(ng)
            for xc in range(solver.mesh.cells_x):
                for e in range(ng):
                    flux_sum[e] += solver.phi[xc*ng+e]
          
            for xc in range(solver.mesh.cells_x):
                for e in range(ng):
                    solver.phi[xc*ng+e] = solver.phi[xc*ng+e] * solver.mesh.cells[0].width / 10.0
    
        for xc in range(solver.mesh.cells_x):
            cell = solver.mesh.cells[xc]
            x_val = sum(solver.mesh.widths[0:xc])
        
            for i in range (100):
            
                x[xc * 100 + i] = x[xc * 100 + i - 1] + cell.width / 100.0
            
                for e in range(ng):
                    if solver.method == 'NEM4':
                        flux[e,j, xc * 100 + i] = solver.phi[xc * ng + e] + solver.coeffs[2 * 4 * xc * ng + 2 * 4 * e] * solver.P1(i / 100.0) + solver.coeffs[2 * 4 * xc * ng + 2 * 4 * e + 1] * solver.P2(i / 100.0) + solver.coeffs[2 * 4 * xc * ng + 2 * 4 * e + 2] * solver.P3(i / 100.0) + solver.coeffs[2 * 4 * xc * ng + 2 * 4 * e + 3] * solver.P4(i / 100.0)
                    elif solver.method == 'NEM2':
                        flux[e,j, xc * 100 + i] = solver.phi[xc * ng + e] + solver.coeffs[2 * 2 * xc * ng + 2 * 2 * e] * solver.P1(i / 100.0) + solver.coeffs[2 * 2 * xc * ng + 2 * 2 * e + 1] * solver.P2(i / 100.0) 
                    else:
                        if xc == 0:
                            if i <= 49:
                                flux[e,j, xc * 100 + i] = solver.phi[xc * ng + e]
                            else:
                                flux[e,j, xc * 100 + i] = solver.phi[xc * ng + e] * (150.0 - i) / 100.0 + solver.phi[(xc + 1) * ng + e] * (i - 50.0) / 100.0
                        elif xc == solver.mesh.cells_x - 1:
                            if i >= 50:
                                flux[e,j, xc * 100 + i] = solver.phi[xc * ng + e]
                            else:
                                flux[e,j, xc * 100 + i] = solver.phi[(xc - 1) * ng + e] * (50.0 - i) / 100.0 + solver.phi[xc * ng + e] * (i + 50.0) / 100.0
                        else:
                            if i <= 49:
                                flux[e,j, xc * 100 + i] = solver.phi[(xc - 1) * ng + e] * (50.0 - i) / 100.0 + solver.phi[xc * ng + e] * (i + 50.0) / 100.0
                            else:
                                flux[e,j, xc * 100 + i] = solver.phi[xc * ng + e] * (150.0 - i) / 100.0 + solver.phi[(xc + 1) * ng + e] * (i - 50.0) / 100.0
                            
    
        for e in range(ng):
            plt.figure(plots[e].number)
            plt.plot(x, flux[e,j, :], label=solver.method)

    
    plt.figure(plots[0].number)
    plt.legend()
    plt.savefig('flux_0.png')
    plt.figure(plots[1].number)
    plt.legend(loc=2)
    plt.savefig('flux_1.png')
    
def plotFlux(solver):
    
    plots = []
    for i in range(solver.ng):
        plots.append(plt.figure())

    x = np.zeros(solver.mesh.cells_x * 100)
    ng = solver.mesh.cells[0].material.num_groups
    flux = np.zeros((ng, solver.mesh.cells_x * 100))

    for xc in range(solver.mesh.cells_x):
        cell = solver.mesh.cells[xc]
        x_val = sum(solver.mesh.widths[0:xc])
    
        for i in range (100):
        
            x[xc * 100 + i] = x[xc * 100 + i - 1] + cell.width / 100.0
        
            for e in range(ng):
                if solver.method == 'NEM4':
                    flux[e, xc * 100 + i] = solver.phi[xc * ng + e] + solver.coeffs[2 * 4 * xc * ng + 2 * 4 * e] * solver.P1(i / 100.0) + solver.coeffs[2 * 4 * xc * ng + 2 * 4 * e + 1] * solver.P2(i / 100.0) + solver.coeffs[2 * 4 * xc * ng + 2 * 4 * e + 2] * solver.P3(i / 100.0) + solver.coeffs[2 * 4 * xc * ng + 2 * 4 * e + 3] * solver.P4(i / 100.0)
                elif solver.method == 'NEM2':
                    flux[e, xc * 100 + i] = solver.phi[xc * ng + e] + solver.coeffs[2 * 2 * xc * ng + 2 * 2 * e] * solver.P1(i / 100.0) + solver.coeffs[2 * 2 * xc * ng + 2 * 2 * e + 1] * solver.P2(i / 100.0) 
                else:
                    if xc == 0:
                        if i <= 49:
                            flux[e, xc * 100 + i] = solver.phi[xc * ng + e]
                        else:
                            flux[e, xc * 100 + i] = solver.phi[xc * ng + e] * (150.0 - i) / 100.0 + solver.phi[(xc + 1) * ng + e] * (i - 50.0) / 100.0
                    elif xc == solver.mesh.cells_x - 1:
                        if i >= 50:
                            flux[e, xc * 100 + i] = solver.phi[xc * ng + e]
                        else:
                            flux[e, xc * 100 + i] = solver.phi[(xc - 1) * ng + e] * (50.0 - i) / 100.0 + solver.phi[xc * ng + e] * (i + 50.0) / 100.0
                    else:
                        if i <= 49:
                            flux[e, xc * 100 + i] = solver.phi[(xc - 1) * ng + e] * (50.0 - i) / 100.0 + solver.phi[xc * ng + e] * (i + 50.0) / 100.0
                        else:
                            flux[e, xc * 100 + i] = solver.phi[xc * ng + e] * (150.0 - i) / 100.0 + solver.phi[(xc + 1) * ng + e] * (i - 50.0) / 100.0
                        

    for e in range(ng):
        plt.figure(plots[e].number)
        plt.plot(x, flux[e, :])
        plt.savefig('flux_' + str(e) + '.png')
    
    
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
    
    
def plotCellFlux(solver, grid = 'on', bit_size = 500):

    # create image
    bit_x = float(bit_size) / solver.mesh.width
    bit_y = float(bit_size) / solver.mesh.height

    ng = solver.ng

    imgs = []
    draws = []
    for i in range(solver.ng):
        img = Image.new('RGB', (bit_size,bit_size), 'white')
        imgs.append(img)
        draw = ImageDraw.Draw(img)
        draws.append(draw)
        
    flux_max = np.zeros(ng)
    for y in range(solver.mesh.cells_y):
        for x in range(solver.mesh.cells_x):
            for e in range(ng):
                flux_max[e]= max([flux_max[e], solver.phi[(y*solver.mesh.cells_x+x)*ng+e]])
    
    # draw cells
    for e in range(ng):
        img = imgs[e]
        draw = draws[e]
        y_val = 0.0
        for y in range(solver.mesh.cells_y):
            x_val = 0.0
            for x in range(solver.mesh.cells_x):
            
                red = 0.0
                blue = 0.0
                green = 0.0
            
                val = solver.phi[(y*solver.mesh.cells_x+x)*ng+e]
            
                # get color
                if (float(val) / flux_max[e] <= 1.0/3.0):
                    red = 0.0
                    green = 3.0 * val / flux_max[e]
                    blue = 1.0
                elif (float(val) / flux_max[e] <= 2.0/3.0):
                    red = 3.0 * val / flux_max[e] - 1.0
                    green = 1.0
                    blue = -3.0 * val / flux_max[e] + 2.0
                else:
                    red = 1.0
                    green = -3.0 * val / flux_max[e] + 3.0
                    blue = 0.0
            
                # convert color to RGB triplet
                red = int(255*red)
                green = int(255*green)
                blue = int(255*blue)

                draw.rectangle([x_val, y_val, x_val + solver.mesh.widths[x] * bit_x, y_val + solver.mesh.heights[y] * bit_y], (red,green,blue))
                
                x_val += solver.mesh.widths[x] * bit_x
        
            y_val += solver.mesh.heights[y] * bit_y
                

    # draw vertical grid lines
    for e in range(ng):
        img = imgs[e]
        draw = draws[e]
        if grid == 'on':
            x_val = 0.0
            for x in range(solver.mesh.cells_x-1):
                x_val += solver.mesh.widths[x]*bit_x
                draw.line((x_val, 0, x_val, bit_size), fill=(0,0,0), width=1)
            
            # draw horizontal grid lines
            y_val = 0.0
            for y in range(solver.mesh.cells_y-1):
                y_val += solver.mesh.heights[y]*bit_y
                draw.line((0, y_val, bit_size, y_val), fill=(0,0,0), width=1)
        

    for e in range(ng):
        img = imgs[e]
        img.save('scalar_flux_' + str(e) + '.png')
    
    
    
