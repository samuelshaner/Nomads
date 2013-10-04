
import Image
import ImageDraw
import ImageFont
from math import *
import matplotlib.pyplot as plt
import numpy as np

def plotMesh(mesh, bit_size = 500):

    cx = mesh.cells_x
    cy = mesh.cells_y
    cz = mesh.cells_z

    for z in range(cz):

        # create image
        bit_x = float(bit_size) / mesh.length_x
        bit_y = float(bit_size) / mesh.length_y
        img = Image.new('RGB', (bit_size,bit_size), 'white')
        draw = ImageDraw.Draw(img)
    
        y_val = 0.0

        # find number of cell types
        num_cell_types = 0
        cell_types = np.array([])
        for y in range(cy):
            for x in range(cx):
                cell = mesh.cells[z*cx*cy+y*cx+x]
            
                if cell.material.id not in cell_types:
                    cell_types = np.append(cell_types, cell.material.id)
                    num_cell_types += 1
    
        # draw cells
        val_max = len(cell_types)
        for y in range(cy):
            x_val = 0.0
            for x in range(cx):
            
                red = 0.0
                blue = 0.0
                green = 0.0
            
                cell = mesh.cells[z*cx*cy+y*cx+x]
            
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

                draw.rectangle([x_val, y_val, x_val + mesh.lengths_x[x] * bit_x, y_val + mesh.lengths_y[y] * bit_y], (red,green,blue))
                
                x_val += mesh.lengths_x[x] * bit_x
        
            y_val += mesh.lengths_y[y] * bit_y
                

        if mesh.cells_x <= 100 and mesh.cells_y <= 100:
            # draw vertical grid lines
            x_val = 0.0
            for x in range(cx-1):
                x_val += mesh.lengths_x[x]*bit_x
                draw.line((x_val, 0, x_val, bit_size), fill=(0,0,0), width=1)

            # draw horizontal grid lines
            y_val = 0.0
            for y in range(cy-1):
                y_val += mesh.lengths_y[y]*bit_y
                draw.line((0, y_val, bit_size, y_val), fill=(0,0,0), width=1)
        
        # save image
        img.save('material_zplane_%04d.png' % z)
    
    
def plotFluxProb(solvers):
    
    plots = [plt.figure(), plt.figure()]
    
    for j in range(len(solvers)):
    
        solver = solvers[j]
        x = np.zeros(solver.mesh.cells_x * 100)
        ng = solver.mesh.cells[0].material.num_groups
        flux = np.zeros((ng, len(solvers), solver.mesh.cells_x * 100))
    
    
        if j == 3:
            flux_sum = np.zeros(ng)
            for xc in range(solver.mesh.cells_x):
                for e in range(ng):
                    flux_sum[e] += solver.phi[xc*ng+e]
          
            for xc in range(solver.mesh.cells_x):
                for e in range(ng):
                    solver.phi[xc*ng+e] = solver.phi[xc*ng+e] * solver.mesh.cells[0].length_x / 10.0
    
        for xc in range(solver.mesh.cells_x):
            cell = solver.mesh.cells[xc]
            x_val = sum(solver.mesh.length_x[0:xc])
        
            for i in range (100):
            
                x[xc * 100 + i] = x[xc * 100 + i - 1] + cell.length_x / 100.0
            
                for e in range(ng):
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
            plt.plot(x, flux[e,j, :])

    
    plt.figure(plots[0].number)
    plt.legend()
    plt.savefig('flux_0.png')
    plt.figure(plots[1].number)
    plt.legend(loc=2)
    plt.savefig('flux_1.png')
    
def plotFlux(solver):
    
    plots = []
    for i in range(solver.num_groups):
        plots.append(plt.figure())

    x = np.zeros(solver.mesh.cells_x * 100)
    ng = solver.mesh.cells[0].material.num_groups
    flux = np.zeros((ng, solver.mesh.cells_x * 100))

    for xc in range(solver.mesh.cells_x):
        cell = solver.mesh.cells[xc]
        x_val = sum(solver.mesh.lengths_x[0:xc])
    
        for i in range (100):
        
            x[xc * 100 + i] = x[xc * 100 + i - 1] + cell.length_x / 100.0
        
            for e in range(ng):
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
        plt.savefig('flux_%03d.png' % e)
   
    
def plotCellFlux(solver, axes=['x','y'], grid = 'on', bit_size = 500):

    cx = solver.mesh.cells_x
    cy = solver.mesh.cells_y
    cz = solver.mesh.cells_z
    ng = solver.num_groups

    if 'x' in axes and 'y' in axes:

        # create image
        bit_x = float(bit_size) / solver.mesh.length_x
        bit_y = float(bit_size) / solver.mesh.length_y

        for z in range(cz):

            imgs = []
            draws = []
            for i in range(solver.num_groups):
                img = Image.new('RGB', (bit_size,bit_size), 'white')
                imgs.append(img)
                draw = ImageDraw.Draw(img)
                draws.append(draw)
        
            flux_max = np.zeros(ng)
            for y in range(cy):
                for x in range(cx):
                    for e in range(ng):
                        flux_max[e]= max([flux_max[e], solver.phi[(z*cx*cy+y*cx+x)*ng+e]])
    
            # draw cells
            for e in range(ng):
                img = imgs[e]
                draw = draws[e]
                y_val = 0.0
                for y in range(cy):
                    x_val = 0.0
                    for x in range(cx):
            
                        red = 0.0
                        blue = 0.0
                        green = 0.0
                        
                        val = solver.phi[(z*cx*cy+y*cx+x)*ng+e]
            
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

                        draw.rectangle([x_val, y_val, x_val + solver.mesh.lengths_x[x] * bit_x, y_val + solver.mesh.lengths_y[y] * bit_y], (red,green,blue))
                
                        x_val += solver.mesh.lengths_x[x] * bit_x
        
                    y_val += solver.mesh.lengths_y[y] * bit_y
                

            # draw vertical grid lines
            for e in range(ng):
                img = imgs[e]
                draw = draws[e]
                if grid == 'on':
                    x_val = 0.0
                    for x in range(cx-1):
                        x_val += solver.mesh.lengths_x[x]*bit_x
                        draw.line((x_val, 0, x_val, bit_size), fill=(0,0,0), width=1)
            
                    # draw horizontal grid lines
                    y_val = 0.0
                    for y in range(cy-1):
                        y_val += solver.mesh.lengths_y[y]*bit_y
                        draw.line((0, y_val, bit_size, y_val), fill=(0,0,0), width=1)
        

            for e in range(ng):
                img = imgs[e]
                img.save('scalar_flux_zplane_%04d_group_%03d.png' % (z,e))


    elif 'x' in axes and 'z' in axes:

        # create image
        bit_x = float(bit_size) / solver.mesh.length_x
        bit_z = float(bit_size) / solver.mesh.length_z

        for y in range(cy):

            imgs = []
            draws = []
            for i in range(solver.num_groups):
                img = Image.new('RGB', (bit_size,bit_size), 'white')
                imgs.append(img)
                draw = ImageDraw.Draw(img)
                draws.append(draw)
        
            flux_max = np.zeros(ng)
            for z in range(cz):
                for x in range(cx):
                    for e in range(ng):
                        flux_max[e]= max([flux_max[e], solver.phi[(z*cx*cy+y*cx+x)*ng+e]])
    
            # draw cells
            for e in range(ng):
                img = imgs[e]
                draw = draws[e]
                z_val = 0.0
                for z in range(cz):
                    x_val = 0.0
                    for x in range(cx):
            
                        red = 0.0
                        blue = 0.0
                        green = 0.0
                        
                        val = solver.phi[(z*cx*cy+y*cx+x)*ng+e]
            
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

                        draw.rectangle([x_val, z_val, x_val + solver.mesh.lengths_x[x] * bit_x, z_val + solver.mesh.lengths_z[z] * bit_z], (red,green,blue))
                
                        x_val += solver.mesh.lengths_x[x] * bit_x
        
                    z_val += solver.mesh.lengths_z[z] * bit_z
                

            # draw vertical grid lines
            for e in range(ng):
                img = imgs[e]
                draw = draws[e]
                if grid == 'on':
                    x_val = 0.0
                    for x in range(cx-1):
                        x_val += solver.mesh.lengths_x[x]*bit_x
                        draw.line((x_val, 0, x_val, bit_size), fill=(0,0,0), width=1)
            
                    # draw horizontal grid lines
                    z_val = 0.0
                    for z in range(cz-1):
                        z_val += solver.mesh.lengths_z[z]*bit_z
                        draw.line((0, z_val, bit_size, z_val), fill=(0,0,0), width=1)
        

            for e in range(ng):
                img = imgs[e]
                img.save('scalar_flux_yplane_%04d_group_%03d.png' % (y,e))
    
    