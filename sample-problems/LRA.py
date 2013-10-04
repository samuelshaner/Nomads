
# LRA problem with implicit buckling

import sys
sys.path.insert(0, '/Users/sam/git/Nomads/src/')
from Nomads import *
    

tol = 1.e-7    
cell_size    = 15.0
iterations = 50

# create mesh
mesh = Mesh([cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size], 
            [cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size])

mesh.setBoundaries('REFLECTIVE', 'ZERO_FLUX', 'REFLECTIVE', 'ZERO_FLUX')
      
# create fuel 1 blade in
fuel1bin = Material(2, 'fuel 1 blade in')
fuel1bin.setSigmaA([0.008252, 0.1003])
fuel1bin.setD([1.255, 0.211])
fuel1bin.setNuSigmaF([0.004602, 0.1091])
fuel1bin.setChi([1.0, 0.0])
fuel1bin.setSigmaS(np.array([[0.0, 0.02533],[0.0, 0.0]]))
fuel1bin.setBuckling(1.e-4)

# create fuel 1 blade out
fuel1bo = Material(2, 'fuel 1 blade out')
fuel1bo.setSigmaA([0.007181, 0.07047])
fuel1bo.setD([1.268, 0.1902])
fuel1bo.setNuSigmaF([0.004609, 0.08675])
fuel1bo.setChi([1.0, 0.0])
fuel1bo.setSigmaS(np.array([[0.0, 0.02767],[0.0, 0.0]]))
fuel1bo.setBuckling(1.e-4)

# create fuel 2 blade in
fuel2bin = Material(2, 'fuel 2 blade in')
fuel2bin.setSigmaA([0.008002, 0.08344])
fuel2bin.setD([1.259, 0.2091])
fuel2bin.setNuSigmaF([0.004663, 0.1021])
fuel2bin.setChi([1.0, 0.0])
fuel2bin.setSigmaS(np.array([[0.0, 0.02617],[0.0, 0.0]]))
fuel2bin.setBuckling(1.e-4)

# create fuel 2 blade out
fuel2bo = Material(2, 'fuel 2 blade out')
fuel2bo.setSigmaA([0.008002, 0.073324])
fuel2bo.setD([1.259, 0.2091])
fuel2bo.setNuSigmaF([0.004663, 0.1021])
fuel2bo.setChi([1.0, 0.0])
fuel2bo.setSigmaS(np.array([[0.0, 0.02617],[0.0, 0.0]]))
fuel2bo.setBuckling(1.e-4)

# create reflector
reflector = Material(2, 'reflector')
reflector.setSigmaA([0.0006034, 0.01911])
reflector.setD([1.257, 0.1592])
reflector.setNuSigmaF([0.0, 0.0])
reflector.setChi([1.0, 0.0])
reflector.setSigmaS(np.array([[0.0, 0.04754],[0.0, 0.0]]))
reflector.setBuckling(1.e-4)

# create cells
for cell in mesh.cells[0:22]: cell.setMaterial(reflector)
for cell in mesh.cells[22:29]: cell.setMaterial(fuel2bin)
for cell in mesh.cells[29:33]: cell.setMaterial(reflector)
for cell in mesh.cells[33:40]: cell.setMaterial(fuel2bin)
for cell in mesh.cells[40:41]: cell.setMaterial(fuel2bo)
for cell in mesh.cells[41:44]: cell.setMaterial(reflector)
for cell in mesh.cells[44:45]: cell.setMaterial(fuel1bo)
for cell in mesh.cells[45:49]: cell.setMaterial(fuel1bin)
for cell in mesh.cells[49:51]: cell.setMaterial(fuel1bo)
for cell in mesh.cells[51:53]: cell.setMaterial(fuel2bin)
for cell in mesh.cells[53:55]: cell.setMaterial(reflector)
for cell in mesh.cells[55:56]: cell.setMaterial(fuel1bo)
for cell in mesh.cells[56:60]: cell.setMaterial(fuel1bin)
for cell in mesh.cells[60:62]: cell.setMaterial(fuel1bo)
for cell in mesh.cells[62:64]: cell.setMaterial(fuel2bin)
for cell in mesh.cells[64:66]: cell.setMaterial(reflector)
for cell in mesh.cells[66:73]: cell.setMaterial(fuel1bin)
for cell in mesh.cells[73:75]: cell.setMaterial(fuel2bin)
for cell in mesh.cells[75:77]: cell.setMaterial(reflector)
for cell in mesh.cells[77:84]: cell.setMaterial(fuel1bin)
for cell in mesh.cells[84:86]: cell.setMaterial(fuel2bin)
for cell in mesh.cells[86:88]: cell.setMaterial(reflector)
for cell in mesh.cells[88:95]: cell.setMaterial(fuel1bin)
for cell in mesh.cells[95:97]: cell.setMaterial(fuel2bin)
for cell in mesh.cells[97:99]: cell.setMaterial(reflector)
for cell in mesh.cells[99:106]: cell.setMaterial(fuel1bin)
for cell in mesh.cells[106:108]: cell.setMaterial(fuel2bin)
for cell in mesh.cells[108:110]: cell.setMaterial(reflector)
for cell in mesh.cells[110:111]: cell.setMaterial(fuel1bo)
for cell in mesh.cells[111:115]: cell.setMaterial(fuel1bin)
for cell in mesh.cells[115:117]: cell.setMaterial(fuel1bo)
for cell in mesh.cells[117:119]: cell.setMaterial(fuel2bin)
for cell in mesh.cells[119:121]: cell.setMaterial(reflector)

mesh = mesh.refineMesh(1.5)
mesh.makeSurfaces()

# plot the mesh
plotMesh(mesh)

# create solver
solver = Solver(mesh)   

# solve the matrix problem to get flux profile and keff
solver.solve(tol)

plotCellFlux(solver)
