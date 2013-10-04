
# LRA problem with implicit buckling

import sys
sys.path.insert(0, '/Users/sam/git/Nomads/src/')
from Nomads import *
    

tol = 1.e-6    
cell_size    = 15.0
fuel_width   = 10.0
iterations = 50

# create mesh
mesh = Mesh([cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size], 
            [cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size],
            [30.0, 300.0, 30.0])

mesh.setBoundaries('REFLECTIVE', 'ZERO_FLUX', 'REFLECTIVE', 'ZERO_FLUX', 'VACUUM', 'VACUUM')
      
# create fuel 1 blade in
fuel1bin = Material(2, 'fuel 1 blade in')
fuel1bin.setSigmaA([0.008252, 0.1003])
fuel1bin.setD([1.255, 0.211])
fuel1bin.setNuSigmaF([0.004602, 0.1091])
fuel1bin.setChi([1.0, 0.0])
fuel1bin.setSigmaS(np.array([[0.0, 0.02533],[0.0, 0.0]]))

# create fuel 1 blade out
fuel1bo = Material(2, 'fuel 1 blade out')
fuel1bo.setSigmaA([0.007181, 0.07047])
fuel1bo.setD([1.268, 0.1902])
fuel1bo.setNuSigmaF([0.004609, 0.08675])
fuel1bo.setChi([1.0, 0.0])
fuel1bo.setSigmaS(np.array([[0.0, 0.02767],[0.0, 0.0]]))

# create fuel 2 blade in
fuel2bin = Material(2, 'fuel 2 blade in')
fuel2bin.setSigmaA([0.008002, 0.08344])
fuel2bin.setD([1.259, 0.2091])
fuel2bin.setNuSigmaF([0.004663, 0.1021])
fuel2bin.setChi([1.0, 0.0])
fuel2bin.setSigmaS(np.array([[0.0, 0.02617],[0.0, 0.0]]))

# create fuel 2 blade out
fuel2bo = Material(2, 'fuel 2 blade out')
fuel2bo.setSigmaA([0.008002, 0.073324])
fuel2bo.setD([1.259, 0.2091])
fuel2bo.setNuSigmaF([0.004663, 0.1021])
fuel2bo.setChi([1.0, 0.0])
fuel2bo.setSigmaS(np.array([[0.0, 0.02617],[0.0, 0.0]]))

# create reflector
reflector = Material(2, 'reflector')
reflector.setSigmaA([0.0006034, 0.01911])
reflector.setD([1.257, 0.1592])
reflector.setNuSigmaF([0.0, 0.0])
reflector.setChi([1.0, 0.0])
reflector.setSigmaS(np.array([[0.0, 0.04754],[0.0, 0.0]]))

# create cells in layer 0
for cell in mesh.cells[0:121]: cell.setMaterial(reflector)

# create cells in layer 1
for cell in mesh.cells[121:143]: cell.setMaterial(reflector)
for cell in mesh.cells[143:150]: cell.setMaterial(fuel2bin)
for cell in mesh.cells[150:154]: cell.setMaterial(reflector)
for cell in mesh.cells[154:161]: cell.setMaterial(fuel2bin)
mesh.cells[161].setMaterial(fuel2bo)
for cell in mesh.cells[162:165]: cell.setMaterial(reflector)
mesh.cells[165].setMaterial(fuel1bo)
for cell in mesh.cells[166:170]: cell.setMaterial(fuel1bin)
for cell in mesh.cells[170:172]: cell.setMaterial(fuel1bo)
for cell in mesh.cells[172:174]: cell.setMaterial(fuel2bin)
for cell in mesh.cells[174:176]: cell.setMaterial(reflector)
mesh.cells[176].setMaterial(fuel1bo)
for cell in mesh.cells[177:181]: cell.setMaterial(fuel1bin)
for cell in mesh.cells[181:183]: cell.setMaterial(fuel1bo)
for cell in mesh.cells[183:185]: cell.setMaterial(fuel2bin)
for cell in mesh.cells[185:187]: cell.setMaterial(reflector)
for cell in mesh.cells[187:194]: cell.setMaterial(fuel1bin)
for cell in mesh.cells[194:196]: cell.setMaterial(fuel2bin)
for cell in mesh.cells[196:198]: cell.setMaterial(reflector)
for cell in mesh.cells[198:205]: cell.setMaterial(fuel1bin)
for cell in mesh.cells[205:207]: cell.setMaterial(fuel2bin)
for cell in mesh.cells[207:209]: cell.setMaterial(reflector)
for cell in mesh.cells[209:216]: cell.setMaterial(fuel1bin)
for cell in mesh.cells[216:218]: cell.setMaterial(fuel2bin)
for cell in mesh.cells[218:220]: cell.setMaterial(reflector)
for cell in mesh.cells[220:227]: cell.setMaterial(fuel1bin)
for cell in mesh.cells[227:229]: cell.setMaterial(fuel2bin)
for cell in mesh.cells[229:231]: cell.setMaterial(reflector)
mesh.cells[231].setMaterial(fuel1bo)
for cell in mesh.cells[231:236]: cell.setMaterial(fuel1bin)
for cell in mesh.cells[236:238]: cell.setMaterial(fuel1bo)
for cell in mesh.cells[238:240]: cell.setMaterial(fuel2bin)
for cell in mesh.cells[240:242]: cell.setMaterial(reflector)

# create cells in layer 2
for cell in mesh.cells[242:363]: cell.setMaterial(reflector)

mesh = mesh.refineMesh(7.5, ['x','y'])
mesh = mesh.refineMesh(30.0, ['z'])
mesh.makeSurfaces()

# plot the mesh
plotMesh(mesh)

# create solver
solver = Solver(mesh)   

# solve the matrix problem to get flux profile and keff
solver.solve(tol)

plotCellFlux(solver, ['x','y'])
