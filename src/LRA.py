
# LRA problem with implicit buckling

from Nomads import *
from Cell import *
from Material import *
from Mesh import *
import plotter as pttr
import sys
import getopt
import time

start = time.time()
    
# parse command line options
try:
    opts, args = getopt.getopt(sys.argv[1:], "t:i:m:", ["tolerance","iterations","method"])
except getopt.GetoptError, err:
    print str(err)
    sys.exit(2)

# set default values
tol = 1.e-8
    
cell_size    = 15.0
fuel_width   = 10.0
solve_method = 'NEM4'
iterations = 50

for o, a in opts:
    if o in ("-t", "--tolerance"):
        tol = float(a)
    elif o in ("-m", "--method"):
        solve_method = a
    elif o in ("-i", "--iterations"):
        iterations = int(a)

    else:
        assert False, "unhandled option"

# create mesh
mesh = Mesh([cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size], 
            [cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size])
      
# create fuel 1 blade in
fuel1bin = Material(2, 'fuel 1 blade in')
fuel1bin.setSigmaA([0.0083775, 0.1003211])
fuel1bin.setD([1.255, 0.211])
fuel1bin.setNuSigmaF([0.004602, 0.1091])
fuel1bin.setChi([1.0, 0.0])
fuel1bin.setSigmaS(np.array([[0.0, 0.02533],[0.0, 0.0]]))

# create fuel 1 blade out
fuel1bo = Material(2, 'fuel 1 blade out')
fuel1bo.setSigmaA([0.0073078, 0.07048902])
fuel1bo.setD([1.268, 0.1902])
fuel1bo.setNuSigmaF([0.004609, 0.08675])
fuel1bo.setChi([1.0, 0.0])
fuel1bo.setSigmaS(np.array([[0.0, 0.02767],[0.0, 0.0]]))

# create fuel 2 blade in
fuel2bin = Material(2, 'fuel 2 blade in')
fuel2bin.setSigmaA([0.0081279, 0.08346091])
fuel2bin.setD([1.259, 0.2091])
fuel2bin.setNuSigmaF([0.004663, 0.1021])
fuel2bin.setChi([1.0, 0.0])
fuel2bin.setSigmaS(np.array([[0.0, 0.02617],[0.0, 0.0]]))

# create fuel 2 blade out
fuel2bo = Material(2, 'fuel 2 blade out')
fuel2bo.setSigmaA([0.0081279, 0.07334491])
fuel2bo.setD([1.259, 0.2091])
fuel2bo.setNuSigmaF([0.004663, 0.1021])
fuel2bo.setChi([1.0, 0.0])
fuel2bo.setSigmaS(np.array([[0.0, 0.02617],[0.0, 0.0]]))

# create reflector
reflector = Material(2, 'reflector')
reflector.setSigmaA([0.0007291, 0.01912592])
reflector.setD([1.257, 0.1592])
reflector.setNuSigmaF([0.0, 0.0])
reflector.setChi([1.0, 0.0])
reflector.setSigmaS(np.array([[0.0, 0.04754],[0.0, 0.0]]))

if solve_method == 'NEM4':
    order = 4
else:
    order = 2

# create cells
for cell in mesh.cells[0:22]: cell.setMaterial(reflector,order)
for cell in mesh.cells[22:29]: cell.setMaterial(fuel2bin,order)
for cell in mesh.cells[29:33]: cell.setMaterial(reflector,order)
for cell in mesh.cells[33:40]: cell.setMaterial(fuel2bin,order)
mesh.cells[40].setMaterial(fuel2bo,order)
for cell in mesh.cells[41:44]: cell.setMaterial(reflector,order)
mesh.cells[44].setMaterial(fuel1bo,order)
for cell in mesh.cells[45:49]: cell.setMaterial(fuel1bin,order)
for cell in mesh.cells[49:51]: cell.setMaterial(fuel1bo,order)
for cell in mesh.cells[51:53]: cell.setMaterial(fuel2bin,order)
for cell in mesh.cells[53:55]: cell.setMaterial(reflector,order)
mesh.cells[55].setMaterial(fuel1bo,order)
for cell in mesh.cells[56:60]: cell.setMaterial(fuel1bin,order)
for cell in mesh.cells[60:62]: cell.setMaterial(fuel1bo,order)
for cell in mesh.cells[62:64]: cell.setMaterial(fuel2bin,order)
for cell in mesh.cells[64:66]: cell.setMaterial(reflector,order)
for cell in mesh.cells[66:73]: cell.setMaterial(fuel1bin,order)
for cell in mesh.cells[73:75]: cell.setMaterial(fuel2bin,order)
for cell in mesh.cells[75:77]: cell.setMaterial(reflector,order)
for cell in mesh.cells[77:84]: cell.setMaterial(fuel1bin,order)
for cell in mesh.cells[84:86]: cell.setMaterial(fuel2bin,order)
for cell in mesh.cells[86:88]: cell.setMaterial(reflector,order)
for cell in mesh.cells[88:95]: cell.setMaterial(fuel1bin,order)
for cell in mesh.cells[95:97]: cell.setMaterial(fuel2bin,order)
for cell in mesh.cells[97:99]: cell.setMaterial(reflector,order)
for cell in mesh.cells[99:106]: cell.setMaterial(fuel1bin,order)
for cell in mesh.cells[106:108]: cell.setMaterial(fuel2bin,order)
for cell in mesh.cells[108:110]: cell.setMaterial(reflector,order)
mesh.cells[110].setMaterial(fuel1bo,order)
for cell in mesh.cells[111:115]: cell.setMaterial(fuel1bin,order)
for cell in mesh.cells[115:117]: cell.setMaterial(fuel1bo,order)
for cell in mesh.cells[117:119]: cell.setMaterial(fuel2bin,order)
for cell in mesh.cells[119:121]: cell.setMaterial(reflector,order)

mesh.makeSurfaces()

# plot the mesh
pttr.plotMesh(mesh)

# create solver
solver = Solver(mesh, solve_method)   

# solve the matrix problem to get flux profile and keff

if solve_method == 'NEM4' or solve_method == 'NEM2':
    for iteration in range(iterations):
    
        print 'CMFD outer iteration ' + str(iteration)
    
        solver.computeDs()
        solver.makeAM()
        solver.solve(tol)
        solver.makeN()
        solver.computeCoeffs()
        solver.computeCurrents()
    
        if abs(solver.keff_old - solver.keff) < 1.e-8:
            print solve_method + ': Converged in ' + str(iteration) + ' iterations --- k_eff = ' + str(solver.keff)[0:10]
            break
        
        
elif solve_method == 'diffusion':
    solver.computeDs()
    solver.makeAM()
    solver.solve(tol)
    print 'DIFFUSION: --- k_eff = ' + str(solver.keff)[0:10]



pttr.plotCellFlux(solver)
if solve_method == 'NEM4' or solve_method == 'NEM2':
    pttr.plotFlux(solver)
#         pttr.plotCurrent(solver)

stop = time.time()

print 'Ran time ' + str(stop-start)[0:5] + ' seconds'

print '----------------------------------------------------------------------'
