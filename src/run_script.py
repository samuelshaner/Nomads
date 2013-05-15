
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
    opts, args = getopt.getopt(sys.argv[1:], "t:i:f:m:w:h:c:", ["tolerance","iterations","flux_plot", "method", "mesh_width", "mesh_height", "cell_size"])
except getopt.GetoptError, err:
    print str(err)
    sys.exit(2)

# set default values
tol = 1.e-8
plot_flux = False
    
mesh_width   = 2
mesh_height  = 1
cell_size    = 10.0
fuel_width   = 10.0
solve_method = 'NEM4'
iterations = 50

for o, a in opts:
    if o in ("-t", "--tolerance"):
        tol = float(a)
    elif o in ("-f", "--flux_plot"):
        plot_flux = True
    elif o in ("-m", "--method"):
        solve_method = a
    elif o in ("-w", "--mesh_width"):
        mesh_width = int(a)
    elif o in ("-h", "--mesh_height"):
        mesh_height = int(a)
    elif o in ("-c", "--cell_size"):
        cell_size = int(a)
    elif o in ("-i", "--iterations"):
        iterations = int(a)

    else:
        assert False, "unhandled option"

# create mesh
mesh = Mesh([cell_size,cell_size], [cell_size])
      
# create fuel
fuel = Material(2, 'fuel')
fuel.setSigmaA([0.005, 0.10])
fuel.setD([1.5, 0.40])
fuel.setNuSigmaF([0.005, 0.15])
fuel.setChi([1.0, 0.0])
fuel.setSigmaS(np.array([[0.0, 0.02],[0.0, 0.0]]))

# create fuel
moderator = Material(2, 'moderator')
moderator.setSigmaA([0.0, 0.01])
moderator.setD([1.5, 0.20])
moderator.setSigmaS(np.array([[0.0, 0.025],[0.0, 0.0]]))

if solve_method == 'NEM4':
    order = 4
else:
    order = 2
    
mesh.cells[0].setMaterial(fuel, order)
mesh.cells[1].setMaterial(moderator, order)
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



if solve_method == 'NEM4' or solve_method == 'NEM2':
    pttr.plotFlux(solver)
#         pttr.plotCurrent(solver)

stop = time.time()

print 'Ran time ' + str(stop-start)[0:5] + ' seconds'

print '----------------------------------------------------------------------'
