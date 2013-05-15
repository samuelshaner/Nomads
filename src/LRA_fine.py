
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
    
cell_size    = 5.0
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
mesh = Mesh([cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,
             cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,
             cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size], 
            [cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,
             cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,
             cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size,cell_size])
      
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
# row 1 and 2
for cell in mesh.cells[0:198]: cell.setMaterial(reflector,order)
# row 3
for cell in mesh.cells[198:219]: cell.setMaterial(fuel2bin,order)
for cell in mesh.cells[219:231]: cell.setMaterial(reflector,order)
for cell in mesh.cells[231:252]: cell.setMaterial(fuel2bin,order)
for cell in mesh.cells[252:264]: cell.setMaterial(reflector,order)
for cell in mesh.cells[264:285]: cell.setMaterial(fuel2bin,order)
for cell in mesh.cells[285:297]: cell.setMaterial(reflector,order)
# row 4
for cell in mesh.cells[297:318]: cell.setMaterial(fuel2bin,order)
for cell in mesh.cells[318:321]: cell.setMaterial(fuel2bo,order)
for cell in mesh.cells[321:330]: cell.setMaterial(reflector,order)
for cell in mesh.cells[330:351]: cell.setMaterial(fuel2bin,order)
for cell in mesh.cells[351:354]: cell.setMaterial(fuel2bo,order)
for cell in mesh.cells[354:363]: cell.setMaterial(reflector,order)
for cell in mesh.cells[363:384]: cell.setMaterial(fuel2bin,order)
for cell in mesh.cells[384:387]: cell.setMaterial(fuel2bo,order)
for cell in mesh.cells[387:396]: cell.setMaterial(reflector,order)
# row 5
for cell in mesh.cells[396:399]: cell.setMaterial(fuel1bo,order)
for cell in mesh.cells[399:411]: cell.setMaterial(fuel1bin,order)
for cell in mesh.cells[411:417]: cell.setMaterial(fuel1bo,order)
for cell in mesh.cells[417:423]: cell.setMaterial(fuel2bin,order)
for cell in mesh.cells[423:429]: cell.setMaterial(reflector,order)
for cell in mesh.cells[429:432]: cell.setMaterial(fuel1bo,order)
for cell in mesh.cells[432:444]: cell.setMaterial(fuel1bin,order)
for cell in mesh.cells[444:450]: cell.setMaterial(fuel1bo,order)
for cell in mesh.cells[450:456]: cell.setMaterial(fuel2bin,order)
for cell in mesh.cells[456:462]: cell.setMaterial(reflector,order)
for cell in mesh.cells[462:465]: cell.setMaterial(fuel1bo,order)
for cell in mesh.cells[465:477]: cell.setMaterial(fuel1bin,order)
for cell in mesh.cells[477:483]: cell.setMaterial(fuel1bo,order)
for cell in mesh.cells[483:489]: cell.setMaterial(fuel2bin,order)
for cell in mesh.cells[489:495]: cell.setMaterial(reflector,order)

# row 6
for cell in mesh.cells[495:498]: cell.setMaterial(fuel1bo,order)
for cell in mesh.cells[498:510]: cell.setMaterial(fuel1bin,order)
for cell in mesh.cells[510:516]: cell.setMaterial(fuel1bo,order)
for cell in mesh.cells[516:522]: cell.setMaterial(fuel2bin,order)
for cell in mesh.cells[522:528]: cell.setMaterial(reflector,order)
for cell in mesh.cells[528:531]: cell.setMaterial(fuel1bo,order)
for cell in mesh.cells[531:543]: cell.setMaterial(fuel1bin,order)
for cell in mesh.cells[543:549]: cell.setMaterial(fuel1bo,order)
for cell in mesh.cells[549:555]: cell.setMaterial(fuel2bin,order)
for cell in mesh.cells[555:561]: cell.setMaterial(reflector,order)
for cell in mesh.cells[561:564]: cell.setMaterial(fuel1bo,order)
for cell in mesh.cells[564:576]: cell.setMaterial(fuel1bin,order)
for cell in mesh.cells[576:582]: cell.setMaterial(fuel1bo,order)
for cell in mesh.cells[582:588]: cell.setMaterial(fuel2bin,order)
for cell in mesh.cells[588:594]: cell.setMaterial(reflector,order)


# row 7
for cell in mesh.cells[594:615]: cell.setMaterial(fuel1bin,order)
for cell in mesh.cells[615:621]: cell.setMaterial(fuel2bin,order)
for cell in mesh.cells[621:627]: cell.setMaterial(reflector,order)
for cell in mesh.cells[627:648]: cell.setMaterial(fuel1bin,order)
for cell in mesh.cells[648:654]: cell.setMaterial(fuel2bin,order)
for cell in mesh.cells[654:660]: cell.setMaterial(reflector,order)
for cell in mesh.cells[660:681]: cell.setMaterial(fuel1bin,order)
for cell in mesh.cells[681:687]: cell.setMaterial(fuel2bin,order)
for cell in mesh.cells[687:693]: cell.setMaterial(reflector,order)

# row 8
for cell in mesh.cells[693:714]: cell.setMaterial(fuel1bin,order)
for cell in mesh.cells[714:720]: cell.setMaterial(fuel2bin,order)
for cell in mesh.cells[720:726]: cell.setMaterial(reflector,order)
for cell in mesh.cells[726:747]: cell.setMaterial(fuel1bin,order)
for cell in mesh.cells[747:753]: cell.setMaterial(fuel2bin,order)
for cell in mesh.cells[753:759]: cell.setMaterial(reflector,order)
for cell in mesh.cells[759:780]: cell.setMaterial(fuel1bin,order)
for cell in mesh.cells[780:786]: cell.setMaterial(fuel2bin,order)
for cell in mesh.cells[786:792]: cell.setMaterial(reflector,order)

# row 9
for cell in mesh.cells[792:813]: cell.setMaterial(fuel1bin,order)
for cell in mesh.cells[813:819]: cell.setMaterial(fuel2bin,order)
for cell in mesh.cells[819:825]: cell.setMaterial(reflector,order)
for cell in mesh.cells[825:846]: cell.setMaterial(fuel1bin,order)
for cell in mesh.cells[846:852]: cell.setMaterial(fuel2bin,order)
for cell in mesh.cells[852:858]: cell.setMaterial(reflector,order)
for cell in mesh.cells[858:879]: cell.setMaterial(fuel1bin,order)
for cell in mesh.cells[879:885]: cell.setMaterial(fuel2bin,order)
for cell in mesh.cells[885:891]: cell.setMaterial(reflector,order)

# row 10
for cell in mesh.cells[891:912]: cell.setMaterial(fuel1bin,order)
for cell in mesh.cells[912:918]: cell.setMaterial(fuel2bin,order)
for cell in mesh.cells[918:924]: cell.setMaterial(reflector,order)
for cell in mesh.cells[924:945]: cell.setMaterial(fuel1bin,order)
for cell in mesh.cells[945:951]: cell.setMaterial(fuel2bin,order)
for cell in mesh.cells[951:957]: cell.setMaterial(reflector,order)
for cell in mesh.cells[957:978]: cell.setMaterial(fuel1bin,order)
for cell in mesh.cells[978:984]: cell.setMaterial(fuel2bin,order)
for cell in mesh.cells[984:990]: cell.setMaterial(reflector,order)

# row 11
for cell in mesh.cells[990:993]: cell.setMaterial(fuel1bo,order)
for cell in mesh.cells[993:1005]: cell.setMaterial(fuel1bin,order)
for cell in mesh.cells[1005:1011]: cell.setMaterial(fuel1bo,order)
for cell in mesh.cells[1011:1017]: cell.setMaterial(fuel2bin,order)
for cell in mesh.cells[1017:1023]: cell.setMaterial(reflector,order)
for cell in mesh.cells[1023:1026]: cell.setMaterial(fuel1bo,order)
for cell in mesh.cells[1026:1038]: cell.setMaterial(fuel1bin,order)
for cell in mesh.cells[1038:1044]: cell.setMaterial(fuel1bo,order)
for cell in mesh.cells[1044:1050]: cell.setMaterial(fuel2bin,order)
for cell in mesh.cells[1050:1056]: cell.setMaterial(reflector,order)
for cell in mesh.cells[1056:1059]: cell.setMaterial(fuel1bo,order)
for cell in mesh.cells[1059:1071]: cell.setMaterial(fuel1bin,order)
for cell in mesh.cells[1071:1077]: cell.setMaterial(fuel1bo,order)
for cell in mesh.cells[1077:1083]: cell.setMaterial(fuel2bin,order)
for cell in mesh.cells[1083:1089]: cell.setMaterial(reflector,order)

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
