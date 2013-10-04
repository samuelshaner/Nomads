
# IAEA problem in 3D

import sys
sys.path.insert(0, '/Users/sam/git/Nomads/src/')
from Nomads import *
    
tol = 1.e-5
dx = 10

# create mesh
mesh = Mesh([dx, dx, dx, dx, dx, dx, dx, dx, dx, dx, dx, dx, dx, dx, dx, dx, dx],
            [dx, dx, dx, dx, dx, dx, dx, dx, dx, dx, dx, dx, dx, dx, dx, dx, dx],
            [20.0, 260.0, 80.0, 20])

mesh.setBoundaries('REFLECTIVE', 'VACUUM', 'REFLECTIVE', 'VACUUM', 'VACUUM', 'VACUUM')
      
# create fuel 1 blade in
mat1 = Material(2, 'material 1')
mat1.setSigmaA([0.01, 0.08])
mat1.setD([1.5, 0.4])
mat1.setNuSigmaF([0.0, 0.135])
mat1.setChi([1.0, 0.0])
mat1.setSigmaS(np.array([[0.0, 0.02],[0.0, 0.0]]))

# create fuel 1 blade out
mat2 = Material(2, 'material 2')
mat2.setSigmaA([0.01, 0.085])
mat2.setD([1.5, 0.4])
mat2.setNuSigmaF([0.0, 0.135])
mat2.setChi([1.0, 0.0])
mat2.setSigmaS(np.array([[0.0, 0.02],[0.0, 0.0]]))

# create fuel 2 blade in
mat3 = Material(2, 'material 3')
mat3.setSigmaA([0.01, 0.13])
mat3.setD([1.5, 0.4])
mat3.setNuSigmaF([0.0, 0.135])
mat3.setChi([1.0, 0.0])
mat3.setSigmaS(np.array([[0.0, 0.02],[0.0, 0.0]]))

# create fuel 2 blade out
mat4 = Material(2, 'material 4')
mat4.setSigmaA([0.0, 0.01])
mat4.setD([2.0, 0.3])
mat4.setNuSigmaF([0.0, 0.0])
mat4.setChi([1.0, 0.0])
mat4.setSigmaS(np.array([[0.0, 0.04],[0.0, 0.0]]))

# create reflector
mat5 = Material(2, 'material 5')
mat5.setSigmaA([0.0, 0.055])
mat5.setD([2.0, 0.3])
mat5.setNuSigmaF([0.0, 0.0])
mat5.setChi([1.0, 0.0])
mat5.setSigmaS(np.array([[0.0, 0.04],[0.0, 0.0]]))

# create cells in layer 0
for cell in mesh.cells[0:289]: cell.setMaterial(mat4)

# create cells in layer 1
for cell in mesh.cells[289:306]: cell.setMaterial(mat4)
for cell in mesh.cells[306:323]: cell.setMaterial(mat4)

for cell in mesh.cells[323:328]: cell.setMaterial(mat1)
for cell in mesh.cells[328:340]: cell.setMaterial(mat4)
for cell in mesh.cells[340:345]: cell.setMaterial(mat1)
for cell in mesh.cells[345:357]: cell.setMaterial(mat4)

for cell in mesh.cells[357:360]: cell.setMaterial(mat2)
for cell in mesh.cells[360:366]: cell.setMaterial(mat1)
for cell in mesh.cells[366:374]: cell.setMaterial(mat4)
for cell in mesh.cells[374:377]: cell.setMaterial(mat2)
for cell in mesh.cells[377:383]: cell.setMaterial(mat1)
for cell in mesh.cells[383:391]: cell.setMaterial(mat4)

for cell in mesh.cells[391:398]: cell.setMaterial(mat2)
for cell in mesh.cells[398:402]: cell.setMaterial(mat1)
for cell in mesh.cells[402:408]: cell.setMaterial(mat4)
for cell in mesh.cells[408:415]: cell.setMaterial(mat2)
for cell in mesh.cells[415:419]: cell.setMaterial(mat1)
for cell in mesh.cells[419:425]: cell.setMaterial(mat4)

for cell in mesh.cells[425:426]: cell.setMaterial(mat3)
for cell in mesh.cells[426:432]: cell.setMaterial(mat2)
for cell in mesh.cells[432:434]: cell.setMaterial(mat3)
for cell in mesh.cells[434:438]: cell.setMaterial(mat1)
for cell in mesh.cells[438:442]: cell.setMaterial(mat4)
for cell in mesh.cells[442:443]: cell.setMaterial(mat3)
for cell in mesh.cells[443:449]: cell.setMaterial(mat2)
for cell in mesh.cells[449:451]: cell.setMaterial(mat3)
for cell in mesh.cells[451:455]: cell.setMaterial(mat1)
for cell in mesh.cells[455:459]: cell.setMaterial(mat4)

for cell in mesh.cells[459:470]: cell.setMaterial(mat2)
for cell in mesh.cells[470:472]: cell.setMaterial(mat1)
for cell in mesh.cells[472:476]: cell.setMaterial(mat4)
for cell in mesh.cells[476:487]: cell.setMaterial(mat2)
for cell in mesh.cells[487:489]: cell.setMaterial(mat1)
for cell in mesh.cells[489:493]: cell.setMaterial(mat4)

for cell in mesh.cells[493:504]: cell.setMaterial(mat2)
for cell in mesh.cells[504:508]: cell.setMaterial(mat1)
for cell in mesh.cells[508:510]: cell.setMaterial(mat4)
for cell in mesh.cells[510:521]: cell.setMaterial(mat2)
for cell in mesh.cells[521:525]: cell.setMaterial(mat1)
for cell in mesh.cells[525:527]: cell.setMaterial(mat4)

for cell in mesh.cells[527:540]: cell.setMaterial(mat2)
for cell in mesh.cells[540:542]: cell.setMaterial(mat1)
for cell in mesh.cells[542:544]: cell.setMaterial(mat4)
for cell in mesh.cells[544:557]: cell.setMaterial(mat2)
for cell in mesh.cells[557:559]: cell.setMaterial(mat1)
for cell in mesh.cells[559:561]: cell.setMaterial(mat4)

for cell in mesh.cells[561:562]: cell.setMaterial(mat3)
for cell in mesh.cells[562:568]: cell.setMaterial(mat2)
for cell in mesh.cells[568:570]: cell.setMaterial(mat3)
for cell in mesh.cells[570:574]: cell.setMaterial(mat2)
for cell in mesh.cells[574:576]: cell.setMaterial(mat1)
for cell in mesh.cells[576:578]: cell.setMaterial(mat4)

# create cells in layer 2
for cell in mesh.cells[578:595]: cell.setMaterial(mat4)
for cell in mesh.cells[595:612]: cell.setMaterial(mat4)

for cell in mesh.cells[612:617]: cell.setMaterial(mat1)
for cell in mesh.cells[617:629]: cell.setMaterial(mat4)
for cell in mesh.cells[629:634]: cell.setMaterial(mat1)
for cell in mesh.cells[634:646]: cell.setMaterial(mat4)

for cell in mesh.cells[646:649]: cell.setMaterial(mat2)
for cell in mesh.cells[649:655]: cell.setMaterial(mat1)
for cell in mesh.cells[655:663]: cell.setMaterial(mat4)
for cell in mesh.cells[663:666]: cell.setMaterial(mat2)
for cell in mesh.cells[666:672]: cell.setMaterial(mat1)
for cell in mesh.cells[672:680]: cell.setMaterial(mat4)

for cell in mesh.cells[680:687]: cell.setMaterial(mat2)
for cell in mesh.cells[687:691]: cell.setMaterial(mat1)
for cell in mesh.cells[691:697]: cell.setMaterial(mat4)
for cell in mesh.cells[697:704]: cell.setMaterial(mat2)
for cell in mesh.cells[704:708]: cell.setMaterial(mat1)
for cell in mesh.cells[708:714]: cell.setMaterial(mat4)

for cell in mesh.cells[714:715]: cell.setMaterial(mat3)
for cell in mesh.cells[715:721]: cell.setMaterial(mat2)
for cell in mesh.cells[721:723]: cell.setMaterial(mat3)
for cell in mesh.cells[723:727]: cell.setMaterial(mat1)
for cell in mesh.cells[727:731]: cell.setMaterial(mat4)
for cell in mesh.cells[731:732]: cell.setMaterial(mat3)
for cell in mesh.cells[732:738]: cell.setMaterial(mat2)
for cell in mesh.cells[738:740]: cell.setMaterial(mat3)
for cell in mesh.cells[740:744]: cell.setMaterial(mat1)
for cell in mesh.cells[744:748]: cell.setMaterial(mat4)

for cell in mesh.cells[748:759]: cell.setMaterial(mat2)
for cell in mesh.cells[759:761]: cell.setMaterial(mat1)
for cell in mesh.cells[761:765]: cell.setMaterial(mat4)
for cell in mesh.cells[765:776]: cell.setMaterial(mat2)
for cell in mesh.cells[776:778]: cell.setMaterial(mat1)
for cell in mesh.cells[778:782]: cell.setMaterial(mat4)

for cell in mesh.cells[782:785]: cell.setMaterial(mat2)
for cell in mesh.cells[785:787]: cell.setMaterial(mat3)
for cell in mesh.cells[787:793]: cell.setMaterial(mat2)
for cell in mesh.cells[793:797]: cell.setMaterial(mat1)
for cell in mesh.cells[797:799]: cell.setMaterial(mat4)
for cell in mesh.cells[799:802]: cell.setMaterial(mat2)
for cell in mesh.cells[802:804]: cell.setMaterial(mat3)
for cell in mesh.cells[804:810]: cell.setMaterial(mat2)
for cell in mesh.cells[810:814]: cell.setMaterial(mat1)
for cell in mesh.cells[814:816]: cell.setMaterial(mat4)

for cell in mesh.cells[816:829]: cell.setMaterial(mat2)
for cell in mesh.cells[829:831]: cell.setMaterial(mat1)
for cell in mesh.cells[831:833]: cell.setMaterial(mat4)
for cell in mesh.cells[833:846]: cell.setMaterial(mat2)
for cell in mesh.cells[846:848]: cell.setMaterial(mat1)
for cell in mesh.cells[848:850]: cell.setMaterial(mat4)

for cell in mesh.cells[850:851]: cell.setMaterial(mat3)
for cell in mesh.cells[851:857]: cell.setMaterial(mat2)
for cell in mesh.cells[857:859]: cell.setMaterial(mat3)
for cell in mesh.cells[859:863]: cell.setMaterial(mat2)
for cell in mesh.cells[863:865]: cell.setMaterial(mat1)
for cell in mesh.cells[865:867]: cell.setMaterial(mat4)


# create cells in layer 3
for cell in mesh.cells[867:1003]: cell.setMaterial(mat4)

for cell in mesh.cells[1003:1004]: cell.setMaterial(mat5)
for cell in mesh.cells[1004:1010]: cell.setMaterial(mat4)
for cell in mesh.cells[1010:1012]: cell.setMaterial(mat5)
for cell in mesh.cells[1012:1020]: cell.setMaterial(mat4)
for cell in mesh.cells[1020:1021]: cell.setMaterial(mat5)
for cell in mesh.cells[1021:1027]: cell.setMaterial(mat4)
for cell in mesh.cells[1027:1029]: cell.setMaterial(mat5)
for cell in mesh.cells[1029:1037]: cell.setMaterial(mat4)

for cell in mesh.cells[1037:1071]: cell.setMaterial(mat4)

for cell in mesh.cells[1071:1074]: cell.setMaterial(mat4)
for cell in mesh.cells[1074:1076]: cell.setMaterial(mat5)
for cell in mesh.cells[1076:1088]: cell.setMaterial(mat4)
for cell in mesh.cells[1088:1091]: cell.setMaterial(mat4)
for cell in mesh.cells[1091:1093]: cell.setMaterial(mat5)
for cell in mesh.cells[1093:1115]: cell.setMaterial(mat4)

for cell in mesh.cells[1105:1139]: cell.setMaterial(mat4)

for cell in mesh.cells[1139:1140]: cell.setMaterial(mat5)
for cell in mesh.cells[1140:1146]: cell.setMaterial(mat4)
for cell in mesh.cells[1146:1148]: cell.setMaterial(mat5)
for cell in mesh.cells[1148:1156]: cell.setMaterial(mat4)


#mesh = mesh.refineMesh(5.0, ['x','y'])
mesh = mesh.refineMesh(20.0, ['z'])
mesh.makeSurfaces()

# plot the mesh
plotMesh(mesh)

# create solver
solver = Solver(mesh)   

# solve the matrix problem to get flux profile and keff
solver.solve(tol)

plotCellFlux(solver, ['x','z'])
