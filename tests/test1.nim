import unittest
import ../src/ModernRobotics
import neo
import neo/statics
import math


## TODO: split test1 into multiple files for each folder in ModernRobotics



proc `=~`(a,b:SomeFloat): bool =
  abs(a-b) < 1e-9


let so3test = matrix([[ 0.0, -3,  2],
                      [ 3.0,  0, -1],
                      [-2.0,  1,  0]])
let omgtest = vector([1.0, 2, 3])

let vtest = vector([1.0,2,3])
let normtest = vector([0.2672612419124244, 0.5345224838248488, 0.8017837257372732])

let rottest = matrix([[0.0,0,1],
                      [1.0,0,0],
                      [0.0,1,0]])

let se3test = matrix( [ [ 0.0, -3,  2, 4],
                        [ 3.0,  0, -1, 5],
                        [-2.0,  1,  0, 6],
                        [ 0.0,  0,  0, 0]])

let vec6test = vector([1.0,2,3,4,5,6])

let homogtranstest = matrix([   [1.0, 0,  0, 0],
                                [0.0, 0, -1, 0],
                                [0.0, 1,  0, 3],
                                [0.0, 0,  0, 1]])

let jointscrewaxis = matrix([   [0.0, 0, 1,   0, 0.2, 0.2],
                                [1.0, 0, 0,   2,   0,   3],
                                [0.0, 1, 0,   0,   2,   1],
                                [1.0, 0, 0, 0.2, 0.3, 0.4]]).T

let thetalisttest = vector([0.2, 1.1, 0.1, 1.2])

let 
  ikinM = matrix([  [-1.0, 0,  0, 0],
                      [ 0.0, 1,  0, 6],
                      [ 0.0, 0, -1, 2],
                      [ 0.0, 0,  0, 1]])
  ikinT = matrix([  [0.0, 1,  0,     -5],
                      [1.0, 0,  0,      4],
                      [0.0, 0, -1, 1.6858],
                      [0.0, 0,  0,      1]])
  ikinthetalist0 = vector([1.5, 2.5, 3])
  ikineomg = 0.01
  ikinev = 0.001

let 
  Xstart = matrix([[1.0, 0, 0, 1],
                   [0.0, 1, 0, 0],
                   [0.0, 0, 1, 1],
                   [0.0, 0, 0, 1]])
  Xend = matrix([  [0.0, 0, 1, 0.1],
                   [1.0, 0, 0,   0],
                   [0.0, 1, 0, 4.1],
                   [0.0, 0, 0,   1]])

let
  dynthetalist = vector([0.1, 0.1, 0.1])
  dyndthetalist = vector([0.1, 0.2, 0.3])
  dynddthetalist = vector([2.0, 1.5, 1])
  g = vector([0.0, 0, -9.8])
  Ftip = vector([1.0, 1, 1, 1, 1, 1])
  M01 = matrix([  [1.0, 0, 0,        0],
                  [0.0, 1, 0,        0],
                  [0.0, 0, 1, 0.089159],
                  [0.0, 0, 0,        1]])

  M12 = matrix([  [ 0.0, 0, 1,    0.28],
                  [ 0.0, 1, 0, 0.13585],
                  [-1.0, 0, 0,       0],
                  [ 0.0, 0, 0,       1]])

  M23 = matrix([  [1.0, 0, 0,       0],
                  [0.0, 1, 0, -0.1197],
                  [0.0, 0, 1,   0.395],
                  [0.0, 0, 0,       1]])

  M34 = matrix([  [1.0, 0, 0,       0],
                  [0.0, 1, 0,       0],
                  [0.0, 0, 1, 0.14225],
                  [0.0, 0, 0,       1]])

  G1 = diag([0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7])
  G2 = diag([0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393])
  G3 = diag([0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275])
  Glist = @[G1, G2, G3]
  Mlist = @[M01, M12, M23, M34]
  Slist = matrix([  [1.0, 0, 1,      0, 1,   0],
                    [0.0, 1, 0, -0.089, 0,     0],
                    [0.0, 1, 0, -0.089, 0, 0.425]]).T




suite "kinematics":
  test "trace":
    let m = matrix([[1.0,2,3],[4.0,5,6],[7.0,8,9]])
    check: m.trace == 15.0
  test "NearZero":
    check: NearZero(-1e-7)
  test "Normalize":
    check: vtest.Normalize == normtest
  test "cross product":
    check: cross(vector([1.0,0,0]),vector([0.0,1,0])) == vector([0.0,0,1])
  test "outer product":
    let u = vector(1.0, 2.0, 3.0)
    let v = vector(3.0, 2.0, 1.0)
    check outer(u,v) == matrix(@[@[3.0, 2, 1],
                                  @[6.0, 4, 2],
                                  @[9.0, 6, 3]])
  test "RotInv":
    check: RotInv(rottest) == matrix([[0.0,1,0],
                                      [0.0,0,1],
                                      [1.0,0,0]])
  test "VecTos03":
    check: omgtest.VecToso3 == so3test
  test "so3ToVec":
    check: so3test.so3ToVec == omgtest
  test "AxisAng3":
    check: vtest.AxisAng3 == (normtest, 3.7416573867739413)

  test "AxisAngToRotation":
    check: AxisAngToRotation(vector([0.0,0,1]),Pi/2) =~ matrix([[0.0,-1,0],
                                                                [1.0, 0,0],
                                                                [0.0, 0,1]])

  test "MatrixExp3":
      
      check: so3test.MatrixExp3 =~ matrix([ [-0.69492056,  0.71352099, 0.08929286],
                                            [-0.19200697, -0.30378504, 0.93319235],
                                            [ 0.69297817,  0.6313497,  0.34810748]])
  test "MatrixLog3":
    check: rottest.MatrixLog3 =~ matrix([ [0.0,         -1.20919958,  1.20919958],
                                          [1.20919958,            0, -1.20919958],
                                          [-1.20919958,  1.20919958,           0]])
  test "RpToTrans":
    let R = matrix([  [1.0, 0,  0],
                      [0.0, 0, -1],
                      [0.0, 1,  0]])
    let p = vector([1.0,2,5])
    check: RptoTrans(R,p) == matrix([[1.0, 0,  0, 1],
                                     [0.0, 0, -1, 2],
                                     [0.0, 1,  0, 5],
                                     [0.0, 0,  0, 1]])
  test "TransToRp":
    let T = matrix([[1.0, 0,  0, 0],
                    [0.0, 0, -1, 0],
                    [0.0, 1,  0, 3],
                    [0.0, 0,  0, 1]])
    let output = TransToRp(T)
    check: output[0] == matrix([ [1.0, 0.0,  0.0],
                                     [0.0, 0.0, -1.0],
                                     [0.0, 1.0,  0.0]])
    check: output[1] == vector( [0.0, 0.0,  3.0])
  
  test "TransInv":
   
    check: homogtranstest.TransInv == matrix([[1.0,  0, 0,  0],
                                 [0.0,  0, 1, -3],
                                 [0.0, -1, 0,  0],
                                 [0.0,  0, 0,  1]])
  test "VecToSe3":
    
    check: vec6test.VecTose3 == se3test

  test "se3ToVec":
    check:  se3test.se3ToVec == vec6test
  
  test "Adjoint":
    let Q = matrix([  [1.0, 0,  0, 0],
                      [0.0, 0, -1, 0],
                      [0.0, 1,  0, 3],
                      [0.0, 0,  0, 1]])
    check: Q.Adjoint == matrix([[1.0, 0,  0, 0, 0,  0],
                                [0.0, 0, -1, 0, 0,  0],
                                [0.0, 1,  0, 0, 0,  0],
                                [0.0, 0,  3, 1, 0,  0],
                                [3.0, 0,  0, 0, 0, -1],
                                [0.0, 0,  0, 0, 1,  0]])
  test "ScrewToAxis":
    let
      q = vector([3.0, 0, 0])
      s = vector([0.0, 0, 1])
      h = 2.0
    check: ScrewToAxis(q,s,h) == vector([0.0,0,1,0,-3,2])

  test "AxisAng6":
    check: vector([1.0, 0, 0, 1, 2, 3]).AxisAng6 == 
          (vector([1.0, 0, 0, 1, 2, 3]), 1.0)

  test "MatrixExp6":
    let se3mat = matrix( [[0.0,          0,           0,          0],
                          [0.0,          0, -1.57079632, 2.35619449],
                          [0.0, 1.57079632,           0, 2.35619449],
                          [0.0,          0,           0,          0]])
    check: se3mat.MatrixExp6 =~ matrix([  [1.0, 0.0,  0.0, 0.0],
                                          [0.0, 0.0, -1.0, 0.0],
                                          [0.0, 1.0,  0.0, 3.0],
                                          [0.0,   0,    0,   1]])
  
  test "MatrixLog6":
    check homogtranstest.MatrixLog6 =~ matrix([   [0.0,          0,           0,           0],
                                                  [0.0,          0, -1.57079633,  2.35619449],
                                                  [0.0, 1.57079633,           0,  2.35619449],
                                                  [0.0,          0,           0,           0]])

  test "ProjectToSO3":
    let mat = matrix(  [[ 0.675,  0.150,  0.720],
                        [ 0.370,  0.771, -0.511],
                        [-0.630,  0.619,  0.472]])
    check mat.ProjectToSO3 =~ matrix([  [ 0.67901136,  0.14894516,  0.71885945],
                                        [ 0.37320708,  0.77319584, -0.51272279],
                                        [-0.63218672,  0.61642804,  0.46942137]])
  
  test "ProjectToSE3":
    check matrix([[ 0.675,  0.150,  0.720,  1.2],
                  [ 0.370,  0.771, -0.511,  5.4],
                  [-0.630,  0.619,  0.472,  3.6],
                  [ 0.003,  0.002,  0.010,  0.9]]).ProjectToSE3 =~ 
                  matrix([[ 0.67901136,  0.14894516,  0.71885945,  1.2 ],
                          [ 0.37320708,  0.77319584, -0.51272279,  5.4 ],
                          [-0.63218672,  0.61642804,  0.46942137,  3.6 ],
                          [ 0.0        ,  0        ,  0        ,  1    ]])
  test "DistanceToSO3":
    check: matrix([     [ 1.0,  0.0,   0.0 ],
                        [ 0.0,  0.1,  -0.95],
                        [ 0.0,  1.0,   0.1 ]]).DistanceToSO3 =~ 0.088352985

  test "DistanceToSE3":
    check: matrix([     [ 1.0,  0.0,   0.0,   1.2 ],
                        [ 0.0,  0.1,  -0.95,  1.5 ],
                        [ 0.0,  1.0,   0.1,  -0.9 ],
                        [ 0.0,  0.0,   0.1,   0.98 ]]).DistanceToSE3 =~ 0.134930537

  test "TestIfSO3":
    check: matrix([     [1.0, 0.0,  0.0 ],
                        [0.0, 0.1, -0.95],
                        [0.0, 1.0,  0.1 ]]).TestIfSO3 == false
  
  test "TestIfSE3":
    check: matrix([     [1.0, 0.0,   0.0,  1.2],
                        [0.0, 0.1, -0.95,  1.5],
                        [0.0, 1.0,   0.1, -0.9],
                        [0.0, 0.0,   0.1, 0.98]]).TestIfSE3 == false
  test "FKinBody":
    let 
      m = matrix([    [-1.0, 0,  0, 0],
                      [ 0.0, 1,  0, 6],
                      [ 0.0, 0, -1, 2],
                      [ 0.0, 0,  0, 1]])

      blist = matrix([    [0.0, 0, -1, 2, 0,   0],
                          [0.0, 0,  0, 0, 1,   0],
                          [0.0, 0,  1, 0, 0, 0.1]]).T
      thetalist = vector([Pi / 2.0, 3, Pi])

    check: FKinBody(m,blist,thetalist) =~ matrix([  [0.0, 1,  0,         -5],
                                                    [1.0, 0,  0,          4],
                                                    [0.0, 0, -1, 1.68584073],
                                                    [0.0, 0,  0,          1]])

  test "FKinSpace":
    let
      M = matrix([    [-1.0, 0,  0, 0],
                      [ 0.0, 1,  0, 6],
                      [ 0.0, 0, -1, 2],
                      [ 0.0, 0,  0, 1]])
      Slist = matrix([    [0.0, 0,  1,  4, 0,    0],
                          [0.0, 0,  0,  0, 1,    0],
                          [0.0, 0, -1, -6, 0, -0.1]]).T
      thetalist = vector([Pi / 2.0, 3, Pi])


    check: FKinSpace(M,Slist,thetalist) =~ matrix([   [0.0, 1,  0,         -5],
                                                      [1.0, 0,  0,          4],
                                                      [0.0, 0, -1, 1.68584073],
                                                      [0.0, 0,  0,          1]])
  test "JacobianBody":
    
    check: JacobianBody(jointscrewaxis,thetalisttest) =~ matrix([   [-0.04528405, 0.99500417,           0,   1],
                                                                    [ 0.74359313, 0.09304865,  0.36235775,   0],
                                                                    [-0.66709716, 0.03617541, -0.93203909,   0],
                                                                    [ 2.32586047,    1.66809,  0.56410831, 0.2],
                                                                    [-1.44321167, 2.94561275,  1.43306521, 0.3],
                                                                    [-2.06639565, 1.82881722, -1.58868628, 0.4]])

  test "JacobianSpace":
    check: JacobianSpace(jointscrewaxis, thetalisttest) =~ matrix([   [  0.0, 0.98006658, -0.09011564,  0.95749426],
                                                                      [  0.0, 0.19866933,   0.4445544,  0.28487557],
                                                                      [  1.0,          0,  0.89120736, -0.04528405],
                                                                      [  0.0, 1.95218638, -2.21635216, -0.51161537],
                                                                      [  0.2, 0.43654132, -2.43712573,  2.77535713],
                                                                      [  0.2, 2.96026613,  3.23573065,  2.22512443]])
  test "IKinBody":
    let 
        Blist = matrix([  [0.0, 0, -1, 2, 0,   0],
                          [0.0, 0,  0, 0, 1,   0],
                          [0.0, 0,  1, 0, 0, 0.1]]).T
       
        ikinres = IkinBody(Blist,ikinM,ikinT,ikinthetalist0,ikineomg,ikinev)
    check: ikinres[0] =~ vector([1.57073819, 2.999667, 3.14153913])
    check: ikinres[1] == true

  test "IKinSpace":
    let 
      Slist = matrix([    [0.0, 0,  1,  4, 0,    0],
                          [0.0, 0,  0,  0, 1,    0],
                          [0.0, 0, -1, -6, 0, -0.1]]).T
      ikinres = IKinSpace(SList,ikinM,ikinT,ikinthetalist0,ikineomg,ikinev)

    check: ikinres[0] =~ vector([1.57073819, 2.999667, 3.14153913])
    check: ikinres[1] == true

  test "ad":
    check: ad(vec6test) =~ matrix([  [ 0.0, -3,  2,  0,  0,  0],
                                     [ 3.0,  0, -1,  0,  0,  0],
                                     [-2.0,  1,  0,  0,  0,  0],
                                     [ 0.0, -6,  5,  0, -3,  2],
                                     [ 6.0,  0, -4,  3,  0, -1],
                                     [-5.0,  4,  0, -2,  1,  0]])

  test "InverseDynamics":
    
    check: InverseDynamics(dynthetalist,dyndthetalist,dynddthetalist,g,Ftip,Mlist,Glist,Slist) =~ vector([74.69616155, -33.06766016, -3.23057314])

  test "MassMatrix":
    check: MassMAtrix(dynthetalist,Mlist,Glist,Slist) =~ matrix([[ 2.25433380e+01, -3.07146754e-01, -7.18426391e-03],
                                                                 [-3.07146754e-01,  1.96850717e+00,  4.32157368e-01],
                                                                 [-7.18426391e-03,  4.32157368e-01,  1.91630858e-01]])
  
  test "VelQuadraticForces":
    check: VelQuadraticForces(dynthetalist,dyndthetalist,Mlist,Glist,Slist) =~ vector([0.26453118, -0.05505157, -0.00689132])
  
  test "GravityForces":
    check: GravityForces(dynthetalist,g,Mlist,Glist,Slist) =~ vector([28.40331262, -37.64094817, -5.4415892])

  test "EndEffectorForces":
    check: EndEffectorForces(dynthetalist,Ftip,Mlist,Glist,Slist) =~ vector([1.40954608, 1.85771497, 1.392409])

  test "ForwardDynamics":
    check: ForwardDynamics(dynthetalist,dyndthetalist,vector([0.5, 0.6, 0.7]),g,Ftip,Mlist,Glist,Slist) =~ vector([-0.97392907, 25.58466784, -32.91499212])

  test "EulerStep":
    let esres = EulerStep(dynthetalist,dyndthetalist,dynddthetalist,0.1)
    check esres[0] =~ vector([ 0.11,  0.12,  0.13])
    check: esres[1] =~ vector([ 0.3 ,  0.35,  0.4 ])
  #test "InverseDynamicsTrajectory":
  #test "ForwardDynamicsTrajectory"
  test "CubicTimeScaling":
    check: CubicTimeScaling(2.0,0.6) =~ 0.216
  test "QuinticTimeScaling":
    check: QuinticTimeScaling(2.0,0.6) =~ 0.16308
  test "JointTrajectory":
    let
      thetastart = vector([1.0, 0, 0, 1, 1, 0.2, 0,1])
      thetaend = vector([1.2, 0.5, 0.6, 1.1, 2, 2, 0.9, 1])
    check: JointTrajectory(thetastart,thetaend,4,6,cubic) =~ matrix([[1.0,    0,     0,      1,      1,     0.2,    0,        1],
                                                                     [1.0208, 0.052, 0.0624, 1.0104, 1.104, 0.3872, 0.0936,   1],
                                                                     [1.0704, 0.176, 0.2112, 1.0352, 1.352, 0.8336, 0.3168,   1],
                                                                     [1.1296, 0.324, 0.3888, 1.0648, 1.648, 1.3664, 0.5832,   1],
                                                                     [1.1792, 0.448, 0.5376, 1.0896, 1.896, 1.8128, 0.8064,   1],
                                                                     [1.2,    0.5,   0.6,    1.1,    2,     2,      0.9,      1]])
      
  test "ScrewTrajectory":
    let
      stres = ScrewTrajectory(Xstart,Xend,5.0,4,cubic)

    check: stres[0] =~ matrix([[1.0, 0, 0, 1],
                               [0.0, 1, 0, 0],
                               [0.0, 0, 1, 1],
                               [0.0, 0, 0, 1]])
    check: stres[1] =~ matrix([ [0.904111, -0.250372, 0.3462609, 0.44095],
                                [0.3462609, 0.904111, -0.250372, 0.528746],
                                [-0.250372, 0.3462609, 0.904111, 1.600666],
                                [  0.0,     0,     0,     1]])

    check: stres[2] =~ matrix([ [0.3462609, -0.250372, 0.904111, -0.117111],
                                [0.904111, 0.3462609, -0.250372,  0.472742],
                                [-0.250372, 0.904111, 0.3462609,  3.273998],
                                [  0.0,     0,     0,      1]])
    check: stres[3] =~ matrix([ [0.0, 0, 1, 0.1],
                                [1.0, 0, 0,   0],
                                [0.0, 1, 0, 4.1],
                                [0.0, 0, 0,   1]])
  test "CartesianTrajectory":
    let ctres = CartesianTrajectory(Xstart,Xend,5.0,4,quintic)
   
    check: ctres[0] =~ matrix([[1.0, 0, 0, 1],
                               [0.0, 1, 0, 0],
                               [0.0, 0, 1, 1],
                               [0.0, 0, 0, 1]])

    check: ctres[1] =~ matrix([[ 0.9366247, -0.214001075,  0.27737633, 0.811111],
                               [ 0.27737633,  0.9366247, -0.214001075,     0],
                               [-0.214001075,  0.27737633,  0.9366247, 1.650617],
                               [   0.0,      0,      0,     1]])

    check: ctres[2] =~ matrix([[ 0.27737633, -0.214001075,  0.9366247, 0.288888],
                               [ 0.9366247,  0.27737633, -0.214001075,     0],
                               [-0.214001075,  0.9366247,  0.27737633, 3.4493827],
                               [    0.0,      0,      0,    1]])

    check: ctres[3] =~ matrix([ [0.0, 0, 1, 0.1],
                                [1.0, 0, 0,   0],
                                [0.0, 1, 0, 4.1],
                                [0.0, 0, 0,   1]])

  test "ComputedTorque":
    let comtor = ComputedTorque(dynthetalist,
                                dyndthetalist,
                                vector([0.2, 0.2, 0.2]),
                                g, Mlist, Glist, Slist,
                                vector([1.0, 1.0, 1.0]),
                                vector([2.0, 1.2, 2]),
                                dynthetalist,
                                1.3, 1.2, 1.1)
    check: comtor =~ vector([133.00525246, -29.94223324, -3.03276856])
                          
