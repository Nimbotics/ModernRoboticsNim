import neo
import neo/statics
import DataTypes
import RigidBodyMotions


func FKinBody*[I:static[int],V](Mat:HomogenousMatrix[V], Blist:StaticMatrix[6,I,V], thetalist:StaticVector[I,V]): HomogenousMatrix[V] =
    ##[
        Computes forward kinematics in the body frame for an open chain robot
    :param M: The home configuration (position and orientation) of the end-
              effector
    :param Blist: The joint screw axes in the end-effector frame when the
                  manipulator is at the home position, in the format of a
                  matrix with axes as the columns
    :param thetalist: A list of joint coordinates
    :return: A homogeneous transformation matrix representing the end-
             effector frame when the joints are at the specified coordinates
             (i.t.o Body Frame)
    Example Input:
        M = matrix([  [-1, 0,  0, 0],
                      [ 0, 1,  0, 6],
                      [ 0, 0, -1, 2],
                      [ 0, 0,  0, 1]])
        Blist = matrix([  [0, 0, -1, 2, 0,   0],
                          [0, 0,  0, 0, 1,   0],
                          [0, 0,  1, 0, 0, 0.1]]).T
        thetalist = vector([Pi / 2.0, 3, Pi])
    Output:
        matrix([  [0, 1,  0,         -5],
                  [1, 0,  0,          4],
                  [0, 0, -1, 1.68584073],
                  [0, 0,  0,          1]])
    ]##
    var tmp = Mat.asDynamic
 
    var i = 0
    for col in Blist.columns:
        tmp = tmp * MatrixExp6((VecTose3((col * thetalist[i])))).asDynamic
        inc i
    result = tmp.asStatic(4,4)

func FKinSpace*[I:static[int],V](M:HomogenousMatrix[V], Slist:StaticMatrix[6,I,V], thetalist:StaticVector[I,V]): HomogenousMatrix[V] =
    ##[
        Computes forward kinematics in the space frame for an open chain robot
    :param M: The home configuration (position and orientation) of the end-
              effector
    :param Slist: The joint screw axes in the space frame when the
                  manipulator is at the home position, in the format of a
                  matrix with axes as the columns
    :param thetalist: A list of joint coordinates
    :return: A homogeneous transformation matrix representing the end-
             effector frame when the joints are at the specified coordinates
             (i.t.o Space Frame)
    Example Input:
        M = matrix([  [-1, 0,  0, 0],
                      [ 0, 1,  0, 6],
                      [ 0, 0, -1, 2],
                      [ 0, 0,  0, 1]])
        Slist = matrix([  [0, 0,  1,  4, 0,    0],
                          [0, 0,  0,  0, 1,    0],
                          [0, 0, -1, -6, 0, -0.1]]).T
        thetalist = vector([Pi / 2.0, 3, Pi])
    Output:
        matrix([  [0, 1,  0,         -5],
                  [1, 0,  0,          4],
                  [0, 0, -1, 1.68584073],
                  [0, 0,  0,          1]])
    ]##

    var i = I - 1
    #var i = M.N - 2
    var tmp = M.asDynamic
    for col in Slist.asDynamic.reverseColumns:
        tmp = MatrixExp6((VecTose3((col * thetalist[i]).asStatic(6)))).asDynamic * tmp
        dec i
    
    result = tmp.asStatic(4,4)

func JacobianBody*[I:static[int],V](Blist:StaticMatrix[6,I,V], thetalist:StaticVector[I,V]): StaticMatrix[6,I,V] =
    ##[Computes the body Jacobian for an open chain robot
    :param Blist: The joint screw axes in the end-effector frame when the
                  manipulator is at the home position, in the format of a
                  matrix with axes as the columns
    :param thetalist: A list of joint coordinates
    :return: The body Jacobian corresponding to the inputs (6xn real
             numbers)
    Example Input:
        Blist = matrix([  [0, 0, 1,   0, 0.2, 0.2],
                          [1, 0, 0,   2,   0,   3],
                          [0, 1, 0,   0,   2,   1],
                          [1, 0, 0, 0.2, 0.3, 0.4]]).T
        thetalist = vector([0.2, 1.1, 0.1, 1.2])
    Output:
        matrix([  [-0.04528405, 0.99500417,           0,   1],
                  [ 0.74359313, 0.09304865,  0.36235775,   0],
                  [-0.66709716, 0.03617541, -0.93203909,   0],
                  [ 2.32586047,    1.66809,  0.56410831, 0.2],
                  [-1.44321167, 2.94561275,  1.43306521, 0.3],
                  [-2.06639565, 1.82881722, -1.58868628, 0.4]])
    ]##
    var Jb = Blist.asDynamic.clone()
    var T = eye[V](4)
    #for i in range(len(thetalist) - 2, -1, -1):
    #    T = T * MatrixExp6(VecTose3(Blist[:, i + 1] \
    #                                     * -thetalist[i + 1])))
    #    Jb[:, i] = Adjoint(T) * Blist[:, i]
    
    #var i = I - 2
    #for c in Blist.asDynamic.reverseColumns:
    for i in countdown(I - 2,0):
        
        T = T * MatrixExp6(VecTose3(Blist.column(i + 1) * -thetalist[i + 1])).asDynamic

        Jb[0..5, i..i] =  Adjoint(T.asStatic(4,4)).asDynamic * Blist.column(i).asDynamic.asMatrix(6,1)
   
    result = Jb.asStatic(6,I)


func JacobianSpace*[I:static[int],V](Slist:StaticMatrix[6,I,V], thetalist:StaticVector[I,V]): StaticMatrix[6,I,V] =
    ##[Computes the space Jacobian for an open chain robot
    :param Slist: The joint screw axes in the space frame when the
                  manipulator is at the home position, in the format of a
                  matrix with axes as the columns
    :param thetalist: A list of joint coordinates
    :return: The space Jacobian corresponding to the inputs (6xn real
             numbers)
    Example Input:
        Slist = matrix([   [0.0, 0, 1,   0, 0.2, 0.2],
                           [1.0, 0, 0,   2,   0,   3],
                           [0.0, 1, 0,   0,   2,   1],
                           [1.0, 0, 0, 0.2, 0.3, 0.4]]).T
        thetalist = vector([0.2, 1.1, 0.1, 1.2])
    Output:
        matrix([  [  0, 0.98006658, -0.09011564,  0.95749426],
                  [  0, 0.19866933,   0.4445544,  0.28487557],
                  [  1,          0,  0.89120736, -0.04528405],
                  [  0, 1.95218638, -2.21635216, -0.51161537],
                  [0.2, 0.43654132, -2.43712573,  2.77535713],
                  [0.2, 2.96026613,  3.23573065,  2.22512443]])
    ]##
    var Js = Slist.asDynamic.clone
    var T = eye[V](4).asStatic(4,4)
    for i in 1 .. I - 1:
        T = T * MatrixExp6(VecTose3(Slist.column(i - 1) * thetalist[i - 1]))
        Js[0..5, i..i] = (Adjoint(T) * Slist.column(i)).asDynamic.asMatrix(6,1)
    result = Js.asStatic(6,I)

func IKinBody*[I:static[int],V](Blist:StaticMatrix[6,I,V], 
                                M:HomogenousMatrix[V], 
                                T:HomogenousMatrix[V], 
                                thetalist0:StaticVector[I,V], 
                                eomg:V, ev:V): (StaticVector[I,V],bool) =
    ##[Computes inverse kinematics in the body frame for an open chain robot
    :param Blist: The joint screw axes in the end-effector frame when the
                  manipulator is at the home position, in the format of a
                  matrix with axes as the columns
    :param M: The home configuration of the end-effector
    :param T: The desired end-effector configuration Tsd
    :param thetalist0: An initial guess of joint angles that are close to
                       satisfying Tsd
    :param eomg: A small positive tolerance on the end-effector orientation
                 error. The returned joint angles must give an end-effector
                 orientation error less than eomg
    :param ev: A small positive tolerance on the end-effector linear position
               error. The returned joint angles must give an end-effector
               position error less than ev
    :return thetalist: Joint angles that achieve T within the specified
                       tolerances,
    :return success: A logical value where TRUE means that the function found
                     a solution and FALSE means that it ran through the set
                     number of maximum iterations without finding a solution
                     within the tolerances eomg and ev.
    Uses an iterative Newton-Raphson root-finding method.
    The maximum number of iterations before the algorithm is terminated has
    been hardcoded in as a variable called maxiterations. It is set to 20 at
    the start of the function, but can be changed if needed.
    Example Input:
        Blist = matrix([[0, 0, -1, 2, 0,   0],
                          [0, 0,  0, 0, 1,   0],
                          [0, 0,  1, 0, 0, 0.1]]).T
        M = matrix([[-1, 0,  0, 0],
                      [ 0, 1,  0, 6],
                      [ 0, 0, -1, 2],
                      [ 0, 0,  0, 1]])
        T = matrix([[0, 1,  0,     -5],
                      [1, 0,  0,      4],
                      [0, 0, -1, 1.6858],
                      [0, 0,  0,      1]])
        thetalist0 = vector([1.5, 2.5, 3])
        eomg = 0.01
        ev = 0.001
    Output:
        (vector([1.57073819, 2.999667, 3.14153913]), True)
    ]##
    let maxiterations = 20
    var 
        thetalist = thetalist0.clone
        i = 0
    
        Vb = se3ToVec(MatrixLog6((TransInv(FKinBody(M, Blist, thetalist)) * T)))
        err = norm(vector([Vb[0], Vb[1], Vb[2]])) > eomg or norm(vector([Vb[3], Vb[4], Vb[5]])) > ev
    while err and (i < maxiterations):
        thetalist = thetalist + (pinv(JacobianBody(Blist, thetalist)) * Vb.asDynamic).asStatic(I)
        inc i
        Vb = se3ToVec(MatrixLog6((TransInv(FKinBody(M, Blist, thetalist)) * T)))
        err = norm(vector([Vb[0], Vb[1], Vb[2]])) > eomg or norm(vector([Vb[3], Vb[4], Vb[5]])) > ev
    return (thetalist, not err)

func IKinSpace*[I:static[int],V](Slist:StaticMatrix[6,I,V], 
                                M:HomogenousMatrix[V], 
                                T:HomogenousMatrix[V], 
                                thetalist0:StaticVector[I,V], 
                                eomg:V, ev:V): (StaticVector[I,V],bool) =
    ##[Computes inverse kinematics in the space frame for an open chain robot
    :param Slist: The joint screw axes in the space frame when the
                  manipulator is at the home position, in the format of a
                  matrix with axes as the columns
    :param M: The home configuration of the end-effector
    :param T: The desired end-effector configuration Tsd
    :param thetalist0: An initial guess of joint angles that are close to
                       satisfying Tsd
    :param eomg: A small positive tolerance on the end-effector orientation
                 error. The returned joint angles must give an end-effector
                 orientation error less than eomg
    :param ev: A small positive tolerance on the end-effector linear position
               error. The returned joint angles must give an end-effector
               position error less than ev
    :return thetalist: Joint angles that achieve T within the specified
                       tolerances,
    :return success: A logical value where TRUE means that the function found
                     a solution and FALSE means that it ran through the set
                     number of maximum iterations without finding a solution
                     within the tolerances eomg and ev.
    Uses an iterative Newton-Raphson root-finding method.
    The maximum number of iterations before the algorithm is terminated has
    been hardcoded in as a variable called maxiterations. It is set to 20 at
    the start of the function, but can be changed if needed.
    Example Input:
        Slist = matrix([[0, 0,  1,  4, 0,    0],
                          [0, 0,  0,  0, 1,    0],
                          [0, 0, -1, -6, 0, -0.1]]).T
        M = matrix([[-1, 0,  0, 0],
                      [ 0, 1,  0, 6],
                      [ 0, 0, -1, 2],
                      [ 0, 0,  0, 1]])
        T = matrix([[0, 1,  0,     -5],
                      [1, 0,  0,      4],
                      [0, 0, -1, 1.6858],
                      [0, 0,  0,      1]])
        thetalist0 = vector([1.5, 2.5, 3])
        eomg = 0.01
        ev = 0.001
    Output:
        (vector([ 1.57073783,  2.99966384,  3.1415342 ]), True)
    ]##
    let maxiterations = 20
    var
        thetalist = thetalist0.clone
        i = 0
        Tsb = FKinSpace(M,Slist, thetalist)
        Vs = Adjoint(Tsb) * se3ToVec(MatrixLog6((TransInv(Tsb) * T)))
        err = norm(vector([Vs[0], Vs[1], Vs[2]])) > eomg or norm(vector([Vs[3], Vs[4], Vs[5]])) > ev
    while err and i < maxiterations:
        thetalist = thetalist + (pinv(JacobianSpace(Slist, thetalist)) * Vs.asDynamic).asStatic(I)
        inc i
        Tsb = FKinSpace(M, Slist, thetalist)
        Vs = Adjoint(Tsb) * se3ToVec(MatrixLog6((TransInv(Tsb) * T)))
        err = norm(vector([Vs[0], Vs[1], Vs[2]])) > eomg or norm(vector([Vs[3], Vs[4], Vs[5]])) > ev
    return (thetalist, not err)