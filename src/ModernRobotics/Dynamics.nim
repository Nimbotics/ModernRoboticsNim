import neo
import neo/statics
import DataTypes
import RigidBodyMotions


template `++`(a:int): untyped =
  a + 1


func ad*[A](V:HexVector[A]):HexMatrix[A] =
    ##[Calculate the 6x6 matrix [adV] of the given 6-vector
    :param V: A 6-vector spatial velocity
    :return: The corresponding 6x6 matrix [adV]
    Used to calculate the Lie bracket [V1, V2] = [adV1]V2
    Example Input:
        V = vector([1, 2, 3, 4, 5, 6])
    Output:
        matrix([  [ 0, -3,  2,  0,  0,  0],
                  [ 3,  0, -1,  0,  0,  0],
                  [-2,  1,  0,  0,  0,  0],
                  [ 0, -6,  5,  0, -3,  2],
                  [ 6,  0, -4,  3,  0, -1],
                  [-5,  4,  0, -2,  1,  0]])
    ]##
    let omgmat = VecToso3(vector([V[0], V[1], V[2]])).asDynamic
    return vstack(hstack(omgmat, zeros(3, 3).asDynamic),
                 hstack(VecToso3(vector([V[3], V[4], V[5]])).asDynamic, omgmat)).asStatic(6,6)

func InverseDynamics*[N,V](thetalist:StaticVector[N,V], 
                          dthetalist:StaticVector[N,V], 
                          ddthetalist:StaticVector[N,V], 
                          g:TriVector[V], 
                          Ftip:HexVector[V],
                          Mlist:seq[HomogenousMatrix[V]], 
                          Glist:seq[HexMatrix[V]], 
                          Slist:StaticMatrix[6,N,V]): 
                          StaticVector[N,V] = 
    ##[Computes inverse dynamics in the space frame for an open chain robot
    :param thetalist: n-vector of joint variables
    :param dthetalist: n-vector of joint rates
    :param ddthetalist: n-vector of joint accelerations
    :param g: Gravity vector g
    :param Ftip: Spatial force applied by the end-effector expressed in frame
                 {n+1}
    :param Mlist: List of link frames {i} relative to {i-1} at the home
                  position
    :param Glist: Spatial inertia matrices Gi of the links
    :param Slist: Screw axes Si of the joints in a space frame, in the format
                  of a matrix with axes as the columns
    :return: The n-vector of required joint forces/torques
    This function uses forward-backward Newton-Euler iterations to solve the
    equation:
    taulist = Mlist(thetalist)ddthetalist + c(thetalist,dthetalist) \
              + g(thetalist) + Jtr(thetalist)Ftip
    Example Input (3 Link Robot):
        thetalist = vector([0.1, 0.1, 0.1])
        dthetalist = vector([0.1, 0.2, 0.3])
        ddthetalist = vector([2, 1.5, 1])
        g = vector([0, 0, -9.8])
        Ftip = vector([1, 1, 1, 1, 1, 1])
        M01 = matrix([[1, 0, 0,        0],
                        [0, 1, 0,        0],
                        [0, 0, 1, 0.089159],
                        [0, 0, 0,        1]])
        M12 = matrix([[ 0, 0, 1,    0.28],
                        [ 0, 1, 0, 0.13585],
                        [-1, 0, 0,       0],
                        [ 0, 0, 0,       1]])
        M23 = matrix([[1, 0, 0,       0],
                        [0, 1, 0, -0.1197],
                        [0, 0, 1,   0.395],
                        [0, 0, 0,       1]])
        M34 = matrix([[1, 0, 0,       0],
                        [0, 1, 0,       0],
                        [0, 0, 1, 0.14225],
                        [0, 0, 0,       1]])
        G1 = diag([0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7])
        G2 = diag([0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393])
        G3 = diag([0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275])
        Glist = @[G1, G2, G3]
        Mlist = @[M01, M12, M23, M34]
        Slist = matrix([[1.0, 0, 1,      0, 1,     0],
                          [0.0, 1, 0, -0.089, 0,     0],
                          [0.0, 1, 0, -0.089, 0, 0.425]]).T
    Output:
        vector([74.69616155, -33.06766016, -3.23057314])
    ]##
    var
        n = thetalist.N
        Mi = eye[V](4)
        Ai = zeros(6, N).asDynamic
        AdTi : array[N+1, HexMatrix[V]]
        Vi = zeros[V](6, ++N)
        Vdi = zeros[V](6, ++N)
        Fi = Ftip.asDynamic.asMatrix(6,1).clone
        taulist = zeros[V](N)
    Vdi[All, 0..0] = hstack(vector([0.0, 0, 0]).asDynamic, -1.0 * g.asDynamic).asmatrix(6,1)

    AdTi[N] = Adjoint(TransInv(Mlist[n]))
    
    for i in 0..(N - 1):
        Mi = Mi * Mlist[i].asDynamic
        Ai[All, i..i] = (Adjoint(TransInv(Mi.asStatic(4,4))) * Slist.column(i)).asMatrix(6,1).asDynamic
        AdTi[i] = Adjoint(MatrixExp6(VecTose3((Ai[All, i..i] * -thetalist[i]).asVector.asStatic(6))) * TransInv(Mlist[i]))
        Vi[All, ++i .. ++i] = (AdTi[i].asDynamic * Vi[All,i..i]) + Ai[All, i..i] * dthetalist[i]

        Vdi[All, ++i .. ++i] = (AdTi[i].asDynamic * Vdi[All, i..i]) + 
            (Ai[All, i..i] * ddthetalist[i]) + 
            ((ad(Vi[All, ++i .. ++i].asVector.asStatic(6)).asDynamic * (Ai[All, i..i] * dthetalist[i])))
        

    for i in countdown(n - 1, 0):
        Fi = ((AdTi[++i].T.asDynamic * Fi) + 
        (Glist[i].asDynamic * Vdi[All, ++i .. ++i]) - 
        (ad(Vi[All, ++i .. ++i].asVector.asStatic(6)).T.asDynamic * (Glist[i].asDynamic * Vi[All, ++i .. ++i])))

        taulist[i] = (Fi.T * Ai[All, i..i])[0,0]
    return taulist.asStatic(N)

func MassMatrix*[N,V](thetalist:StaticVector[N,V], 
                     Mlist:seq[HomogenousMatrix[V]], 
                     Glist:seq[HexMatrix[V]], 
                     Slist:StaticMatrix[6,N,V]): SquareMatrix[N,V] = 
    ##[Computes the mass matrix of an open chain robot based on the given
    configuration
    :param thetalist: A list of joint variables
    :param Mlist: List of link frames i relative to i-1 at the home position
    :param Glist: Spatial inertia matrices Gi of the links
    :param Slist: Screw axes Si of the joints in a space frame, in the format
                  of a matrix with axes as the columns
    :return: The numerical inertia matrix M(thetalist) of an n-joint serial
             chain at the given configuration thetalist
    This function calls InverseDynamics n times, each time passing a
    ddthetalist vector with a single element equal to one and all other
    inputs set to zero.
    Each call of InverseDynamics generates a single column, and these columns
    are assembled to create the inertia matrix.
    Example Input (3 Link Robot):
        thetalist = matrix([0.1, 0.1, 0.1])
        M01 = matrix([[1, 0, 0,        0],
                        [0, 1, 0,        0],
                        [0, 0, 1, 0.089159],
                        [0, 0, 0,        1]])
        M12 = matrix([[ 0, 0, 1,    0.28],
                        [ 0, 1, 0, 0.13585],
                        [-1, 0, 0,       0],
                        [ 0, 0, 0,       1]])
        M23 = matrix([[1, 0, 0,       0],
                        [0, 1, 0, -0.1197],
                        [0, 0, 1,   0.395],
                        [0, 0, 0,       1]])
        M34 = matrix([[1, 0, 0,       0],
                        [0, 1, 0,       0],
                        [0, 0, 1, 0.14225],
                        [0, 0, 0,       1]])
        G1 = diag([0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7])
        G2 = diag([0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393])
        G3 = diag([0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275])
        Glist = seq([G1, G2, G3])
        Mlist = seq([M01, M12, M23, M34])
        Slist = matrix([[1, 0, 1,      0, 1,     0],
                          [0, 1, 0, -0.089, 0,     0],
                          [0, 1, 0, -0.089, 0, 0.425]]).T
    Output:
        matrix([[ 2.25433380e+01, -3.07146754e-01, -7.18426391e-03],
                  [-3.07146754e-01,  1.96850717e+00,  4.32157368e-01],
                  [-7.18426391e-03,  4.32157368e-01,  1.91630858e-01]])
    ]##
    
    var m = zeros[V](N, N)
    for i in 0..(N - 1):
        var ddthetalist = statics.zeros(N,V)
        ddthetalist[i] = 1.0
        m[All, i..i] = InverseDynamics(thetalist, 
                                            statics.zeros(N,V), 
                                            ddthetalist, 
                                            vector([0.0, 0, 0]), 
                                            vector([0.0, 0, 0, 0, 0, 0]), 
                                            Mlist, 
                                            Glist, 
                                            Slist).asDynamic.asMatrix(N,1)
    result = m.asStatic(N,N)
    
func VelQuadraticForces*[N,V](thetalist:StaticVector[N,V],
                             dthetalist:StaticVector[N,V], 
                             Mlist:seq[HomogenousMatrix[V]], 
                             Glist:seq[HexMatrix[V]], 
                             Slist:StaticMatrix[6,N,V]):StaticVector[N,V] =
    ##[Computes the Coriolis and centripetal terms in the inverse dynamics of
    an open chain robot
    :param thetalist: A list of joint variables,
    :param dthetalist: A list of joint rates,
    :param Mlist: List of link frames i relative to i-1 at the home position,
    :param Glist: Spatial inertia matrices Gi of the links,
    :param Slist: Screw axes Si of the joints in a space frame, in the format
                  of a matrix with axes as the columns.
    :return: The vector c(thetalist,dthetalist) of Coriolis and centripetal
             terms for a given thetalist and dthetalist.
    This function calls InverseDynamics with g = 0, Ftip = 0, and
    ddthetalist = 0.
    Example Input (3 Link Robot):
        thetalist = np.array([0.1, 0.1, 0.1])
        dthetalist = np.array([0.1, 0.2, 0.3])
        M01 = matrix([[1, 0, 0,        0],
                        [0, 1, 0,        0],
                        [0, 0, 1, 0.089159],
                        [0, 0, 0,        1]])
        M12 = matrix([[ 0, 0, 1,    0.28],
                        [ 0, 1, 0, 0.13585],
                        [-1, 0, 0,       0],
                        [ 0, 0, 0,       1]])
        M23 = matrix([[1, 0, 0,       0],
                        [0, 1, 0, -0.1197],
                        [0, 0, 1,   0.395],
                        [0, 0, 0,       1]])
        M34 = matrix([[1, 0, 0,       0],
                        [0, 1, 0,       0],
                        [0, 0, 1, 0.14225],
                        [0, 0, 0,       1]])
        G1 = diag([0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7])
        G2 = diag([0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393])
        G3 = diag([0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275])
        Glist = seq([G1, G2, G3])
        Mlist = seq([M01, M12, M23, M34])
        Slist = seq([[1, 0, 1,      0, 1,     0],
                          [0, 1, 0, -0.089, 0,     0],
                          [0, 1, 0, -0.089, 0, 0.425]]).T
    Output:
        vector([0.26453118, -0.05505157, -0.00689132])
    ]##
    InverseDynamics(thetalist, 
                    dthetalist,
                    statics.zeros(N,V), 
                    statics.zeros(3,V), 
                    statics.zeros(6,V), 
                    Mlist, 
                    Glist, 
                    Slist)

func GravityForces*[N,V](thetalist:StaticVector[N,V], 
                         g:TriVector[V], 
                         Mlist:seq[HomogenousMatrix[V]], 
                         Glist:seq[HexMatrix[V]], 
                         Slist:StaticMatrix[6,N,V]): StaticVector[N,V] =
    ##[Computes the joint forces/torques an open chain robot requires to
    overcome gravity at its configuration
    :param thetalist: A list of joint variables
    :param g: 3-vector for gravitational acceleration
    :param Mlist: List of link frames i relative to i-1 at the home position
    :param Glist: Spatial inertia matrices Gi of the links
    :param Slist: Screw axes Si of the joints in a space frame, in the format
                  of a matrix with axes as the columns
    :return grav: The joint forces/torques required to overcome gravity at
                  thetalist
    This function calls InverseDynamics with Ftip = 0, dthetalist = 0, and
    ddthetalist = 0.
    Example Inputs (3 Link Robot):
        thetalist = np.array([0.1, 0.1, 0.1])
        g = np.array([0, 0, -9.8])
        M01 = np.array([[1, 0, 0,        0],
                        [0, 1, 0,        0],
                        [0, 0, 1, 0.089159],
                        [0, 0, 0,        1]])
        M12 = np.array([[ 0, 0, 1,    0.28],
                        [ 0, 1, 0, 0.13585],
                        [-1, 0, 0,       0],
                        [ 0, 0, 0,       1]])
        M23 = np.array([[1, 0, 0,       0],
                        [0, 1, 0, -0.1197],
                        [0, 0, 1,   0.395],
                        [0, 0, 0,       1]])
        M34 = np.array([[1, 0, 0,       0],
                        [0, 1, 0,       0],
                        [0, 0, 1, 0.14225],
                        [0, 0, 0,       1]])
        G1 = np.diag([0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7])
        G2 = np.diag([0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393])
        G3 = np.diag([0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275])
        Glist = np.array([G1, G2, G3])
        Mlist = np.array([M01, M12, M23, M34])
        Slist = np.array([[1, 0, 1,      0, 1,     0],
                          [0, 1, 0, -0.089, 0,     0],
                          [0, 1, 0, -0.089, 0, 0.425]]).T
    Output:
        vector([28.40331262, -37.64094817, -5.4415892])
    ]##
   
    InverseDynamics(thetalist, 
                    statics.zeros(N,V), 
                    statics.zeros(N,V), 
                    g, 
                    statics.zeros(6,V), 
                    Mlist, 
                    Glist, 
                    Slist)

func EndEffectorForces*[N,V](thetalist:StaticVector[N,V], 
                             Ftip:HexVector[V], 
                             Mlist:seq[HomogenousMatrix[V]], 
                             Glist:seq[HexMatrix[V]], 
                             Slist:StaticMatrix[6,N,V]): StaticVector[N,V] =
    ##[Computes the joint forces/torques an open chain robot requires only to
    create the end-effector force Ftip
    :param thetalist: A list of joint variables
    :param Ftip: Spatial force applied by the end-effector expressed in frame
                 {n+1}
    :param Mlist: List of link frames i relative to i-1 at the home position
    :param Glist: Spatial inertia matrices Gi of the links
    :param Slist: Screw axes Si of the joints in a space frame, in the format
                  of a matrix with axes as the columns
    :return: The joint forces and torques required only to create the
             end-effector force Ftip
    This function calls InverseDynamics with g = 0, dthetalist = 0, and
    ddthetalist = 0.
    Example Input (3 Link Robot):
        thetalist = vector([0.1, 0.1, 0.1])
        Ftip = vector([1, 1, 1, 1, 1, 1])
        M01 = matrix([[1, 0, 0,        0],
                        [0, 1, 0,        0],
                        [0, 0, 1, 0.089159],
                        [0, 0, 0,        1]])
        M12 = matrix([[ 0, 0, 1,    0.28],
                        [ 0, 1, 0, 0.13585],
                        [-1, 0, 0,       0],
                        [ 0, 0, 0,       1]])
        M23 = matrix([[1, 0, 0,       0],
                        [0, 1, 0, -0.1197],
                        [0, 0, 1,   0.395],
                        [0, 0, 0,       1]])
        M34 = matrix([[1, 0, 0,       0],
                        [0, 1, 0,       0],
                        [0, 0, 1, 0.14225],
                        [0, 0, 0,       1]])
        G1 = diag([0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7])
        G2 = diag([0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393])
        G3 = diag([0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275])
        Glist = seq([G1, G2, G3])
        Mlist = seq([M01, M12, M23, M34])
        Slist = matrix([[1, 0, 1,      0, 1,     0],
                          [0, 1, 0, -0.089, 0,     0],
                          [0, 1, 0, -0.089, 0, 0.425]]).T
    Output:
        vector([1.40954608, 1.85771497, 1.392409])
    ]##
   
    InverseDynamics(thetalist, 
                    statics.zeros(N,V), 
                    statics.zeros(N,V), 
                    statics.zeros(3,V), 
                    Ftip, 
                    Mlist, 
                    Glist, 
                    Slist)

func ForwardDynamics*[N,V](thetalist:StaticVector[N,V], 
                           dthetalist:StaticVector[N,V], 
                           taulist:StaticVector[N,V],
                           g:TriVector[V], 
                           Ftip:HexVector[V], 
                           Mlist:seq[HomogenousMatrix[V]], 
                           Glist:seq[HexMatrix[V]], 
                           Slist:StaticMatrix[6,N,V]): StaticVector[N,V] =
    ##[
    Computes forward dynamics in the space frame for an open chain robot
    :param thetalist: A list of joint variables
    :param dthetalist: A list of joint rates
    :param taulist: An n-vector of joint forces/torques
    :param g: Gravity vector g
    :param Ftip: Spatial force applied by the end-effector expressed in frame
                 {n+1}
    :param Mlist: List of link frames i relative to i-1 at the home position
    :param Glist: Spatial inertia matrices Gi of the links
    :param Slist: Screw axes Si of the joints in a space frame, in the format
                  of a matrix with axes as the columns
    :return: The resulting joint accelerations
    This function computes ddthetalist by solving:
    Mlist(thetalist) * ddthetalist = taulist - c(thetalist,dthetalist) \
                                     - g(thetalist) - Jtr(thetalist) * Ftip
    Example Input (3 Link Robot):
        thetalist = vector([0.1, 0.1, 0.1])
        dthetalist = vector([0.1, 0.2, 0.3])
        taulist = vector([0.5, 0.6, 0.7])
        g = vector([0, 0, -9.8])
        Ftip = vector([1, 1, 1, 1, 1, 1])
        M01 = matrix([[1, 0, 0,        0],
                        [0, 1, 0,        0],
                        [0, 0, 1, 0.089159],
                        [0, 0, 0,        1]])
        M12 = matrix([[ 0, 0, 1,    0.28],
                        [ 0, 1, 0, 0.13585],
                        [-1, 0, 0,       0],
                        [ 0, 0, 0,       1]])
        M23 = matrix([[1, 0, 0,       0],
                        [0, 1, 0, -0.1197],
                        [0, 0, 1,   0.395],
                        [0, 0, 0,       1]])
        M34 = matrix([[1, 0, 0,       0],
                      [0, 1, 0,       0],
                      [0, 0, 1, 0.14225],
                      [0, 0, 0,       1]])
        G1 = diag([0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7])
        G2 = diag([0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393])
        G3 = diag([0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275])
        Glist = seq([G1, G2, G3])
        Mlist = seq([M01, M12, M23, M34])
        Slist = matrix([[1, 0, 1,      0, 1,     0],
                        [0, 1, 0, -0.089, 0,     0],
                        [0, 1, 0, -0.089, 0, 0.425]]).T
    Output:
        vector([-0.97392907, 25.58466784, -32.91499212])
    ]##
    return inv(MassMatrix(thetalist, Mlist, Glist, Slist)) * (taulist -
                   VelQuadraticForces(thetalist, dthetalist, Mlist, Glist, Slist) -
                   GravityForces(thetalist, g, Mlist, Glist, Slist) -
                   EndEffectorForces(thetalist, Ftip, Mlist, Glist, 
                                      Slist))


func EulerStep*[N,V](thetalist:StaticVector[N,V], 
                    dthetalist:StaticVector[N,V], 
                    ddthetalist:StaticVector[N,V], dt:V): 
                    (StaticVector[N,V], StaticVector[N,V]) =
    ##[
    Compute the joint angles and velocities at the next timestep using from here
    first order Euler integration
    :param thetalist: n-vector of joint variables
    :param dthetalist: n-vector of joint rates
    :param ddthetalist: n-vector of joint accelerations
    :param dt: The timestep delta t
    :return thetalistNext: Vector of joint variables after dt from first
                           order Euler integration
    :return dthetalistNext: Vector of joint rates after dt from first order
                            Euler integration
    Example Inputs (3 Link Robot):
        thetalist = vector([0.1, 0.1, 0.1])
        dthetalist = vector([0.1, 0.2, 0.3])
        ddthetalist = vector([2, 1.5, 1])
        dt = 0.1
    Output:
        thetalistNext:
        vector([ 0.11,  0.12,  0.13])
        dthetalistNext:
        vector([ 0.3 ,  0.35,  0.4 ])
    ]##
    return (thetalist + dt * dthetalist, dthetalist + dt * ddthetalist)

func InverseDynamicsTrajectory[N,n,V](thetamat:StaticMatrix[N,n,V], 
                                      dthetamat:StaticMatrix[N,n,V], 
                                      ddthetamat:StaticMatrix[N,n,V], 
                                      g:TriVector[V], 
                                      Ftipmat:StaticMatrix[6,N,V], 
                                      Mlist:seq[HomogenousMatrix[V]], 
                                      Glist:seq[HexMatrix[V]], 
                                      Slist:StaticMatrix[6,N,V]) =
    ##[ Calculates the joint forces/torques required to move the serial chain
    along the given trajectory using inverse dynamics
    :param thetamat: An N x n matrix of robot joint variables
    :param dthetamat: An N x n matrix of robot joint velocities
    :param ddthetamat: An N x n matrix of robot joint accelerations
    :param g: Gravity vector g
    :param Ftipmat: An N x 6 matrix of spatial forces applied by the end-
                    effector (If there are no tip forces the user should
                    input a zero and a zero matrix will be used)
    :param Mlist: List of link frames i relative to i-1 at the home position
    :param Glist: Spatial inertia matrices Gi of the links
    :param Slist: Screw axes Si of the joints in a space frame, in the format
                  of a matrix with axes as the columns
    :return: The N x n matrix of joint forces/torques for the specified
             trajectory, where each of the N rows is the vector of joint
             forces/torques at each time step
    Example Inputs (3 Link Robot):
        
        import neo
        import modern_robotics
        # Create a trajectory to follow using functions from Chapter 9
        thetastart =  vector([0, 0, 0])
        thetaend =  vector([Pi / 2, Pi / 2, Pi / 2])
        Tf = 3
        N= 1000
        method = 5
        traj = JointTrajectory(thetastart, thetaend, Tf, N, method)
        thetamat = traj.clone
        dthetamat = zeros(1000, 3)
        ddthetamat = zeros(1000, 3)
        dt = Tf / (N - 1.0)
        for i in range(np.array(traj).shape[0] - 1):
            dthetamat[i + 1, :] = (thetamat[i + 1, :] - thetamat[i, :]) / dt
            ddthetamat[i + 1, :] \
            = (dthetamat[i + 1, :] - dthetamat[i, :]) / dt
        # Initialize robot description (Example with 3 links)
        g =  vector([0, 0, -9.8])
        Ftipmat = ones(N, 6)
        M01 = matrix([  [1, 0, 0,        0],
                        [0, 1, 0,        0],
                        [0, 0, 1, 0.089159],
                        [0, 0, 0,        1]])
        M12 = matrix([  [ 0, 0, 1,    0.28],
                        [ 0, 1, 0, 0.13585],
                        [-1, 0, 0,       0],
                        [ 0, 0, 0,       1]])
        M23 = matrix([[1, 0, 0,       0],
                        [0, 1, 0, -0.1197],
                        [0, 0, 1,   0.395],
                        [0, 0, 0,       1]])
        M34 = matrix([  [1, 0, 0,       0],
                        [0, 1, 0,       0],
                        [0, 0, 1, 0.14225],
                        [0, 0, 0,       1]])
        G1 = diag([0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7])
        G2 = diag([0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393])
        G3 = diag([0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275])
        Glist = seq([G1, G2, G3])
        Mlist = seq([M01, M12, M23, M34])
        Slist = matrix([[1, 0, 1,      0, 1,     0],
                          [0, 1, 0, -0.089, 0,     0],
                          [0, 1, 0, -0.089, 0, 0.425]]).T
        taumat \
        = mr.InverseDynamicsTrajectory(thetamat, dthetamat, ddthetamat, g, \
                                       Ftipmat, Mlist, Glist, Slist)
    
    ]##
    var
        thetamat = thetamat.T
        dthetamat = dthetamat.T
        ddthetamat = ddthetamat.T
        Ftipmat = Ftipmat.T
        taumat = thetamat.clone
        
    for i in 0 .. thetamat.N - 1:
        taumat[All, i..i] = InverseDynamics(thetamat[All, i..i], dthetamat[All, i..i], 
                          ddthetamat[All, i..i], g, Ftipmat[All, i..i], Mlist, 
                          Glist, Slist)

    return taumat.T

#func ForwardDynamicsTrajectory