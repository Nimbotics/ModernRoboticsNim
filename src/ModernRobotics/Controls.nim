import DataTypes
import Dynamics
import neo
import neo/statics



func ComputedTorque*[n,V](thetalist:StaticVector[n,V], 
                         dthetalist:StaticVector[n,V], 
                         eint:StaticVector[n,V], 
                         g:TriVector[V], 
                         Mlist:seq[HomogenousMatrix[V]], 
                         Glist:seq[HexMatrix[V]], 
                         Slist:StaticMatrix[6,n,V], 
                         thetalistd:StaticVector[n,V], 
                         dthetalistd:StaticVector[n,V], 
                         ddthetalistd:StaticVector[n,V], 
                         Kp:V, Ki:V, Kd:V): StaticVector[n,V] =
    ##[
    Computes the joint control torques at a particular time instant
    :param thetalist: n-vector of joint variables
    :param dthetalist: n-vector of joint rates
    :param eint: n-vector of the time-integral of joint errors
    :param g: Gravity vector g
    :param Mlist: List of link frames {i} relative to {i-1} at the home
                  position
    :param Glist: Spatial inertia matrices Gi of the links
    :param Slist: Screw axes Si of the joints in a space frame, in the format
                  of a matrix with axes as the columns
    :param thetalistd: n-vector of reference joint variables
    :param dthetalistd: n-vector of reference joint velocities
    :param ddthetalistd: n-vector of reference joint accelerations
    :param Kp: The feedback proportional gain (identical for each joint)
    :param Ki: The feedback integral gain (identical for each joint)
    :param Kd: The feedback derivative gain (identical for each joint)
    :return: The vector of joint forces/torques computed by the feedback
             linearizing controller at the current instant
    Example Input:
        thetalist = vector([0.1, 0.1, 0.1])
        dthetalist = vector([0.1, 0.2, 0.3])
        eint = vector([0.2, 0.2, 0.2])
        g = vector([0, 0, -9.8])
        M01 = matrix([  [1, 0, 0,        0],
                        [0, 1, 0,        0],
                        [0, 0, 1, 0.089159],
                        [0, 0, 0,        1]])
        M12 = matrix([  [ 0, 0, 1,    0.28],
                        [ 0, 1, 0, 0.13585],
                        [-1, 0, 0,       0],
                        [ 0, 0, 0,       1]])
        M23 = matrix([  [1, 0, 0,       0],
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
        Glist = @[G1, G2, G3]
        Mlist = @[M01, M12, M23, M34]
        Slist = matrix([[1, 0, 1,      0, 1,     0],
                          [0, 1, 0, -0.089, 0,     0],
                          [0, 1, 0, -0.089, 0, 0.425]]).T
        thetalistd = vector([1.0, 1.0, 1.0])
        dthetalistd = vector([2, 1.2, 2])
        ddthetalistd = vector([0.1, 0.1, 0.1])
        Kp = 1.3
        Ki = 1.2
        Kd = 1.1
    Output:
        vector([133.00525246, -29.94223324, -3.03276856])
    ]##
    let e = thetalistd - thetalist
    return (MassMatrix(thetalist, Mlist, Glist, Slist) * 
                  (Kp * e + Ki * (eint + e) +
                  Kd * (dthetalistd - dthetalist))) +
            InverseDynamics(thetalist, dthetalist, ddthetalistd, g, 
                             vector([0.0, 0, 0, 0, 0, 0]), Mlist, Glist, Slist)