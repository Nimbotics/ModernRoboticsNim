import neo
import neo/statics
import math
import DataTypes
import RigidBodyMotions


func CubicTimeScaling*[V:SomeFloat](Tf:V, t:V): V =
    ##[Computes s(t) for a cubic time scaling
    :param Tf: Total time of the motion in seconds from rest to rest
    :param t: The current time t satisfying 0 < t < Tf
    :return: The path parameter s(t) corresponding to a third-order
             polynomial motion that begins and ends at zero velocity
    Example Input:
        Tf = 2
        t = 0.6
    Output:
        0.216
    ]##
    return (3 * (1.0 * t / Tf) ^ 2) - (2 * (1.0 * t / Tf) ^ 3)

func QuinticTimeScaling*[V:SomeFloat](Tf:V, t:V): V =
    ##[Computes s(t) for a quintic time scaling
    :param Tf: Total time of the motion in seconds from rest to rest
    :param t: The current time t satisfying 0 < t < Tf
    :return: The path parameter s(t) corresponding to a fifth-order
             polynomial motion that begins and ends at zero velocity and zero
             acceleration
    Example Input:
        Tf = 2
        t = 0.6
    Output:
        0.16308
    ]##
    return 10 * (1.0 * t / Tf) ^ 3 - 15 * (1.0 * t / Tf) ^ 4 + 6 * (1.0 * t / Tf) ^ 5

func JointTrajectory*[N,V](thetastart:StaticVector[N,V], 
                         thetaend:StaticVector[N,V], Tf:V, 
                         n:static[int], scale:poly): StaticMatrix[n,N,V] =
    ##[
    Computes a straight-line trajectory in joint space
    :param thetastart: The initial joint variables
    :param thetaend: The final joint variables
    :param Tf: Total time of the motion in seconds from rest to rest
    :param n: The number of points n > 1 (Start and stop) in the discrete
              representation of the trajectory
    :param method: The time-scaling method, where 3 indicates cubic (third-
                   order polynomial) time scaling and 5 indicates quintic
                   (fifth-order polynomial) time scaling
    :return: A trajectory as an n x N matrix, where each row is an n-vector
             of joint variables at an instant in time. The first row is
             thetastart and the Nth row is thetaend . The elapsed time
             between each row is Tf / (N - 1)
    Example Input:
        thetastart = vector([1, 0, 0, 1, 1, 0.2, 0,1])
        thetaend = vector([1.2, 0.5, 0.6, 1.1, 2, 2, 0.9, 1])
        Tf = 4
        N = 6
        scale = cubic
    Output:
        matrix([  [     1,     0,      0,      1,     1,    0.2,      0, 1]
                  [1.0208, 0.052, 0.0624, 1.0104, 1.104, 0.3872, 0.0936, 1]
                  [1.0704, 0.176, 0.2112, 1.0352, 1.352, 0.8336, 0.3168, 1]
                  [1.1296, 0.324, 0.3888, 1.0648, 1.648, 1.3664, 0.5832, 1]
                  [1.1792, 0.448, 0.5376, 1.0896, 1.896, 1.8128, 0.8064, 1]
                  [   1.2,   0.5,    0.6,    1.1,     2,      2,    0.9, 1]])
    ]##
    var
        timegap = Tf / (n.V - 1.0)
        traj = zeros[V](N, n)
        s:V
    for i in 0 .. n-1:
        if scale == cubic:
            s = CubicTimeScaling(Tf, timegap * i.V)
        else:
            s = QuinticTimeScaling(Tf, timegap * i.V)
        traj[All, i..i] = (s * thetaend + (1 - s) * thetastart).asDynamic.asMatrix(N,1)
    
    return traj.T.asStatic(n,N)

func ScrewTrajectory*[V](Xstart:HomogenousMatrix[V], 
                        Xend:HomogenousMatrix[V], Tf:V, 
                        n:static[int], 
                        scale:poly): array[n, HomogenousMatrix[V]] =
    ##[
    Computes a trajectory as a list of N SE(3) matrices corresponding to
    the screw motion about a space screw axis
    :param Xstart: The initial end-effector configuration
    :param Xend: The final end-effector configuration
    :param Tf: Total time of the motion in seconds from rest to rest
    :param n: The number of points n > 1 (Start and stop) in the discrete
              representation of the trajectory
    :param method: The time-scaling method, where 3 indicates cubic (third-
                   order polynomial) time scaling and 5 indicates quintic
                   (fifth-order polynomial) time scaling
    :return: The discretized trajectory as a list of n matrices in SE(3)
             separated in time by Tf/(n-1). The first in the list is Xstart
             and the Nth is Xend
    Example Input:
        Xstart = matrix([[1, 0, 0, 1],
                           [0, 1, 0, 0],
                           [0, 0, 1, 1],
                           [0, 0, 0, 1]])
        Xend = matrix([[0, 0, 1, 0.1],
                         [1, 0, 0,   0],
                         [0, 1, 0, 4.1],
                         [0, 0, 0,   1]])
        Tf = 5
        n = 4
        method = cubic
    Output:
        [matrix([  [1.0, 0, 0, 1]
                   [0.0, 1, 0, 0]
                   [0.0, 0, 1, 1]
                   [0.0, 0, 0, 1]]),
         matrix([  [0.904, -0.25, 0.346, 0.441]
                   [0.346, 0.904, -0.25, 0.529]
                   [-0.25, 0.346, 0.904, 1.601]
                   [  0.0,     0,     0,     1]]),
         matrix([[0.346, -0.25, 0.904, -0.117]
                   [0.904, 0.346, -0.25,  0.473]
                   [-0.25, 0.904, 0.346,  3.274]
                   [  0.0,     0,     0,      1]]),
         matrix([  [0.0, 0, 1, 0.1]
                   [1.0, 0, 0,   0]
                   [0.0, 1, 0, 4.1]
                   [0.0, 0, 0,   1]])]
    ]##
    var 
        timegap = Tf / (n.V - 1.0)
        s:V
    
    for i in 0.. n - 1:
        if scale == cubic:
            s = CubicTimeScaling(Tf, timegap * i.V)
        else:
            s = QuinticTimeScaling(Tf, timegap * i.V)
        result[i] = (Xstart * MatrixExp6(MatrixLog6((TransInv(Xstart) * Xend)) * s))
    

func CartesianTrajectory*[V](Xstart:HomogenousMatrix[V], 
                            Xend:HomogenousMatrix[V], 
                            Tf:V, n:static[int], 
                            scale:poly): array[n,HomogenousMatrix[V]] =
    ##[Computes a trajectory as a list of n SE(3) matrices corresponding to
    the origin of the end-effector frame following a straight line
    :param Xstart: The initial end-effector configuration
    :param Xend: The final end-effector configuration
    :param Tf: Total time of the motion in seconds from rest to rest
    :param N: The number of points n > 1 (Start and stop) in the discrete
              representation of the trajectory
    :param method: The time-scaling method, where 3 indicates cubic (third-
                   order polynomial) time scaling and 5 indicates quintic
                   (fifth-order polynomial) time scaling
    :return: The discretized trajectory as a list of N matrices in SE(3)
             separated in time by Tf/(n-1). The first in the list is Xstart
             and the nth is Xend
    This function is similar to ScrewTrajectory, except the origin of the
    end-effector frame follows a straight line, decoupled from the rotational
    motion.
    Example Input:
        Xstart = matrix([[1, 0, 0, 1],
                         [0, 1, 0, 0],
                         [0, 0, 1, 1],
                         [0, 0, 0, 1]])
        Xend = matrix([[0, 0, 1, 0.1],
                       [1, 0, 0,   0],
                       [0, 1, 0, 4.1],
                       [0, 0, 0,   1]])
        Tf = 5
        n = 4
        method = quintic
    Output:
        [matrix([[1, 0, 0, 1]
                   [0, 1, 0, 0]
                   [0, 0, 1, 1]
                   [0, 0, 0, 1]]),
         matrix([[ 0.937, -0.214,  0.277, 0.811]
                   [ 0.277,  0.937, -0.214,     0]
                   [-0.214,  0.277,  0.937, 1.651]
                   [     0,      0,      0,     1]]),
         matrix([[ 0.277, -0.214,  0.937, 0.289]
                   [ 0.937,  0.277, -0.214,     0]
                   [-0.214,  0.937,  0.277, 3.449]
                   [     0,      0,      0,     1]]),
         matrix([[0, 0, 1, 0.1]
                   [1, 0, 0,   0]
                   [0, 1, 0, 4.1]
                   [0, 0, 0,   1]])]
    ]##

    var 
        timegap = Tf / (n.V - 1.0)
        (Rstart, pstart) = TransToRp(Xstart)
        (Rend, pend) = TransToRp(Xend)
        s:V
    for i in 0 .. n - 1:
        if scale == cubic:
            s = CubicTimeScaling(Tf, timegap * i.V)
        else:
            s = QuinticTimeScaling(Tf, timegap * i.V)
        result[i] = (vstack(hstack((Rstart * MatrixExp3(MatrixLog3(Rstart.T * Rend) * s)).asDynamic, 
                    (s * pend + (1 - s) * pstart).asDynamic.asMatrix(3,1)), 
                   matrix(@[@[0.0, 0, 0, 1]]))).asStatic(4,4)