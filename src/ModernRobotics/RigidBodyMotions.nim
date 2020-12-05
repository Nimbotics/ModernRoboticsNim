import neo
import neo/statics
import DataTypes
import math


func RotInv*[V](R:RotationMatrix[V]): RotationMatrix[V] =
    ##[
        Inverts a rotation matrix
    :param R: A rotation matrix
    :return: The inverse of R
    Example Input:
        R = matrix( [  [0, 0, 1],
                       [1, 0, 0],
                       [0, 1, 0]])
    Output:
        matrix( [ [0, 1, 0],
                  [0, 0, 1],
                  [1, 0, 0]])
    ]##
    R.T

func VecToso3*[V](omg:Omega[V]):so3mat[V] =
    ##[
        Converts a 3-vector to an so(3) representation
    :param omg: A 3-vector
    :return: The skew symmetric representation of omg
    Example Input:
        omg = vector([1, 2, 3])
    Output:
        matrix( [ [ 0, -3,  2],
                  [ 3,  0, -1],
                  [-2,  1,  0]])
    ]##
    result =  [[    0.V,  -omg[2],    omg[1]],
               [ omg[2],        0,   -omg[0]],
               [-omg[1],   omg[0],        0]].matrix


func so3ToVec*[V](m:so3mat[V]): Omega[V] =
    ##[
        Converts an so(3) representation to a 3-vector
    :param so3mat: A 3x3 skew-symmetric matrix
    :return: The 3-vector corresponding to so3mat
    Example Input:
        so3mat = matrix( [ [ 0, -3,  2],
                           [ 3,  0, -1],
                           [-2,  1,  0]])
    Output:
        vector(@[1, 2, 3])
    ]##
    vector([m[2, 1], m[0, 2], m[1, 0]])


func VecsToRotation*[V](a,b:TriVector[V]): RotationMatrix[V] =
    ##[
        converts a vector and an origin vector into a rotation matrix 
        between the two vectors
        example Input:
            a = vector([0,0,1])
            b = vector([1,0,0])
        Output:
            matrix[

    ]##
    let v = cross(a,b)
    let c = a * b
    let so3 = VecToso3(v)
    result = (eye(3,V) + so3 + so3 * so3 * (1/(1+c)))


func AxisAng3*(expc3:TriVector):(TriVector, float) =
    ##[
        Converts a 3-vector of exponential coordinates for rotation into
    axis-angle form
    :param expc3: A 3-vector of exponential coordinates for rotation
    :return Trivector: A unit rotation axis
    :return float: The corresponding rotation angle
    Example Input:
        expc3 = vector([1, 2, 3])
    Output:
        (vector([0.26726124, 0.53452248, 0.80178373]), 3.7416573867739413)
    ]##
    if expc3 in [vector [0.0, 0.0, 1.0], vector [1.0, 0.0, 0.0], vector [0.0, 1.0, 0.0]]:
        result = (expc3, 0.0)
    else:
        result = (expc3.Normalize,expc3.norm)

func AxisAngToRotation*[V](u:TriVector[V],theta:V):RotationMatrix[V] =
    ##[
        Converts a 3-vector and and an angle theta into a rotation matrix
        :param u: A 3-vector representing the axis
        :param theta: a float representing the angle around the axis
        :return Rotation Marix: the resulting Rotation Matrix
        Example Input:
            u = vector([0,0,1])
            theta = Pi/2
        Output: matrix([[0.0,-1,0],
                        [1.0, 0,0],
                        [0.0, 0,1]])

    ]##
    cos(theta) * eye(3,V) + sin(theta) * VecToso3(u) + (1 - cos(theta)) * outer(u,u)

func AxisAngToRotation*[V](aa:AxisAng[V]):RotationMatrix[V] =
    AxisAngToRotation(aa[0],aa[1])

func MatrixExp3*[V](m:so3mat[V]): TriMatrix[V] =
    ##[
        Computes the matrix exponential of a matrix in so(3)
    :param so3mat: A 3x3 skew-symmetric matrix
    :return: The matrix exponential of so3mat
    Example Input:
        so3mat = matrix(@[ @[ 0, -3,  2],
                           @[ 3,  0, -1],
                           @[-2,  1,  0]])
    Output:
        matrix(@[ @[-0.69492056,  0.71352099,  0.08929286],
                  @[-0.19200697, -0.30378504,  0.93319235],
                  @[ 0.69297817,  0.6313497 ,  0.34810748]])
    ]##
    let omgtheta = so3ToVec(m)
    if NearZero(omgtheta.norm):

        return statics.eye(3,V)
    else:
        let
            theta = AxisAng3(omgtheta)[1]
            omgmat = m / theta
        return eye(3,V) + sin(theta) * omgmat + (1-cos(theta)) * (omgmat * omgmat)

func MatrixLog3*[V](R:RotationMatrix[V]): TriMatrix[V] =
    ##[
        Computes the matrix logarithm of a rotation matrix
    :param R: A 3x3 rotation matrix
    :return: The matrix logarithm of R
    Example Input:
        R = matrix(@[ @[0, 0, 1],
                      @[1, 0, 0],
                      @[0, 1, 0]])
    Output:
        matrix(@[ @[          0, -1.20919958,  1.20919958],
                  @[ 1.20919958,           0, -1.20919958],
                  @[-1.20919958,  1.20919958,           0]])
    ]##
    let acosinput = (R.trace - 1.0) / 2.0
    var omg:TriVector[R.A]
    if acosinput >= 1.0:
        return statics.zeros(3,3,R.A)
    elif acosinput <= -1.0:
        if not NearZero(1.0 + R[2,2]):
            omg = ((1.0 / sqrt(2.0 * (R[2,2]))) * vector([R[0,2], R[1, 2], 1.0 + R[2,2]]))
        elif not NearZero(1.0 + R[1,1]):
            omg = ((1.0 / sqrt(2.0 * (R[1,1]))) * vector([R[0,1], R[1, 1], 1.0 + R[2,1]]))
        else:
            omg = ((1.0 / sqrt(2.0 * (R[0,0]))) * vector([R[0,0], R[1, 0], 1.0 + R[2,0]]))
        return VecToso3(Pi * omg)
    else:
        let theta = arccos(acosinput)
        return (theta / 2.0 / sin(theta)) * (R - R.T)

func RpToTrans*[V](R:RotationMatrix[V], p:PositionVec[V]): QuadMatrix[V] = 
    ##[
        Converts a rotation matrix and a position vector into homogeneous
    transformation matrix
    :param R: A 3x3 rotation matrix
    :param p: A 3-vector
    :return: A homogeneous transformation matrix corresponding to the inputs
    Example Input:
        R = matrix(@[  @[1, 0,  0],
                       @[0, 0, -1],
                       @[0, 1,  0]])
        p = vector(    @[1, 2,  5])
    Output:
        matrix(@[@[1, 0,  0, 1],
                  [0, 0, -1, 2],
                  [0, 1,  0, 5],
                  [0, 0,  0, 1]])
    ]##
    let R = R.asDynamic
    
    result = R.hstack(p.asMatrix(3,1).asDynamic).vstack(matrix(@[@[0.0,0,0,1]])).asStatic(4,4)

func TransToRp*[V](Q:QuadMatrix[V]): (RotationMatrix[V],PositionVec[V]) =
    ##[
        Converts a homogeneous transformation matrix into a rotation matrix
    and position vector
    :param T: A homogeneous transformation matrix
    :return R: The corresponding rotation matrix,
    :return p: The corresponding position vector.
    Example Input:
        T = matrix(@[@[1, 0,  0, 0],
                      [0, 0, -1, 0],
                      [0, 1,  0, 3],
                      [0, 0,  0, 1]])
    Output:
        (matrix(@[ @[1, 0,  0],
                   @[0, 0, -1],
                   @[0, 1,  0]]),
        vector(   @[0, 0, 3])
    ]##
    (Q[0..2, 0..2], Q[0..2, 3..3].asVector)

func TransInv*[V](Q:QuadMatrix[V]):QuadMatrix[V] = 
    ##[
        Inverts a homogeneous transformation matrix
    :param Q: A homogeneous transformation matrix
    :return: The inverse of Q
    Uses the structure of transformation matrices to avoid taking a matrix
    inverse, for efficiency.
    Example input:
        Q = matrix(@[ @[1, 0,  0, 0],
                      @[0, 0, -1, 0],
                      @[0, 1,  0, 3],
                      @[0, 0,  0, 1]])
    Output:
        matrix(@[  @[1,  0, 0,  0],
                   @[0,  0, 1, -3],
                   @[0, -1, 0,  0],
                   @[0,  0, 0,  1]])
    ]##
    let 
        (R, p) = Q.TransToRp
        Rt = R.T.asDynamic
        pd = p.asDynamic
    result = Rt.hstack((-1.V * (Rt * pd)).asMatrix(3,1)).vstack(matrix(@[@[0.V,0,0,1]])).asStatic(4,4)

func VecTose3*[V](H:HexVector[V]):QuadMatrix[V] =
    ##[
        Converts a spatial velocity vector into a 4x4 matrix in se3
    :param H: A 6-vector representing a spatial velocity
    :return: The 4x4 se3 representation of V
    Example Input:
        V = vector([1, 2, 3, 4, 5, 6])
    Output:
        matrix( [ [ 0, -3,  2, 4],
                  [ 3,  0, -1, 5],
                  [-2,  1,  0, 6],
                  [ 0,  0,  0, 0]])
    ]##
    let H = H.asMatrix(1,6).asDynamic
    (VecToso3(H[0..0,0..2].asVector.asStatic(3)).asDynamic.hstack(H[0..0,3..5].T).vstack(matrix(@[@[0.V,0,0,0]]))).asStatic(4,4)

func se3ToVec*[V](Q:QuadMatrix[V]): HexVector[V] =
    ##[
        Converts an se3 matrix into a spatial velocity vector
    :param se3mat: A 4x4 matrix in se3
    :return: The spatial velocity 6-vector corresponding to se3mat
    Example Input:
        se3mat = matrix(@[ @[ 0, -3,  2, 4],
                           @[ 3,  0, -1, 5],
                           @[-2,  1,  0, 6],
                           @[ 0,  0,  0, 0]])
    Output:
        vector(@[1, 2, 3, 4, 5, 6])
    ]##
    vector([Q[2,1],Q[0,2],Q[1,0],Q[0,3],Q[1,3],Q[2,3]])

func Adjoint*[V](Q:QuadMatrix[V]): HexMatrix[V] =
    ##[
        Computes the adjoint representation of a homogeneous transformation
    matrix
    :param T: A homogeneous transformation matrix
    :return: The 6x6 adjoint representation [AdT] of T
    Example Input:
        Q = matrix([[1, 0,  0, 0],
                      [0, 0, -1, 0],
                      [0, 1,  0, 3],
                      [0, 0,  0, 1]])
    Output:
        matrix([[1, 0,  0, 0, 0,  0],
                  [0, 0, -1, 0, 0,  0],
                  [0, 1,  0, 0, 0,  0],
                  [0, 0,  3, 1, 0,  0],
                  [3, 0,  0, 0, 0, -1],
                  [0, 0,  0, 0, 1,  0]])
    ]##
    let 
        (R, p) = Q.TransToRp
        Rd = R.asDynamic
    (Rd.hstack(zeros(3,3,V).asDynamic).vstack(((VecToso3(p) * R).asDynamic).hstack(Rd))).asStatic(6,6)
    
func ScrewToAxis*[V](q:PositionVec[V], s:TriVector[V], h:V): HexVector[V] =
    ##[
        Takes a parametric description of a screw axis and converts it to a
    normalized screw axis
    :param q: A point lying on the screw axis
    :param s: A unit vector in the direction of the screw axis
    :param h: The pitch of the screw axis
    :return: A normalized screw axis described by the inputs
    Example Input:
        q = vector([3, 0, 0])
        s = vector([0, 0, 1])
        h = 2
    Output:
        vector([0, 0, 1, 0, -3, 2])
    ]##
    (s.asMatrix(3,1).asDynamic.hstack((cross(q,s) + (h * s)).asDynamic.asMatrix(3,1))).asVector.asStatic(6)

func AxisAng6*[V](expc6:HexVector[V]):(HexVector[V], V) =
    ##[
        Converts a 6-vector of exponential coordinates into screw axis-angle
    form
    :param expc6: A 6-vector of exponential coordinates for rigid-body motion
                  S*theta
    :return S: The corresponding normalized screw axis
    :return theta: The distance traveled along/about S
    Example Input:
        expc6 = vector([1, 0, 0, 1, 2, 3])
    Output:
        (vector([1.0, 0.0, 0.0, 1.0, 2.0, 3.0]), 1.0)
    ]##
    var theta : float
    theta = vector([expc6[0], expc6[1], expc6[2]]).norm
    if theta.NearZero:
        theta = vector([expc6[3], expc6[4], expc6[5]]).norm
    result = (expc6/theta, theta)

func MatrixExp6*[V](se3mat:QuadMatrix[V]):QuadMatrix[V] =
    ##[
        Computes the matrix exponential of an se3 representation of
    exponential coordinates
    :param se3mat: A matrix in se3
    :return: The matrix exponential of se3mat
    Example Input:
        se3mat = matrix( [[0,          0,           0,          0],
                          [0,          0, -1.57079632, 2.35619449],
                          [0, 1.57079632,           0, 2.35619449],
                          [0,          0,           0,          0]])
    Output:
        matrix([ [1.0, 0.0,  0.0, 0.0],
                 [0.0, 0.0, -1.0, 0.0],
                 [0.0, 1.0,  0.0, 3.0],
                 [  0,   0,    0,   1]])
    ]##
    let 
        subse3mat = se3mat[0..2,0..2]
        omgtheta = so3ToVec(subse3mat)
        se3matdyn = se3mat.asDynamic
        
    
    if omgtheta.norm.NearZero:
        return (eye(3,V).asDynamic.hstack(se3matdyn[0 .. 2, 3..3]).vstack(matrix(@[@[0.V,0,0,1]]))).asStatic(4,4)
    let 
        theta = AxisAng3(omgtheta)[1]
        omgmat = (subse3mat / theta).asDynamic

    var 
        tmp1 = MatrixExp3(subse3mat).asDynamic

        tmp2 =  theta * eye[V](3) + 
                (1 - cos(theta)) * omgmat + 
                (theta - sin(theta)) * (omgmat * omgmat)
        tmp3 = (se3matdyn[0 .. 2, 3..3] / theta)

    return tmp1.hstack(tmp2 * tmp3).vstack(matrix(@[@[0.V,0,0,1]])).asStatic(4,4)
    

func MatrixLog6*[V](T:QuadMatrix[V]):QuadMatrix[V] = 
    ##[
        Computes the matrix logarithm of a homogeneous transformation matrix
    :param R: A matrix in SE3
    :return: The matrix logarithm of R
    Example Input:
        T = np.array([[1, 0,  0, 0],
                      [0, 0, -1, 0],
                      [0, 1,  0, 3],
                      [0, 0,  0, 1]])
    Output:
        matrix([  [0,          0,           0,           0],
                  [0,          0, -1.57079633,  2.35619449],
                  [0, 1.57079633,           0,  2.35619449],
                  [0,          0,           0,           0]])
    ]##

    let 
        (R, p) = TransToRp(T)
        omgmat = MatrixLog3(R)
    if omgmat == statics.zeros(3,3,V):
        return (zeros[V](3,3).vstack(matrix(@[@[T[0,3],T[1,3],T[2,3]]])).hstack(matrix(@[@[0.0,0,0,0]]))).asStatic(4,4)
    else:
        let 
            theta = arccos((trace(R) - 1.V) / 2.0)
            invtheta = 1/theta
            tmp1 = (eye[V](3) - 
            (omgmat / 2.0).asDynamic) + 
            ((((invtheta - cot(theta / 2.V) / 2.V)) * (omgmat * omgmat * invtheta))).asDynamic
            tmp2 = (tmp1 * p.asDynamic).asMatrix(3,1)
            tmp3 = (omgmat.asDynamic.hstack(tmp2))
     
        result = tmp3.vstack(matrix(@[@[0.0,0,0,0]])).asStatic(4,4)
        

func ProjectToSO3*[V](mat:TriMatrix[V]): TriMatrix[V] =
    ##[
        Returns a projection of mat into SO(3)
    :param mat: A matrix near SO(3) to project to SO(3)
    :return: The closest matrix to R that is in SO(3)
    Projects a matrix mat to the closest matrix in SO(3) using singular-value
    decomposition (see
    http://hades.mech.northwestern.edu/index.php/Modern_Robotics_Linear_Algebra_Review).
    This function is only appropriate for matrices close to SO(3).
    Example Input:
        mat = matrix(  [[ 0.675,  0.150,  0.720],
                        [ 0.370,  0.771, -0.511],
                        [-0.630,  0.619,  0.472]])
    Output:
        matrix([  [ 0.67901136,  0.14894516,  0.71885945],
                  [ 0.37320708,  0.77319584, -0.51272279],
                  [-0.63218672,  0.61642804,  0.46942137]])
    ]##
    let 
        (U,s,Vh) = mat.asDynamic.svd
        R = (U * Vh)
    #if R.det < 0.V:
    #    R[_, s[2,2]] = -R[_, s[2,2]]

    return R.asStatic(3,3)

func ProjectToSE3*[V](mat:QuadMatrix[V]): QuadMatrix[V] =
    ##[Returns a projection of mat into SE(3)
    :param mat: A 4x4 matrix to project to SE(3)
    :return: The closest matrix to T that is in SE(3)
    Projects a matrix mat to the closest matrix in SE(3) using singular-value
    decomposition (see
    http://hades.mech.northwestern.edu/index.php/Modern_Robotics_Linear_Algebra_Review).
    This function is only appropriate for matrices close to SE(3).
    Example Input:
        mat = matrix([[ 0.675,  0.150,  0.720,  1.2],
                        [ 0.370,  0.771, -0.511,  5.4],
                        [-0.630,  0.619,  0.472,  3.6],
                        [ 0.003,  0.002,  0.010,  0.9]])
    Output:
        matrix([[ 0.67901136,  0.14894516,  0.71885945,  1.2 ],
                  [ 0.37320708,  0.77319584, -0.51272279,  5.4 ],
                  [-0.63218672,  0.61642804,  0.46942137,  3.6 ],
                  [ 0.        ,  0.        ,  0.        ,  1.  ]])
    ]##
    return RpToTrans(ProjectToSO3(mat[0..2, 0..2]), mat[0..2, 3..3].T.asVector)

func DistanceToSO3*[V](mat:TriMatrix[V]): V = 
    ##[Returns the Frobenius norm to describe the distance of mat from the
    SO(3) manifold
    :param mat: A 3x3 matrix
    :return: A quantity describing the distance of mat from the SO(3)
             manifold
    Computes the distance from mat to the SO(3) manifold using the following
    method:
    If det(mat) <= 0, return a large number.
    If det(mat) > 0, return norm(mat^T.mat - I).
    Example Input:
        mat = matrix([[ 1.0,  0.0,   0.0 ],
                        [ 0.0,  0.1,  -0.95],
                        [ 0.0,  1.0,   0.1 ]])
    Output:
        0.08835
    ]##
    let mat = mat.asDynamic
    if mat.det > 0.V:
        
        return norm(((mat.T * mat) - eye[V](3)))
    else:
        return 1e+9

func DistanceToSE3*[V](mat:QuadMatrix[V]): V =
    ##[Returns the Frobenius norm to describe the distance of mat from the
    SE(3) manifold
    :param mat: A 4x4 matrix
    :return: A quantity describing the distance of mat from the SE(3)
              manifold
    Computes the distance from mat to the SE(3) manifold using the following
    method:
    Compute the determinant of matR, the top 3x3 submatrix of mat.
    If det(matR) <= 0, return a large number.
    If det(matR) > 0, replace the top 3x3 submatrix of mat with matR^T.matR,
    and set the first three entries of the fourth column of mat to zero. Then
    return norm(mat - I).
    Example Input:
        mat = matrix([[ 1.0,  0.0,   0.0,   1.2 ],
                        [ 0.0,  0.1,  -0.95,  1.5 ],
                        [ 0.0,  1.0,   0.1,  -0.9 ],
                        [ 0.0,  0.0,   0.1,   0.98 ]])
    Output:
        0.134931
    ]##
    let mat = mat.asDynamic
    let matR = mat[0..2, 0..2]
    if det(matR) > 0.V:
        return norm((vstack(hstack((matR.T * matR),
                                          zeros[V](3, 1)),
                              mat[3..3, 0..3]) - eye[V](4)))
    else:
        return 1e+9

func TestIfSO3*(mat:TriMatrix): bool =
    ##[Returns true if mat is close to or on the manifold SO(3)
    :param mat: A 3x3 matrix
    :return: True if mat is very close to or in SO(3), false otherwise
    Computes the distance d from mat to the SO(3) manifold using the
    following method:
    If det(mat) <= 0, d = a large number.
    If det(mat) > 0, d = norm(mat^T.mat - I).
    If d is close to zero, return true. Otherwise, return false.
    Example Input:
        mat = matrix([  [1.0, 0.0,  0.0 ],
                        [0.0, 0.1, -0.95],
                        [0.0, 1.0,  0.1 ]])
    Output:
        false
    ]##
    return abs(DistanceToSO3(mat)) < 1e-3

func TestIfSE3*(mat:QuadMatrix): bool =
    ##[Returns true if mat is close to or on the manifold SE(3)
    :param mat: A 4x4 matrix
    :return: True if mat is very close to or in SE(3), false otherwise
    Computes the distance d from mat to the SE(3) manifold using the
    following method:
    Compute the determinant of the top 3x3 submatrix of mat.
    If det(mat) <= 0, d = a large number.
    If det(mat) > 0, replace the top 3x3 submatrix of mat with mat^T.mat, and
    set the first three entries of the fourth column of mat to zero.
    Then d = norm(T - I).
    If d is close to zero, return true. Otherwise, return false.
    Example Input:
        mat = matrix([[1.0, 0.0,   0.0,  1.2],
                        [0.0, 0.1, -0.95,  1.5],
                        [0.0, 1.0,   0.1, -0.9],
                        [0.0, 0.0,   0.1, 0.98]])
    Output:
        false
    ]##
    return abs(DistanceToSE3(mat)) < 1e-3
