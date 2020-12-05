import neo, neo/statics
import sequtils
import strutils

type 
    SquareMatrix*[M:static[int];A] = StaticMatrix[M,M,A]
    TriMatrix*[A] = SquareMatrix[3,A]
    QuadMatrix*[A] = SquareMatrix[4,A]
    HexMatrix*[A] = SquareMatrix[6,A]
    RotationMatrix* = TriMatrix
    TransformationMatrix* = QuadMatrix
    HomogenousMatrix* = QuadMatrix
    so3mat* = TriMatrix

    TriVector*[V] = StaticVector[3,V]
    HexVector*[V] = StaticVector[6,V]
    Omega* = TriVector
    PositionVec* = TriVector

    AxisAng*[V] = (TriVector[V], V)

    Matrices* = Matrix or StaticMatrix
    Vectors* = Vector or StaticVector


proc truncate(x:float):float = 
  if (x > 1e-7) or (x < -1e-7):
    result = x

proc numericalformatter*[T](mat:T, precision=2,trim=false,center=10):string =
  mat.mapit( (it.truncate.formateng(precision, trim=trim)).center(center) ).join(" ")

proc pretty*[N,M,V](mat:StaticMatrix[N,M,V],precision=2,trim=false,center=10):string =

  for val in mat.rows:
    result.add "|"
    
    result.add val.numericalformatter(precision,trim,center)
    result.add "|\n"

proc pretty*[N,V](mat:StaticVector[N,V],precision=2,trim=false,center=10):string =

    result.add "|"    
    result.add mat.numericalformatter(precision,trim,center)
    result.add "|\n"

proc pretty*[V](mat:(RotationMatrix[V],PositionVec[V]),precision=2,trim=false,center=10):string =
  result.add "Rotation Matrix:\n"
  result.add mat[0].pretty
  result.add "\nPosition Vector\n"
  result.add mat[1].pretty

type poly* {.pure.} = enum
  cubic, quintic 

template norm*(mat:Matrices or Vectors):untyped =
    mat.l_2

template norm*[V](arr:array[3,V]):untyped =
  Vector(arr).l_2

template trace*(mat:Matrices):untyped =
    mat.tr

func cross*[V](v1,v2:TriVector[V]): TriVector[V] =
    ##[
      implements the cross product for 3-d vectors
    ]##
    #let (v1,v2) = (v1.asDynamic,v2.asDynamic)
 
    result = [ v1[1]*v2[2] - v1[2]*v2[1],
               v1[2]*v2[0] - v1[0]*v2[2],
               v1[0]*v2[1] - v1[1]*v2[0]].vector

template `><`*[V](v1,v2:Trivector[V]):untyped =
  ##[
    Convenient template for making a cross product operator that actually 
    looks like a cross product
  ]##
  cross(v1,v2)

proc outer*[V](u,v:Vector[V]):Matrix[V] =
  
  u.asMatrix(u.len,1) * v.asMatrix(1,v.len)

proc outer*[N,M,V](u:StaticVector[N,V],v:StaticVector[M,V]):StaticMatrix[N,M,V] =
  outer(u.asDynamic,v.asDynamic).asStatic(N,M)

func NearZero*(f:float):bool =
    ##[
        Determines whether a scalar is small enough to be treated as zero
    :param f: A scalar input to check
    :return: True if z is close to zero, false otherwise
    Example Input:
        f = -1e-7
    Output:
        True
    ]##
    abs(f) < 1e-6

func Normalize*[V](v:TriVector[V]): TriVector[V] =
    ##[
        Normalizes a vector
    :param v: 3 vector
    :return: A unit vector pointing in the same direction as v
    Example Input:
        v = vector(@[1, 2, 3])
    Output:
        vector([0.26726124, 0.53452248, 0.80178373])
    
    ]##
    v/v.norm


func pinv*[V](G:Matrix[V]):Matrix[V] =
  ##[
    calculates the Moore-Penrose Inverse of a matrix based on the research paper
    obtained at https://arxiv.org/ftp/arxiv/papers/0804/0804.4809.pdf

  ]##
  var 
    m = G.M
    n = G.N
    transpose=false
    A:Matrix[V]
  if m<n:
    transpose=true
    A=G * G.T
    n=m
  else:
    A=G.T * G
  var 
    dA=diag(A.data)
    tol = min(da) * 1e-9
    L = zeros(A.M,A.N)
    r = 0

  for k in 0 .. n - 1:

    if r == 0:
      L[k..n - 1,r..r] = A[k .. n-1,k..k]
    else:
      
      L[k..n - 1,r..r] = A[k .. n-1,k..k] - L[k .. n - 1,0 .. (r-1)] * L[k..k, 0 .. (r-1)].T
    
    if L[k,r]>tol:
      L[k,r]=sqrt(L[k,r])
      if k < n - 1:
        L[(k+1)..n - 1,r..r] = L[(k+1)..n - 1,r..r]/L[k,r]
    else:
      dec r
    inc r

  L=L[ALL,0..r-1]
  var M=inv(L.T*L);
  if transpose:
    result=G.T*L*M*M*L.T
  else:
   result=L*M*M*L.T*G.T
  
proc pinv*[M,N,V](A:StaticMatrix[M,N,V]):Matrix[V] =
  A.asDynamic.pinv


iterator reversecolumns*[A](m: Matrix[A]): auto {. inline .} =
  let
    mp = cast[CPointer[A]](m.fp)
    step = if m.order == colMajor: m.ld else: 1
  var v = m.column(m.N - 1)
  yield v
  for j in countdown(m.N-2, 0):
    v.fp = addr(mp[j * step])
    yield v


when isMainModule:

  let
    a = statics.randomMatrix(3, 3)
    b = eye(5)
    c = matrix(@[@[1.0,2,3],
                 @[4.0,5,6],
                 @[7.0,8,9]])
    d = matrix(@[@[7.0, 2], @[3.0, 4], @[5.0, 3]])
    e = vector(@[1.0, 2, 3])#.asDynamic
    f = vector(@[3.0, 2, 1])#.asDynamic
  #echo e.data.diag
  echo pinv(d)
  echo outer(e,f)