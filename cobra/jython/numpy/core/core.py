# (c) Simon J Galbraith 2007.
# Daniel Hyduke 2010
import java, jarray
from copy import deepcopy
from cern.colt.list import IntArrayList, DoubleArrayList
from cern.colt.matrix import DoubleFactory2D
from cern.colt.matrix.impl import DenseDoubleMatrix2D
from cern.colt.matrix.impl import SparseDoubleMatrix2D
from cern.colt.matrix.linalg import Algebra, EigenvalueDecomposition, LUDecomposition;
from org.python.core.exceptions import ValueError as PyValueException;
from cern.colt.matrix.doublealgo.Statistic import covariance, correlation;
from cern.jet.math.Functions import abs as cern_abs
from multiarray import ndarray, array
def mean(A):
    """Calculate the mean of Matrix object A

    """
    [r,c]=size(A)
    s=0
    for j in range(1,c):
        s = s+A[r,j]
    return s/c

def sum(A):
    """Calculate the sum of Matrix object A

    """
    return A._M.zSum()


def norm(A,ntype=None):
    F = Algebra();
    if ntype=='fro':
        r=F.normF(A._M);
    elif ntype == 2:
        r=F.norm2(A._M);
    else:
        r=F.norm2(A._M);
    return r;

def rank(A):
    
    if isinstance(A,ndarray):
        F = Algebra();
        r=F.rank(A._M);
        return int(r);
    else:
        raise PyValueException, "Rank function can only be called on matrix objects"

def cond(A):
    F = Algebra();
    return F.cond(A._M);

def size(A):
    return A.__sz__();

def transpose(A):
    F = Algebra();
    if isinstance(A,float):
        return A;
    else:
        return ndarray(F.transpose(A._M));

def inverse(A):
    F = Algebra();
    x=F.inverse(A._M);
    return ndarray(x)


def eig(A):
        # check that _M is square
    try:
        E = EigenvalueDecomposition(A._M);
        U = E.getV();
        eigs = E.getD();
    except PyValueException, e:
        print e;
        raise PyValueException,"Error in eig, check warnings()";  
    return [ndarray(eigs),ndarray(U)];

def solve(A, B):
    F = Algebra();
    if isinstance(A, ndarray) and isinstance(B, float):
            return F.solve(A._M, B);
    elif isinstance(B, ndarray) and  isinstance(A, float):
            return F.solve(A, B._M);
    elif isinstance(A,ndarray) and isinstance(B, ndarray):
            return ndarray(F.solve(A._M, B._M))
    else:
        return A / B

def solve_transpose(A, B):
    F = Algebra();
    if isinstance(A, ndarray) and isinstance(B, float):
            return F.solveTranspose(A._M,B);
    elif isinstance(B, ndarray) and  isinstance(A, float):
            return F.solveTranspose(A, B._M);
    elif isinstance(A, ndarray) and isinstance(B, ndarray):
            return ndarray(F.solveTranspose(A._M, B._M));
    else:
        return A / B

def solve_LR(A, B):
    T = A._M.copy()
    F = LUDecomposition(T);
    if isinstance(A, ndarray) and isinstance(B, float):
            return F.solve(B);
    elif isinstance(A, ndarray) and isinstance(B, ndarray):
        C = F.solve(B._M);
        return ndarray(C)

def cov(A):
    return ndarray(covariance(A._M))

def cor(A):  # handle multidimensional matrix case
    B = cov(A);
    return  ndarray(correlation(B._M))


def abs(A):
    F = cern_abs
    if isinstance(A,float):
        return java.lang.Math.abs(A)
    else:
        X = A._M.assign(F)
        return ndarray(X);

def svd(A):
    X = SingularValueDecomposition(A._M)
    u = X.getU()
    v = X.getV()
    e = X.getS()
    return [ndarray(u), ndarray(e), ndarray(v)]

#TODO:
#Make sure all of the functions below are are defined.  To know
#how they should operate, look at the docstring:
# pydoc numpy.the_function
#
# or in ipython:
# from numpy import *
# ?the_function
#
#
#TODO: These are the java mappings to the type.  If there's any
#difficulties then check the cern.colt data types
int32 = int
int64 = long

def ones(shape, dtype=float, order='C'):
    """Return a new array of given shape and type, filled with ones.
    
     
    See Also
    --------
    zeros
    
    Examples
    --------
    >>> numjy.ones(5)
    array([ 1.,  1.,  1.,  1.,  1.])
    
    >>> numjy.ones((5,), dtype=numjy.int)
    array([1, 1, 1, 1, 1])
    
    >>> numjy.ones((2, 1))
    array([[ 1.],
           [ 1.]])
    
    >>> s = (2,2)
    >>> numjy.ones(s)
    array([[ 1.,  1.],
           [ 1.,  1.]])"""
    return(ndarray(shape[0], shape[1], 1))

def sign(x):
    """ sign(x[, out])
    
    Returns an element-wise indication of the sign of a number.
    
    The `sign` function returns ``-1 if x < 0, 0 if x==0, 1 if x > 0``.
    
    Parameters
    ----------
    x : array_like
      Input values.
    
    Returns
    -------
    y : ndarray
      The sign of `x`.
    
    Examples
    --------
    >>> numjy.sign([-5., 4.5])
    array([-1.,  1.])
    >>> numjy.sign(0)
    0
 
    """
    def sign_int(a_number):
        if a_number < 0:
            return(-1)
        elif a_number == 0:
            return(0)
        else:
            return(1)

    if hasattr(x, '__iter__' ) or isinstance(x, ndarray):
        return_value = array([map(sign_int, list(array(x)._M.toArray()))])
    else:
        #In the case the input is just a number return the sign of the number.
        return_value = sign_int(x)
    return(return_value)       


def vstack(tup):
    """     Stack arrays in sequence vertically (row wise).
    
    Take a sequence of arrays and stack them vertically to make a single
    array. Rebuild arrays divided by `vsplit`.
    Parameters
    ----------
    tup : sequence of ndarrays
        Tuple containing arrays to be stacked. The arrays must have the same
        shape along all but the first axis.
    Returns
    -------
    stacked : ndarray
        The array formed by stacking the given arrays.
    See Also
    --------
    hstack : Stack arrays in sequence horizontally (column wise).
    dstack : Stack arrays in sequence depth wise (along third dimension).
    concatenate : Join a sequence of arrays together.
    vsplit : Split array into a list of multiple sub-arrays vertically.
    
    
    Notes
    -----
    Equivalent to ``np.concatenate(tup, axis=0)``
    
    Examples
    --------
    >>> a = np.array([1, 2, 3])
    >>> b = np.array([2, 3, 4])
    >>> np.vstack((a,b))
    array([[1, 2, 3],
           [2, 3, 4]])
    
    >>> a = np.array([[1], [2], [3]])
    >>> b = np.array([[2], [3], [4]])
    >>> np.vstack((a,b))
    array([[1],
           [2],
           [3],
           [2],
           [3],
           [4]])

           """
    if isinstance(tup[0], sdarray):
        matrix_factory = DoubleFactory2D.sparse
    else:
        #Allow for the case that python lists or tuples are fed to the function
        tup = map(array, tup)
        matrix_factory = DoubleFactory2D.dense

    stacked_matrix = tup[0]._M
    for the_array in tup[1:]:
        stacked_matrix = matrix_factory.appendRows(stacked_matrix, the_array._M)
    return(ndarray(stacked_matrix))


def hstack(tup):
    """
    Stack arrays in sequence horizontally (column wise).

    Take a sequence of arrays and stack them horizontally to make a single array. Rebuild arrays divided by hsplit.

    Parameters:
    tup : sequence of ndarrays
    All arrays must have the same shape along all but the second axis.
    Returns:
    stacked : ndarray
    The array formed by stacking the given arrays.
    See also
    vstack
    Stack arrays in sequence vertically (row wise).
    dstack
    Stack arrays in sequence depth wise (along third axis).
    concatenate
    Join a sequence of arrays together.
    hsplit
    Split array along second axis.
    Notes

    Equivalent to np.concatenate(tup, axis=1)

    Examples

    >>> a = np.array((1,2,3))
    >>> b = np.array((2,3,4))
    >>> np.hstack((a,b))
    array([1, 2, 3, 2, 3, 4])
    >>> a = np.array([[1],[2],[3]])
    >>> b = np.array([[2],[3],[4]])
    >>> np.hstack((a,b))
    array([[1, 2],
           [2, 3],
           [3, 4]])
    """

    if isinstance(tup[0], sdarray):
        matrix_factory = DoubleFactory2D.sparse
    else:
        tup = map(array, tup)
        matrix_factory = DoubleFactory2D.dense

    hstacked_matrix = tup[0]._M
    for the_array in tup[1:]:
        hstacked_matrix = matrix_factory.appendColumns(hstacked_matrix, the_array._M)
    return(ndarray(hstacked_matrix))


def where (condition, x=None, y=None):
    """
    Return a masked array with elements from x or y, depending on condition.

    Returns a masked array, shaped like condition, where the elements
    are from `x` when `condition` is True, and from `y` otherwise.
    If neither `x` nor `y` are given, the function returns a tuple of
    indices where `condition` is True (the result of
    ``condition.nonzero()``).

    Parameters
    ----------
    condition : array_like, bool
        The condition to meet. For each True element, yield the corresponding
        element from `x`, otherwise from `y`.
    x, y : array_like, optional
        Values from which to choose. `x` and `y` need to have the same shape
        as condition, or be broadcast-able to that shape.

    Returns
    -------
    out : MaskedArray or tuple of ndarrays
        The resulting masked array if `x` and `y` were given, otherwise
        the result of ``condition.nonzero()``.

    See Also
    --------
    numpy.where : Equivalent function in the top-level NumPy module.

    Examples
    --------
    >>> x = np.ma.array(np.arange(9.).reshape(3, 3), mask=[[0, 1, 0],
    ...                                                    [1, 0, 1],
    ...                                                    [0, 1, 0]])
    >>> print x
    [[0.0 -- 2.0]
     [-- 4.0 --]
     [6.0 -- 8.0]]
    >>> np.ma.where(x > 5)    # return the indices where x > 5
    (array([2, 2]), array([0, 2]))

    >>> print np.ma.where(x > 5, x, -3.1416)
    [[-3.1416 -- -3.1416]
     [-- -3.1416 --]
     [6.0 -- 8.0]]

    """

##     if x is None and y is None:
##         return filled(condition, 0).nonzero()
##     elif x is None or y is None:
##         raise ValueError, "Either both or neither x and y should be given."
##     # Get the condition ...............
##     fc = filled(condition, 0).astype(MaskType)
##     notfc = np.logical_not(fc)
##     # Get the data ......................................
##     xv = getdata(x)
##     yv = getdata(y)
##     if x is masked:
##         ndtype = yv.dtype
##     elif y is masked:
##         ndtype = xv.dtype
##     else:
##         ndtype = np.find_common_type([xv.dtype, yv.dtype], [])
##     # Construct an empty array and fill it
##     d = np.empty(fc.shape, dtype=ndtype).view(MaskedArray)
##     _data = d._data
##     np.putmask(_data, fc, xv.astype(ndtype))
##     np.putmask(_data, notfc, yv.astype(ndtype))
##     # Create an empty mask and fill it
##     _mask = d._mask = np.zeros(fc.shape, dtype=MaskType)
##     np.putmask(_mask, fc, getmask(x))
##     np.putmask(_mask, notfc, getmask(y))
##     _mask |= getmaskarray(condition)
##     if not _mask.any():
##         d._mask = nomask
##     return d








## def where(condition  , tup = 'blank', arg3 = 'blank'):








##     def conditioncheck(x):
##         if x == 'True':
##             return(1)
##         else:
##             return(0)

##     if isinstance(condition, bool):
##         truth_matrix = zeros(z.columns(), z.rows())
##         for i in xrange(truth_matrix.rows):
##             for j in xrange(truth_matrix.columns):
##                 truth_matrix.set(i,j, conditioncheck(z.get(i,j)))
##     else:
##         truth_matrix = zeros((len(condition), len(condition[0])))._M
##         for i in xrange(truth_matrix.rows()):
##             for j in xrange(truth_matrix.columns()):
##                 truth_matrix.set(i,j, conditioncheck(condition[i][j]))

##     if arg3 == 'blank':
##         if tup == 'blank':
##                 row_list = IntArrayList()
##                 column_list = IntArrayList()
##                 coordinate_list = DoubleArrayList()
##                 truth_matrix.getNonZeros(row_list, column_list, coordinate_list)
##                 return( (array(row_list), array(column_list)) )
##         else:
##             value_matrix = list(tup)
##             for i, x in enumerate(value_matrix):
##                 if isinstance(x, int):
##                     matfac = DoubleFactory2D.dense
##                     matrix = matfac.make(truth_matrix.rows(), truth_matrix.columns(), x)
##                     value_matrix[i] = matrix
##                 else:
##                     value_matrix[i] = x._M
##             return_array = zeros((truth_matrix.rows(), truth_matrix.columns()))._M
##             for i in xrange(truth_matrix.rows()):
##                 for j in xrange(truth_matrix.columns()):
##                     if truth_matrix.get(i,j) == 1:
##                         return_array.set(i,j, value_matrix[0].get(i,j))
##                     else:
##                         return_array.set(i,j, value_matrix[1].get(i,j))
##     else:
##         value_matrix = [ tup , arg3 ]
##         for i, x in enumerate(value_matrix):
##             if isinstance(x, int):
##                 matfac = DoubleFactory2D.dense
##                 matrix = matfac.make(truth_matrix.rows(), truth_matrix.columns(), x)
##                 value_matrix[i] = matrix
##             else:
##                 value_matrix[i] = x
##         return_array = zeros((truth_matrix.rows(), truth_matrix.columns()))._M
##         for i in xrange(truth_matrix.rows()):
##             for j in xrange(truth_matrix.columns()):
##                 if truth_matrix.get(i,j) == 1:
##                     print(type(value_matrix[0]))
##                     return_array.set(i,j, value_matrix[0].get(i,j))
##                 else:
##                     return_array.set(i,j, value_matrix[1].get(i,j))


##     return( return_array )

def nonzero(array):
    """Return the indices of the elements that are non-zero.
    Returns a tuple of arrays, one for each dimension of a, containing the indices of the non-zero elements
    in that dimension. The corresponding non-zero values can be obtained with:
    Parameters:	
    a : array_like
    Input array.
    Returns:	
    tuple_of_arrays : tuple
    Indices of elements that are non-zero.

    >>> x = np.eye(3)

    >>> x
    array([[ 1.,  0.,  0.],
    [ 0.,  1.,  0.],
    [ 0.,  0.,  1.]])

    >>> np.nonzero(x)
    (array([0, 1, 2]), array([0, 1, 2]))

    >>> x[np.nonzero(x)]
    array([ 1.,  1.,  1.])

    >>> np.transpose(np.nonzero(x))
    array([[0, 0],
       [1, 1],
       [2, 2]])"""
    
    rowList = IntArrayList()
    columnList = IntArrayList()
    coordinateList = DoubleArrayList()
    array._M.getNonZeros(rowList, columnList, coordinateList)
#TODO update array function to deal with Int....
    return(array(rowList), array(columnList))

def repeat(array, repeat_tup, axis = 2):
    """
    Repeat elements of an array.

    Parameters:
    a : array_like
    Input array.
    repeats : {int, array of ints}
    The number of repetitions for each element. repeats is broadcasted to fit the shape of the given axis.
    axis : int, optional
    The axis along which to repeat values. By default, use the flattened input array, and return a flat output array.
    Returns:
    repeated_array : ndarray
    Output array which has the same shape as a, except along the given axis.
    See also
    tile
    Tile an array.
    Examples

    >>> x = np.array([[1,2],[3,4]])
    >>> np.repeat(x, 2)
    array([1, 1, 2, 2, 3, 3, 4, 4])
    >>> np.repeat(x, 3, axis=1)
    array([[1, 1, 1, 2, 2, 2],
           [3, 3, 3, 4, 4, 4]])
    >>> np.repeat(x, [1, 2], axis=0)
    array([[1, 2],
           [3, 4],
           [3, 4]])
        """
    repeat = list(repeat_tup)
    #TODO: Make sure this can handle sparse matrices as well.  See the vstack function for ideas.
    matrix_factory = DoubleFactory2D.dense
    if axis == 1:
        repeat_matrix = matrix_factory.repeat(array._M.viewPart(0,0,array._M.rows(),1),1,repeat[0])
        for i in range( array._M.columns())[1:]:
            repeated_slice =  matrix_factory.repeat(array._M.viewPart(0, i, array._M.rows(), 1) , 1, repeat[i])
            repeat_matrix = matrix_factory.appendColumns(repeat_matrix, repeated_slice )
        return(repeat_matrix)
    elif axis == 0:
        repeat_matrix = matrix_factory.repeat(array._M.viewPart(0,0,1,array._M.columns()),repeat[0],1)
        for i in range(array._M.columns())[1:]:   
            repeat_matrix = matrix_factory.appendRows(repeat_matrix,
                                                      matrix_factory.repeat(array._M.viewPart(i,0,1,array._M.columns()),
                                                                            repeat[i], 1))
            return(repeat_matrix)
    else:
        pass





#One
#From scipy sparse.lil_matrix, sparse.hstack, sparse.vstack

def matrix():
    #This is not pressing.  I updated my cobra modules to use array instead of
    #matrix so this can be dealt with later.
    pass



    


