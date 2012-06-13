#cobra.numjy.Matrix.py
#Basic matrix class that is going to be used to mimic numpy.ndarray
#capabilities.
#
#Derived from Simon Galbraith's ncajava Matrix.py file
import java, javax, jarray
from copy import deepcopy
from cern.colt.matrix import DoubleMatrix2D
from cern.colt.matrix.DoubleFactory2D import dense,sparse;
from cern.colt.matrix.impl import DenseDoubleMatrix2D as ndarray
from cern.colt.matrix.impl import SparseDoubleMatrix2D as sdarray
from cern.colt.matrix.linalg import Algebra;
from org.python.core.exceptions import ValueError as PyValueException;
from org.python.core import PyString,PySlice,PySequence,PyList;

class Matrix(javax.swing.table.AbstractTableModel):
    _M = None;   
    _name = 'data matrix'  
    varname = ''
    column_names = []
    row_names = []
           
    def __init__(self,M=None,N=None,v=None,sparse=None):

        """
        M is the number of rows
        N is the number of columns
        v is the default value
        """
        if isinstance(M,DoubleMatrix2D):
            self._M = M.copy();  
        elif (isinstance(M,int) and isinstance(N,int)):
            if sparse:
                F = sparse
            else:
                F = dense;
            if v is None:
                self._M = F.make(M,N,0);
            elif isinstance(v,int): 
                self._M = F.make(M,N,v);
            elif isinstance(v,PyList):
                self._M = F.make(jarray.array(v,'d'),1)
            elif isinstance(v,PyString):
                self._M = F.random(M,N);
            else:
                if sparse:
                    self._M = SparseDoubleMatrix2D(v)
                else:
                    self._M = ndarray(v)
        self.shape = (self._M.rows(),self._M.columns())


    def __copy__(self):
        r = new.instance(self.__class__, self.__dict__.copy() )
        r._M = self._M.copy();
        print "in copy"
        return r
    
    def __sub__(A,B):
        
        [ar,ac]=size(A);        
        C = Matrix(ar,ac,0);
        for i in range(1,ar):
            for j in range(1,ac):
                C[i,j]=A[i,j]-B[i,j]
        return C;
    
    def __mul__(A,B):  
        
        # need to check types and multiple based on them..     
        try:
            F = Algebra();
            C=(F.mult(A._M,B._M));
        except:
            raise PyValueException, "Inner dimension mismatch in matrix multiply.";
            return None;
        return Matrix(C) 
    def __div__(A,B):
        #print size(A)

        try:
            F = Algebra();
            R = F.solve(A._M,B._M);
            return R;
        except (java.lang.IllegalArgumentException) , e :
            # check the error class types according to the matrix class so we can intelligently report the error.
            print e.getMessage();
            return None;        
        
    
    def __repr__(self):
        return self._M.toString();
            
    def __str__(self):
        return self._M.toString();

    def __sz__(self):        
        if isinstance(self,Matrix):
            x=self._M.rows();
            y=self._M.columns();
            return (x,y);
        else:
            raise PyValueException, "Argument must be a matrix.";   
           
    def __setitem__(self,idx,v):

        if v is None:
            print idx
            raise PyValueException, "v is none"
        
        if isinstance(v,float): 
            self._M.set(idx[0],idx[1],v);
            return;
        
        Y=idx[1];
        X=idx[0];
            
        if isinstance(X,PyList):
            X=map(lambda x: x, X)
        elif isinstance(X,PySlice):
            if X.start == None:
               X=range(0,self._M.rows())
        elif isinstance(X,int):
            X=[X];
        
        if isinstance(Y,PyList):
            Y=map(lambda x: x,Y);
        elif isinstance(Y,PySlice):
            if Y.start == None:
               Y=range(0,self._M.cols())        
        elif isinstance(Y,int):
            Y=[Y];

        order=0;
        if len(X)>len(Y):
            order=1;
            
        #print "the order is " , order    
        if order==0:
            y=1
            for q in Y:
                    x=1
                    for z in X: 
         #               print z,q,x,y,v
                        self._M.set(z,q,v[x,y]);
                        x+=1;
                    y+=1;
        else:
            x=1
            for z in X:                    
                    y=1;
                    for q in Y:             
                        self._M.set(z,q,v[x,y]);
                        y+=1;
                    x+=1;
            
                
    def __getslice__(self, i, j):

        if i.start != None:
            x=range(i.start,i.stop);
        else:
            x=range(0,self._M.rows())
        if j.start !=None:    
            y=range(j.start,j.stop)
        else:
            y=range(0,self._M.columns())
        
        return Matrix(self._M.viewSelection(x,y));

         
    def __getitem__(self,idx):
        x=idx[0];
        y=idx[1];
        if x<0 or y<0:
            raise PyValueException, "Index must be positive number";     
          # this will fail on pyslice
        
        if isinstance(x,PySlice):
            if x.start != None:
               x=range(x.start,x.stop);
            else:
               x=range(0,self._M.rows())                   
        elif isinstance(x,int):
            x=x;
            x=[x]; 
        elif isinstance(x,PyList):
            x=map(lambda x: x, x)
       
        if isinstance(y,int):
               y=y;
               y=[y];
        elif isinstance(y,PySlice):
            if y.start != None:
                y=range(y.start,y.stop);
            else:
                y=range(0,self._M.columns())         
        elif isinstance(y,PySlice):
            if y.start !=None:    
                y=range(y.start,y.stop)
            else:
                y=range(0,self._M.columns())
        elif isinstance(y,PyList):
            y=map(lambda x: x, y)

        
        if len(x)<2 and len(y)<2:
            r = self._M.getQuick(x[0],y[0]);
            return float(r)  # this is a specific element
        else:
            return Matrix(self._M.viewSelection(x,y));
