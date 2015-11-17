import numpy as np 
import scipy as sp
import sys
from sparsesvd import sparsesvd

def ER1MP(Y_sp, rank):

  nrows, ncols = obsSpMat.shape
  X = np.zeros((nrows, ncols))
  
  for k in range(rank):
    #find top singular vector of Y_sp - X
    (u, v) = topSVD(Y_sp - X)
    M_k = np.outer(u,v)
    
    #solve least square to find optimal alphas, \alpha_k
    #TODO:
    np.linalg.lstsq(A,b)




def main():
  
  #observed sparse matrix
  obsMatName = sys.argv[1]
  rank = int(sys.argv[2])
  




if __name__ == '__main__':
  main()


