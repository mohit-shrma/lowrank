import numpy as np 
import scipy as sp
import sys
from sparsesvd import sparsesvd

def ER1MP(Y_sp, rank):

  nrows, ncols = obsSpMat.shape
  X = np.zeros((nrows, ncols))
  
  #TODO
  (row_inds, col_inds) = Y_sp.nnzind
  
  #TODO: impleament proj operator
  Y_proj
  b = Y_sp[row_inds, col_inds]
  theta = np.zeros(rank)
  us = []
  vs = []
  for k in range(rank):
    #find top singular vector of Y_sp - X
    u, s, v = sp.linalg.svds(Y_proj - X, 1)
    us.append(u)
    vs.append(v)
    M_k = np.outer(u,v)
    
    #solve least square to find optimal alphas, \alpha_k
    #TODO:
    a1 = X.reshape(nrows*ncols, 1)
    a2 = M_k_proj.reshape(nrows*ncols, 1)
    b = Y_proj.reshape(nrows*ncols, 1)
    A = np.vstack((a1, a2))
    alpha1, alpha2 = np.linalg.lstsq(A.T, b)
    
    X = alpha1*X + alpha2*M_k_proj
    theta[k] = alpha2
    for i in range(k-1):
      theta[i] = theta[i]*alpha1

  #return output matrix
  opMat = np.zeros(nrows, ncols)
  for k in range(rank):
    opMat += theta[k]*np.outer(us[k],vs[k])
  return opMat


def main():
  
  #observed sparse matrix
  obsMatName = sys.argv[1]
  rank = int(sys.argv[2])
  




if __name__ == '__main__':
  main()


