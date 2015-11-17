import numpy as np 
import scipy as sp
import scipy.sparse.linalg as splinalg
import sys
from sparsesvd import sparsesvd


def getProjMat(mat, rowInds, colInds):
  nrows, ncols = mat.shape 
  mMat = np.zeros((nrows, ncols))
  mMat[rowInds, colInds] = 1 
  projMat = np.multiply(mat, mMat)
  return projMat


def ER1MP(Y, rank):

  nrows, ncols = Y.shape
  X = np.zeros((nrows, ncols))
  row_inds, col_inds = np.nonzero(Y)
  
  #get non-zero indices projection of Y
  Y_proj = getProjMat(Y, row_inds, col_inds)
  b = Y_proj.reshape(nrows*ncols, 1)
  theta = np.zeros(rank)
  us = []
  vs = []

  for k in range(rank):
    
    #find top singular vector of Y_sp - X
    u, s, v = splinalg.svds(Y_proj - X, 1)
    us.append(u)
    vs.append(v)
    M_k = np.outer(u,v)
    M_k_proj = getProjMat(M_k, row_inds, col_inds) 
    
    #solve least square to find optimal alphas, \alpha_k
    a1 = X.reshape(nrows*ncols, 1)
    a2 = M_k_proj.reshape(nrows*ncols, 1)
    A = np.hstack((a1, a2))
    sol = np.linalg.lstsq(A, b)
    alphas = sol[0]
    alpha1 = alphas[0]
    alpha2 = alphas[1]

    X = alpha1*X + alpha2*M_k_proj
    theta[k] = alpha2
    for i in range(k-1):
      theta[i] = theta[i]*alpha1

  #return output matrix
  opMat = np.zeros((nrows, ncols))
  for k in range(rank):
    opMat += theta[k]*np.outer(us[k],vs[k])

  return opMat


def main():
  
  #observed sparse matrix
  obsMatName = sys.argv[1]
  rank = int(sys.argv[2])
  

if __name__ == '__main__':
  main()


