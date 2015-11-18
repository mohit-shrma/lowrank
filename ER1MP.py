import numpy as np 
import scipy as sp
import scipy.sparse as sparse
import scipy.sparse.linalg as splinalg
import sys


def rmse(X, Y, rowInds, colInds):
  x = X[rowInds, colInds]
  y = Y[rowInds, colInds]
  nElems = len(rowInds)
  x = x.reshape(nElems, 1)
  y = y.reshape(nElems, 1)
  rmse = np.sqrt(np.sum(np.square(x-y))/nElems)
  return rmse


def getProjMat(mat, rowInds, colInds):
  nrows, ncols = mat.shape 
  mMat = np.zeros((nrows, ncols))
  mMat[rowInds, colInds] = 1 
  projMat = np.multiply(mat, mMat)
  return projMat


def SVP(Y, rank, delta, maxIter, trainRowInds, trainColInds, testRowInds,
    testColInds):
  
  nrows, ncols = Y.shape
  X_star       = getProjMat(Y, trainRowInds, trainColInds)
  nElems       = len(trainRowInds)
  
  print 'Size of indices: ', nElems

  X            = np.zeros((nrows, ncols))
  p            = (nElems*1.0)/(nrows*ncols)
 
  rate = 1.0/((1.0 + delta)*p)
  print 'p: ', p, 'delta: ', delta, 'Rate: ', rate

  for t in range(maxIter):
    X_proj = getProjMat(X, trainRowInds, trainColInds)
    diffX = X -  ((X_proj - X_star)*(rate))
    
    #perform orthogonal projection of above
    U, s, Vh = splinalg.svds(diffX, rank)
    X = np.dot(U, np.dot(np.diag(s), Vh))

    trainErr = rmse(X, Y, trainRowInds, trainColInds)
    
    testErr = -1
    if len(testRowInds) > 0:
      testErr  = rmse(X, Y, testRowInds, testColInds)
    
    if t % 10 == 0:
      print 't: ', t, 'Train Error: ', trainErr, 'Test Err: ', testErr
    
    if trainErr < 0.001:
      break

  return X


def ER1MP(Y, rank, row_inds, col_inds):

  nrows, ncols = Y.shape
  X = np.zeros((nrows, ncols))
  
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
    alpha1 = alphas[0,0]
    alpha2 = alphas[1,0]

    X = alpha1*X + alpha2*M_k_proj
    theta[k] = alpha2
    for i in range(k-1):
      theta[i] = theta[i]*alpha1

  #return output matrix
  opMat = np.zeros((nrows, ncols))
  for k in range(rank):
    opMat += theta[k]*np.outer(us[k],vs[k])
  
  return opMat


def genLowRankMat(nrows, ncols, rank):
  A = np.random.rand(nrows, ncols)
  U,s,Vh = np.linalg.svd(A)

  A_loRank = np.dot(U[:,0:rank], np.dot(np.diag(s[:rank]), Vh[0:rank, :]))
  return A_loRank


def readTriplet0IndMat(fileName):
  ipMat = np.loadtxt(fileName)
  rows, cols, data = ipMat.T
  spMat = sparse.coo_matrix((data, (rows, cols)))
  return spMat


def main():
  
  #observed sparse matrix
  obsMatName = sys.argv[1]
  rank = int(sys.argv[2])
  

if __name__ == '__main__':
  main()


