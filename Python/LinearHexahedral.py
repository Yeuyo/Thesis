import numpy
import LinearHexahedralLoading
LinearHexahedralLoading.initialize()
coord=LinearHexahedralLoading.coord; etpl=LinearHexahedralLoading.etpl.astype(int); bc=LinearHexahedralLoading.bc.astype(int); f=LinearHexahedralLoading.f; E=LinearHexahedralLoading.E; v=LinearHexahedralLoading.v
#coord=numpy.array([numpy.array([0, 0, 0]),numpy.array([0, 0, 1]),numpy.array([1,0,1]),numpy.array([1,0,0]),numpy.array([0,1,0]),numpy.array([0,1,1]),numpy.array([1,1,1]),numpy.array([1,1,0])])
#etpl=numpy.array([1,2,3,4,5,6,7,8])
#bc=numpy.array([numpy.array([1,0]),numpy.array([2,0]),numpy.array([3,0]),numpy.array([4,0]),numpy.array([5,0]),numpy.array([8,0]),numpy.array([11,0]),numpy.array([12,0]),numpy.array([13,0]),numpy.array([15,0]),numpy.array([16,0])])
#f=numpy.zeros([24,1]); f[23]=1; E=1.2e9; v=0
if etpl.ndim==1:
    etpl=numpy.expand_dims(etpl,0)
nodes=numpy.shape(coord)[0]; nDoF=3*nodes; nels=numpy.shape(etpl)[0]
eDoF=numpy.zeros([24,1]); DoF=numpy.reshape(list(range(1,nDoF+1)),[3,nodes],order='F').T
K=numpy.zeros([nDoF,nDoF]); d=numpy.zeros([nDoF,1])
Iz=numpy.zeros([6,6]); Iz[0:3,0:3]=numpy.ones([3,3])
I2=numpy.eye(6); I2[3:6,3:6]=0.5*numpy.eye(3)
D=E/((1+v)*(1-2*v))*(v*Iz+((1-2*v)*I2))
wp=numpy.ones([8,1]); dN=numpy.zeros([24,8])
gpcoord=numpy.zeros([nels,8,3]); BeGp=numpy.zeros([nels,8,6,24])
xsi=numpy.divide([-1,-1, 1, 1,-1,-1, 1, 1],numpy.sqrt(3))
eta=numpy.divide([-1,-1,-1,-1, 1, 1, 1, 1],numpy.sqrt(3))
zet=numpy.divide([-1, 1, 1,-1,-1, 1, 1,-1],numpy.sqrt(3))
N=numpy.zeros([8,8])
for n in range(8):
    sxsi=numpy.sign(xsi[n]); seta=numpy.sign(eta[n]); szet=numpy.sign(zet[n])
    for gp in range(8):
        N[gp,n]=(1+sxsi*xsi[gp])*(1+seta*eta[gp])*(1+szet*zet[gp])
        dN[3*(gp+1)-3,n]=sxsi*(1+seta*eta[gp])*(1+szet*zet[gp])
        dN[3*(gp+1)-2,n]=seta*(1+sxsi*xsi[gp])*(1+szet*zet[gp])
        dN[3*(gp+1)-1,n]=szet*(1+sxsi*xsi[gp])*(1+seta*eta[gp])
N=N/8; dN=dN/8
for nel in range (nels):
    for n in range (8):
        eDoF[3*(n+1)-3:3*(n+1)]=DoF[etpl[nel,n]-1,:][numpy.newaxis].T
    JT=numpy.dot(dN,coord[etpl[nel,:]-1,:]); ke=numpy.zeros([24,24])
    for n in range(8):
        indx=numpy.array(range(3*(n+1)-3,3*(n+1)))
        gpcoord[nel,n,0:3]=numpy.dot(N[n,:],coord[etpl[nel,:]-1,0:3])
        detJ=numpy.linalg.det(JT[indx,:]); dNx=numpy.linalg.lstsq(JT[indx,:],dN[indx,:],rcond=None)[0]; B=numpy.zeros([6,24])
        B[0,range(0,24,3)]=dNx[0,:]; B[1,range(1,24,3)]=dNx[1,:]; B[2,range(2,24,3)]=dNx[2,:]
        B[3,range(0,24,3)]=dNx[1,:]; B[3,range(1,24,3)]=dNx[0,:]; B[4,range(1,24,3)]=dNx[2,:]
        B[4,range(2,24,3)]=dNx[1,:]; B[5,range(0,24,3)]=dNx[2,:]; B[5,range(2,24,3)]=dNx[0,:]
        BeGp[nel,n,:,:]=B
        ke=ke+(B.T@D@B*detJ*wp[n])
    K[numpy.ix_(numpy.ravel((eDoF-1).astype(int)),numpy.ravel((eDoF-1).astype(int)))]=K[numpy.ix_(numpy.ravel((eDoF-1).astype(int)),numpy.ravel((eDoF-1).astype(int)))]+ke
fDoF=list(range(0,nDoF)); pDoF=bc[:,0]; d[pDoF]=numpy.expand_dims(bc[:,1],0).T
fDoF=[fDoF for fDoF in fDoF if fDoF not in pDoF-1]
d[fDoF]=numpy.linalg.lstsq(K[numpy.ix_(fDoF,fDoF)],f[fDoF]-(K[numpy.ix_(fDoF,pDoF)]@d[pDoF]),rcond=None)[0]
print(d)
