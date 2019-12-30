import numpy
import scipy.sparse
import scipy.sparse.linalg
import copy
def GNLcantilever_endload():
    NRitmax=50; NRtol=1.0e-6; ngp=8; lstps=20
    E=1.2e6; v=0; nels=300
    len=10; wid=1; thk=0.1
    coord=numpy.zeros([(nels-1)*4+8,3]); etpl=numpy.zeros([nels+1,8])
    Iz=numpy.zeros([6,6]); Iz[0:3,0:3]=numpy.ones([3,3])
    I2=numpy.eye(6); I2[3:6,3:6]=0.5*numpy.eye(3)
    D=numpy.zeros([6,6,ngp,nels+1]); bc=numpy.zeros([12,2])
    for nel in range(nels+1):
        coord[nel*4:(nel+1)*4,:]=numpy.array([numpy.array([nel*(len/nels),0,0]),numpy.array([nel*(len/nels),0,thk]),numpy.array([nel*(len/nels),wid,thk]),numpy.array([nel*(len/nels),wid,0])])
        n1=nel*4+1; n2=n1+1; n3=(nel+1)*4+2; n4=n3-1; n5=n2+2; n6=n5-1; n7=n3+1; n8=n7+1
        etpl[nel,:]=[n1,n2,n3,n4,n5,n6,n7,n8]
        for gp in range (ngp):
            D[:,:,gp,nel]=E/((1+v)*(1-2*v))*(v*Iz+((1-2*v)*I2))
    etpl=etpl[:-1,:]; D=D[:,:,:,:-1]
    bc[:,0]=numpy.array(range(0,12))
    fext=numpy.zeros([coord.shape[0]*3,1])
    fext[range(3*coord.shape[0]-10,3*coord.shape[0],3)]=1
    return coord, etpl, fext, bc, D, ngp, lstps, NRitmax, NRtol

def GpPos(ngp):
    if ngp==8: #2x2x2
        wp=numpy.ones([8,1]); g2=1/numpy.sqrt(3)
        xsi=numpy.array([-1,-1, 1, 1,-1,-1, 1, 1])*g2
        eta=numpy.array([-1,-1,-1,-1, 1, 1, 1, 1])*g2
        zet=numpy.array([-1, 1, 1,-1,-1, 1, 1,-1])*g2
        GpLoc=numpy.array([xsi,eta,zet]).T
    elif ngp==12: #2x2x3
        g2=1/numpy.sqrt(3); g3=numpy.sqrt(3/5); w1=8/9; w2=5/9
        xsi=numpy.array([-1,-1,-1, 1, 1, 1,-1,-1,-1, 1, 1, 1])*g2
        eta=numpy.array([-1,-1,-1,-1,-1,-1, 1, 1, 1, 1, 1, 1])*g2
        xsi=numpy.array([-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1])*g3
        wp=numpy.array([w2,w1,w2,w2,w1,w2,w2,w1,w2,w2,w1,w2]); GpLoc=numpy.array([xsi,eta,zet]).T
    elif ngp==18: #3x3x2
        w=numpy.zeros([18,2]); g3=numpy.sqrt(3/5); g2=1/numpy.sqrt(3); w1=8/9; w2=5/9
        xsi=   numpy.array([-1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 1, 1])*g3
        eta=   numpy.array([-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1])*g3
        zet=   numpy.array([-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1])*g3
        w[:,0]=numpy.array([w2,w2,w1,w1,w2,w2,w2,w2,w1,w1,w2,w2,w2,w2,w1,w1,w2,w2]).T
        w[:,1]=numpy.array([w2,w2,w2,w2,w2,w2,w1,w1,w1,w1,w1,w1,w2,w2,w2,w2,w2,w2]).T
        wp=w[:,0]*w[:,1]; GpLoc=numpy.array([xsi,eta,zet]).T
    elif ngp==27: #3x3x3
        w=numpy.zeros([27,3]); g3=numpy.sqrt(3/5); w1=8/9; w2=5/9
        xsi=numpy.array([-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1])*g3
        eta=numpy.array([-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1])*g3
        zet=numpy.array([-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1])*g3
        w[:,0]=numpy.array([w2,w1,w2,w2,w1,w2,w2,w1,w2, w2,w1,w2,w2,w1,w2,w2,w1,w2, w2,w1,w2,w2,w1,w2,w2,w1,w2])
        w[:,1]=numpy.array([w2,w2,w2,w2,w2,w2,w2,w2,w2, w1,w1,w1,w1,w1,w1,w1,w1,w1, w2,w2,w2,w2,w2,w2,w2,w2,w2])
        w[:,2]=numpy.array([w2,w2,w2,w1,w1,w1,w2,w2,w2, w2,w2,w2,w1,w1,w1,w2,w2,w2, w2,w2,w2,w1,w1,w1,w2,w2,w2])
        wp=w[:,0]*w[:,1]*w[:,2]; GpLoc=numpy.array([xsi,eta,zet]).T
    return wp,GpLoc

def formB(dNx,nen):
    B=numpy.zeros([6,nen*3])
    B[0,range(0,nen*3,3)]=dNx[0,:]; B[1,range(1,nen*3,3)]=dNx[1,:]; B[2,range(2,nen*3,3)]=dNx[2,:]
    B[3,range(0,nen*3,3)]=dNx[1,:]; B[3,range(1,nen*3,3)]=dNx[0,:]; B[4,range(1,nen*3,3)]=dNx[2,:]
    B[4,range(2,nen*3,3)]=dNx[1,:]; B[5,range(0,nen*3,3)]=dNx[2,:]; B[5,range(2,nen*3,3)]=dNx[0,:]
    return B

def dershapefunc2D(xsi,eta,zet):
    dNr=numpy.zeros([3,8]); r2=3
    dNr[range(0,r2,2),0]=-0.125*(1-eta)*(1-zet)
    dNr[range(0,r2,2),1]=-0.125*(1-eta)*(1+zet)
    dNr[range(0,r2,2),2]= 0.125*(1-eta)*(1+zet)
    dNr[range(0,r2,2),3]= 0.125*(1-eta)*(1-zet)
    dNr[range(0,r2,2),4]=-0.125*(1+eta)*(1-zet)
    dNr[range(0,r2,2),5]=-0.125*(1+eta)*(1+zet)
    dNr[range(0,r2,2),6]= 0.125*(1+eta)*(1+zet)
    dNr[range(0,r2,2),7]= 0.125*(1+eta)*(1-zet)
    dNr[range(1,r2,2),0]=-0.125*(1-xsi)*(1-zet)
    dNr[range(1,r2,2),1]=-0.125*(1-xsi)*(1+zet)
    dNr[range(1,r2,2),2]=-0.125*(1+xsi)*(1+zet)
    dNr[range(1,r2,2),3]=-0.125*(1+xsi)*(1-zet)
    dNr[range(1,r2,2),4]= 0.125*(1-xsi)*(1-zet)
    dNr[range(1,r2,2),5]= 0.125*(1-xsi)*(1+zet)
    dNr[range(1,r2,2),6]= 0.125*(1+xsi)*(1+zet)
    dNr[range(1,r2,2),7]= 0.125*(1+xsi)*(1-zet)
    dNr[range(2,r2,2),0]=-0.125*(1-xsi)*(1-eta)
    dNr[range(2,r2,2),1]= 0.125*(1-xsi)*(1-eta)
    dNr[range(2,r2,2),2]= 0.125*(1+xsi)*(1-eta)
    dNr[range(2,r2,2),3]=-0.125*(1+xsi)*(1-eta)
    dNr[range(2,r2,2),4]=-0.125*(1-xsi)*(1+eta)
    dNr[range(2,r2,2),5]= 0.125*(1-xsi)*(1+eta)
    dNr[range(2,r2,2),6]= 0.125*(1+xsi)*(1+eta)
    dNr[range(2,r2,2),7]=-0.125*(1+xsi)*(1+eta)
    return dNr

def TLFE(coord, D, duvw, ngp, epsEn, sigN, oduX):
    nen=numpy.shape(coord)[0]; epsE=numpy.zeros([6,ngp]); sig=numpy.zeros([6,ngp])
    wp,GpLoc=GpPos(ngp); ke=numpy.zeros([nen*3,nen*3]); Ks=numpy.zeros([nen*3,nen*3])
    fint=numpy.zeros([nen*3,1]); duX=numpy.zeros([3,3,ngp])
    for gp in range(ngp):
        xsi=GpLoc[gp,0]; eta=GpLoc[gp,1]; zet=GpLoc[gp,2]
        dNr=dershapefunc2D(xsi,eta,zet)
        dXr=dNr@coord; detJ=numpy.linalg.det(dXr); dNX=numpy.linalg.lstsq(dXr,dNr,rcond=None)[0]
        duX[:,:,gp]=dNX@numpy.reshape(duvw,[3,nen],order='F').T
        epsEt=numpy.array([0.5*(duX[0,0,gp] + duX[0,0,gp] + \
                                oduX[0,0,gp]* duX[0,0,gp] + oduX[0,1,gp]* duX[0,1,gp] + oduX[0,2,gp]* duX[0,2,gp] + \
                                 duX[0,0,gp]*oduX[0,0,gp] +  duX[0,1,gp]*oduX[0,1,gp] +  duX[0,2,gp]*oduX[0,2,gp] + \
                                 duX[0,0,gp]* duX[0,0,gp] +  duX[0,1,gp]* duX[0,1,gp] +  duX[0,2,gp]* duX[0,2,gp]), \
                           0.5*(duX[1,1,gp] + duX[1,1,gp] + \
                                oduX[1,0,gp]* duX[1,0,gp] + oduX[1,1,gp]* duX[1,1,gp] + oduX[1,2,gp]* duX[1,2,gp] + \
                                 duX[1,0,gp]*oduX[1,0,gp] +  duX[1,1,gp]*oduX[1,1,gp] +  duX[1,2,gp]*oduX[1,2,gp] + \
                                 duX[1,0,gp]* duX[1,0,gp] +  duX[1,1,gp]* duX[1,1,gp] +  duX[1,2,gp]* duX[1,2,gp]), \
                           0.5*(duX[2,2,gp] + duX[2,2,gp] + \
                                oduX[2,0,gp]* duX[2,0,gp] + oduX[2,1,gp]* duX[2,1,gp] + oduX[2,2,gp]* duX[2,2,gp] + \
                                 duX[2,0,gp]*oduX[2,0,gp] +  duX[2,1,gp]*oduX[2,1,gp] +  duX[2,2,gp]*oduX[2,2,gp] + \
                                 duX[2,0,gp]* duX[2,0,gp] +  duX[2,1,gp]* duX[2,1,gp] +  duX[2,2,gp]* duX[2,2,gp]), \
                           1.0*(duX[1,0,gp] + duX[0,1,gp] + \
                                oduX[0,0,gp]* duX[1,0,gp] + oduX[0,1,gp]* duX[1,1,gp] + oduX[0,2,gp]* duX[1,2,gp] + \
                                 duX[0,0,gp]*oduX[1,0,gp] +  duX[0,1,gp]*oduX[1,1,gp] +  duX[0,2,gp]*oduX[1,2,gp] + \
                                 duX[0,0,gp]* duX[1,0,gp] +  duX[0,1,gp]* duX[1,1,gp] +  duX[0,2,gp]* duX[1,2,gp]), \
                           1.0*(duX[1,2,gp] + duX[2,1,gp] + \
                                oduX[2,0,gp]* duX[1,0,gp] + oduX[2,1,gp]* duX[1,1,gp] + oduX[2,2,gp]* duX[1,2,gp] + \
                                 duX[2,0,gp]*oduX[1,0,gp] +  duX[2,1,gp]*oduX[1,1,gp] +  duX[2,2,gp]*oduX[1,2,gp] + \
                                 duX[2,0,gp]* duX[1,0,gp] +  duX[2,1,gp]* duX[1,1,gp] +  duX[2,2,gp]* duX[1,2,gp]), \
                           1.0*(duX[2,0,gp] + duX[0,2,gp] + \
                                oduX[0,0,gp]* duX[2,0,gp] + oduX[0,1,gp]* duX[2,1,gp] + oduX[0,2,gp]* duX[2,2,gp] + \
                                 duX[0,0,gp]*oduX[2,0,gp] +  duX[0,1,gp]*oduX[2,1,gp] +  duX[0,2,gp]*oduX[2,2,gp] + \
                                 duX[0,0,gp]* duX[2,0,gp] +  duX[0,1,gp]* duX[2,1,gp] +  duX[0,2,gp]* duX[2,2,gp])])
        epsE[:,gp]=epsEt+epsEn[:,gp]
        sig[:,gp]=(D[:,:,gp]@epsEt)+sigN[:,gp]
        s=numpy.zeros([3,3])
        s[0,:]=numpy.array([sig[0,gp],sig[3,gp],sig[5,gp]])
        s[1,:]=numpy.array([sig[3,gp],sig[1,gp],sig[4,gp]])
        s[2,:]=numpy.array([sig[5,gp],sig[4,gp],sig[2,gp]])
        BL=formB(dNX,nen); BNL=numpy.zeros([6,nen*3])
        tduX=oduX[:,:,gp]+duX[:,:,gp]
        A=numpy.zeros([6,9])
        A[0,:]=numpy.block([         tduX[0,:],numpy.zeros([1,3]),numpy.zeros([1,3])])
        A[1,:]=numpy.block([numpy.zeros([1,3]),         tduX[1,:],numpy.zeros([1,3])])
        A[2,:]=numpy.block([numpy.zeros([1,3]),numpy.zeros([1,3]),         tduX[2,:]])
        A[3,:]=numpy.block([         tduX[1,:],         tduX[0,:],numpy.zeros([1,3])])
        A[4,:]=numpy.block([numpy.zeros([1,3]),         tduX[2,:],         tduX[1,:]])
        A[5,:]=numpy.block([         tduX[2,:],numpy.zeros([1,3]),         tduX[0,:]])
        for n in range(nen):
            for m in range(nen):
                Ks[n*3:(n+1)*3,m*3:(m+1)*3]=dNX[:,n].T@s@dNX[:,m]*numpy.eye(3)
            BNL[:,n*3:(n+1)*3]=A@numpy.block([[dNX[0,n]*numpy.eye(3)],[dNX[1,n]*numpy.eye(3)],[dNX[2,n]*numpy.eye(3)]])
        Bt=BL+BNL; ke=ke+(Bt.T@D[:,:,gp]@Bt+Ks)*detJ*wp[gp]; fint=fint+(Bt.T@sig[:,gp]*detJ*wp[gp])[numpy.newaxis].T
    return ke, fint, epsE, sig, duX

coord,etpl,fext0,bc,D,ngp,lstps,NRitmax,NRtol=GNLcantilever_endload()
nels,nen=etpl.shape; nodes=coord.shape[0]
ndim=3; nDoF=nodes*ndim; neDoF=(nen*ndim)**2
krow=numpy.zeros([neDoF*nels,1]); kcol=numpy.zeros([neDoF*nels,1]); kval=numpy.zeros([neDoF*nels,1])
uvw=numpy.zeros([nDoF,1]); uvwold=numpy.zeros([nDoF,1]); fint=numpy.zeros([nDoF,1]); react=numpy.zeros([nDoF,1])
fd=list(range(0,nDoF)); fd=[fd for fd in fd if fd not in bc[:,0]]
epsE=numpy.zeros([6,ngp,nels]); epsEn=numpy.zeros([6,ngp,nels]); sig=numpy.zeros([6,ngp,nels]); sigN=numpy.zeros([6,ngp,nels])
oduX=numpy.zeros([3,3,ngp,nels]); duX=numpy.zeros([3,3,ngp,nels])
for lstp in range(lstps):
    fext=(lstp/lstps)*fext0; oobf=react+fext-fint; oobfnorm=2*NRtol; NRit=0
    while NRit<NRitmax and oobfnorm>NRtol:
        NRit=NRit+1; fint=numpy.zeros([nDoF,1]); dreact=numpy.zeros([nDoF,1]); dduvw=numpy.zeros([nDoF,1])
        if lstp>=1:
            Kt=scipy.sparse.csc_matrix((numpy.ravel(kval), (numpy.ravel(krow), numpy.ravel(kcol))), shape=(nDoF, nDoF), dtype=float)
            dduvw[bc[:,0].astype(int)]=((1+numpy.sign(1-NRit))*bc[:,1]/lstps)[numpy.newaxis].T #so the boundary displacement only apply a boundary once on the next line once, also won't add on later
            #dduvw[fd]=numpy.linalg.lstsq(Kt[numpy.ix_(fd,fd)].todense(),oobf[fd]-(Kt[numpy.ix_(fd,bc[:,0].astype(int))]@dduvw[bc[:,0].astype(int)]),rcond=None)[0]
            #dduvw[fd]=scipy.sparse.linalg.lsqr(Kt[numpy.ix_(fd,fd)],oobf[fd]-(Kt[numpy.ix_(fd,bc[:,0].astype(int))]@dduvw[bc[:,0].astype(int)]),iter_lim=None,atol=1e-15,btol=1e-15)
            dduvw[fd]=scipy.sparse.linalg.spsolve(Kt[numpy.ix_(fd,fd)],oobf[fd]-(Kt[numpy.ix_(fd,bc[:,0].astype(int))]@dduvw[bc[:,0].astype(int)]))[numpy.newaxis].T
            dreact[bc[:,0].astype(int)]=Kt[bc[:,0].astype(int),:]@dduvw-oobf[bc[:,0].astype(int)]
        uvw=uvw+dduvw; react=react+dreact; duvw=uvw-uvwold
        for nel in range(nels):
            ed=numpy.reshape((numpy.ones([ndim,1])*etpl[nel,:]*ndim-(numpy.array(range(ndim-1,-1,-1))[numpy.newaxis].T@numpy.ones([1,nen])))-1,[1,nen*ndim],order='F')
            ke,felem,epsE[:,:,nel],sig[:,:,nel],duX[:,:,:,nel]=TLFE(coord[etpl[nel,:].astype(int)-1,:],D[:,:,:,nel],duvw[ed.astype(int)],ngp,epsEn[:,:,nel],sigN[:,:,nel],oduX[:,:,:,nel])
            if lstp==0:
                krow[nel*neDoF:(nel+1)*neDoF]=numpy.reshape(ed.T*numpy.ones([1,nen*ndim]),[neDoF,1],order='F')
                kcol[nel*neDoF:(nel+1)*neDoF]=numpy.reshape(numpy.ones([nen*ndim,1])*ed,[neDoF,1],order='F')
            kval[nel*neDoF:(nel+1)*neDoF]=numpy.reshape(ke,[neDoF,1],order='F')
            fint[ed.astype(int)]=fint[ed.astype(int)]+felem
        oobf=fext+react-fint; oobfnorm=numpy.linalg.norm(oobf)/numpy.linalg.norm(fext+react+numpy.spacing(1))
        print("%.0f %.0f %.3e" % (lstp,NRit,oobfnorm))
    uvwold=copy.deepcopy(uvw); epsEn=copy.deepcopy(epsE); sigN=copy.deepcopy(sig); oduX=oduX+duX