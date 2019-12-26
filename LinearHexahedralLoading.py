import numpy
def initialize():
    global E,v,coord,etpl,bc,f
    E=1e9; v=0; nels=50; coord=numpy.zeros([(nels-1)*4+8,3])
    len=1; wid=0.1; thk=0.1; etpl=numpy.zeros([nels+1,8]); bc=numpy.zeros([12,2])
    for nel in range(nels+1):
        coord[(nel)*4:(nel+1)*4,:]=[[nel*(len/nels),0,0],[nel*(len/nels),0,thk],[nel*(len/nels),wid,thk],[nel*(len/nels),wid,0]]
        n1=nel*4+1; n2=n1+1; n3=(nel+1)*4+2; n4=n3-1; n5=n2+2; n6=n5-1; n7=n3+1; n8=n7+1
        etpl[nel,:]=[n1,n2,n3,n4,n5,n6,n7,n8]
    etpl=etpl[:-1,:]
    bc[:,0]=numpy.array(range(1,13))
    f=numpy.zeros([coord.shape[0]*3,1])
    f[range(3*coord.shape[0]-10,3*coord.shape[0],3)]=-0.25e+3