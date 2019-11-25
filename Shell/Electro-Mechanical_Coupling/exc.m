function [coord,etpl,fext,bc,ngp,lstps,NRitmax,NRtol,ndim,E,nu,ka,oVn]=exc
ngp=18; NRitmax=50; NRtol=1e-6; ndim=5; lstps=10;
len=10; dep=1; wid=1; nen=9;
E=1.2e6; nu=0; ka=5/6;
lay=1; E=E*ones(1,lay);
nelx=10; nely=2; nels=nelx*nely;
coord=zeros((nelx*2+1)*(nely*2+1),3+lay); coord(:,4:4+lay-1)=dep/lay;
etpl=zeros(nels,9); k=1;
for i=1:(nelx*2+1)
  for j=1:(nely*2+1)
    coord(k,1:2)=[(i-1)*len/(2*nelx) (j-1)*wid/(2*nely)]; k=k+1;
  end
end
k=1;
for i=1:nelx
  for j=1:nely
    n1=(i-1)*2*(nely*2+1)+(j-1)*2+1;
    n2=n1+1; n3=n2+1; n4=n3+(nely*2+1); n5=n4+(nely*2+1);
    n6=n5-1; n7=n6-1; n8=n4-2; n9=n8+1;
    etpl(k,:)=[n1 n2 n3 n4 n5 n6 n7 n8 n9];
    oVn(:,:,k)=[zeros(nen,2) ones(nen,1)]; k=k+1;
  end
end
nodes=size(coord,1);
fext=zeros(nodes*ndim,1);
bc=zeros(nodes*5,2); k=0;
for i=1:nodes
  if coord(i,1)==0
    bc(i*5-4,:)=[i*5-4 0]; k=k+1;
    bc(i*5-3,:)=[i*5-3 0]; k=k+1; 
    bc(i*5-2,:)=[i*5-2 0]; k=k+1;
    bc(i*5-1,:)=[i*5-1 0]; k=k+1;
    bc(i*5  ,:)=[i*5   0]; k=k+1;
  end
end
bc=sortrows(bc,1);
bc(1:nodes*5-k,:)=[];