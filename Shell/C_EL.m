function [coord,etpl,fext,bc,ngp,lstps,NRitmax,NRtol,ndim,D,dxr]=C_EL
ngp=8; NRitmax=50; NRtol=1e-6; ndim=5; lstps=20;
nels=8; len=10; dep=0.1; wid=1; lay=1;
E=1.2e6; nu=0; k=5/6; nen=9;
E=E*ones(1,lay); oVn=zeros(nen,3,nels);
coord=zeros((nels*2+1)*3,3+lay); D=zeros(6,6,ngp,lay,nels);
etpl=zeros(nels,nen); dxr=zeros(3,3,ngp,lay,nels);
coord(:,4:4+lay-1)=dep/lay;
for i=1:(nels*2+1)
  coord((i-1)*3+1,1:2)=[(i-1)*len/(2*nels) 0];  
  coord((i-1)*3+2,1:2)=[(i-1)*len/(2*nels) wid/2];
  coord((i-1)*3+3,1:2)=[(i-1)*len/(2*nels) wid];
end
for i=1:nels
  n1=(i-1)*6+1; n2=n1+1; n3=n1+2; n4=n1+5; n5=n1+8;
  n6=n1+7; n7=n1+6; n8=n1+3; n9=n1+4;
  etpl(i,:)=[n1 n2 n3 n4 n5 n6 n7 n8 n9];
  oVn(:,:,i)=[zeros(nen,2) ones(nen,1)];
end
nodes=size(coord,1);
fext=zeros(nodes*ndim,1);
fext( nodes*5   -2)=04*1/6;
fext((nodes-1)*5-2)=04*2/3;
fext((nodes-2)*5-2)=04*1/6;
bc=zeros(nodes*5,2); n=0;
for i=1:nodes
  if coord(i,1)==0
    bc(i*5-4,:)=[i*5-4 0]; n=n+1;
    bc(i*5-3,:)=[i*5-3 0]; n=n+1; 
    bc(i*5-2,:)=[i*5-2 0]; n=n+1;
    bc(i*5-1,:)=[i*5-1 0]; n=n+1;
    bc(i*5  ,:)=[i*5   0]; n=n+1;
  end
end
bc=sortrows(bc,1);
bc(1:nodes*5-n,:)=[];

ex=[1 0 0].'; ey=[0 1 0].'; ez=[0 0 1].';
for nel=1:nels
  elcoord=coord(etpl(nel,:),1:3);
  t=coord(etpl(nel,:),4:4+lay-1); T=sum(t,2);
  [~,GpLoc]=GpPos(ngp);
  for n=1:lay
    Dn=(E(n)/(1-nu^2))*[[1 nu; nu 1; 0 0] zeros(3,4);...
      zeros(3,3) eye(3)*k*(1-nu)*0.5];
    Dn(4,4)=Dn(4,4)/k;
    for gp=1:ngp
      xsi=GpLoc(gp,1); eta=GpLoc(gp,2); zeta=GpLoc(gp,3);
      N  =shapefunc(xsi,eta);
      zet=-1+(2*sum(N*t(:,1:n))-N*(t(:,n).*(1-zeta)))/(N*T);
      dNr=dershapefunc2D(xsi,eta);
      dxr(1:2,:,gp,n,nel)=(dNr*(elcoord+0.5*zet*(oVn(:,:,nel).*(T*ones(1,3)))));
      dxr(3,:,gp,n,nel)=(0.5*(N.*T')*oVn(:,:,nel));
      VnGp=(N*oVn(:,:,nel)).';
      sdir=dxr(2,:,gp,n,nel)'/norm(dxr(2,:,gp,n,nel)');
      er=cross(sdir,VnGp); er=er/norm(er);
      es=cross(VnGp,er);   es=es/norm(es);
      et=VnGp/norm(VnGp);
      l1=(ex.'*er); m1=(ey.'*er); n1=(ez.'*er);
      l2=(ex.'*es); m2=(ey.'*es); n2=(ez.'*es);
      l3=(ex.'*et); m3=(ey.'*et); n3=(ez.'*et);
      Q=[l1*l1     m1*m1   n1*n1 l1*m1       m1*n1       n1*l1       ;
        l2*l2     m2*m2   n2*n2 l2*m2       m2*n2       n2*l2       ;
        l3*l3     m3*m3   n3*n3 l3*m3       m3*n3       n3*l3       ;
        2*l1*l2 2*m1*m2 2*n1*n2 l1*m2+l2*m1 m1*n2+m2*n1 n1*l2+n2*l1 ;
        2*l2*l3 2*m2*m3 2*n2*n3 l2*m3+l3*m2 m2*n3+m3*n2 n2*l3+n3*l2 ;
        2*l3*l1 2*m3*m1 2*n3*n1 l3*m1+l1*m3 m3*n1+m1*n3 n3*l1+n1*l3];
      D(:,:,gp,n,nel)=Q'*Dn*Q;
    end
  end
end

function [wp,GpLoc]=GpPos(ngp)
if ngp==8 %2x2x2
  wp=ones(8,1); g2=1/sqrt(3); 
  xsi=[-1 -1 1 1 -1 -1 1 1].'*g2;
  eta=[-1 -1 -1 -1 1 1 1 1].'*g2;
  zet=[-1 1 1 -1 -1 1 1 -1].'*g2;
  GpLoc=[xsi eta zet];
elseif ngp==4 %2x2x1
  wp=ones(4,1); g2=1/sqrt(3); 
  xsi=[-1  1 -1  1].'*g2;
  eta=[-1 -1  1  1].'*g2;
  zet=[ 0  0  0  0].'*g2;
  GpLoc=[xsi eta zet];
elseif ngp==12 %2x2x3
  g2=1/sqrt(3); g3=sqrt(3/5); w1=8/9; w2=5/9;
  xsi=[-1 -1 -1  1  1  1 -1 -1 -1  1  1  1].'*g2;
  eta=[-1 -1 -1 -1 -1 -1  1  1  1  1  1  1].'*g2;
  zet=[-1  0  1 -1  0  1 -1  0  1 -1  0  1].'*g3;
  wp=[w2;w1;w2;w2;w1;w2;w2;w1;w2;w2;w1;w2];
  GpLoc=[xsi eta zet];
elseif ngp==18 %3x3x2
  g3=sqrt(3/5); g2=1/sqrt(3); w1=8/9; w2=5/9;
  xsi=   [-1;-1; 0; 0; 1; 1;-1;-1; 0; 0; 1; 1;-1;-1; 0; 0; 1; 1]*g3;
  eta=   [-1;-1;-1;-1;-1;-1; 0; 0; 0; 0; 0; 0; 1; 1; 1; 1; 1; 1]*g3;
  zet=   [-1; 1;-1; 1;-1; 1;-1; 1;-1; 1;-1; 1;-1; 1;-1; 1;-1; 1]*g2;
  w(:,1)=[w2;w2;w1;w1;w2;w2;w2;w2;w1;w1;w2;w2;w2;w2;w1;w1;w2;w2];
  w(:,2)=[w2;w2;w2;w2;w2;w2;w1;w1;w1;w1;w1;w1;w2;w2;w2;w2;w2;w2];
  wp=w(:,1).*w(:,2);
  GpLoc=[xsi eta zet];
end

function [N]=shapefunc(xsi,eta)
N(:,1)= 0.25*xsi.*(xsi-1).*eta.*(eta-1);
N(:,2)=-0.50*xsi.*(xsi-1).*(eta+1).*(eta-1);
N(:,3)= 0.25*xsi.*(xsi-1).*eta.*(eta+1);
N(:,4)=-0.50*(xsi+1).*(xsi-1).*eta.*(eta+1);
N(:,5)= 0.25*xsi.*(xsi+1).*eta.*(eta+1);
N(:,6)=-0.50*xsi.*(xsi+1).*(eta+1).*(eta-1);
N(:,7)= 0.25*xsi.*(xsi+1).*eta.*(eta-1);
N(:,8)=-0.50*(xsi+1).*(xsi-1).*eta.*(eta-1);
N(:,9)=(xsi+1).*(xsi-1).*(eta+1).*(eta-1);

function [dNr]=dershapefunc2D(xsi,eta)
ngp=size(xsi,1); r2=ngp*2;
dNr(1:2:r2  ,1)= 1/4*eta.*(eta-1).*(2*xsi-1);
dNr(1:2:r2  ,2)=-1/2*(eta+1).*(eta-1).*(2*xsi-1);
dNr(1:2:r2  ,3)= 1/4*eta.*(eta+1).*(2*xsi-1);
dNr(1:2:r2  ,4)=-1*eta.*(eta+1).*xsi;
dNr(1:2:r2  ,5)= 1/4*eta.*(eta+1).*(2*xsi+1);
dNr(1:2:r2  ,6)=-1/2*(eta+1).*(eta-1).*(2*xsi+1);
dNr(1:2:r2  ,7)= 1/4*eta.*(eta-1).*(2*xsi+1);
dNr(1:2:r2  ,8)=-1*eta.*(eta-1).*xsi;
dNr(1:2:r2  ,9)= 2*(eta+1).*(eta-1).*xsi;
dNr(2:2:r2+1,1)= 1/4*xsi.*(xsi-1).*(2*eta-1);
dNr(2:2:r2+1,2)=-1*xsi.*(xsi-1).*eta;
dNr(2:2:r2+1,3)= 1/4*xsi.*(xsi-1).*(2*eta+1);
dNr(2:2:r2+1,4)=-1/2*(xsi+1).*(xsi-1).*(2*eta+1);
dNr(2:2:r2+1,5)= 1/4*xsi.*(xsi+1).*(2*eta+1);
dNr(2:2:r2+1,6)=-1*xsi.*(xsi+1).*eta;
dNr(2:2:r2+1,7)= 1/4*xsi.*(xsi+1).*(2*eta-1);
dNr(2:2:r2+1,8)=-1/2*(xsi+1).*(xsi-1).*(2*eta-1); 
dNr(2:2:r2+1,9)= 2*(xsi+1).*(xsi-1).*eta;