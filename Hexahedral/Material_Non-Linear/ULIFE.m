function[ke,fint,epsE,sig,F]=ULIFE(nodeData,K,miu,uvw,duvw,ngp,Fold)
nen=size(nodeData,1); coord=nodeData(:,1:3); epsE=zeros(6,ngp);
[wp,GpLoc]=GpPos(ngp); ke=zeros(nen*3); Ks=ke; Kp=ke;
fint=zeros(nen*3,1); F=zeros(3,3,ngp); sig=zeros(6,ngp);
I39=[ones(3,3) zeros(3,6); zeros(6,9)];
I9=eye(9);
I93=zeros(9); I93([1, 11, 21, 32, 40, 52, 60, 72, 80])=1;
ve=0; Ve=0; dNxm=zeros(3,nen);
for gp=1:ngp
  xsi=GpLoc(gp,1); eta=GpLoc(gp,2); zet=GpLoc(gp,3);
  dNr=dershapefunc2D(xsi,eta,zet);
  dXr=dNr*coord; dxr=dNr*(coord+uvw(reshape(1:3*nen,3,nen))'); 
  dNx=dxr\dNr; detJ0=det(dXr); detJ=det(dxr);
  JW=detJ*wp(gp); ve=ve+JW;
  dNxm=dNxm+dNx*JW;
  Ve=Ve+detJ0*wp(gp);
end
dNxm=dNxm/ve; Jbar=ve/Ve;
p=K*(Jbar-1); Kbar=K*Jbar;
for gp=1:ngp
  xsi=GpLoc(gp,1); eta=GpLoc(gp,2); zet=GpLoc(gp,3);
  dNr=dershapefunc2D(xsi,eta,zet);
  dxr=dNr*(coord+uvw(reshape(1:3*nen,3,nen))');
  detJ=det(dxr); dNx=dxr\dNr;
  dF=inv(eye(3)-duvw(reshape(1:3*nen,3,nen))*dNx.');
  F(:,:,gp)=dF*Fold(:,:,gp); 
  J=det(F(:,:,gp)); b=F(:,:,gp)*F(:,:,gp)';
  Ib=trace(b); b9=b([1,5,9,4,2,8,6,3,7])';
  I9b1=[b9 b9 b9 zeros(9,6)]; I9b2=[b9'; b9'; b9'; zeros(6,9)];
  [B,G]=formBG(dNx,nen);
  s=(miu*J^(-5/3))*(b - 1/3*Ib*eye(3)); s=s+p*eye(3);
  for n=1:nen
    for m=1:nen
      Ks(n*3-2:n*3,m*3-2:m*3)=dNx(:,n)'*s*dNx(:,m)*eye(3);
    end
  end
  sig(:,gp)=s([1,5,9,2,6,3])';
  D=2*(miu*J^(-5/3))*...
      ((1/6)*Ib*(I9+I93)-(1/3)*I9b1-(1/3)*I9b2+(1/9)*Ib*I39);
  D=D+p*(I39-I9-I93);
  ke=ke+(G'*D*G+Ks)*detJ*wp(gp);
  fint=fint+B.'*sig(:,gp)*detJ*wp(gp);
end
for n=1:nen
  for m=1:nen
    Kp(n*3-2:n*3,m*3-2:m*3)=Kbar*ve*(dNxm(:,n)*dNxm(:,m)');
  end
end
ke=ke+Kp;

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
end

function [dNr]=dershapefunc2D(xsi,eta,zet)
ngp=size(xsi,1); r2=ngp*3;
dNr(1:2:r2,1)=-0.125*(1-eta).*(1-zet);
dNr(1:2:r2,2)=-0.125*(1-eta).*(1+zet);
dNr(1:2:r2,3)= 0.125*(1-eta).*(1+zet);
dNr(1:2:r2,4)= 0.125*(1-eta).*(1-zet);
dNr(1:2:r2,5)=-0.125*(1+eta).*(1-zet);
dNr(1:2:r2,6)=-0.125*(1+eta).*(1+zet);
dNr(1:2:r2,7)= 0.125*(1+eta).*(1+zet);
dNr(1:2:r2,8)= 0.125*(1+eta).*(1-zet);
dNr(2:2:r2,1)=-0.125*(1-xsi).*(1-zet);
dNr(2:2:r2,2)=-0.125*(1-xsi).*(1+zet);
dNr(2:2:r2,3)=-0.125*(1+xsi).*(1+zet);
dNr(2:2:r2,4)=-0.125*(1+xsi).*(1-zet);
dNr(2:2:r2,5)= 0.125*(1-xsi).*(1-zet);
dNr(2:2:r2,6)= 0.125*(1-xsi).*(1+zet);
dNr(2:2:r2,7)= 0.125*(1+xsi).*(1+zet);
dNr(2:2:r2,8)= 0.125*(1+xsi).*(1-zet);
dNr(3:2:r2,1)=-0.125*(1-xsi).*(1-eta);
dNr(3:2:r2,2)= 0.125*(1-xsi).*(1-eta);
dNr(3:2:r2,3)= 0.125*(1+xsi).*(1-eta);
dNr(3:2:r2,4)=-0.125*(1+xsi).*(1-eta);
dNr(3:2:r2,5)=-0.125*(1-xsi).*(1+eta);
dNr(3:2:r2,6)= 0.125*(1-xsi).*(1+eta);
dNr(3:2:r2,7)= 0.125*(1+xsi).*(1+eta);
dNr(3:2:r2,8)=-0.125*(1+xsi).*(1+eta);

function [B,G]=formBG(dNx,nen)
B=zeros(6,nen*3);
B(1,1:3:end)=dNx(1,:); B(2,2:3:end)=dNx(2,:);
B(3,3:3:end)=dNx(3,:); B(4,1:3:end)=dNx(2,:);
B(4,2:3:end)=dNx(1,:); B(5,2:3:end)=dNx(3,:);
B(5,3:3:end)=dNx(2,:); B(6,1:3:end)=dNx(3,:);
B(6,3:3:end)=dNx(1,:);
G=zeros(9,nen*3); G(1:3,:)=B(1:3,:);
G(4,1:3:end)=dNx(2,:); G(5,2:3:end)=dNx(1,:);
G(6,2:3:end)=dNx(3,:); G(7,3:3:end)=dNx(2,:);
G(8,3:3:end)=dNx(1,:); G(9,1:3:end)=dNx(3,:);
