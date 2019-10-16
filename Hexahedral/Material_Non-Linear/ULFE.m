function[ke,fint,epsE,sig,F]=ULFE(nodeData,lam,miu,uvw,duvw,ngp,Fold)
nen=size(nodeData,1); coord=nodeData(:,1:3); epsE=zeros(6,ngp);
[wp,GpLoc]=GpPos(ngp); ke=zeros(nen*3); Ks=ke;
fint=zeros(nen*3,1); F=zeros(3,3,ngp); sig=zeros(6,ngp);
for gp=1:ngp
  xsi=GpLoc(gp,1); eta=GpLoc(gp,2); zet=GpLoc(gp,3);
  dNr=dershapefunc2D(xsi,eta,zet);
  dxr=dNr*(coord+uvw(reshape(1:3*nen,3,nen))');
  detJ=det(dxr); dNx=dxr\dNr;
  dF=inv(eye(3)-duvw(reshape(1:3*nen,3,nen))*dNx.');
  F(:,:,gp)=dF*Fold(:,:,gp);
  J=det(F(:,:,gp));
  lamb=lam/J; miub=(miu-lam*log(J))/J;
  B=formB(dNx,nen);
  D=[lamb*ones(3)+2*miub*eye(3) zeros(3); zeros(3) miub*eye(3)];
  s=(miu/J)*(F(:,:,gp)*F(:,:,gp)'-eye(3)) + (lam/J)*log(J)*eye(3);
  for n=1:nen
    for m=1:nen
      Ks(n*3-2:n*3,m*3-2:m*3)=dNx(:,n)'*s*dNx(:,m)*eye(3);
    end
  end
  sig(:,gp)=s([1,5,9,2,6,3])';
  ke=ke+(B'*D*B+Ks)*detJ*wp(gp);
  fint=fint+B.'*sig(:,gp)*detJ*wp(gp);
end

function [wp,GpLoc]=GpPos(ngp)
if ngp==8 %2x2x2
  wp=ones(8,1); g2=1/sqrt(3); 
  xsi=[-1 -1 1 1 -1 -1 1 1].'*g2;
  eta=[-1 -1 -1 -1 1 1 1 1].'*g2;
  zet=[-1 1 1 -1 -1 1 1 -1].'*g2;
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
elseif ngp==27 %3x3x3
  g3=sqrt(3/5); w1=8/9; w2=5/9;
  xsi=[-1; 0; 1;-1; 0; 1;-1; 0; 1; -1; 0; 1;-1; 0; 1;-1; 0; 1; -1; 0; 1;-1; 0; 1;-1; 0; 1]*g3;
  eta=[-1;-1;-1;-1;-1;-1;-1;-1;-1;  0; 0; 0; 0; 0; 0; 0; 0; 0;  1; 1; 1; 1; 1; 1; 1; 1; 1]*g3;
  zet=[-1;-1;-1; 0; 0; 0; 1; 1; 1; -1;-1;-1; 0; 0; 0; 1; 1; 1; -1;-1;-1; 0; 0; 0; 1; 1; 1]*g3;
  w(:,1)=[w2;w1;w2;w2;w1;w2;w2;w1;w2; w2;w1;w2;w2;w1;w2;w2;w1;w2; w2;w1;w2;w2;w1;w2;w2;w1;w2];
  w(:,2)=[w2;w2;w2;w2;w2;w2;w2;w2;w2; w1;w1;w1;w1;w1;w1;w1;w1;w1; w2;w2;w2;w2;w2;w2;w2;w2;w2];
  w(:,3)=[w2;w2;w2;w1;w1;w1;w2;w2;w2; w2;w2;w2;w1;w1;w1;w2;w2;w2; w2;w2;w2;w1;w1;w1;w2;w2;w2];
  wp=w(:,1).*w(:,2).*w(:,3);
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

function B=formB(dNx,nen)
B=zeros(6,nen*3);
B(1,1:3:end)=dNx(1,:); B(2,2:3:end)=dNx(2,:); B(3,3:3:end)=dNx(3,:);
B(4,1:3:end)=dNx(2,:); B(4,2:3:end)=dNx(1,:); B(5,2:3:end)=dNx(3,:);
B(5,3:3:end)=dNx(2,:); B(6,1:3:end)=dNx(3,:); B(6,3:3:end)=dNx(1,:);
