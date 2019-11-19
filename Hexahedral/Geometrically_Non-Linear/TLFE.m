function[ke,fint,epsE,sig,duX]=TLFE(coord,D,duvw,ngp,epsEn,sigN,oduX)
nen=size(coord,1); epsE=zeros(6,ngp); sig=epsE;
[wp,GpLoc]=GpPos(ngp); ke=zeros(nen*3); Ks=ke;
fint=zeros(nen*3,1); duX=zeros(3,3,ngp);
for gp=1:ngp
  xsi=GpLoc(gp,1); eta=GpLoc(gp,2); zet=GpLoc(gp,3);
  dNr=dershapefunc2D(xsi,eta,zet);
  dXr=dNr*coord; detJ=det(dXr); dNX=dXr\dNr;
  duX(:,:,gp)=dNX*reshape(duvw,3,nen)';
  epsEt=[0.5*(duX(1,1,gp) + duX(1,1,gp) + ...
         oduX(1,1,gp)*duX(1,1,gp) + oduX(1,2,gp)*duX(1,2,gp) + ...
         oduX(1,3,gp)*duX(1,3,gp) + duX(1,1,gp)*oduX(1,1,gp) + ...
         duX(1,2,gp)*oduX(1,2,gp) + duX(1,3,gp)*oduX(1,3,gp) + ...
         duX(1,1,gp)*duX(1,1,gp) + duX(1,2,gp)*duX(1,2,gp) + ...
         duX(1,3,gp)*duX(1,3,gp));
         0.5*(duX(2,2,gp) + duX(2,2,gp) + ...
         oduX(2,1,gp)*duX(2,1,gp) + oduX(2,2,gp)*duX(2,2,gp) + ...
         oduX(2,3,gp)*duX(2,3,gp) + duX(2,1,gp)*oduX(2,1,gp) + ...
         duX(2,2,gp)*oduX(2,2,gp) + duX(2,3,gp)*oduX(2,3,gp) + ...
         duX(2,1,gp)*duX(2,1,gp) + duX(2,2,gp)*duX(2,2,gp) + ...
         duX(2,3,gp)*duX(2,3,gp));
         0.5*(duX(3,3,gp) + duX(3,3,gp) + ...
         oduX(3,1,gp)*duX(3,1,gp) + oduX(3,2,gp)*duX(3,2,gp) + ...
         oduX(3,3,gp)*duX(3,3,gp) + duX(3,1,gp)*oduX(3,1,gp) + ...
         duX(3,2,gp)*oduX(3,2,gp) + duX(3,3,gp)*oduX(3,3,gp) + ...
         duX(3,1,gp)*duX(3,1,gp) + duX(3,2,gp)*duX(3,2,gp) + ...
         duX(3,3,gp)*duX(3,3,gp));
         1.0*(duX(2,1,gp) + duX(1,2,gp) + ...
         oduX(1,1,gp)*duX(2,1,gp) + oduX(1,2,gp)*duX(2,2,gp) + ...
         oduX(1,3,gp)*duX(2,3,gp) + duX(1,1,gp)*oduX(2,1,gp) + ...
         duX(1,2,gp)*oduX(2,2,gp) + duX(1,3,gp)*oduX(2,3,gp) + ...
         duX(1,1,gp)*duX(2,1,gp) + duX(1,2,gp)*duX(2,2,gp) + ...
         duX(1,3,gp)*duX(2,3,gp));
         1.0*(duX(2,3,gp) + duX(3,2,gp) + ...
         oduX(3,1,gp)*duX(2,1,gp) + oduX(3,2,gp)*duX(2,2,gp) + ...
         oduX(3,3,gp)*duX(2,3,gp) + duX(3,1,gp)*oduX(2,1,gp) + ...
         duX(3,2,gp)*oduX(2,2,gp) + duX(3,3,gp)*oduX(2,3,gp) + ...
         duX(3,1,gp)*duX(2,1,gp) + duX(3,2,gp)*duX(2,2,gp) + ...
         duX(3,3,gp)*duX(2,3,gp));
         1.0*(duX(3,1,gp) + duX(1,3,gp) + ...
         oduX(1,1,gp)*duX(3,1,gp) + oduX(1,2,gp)*duX(3,2,gp) + ...
         oduX(1,3,gp)*duX(3,3,gp) + duX(1,1,gp)*oduX(3,1,gp) + ...
         duX(1,2,gp)*oduX(3,2,gp) + duX(1,3,gp)*oduX(3,3,gp) + ...
         duX(1,1,gp)*duX(3,1,gp) + duX(1,2,gp)*duX(3,2,gp) + ...
         duX(1,3,gp)*duX(3,3,gp))];
  epsE(:,gp)=epsEt+epsEn(:,gp);
  sig(:,gp)=(D(:,:,gp)*epsEt)+sigN(:,gp);
  s=[sig(1,gp) sig(4,gp) sig(6,gp);
     sig(4,gp) sig(2,gp) sig(5,gp);
     sig(6,gp) sig(5,gp) sig(3,gp)];
  BL=formB(dNX,nen); BNL=zeros(6,nen*3);
  tduX=oduX(:,:,gp)+duX(:,:,gp);
  A=[ tduX(1,:) zeros(1,3) zeros(1,3);
     zeros(1,3)  tduX(2,:) zeros(1,3);
     zeros(1,3) zeros(1,3)  tduX(3,:);
      tduX(2,:)  tduX(1,:) zeros(1,3);
     zeros(1,3)  tduX(3,:)  tduX(2,:);
      tduX(3,:) zeros(1,3)  tduX(1,:)];
  for n=1:nen
    for m=1:nen
      Ks(n*3-2:n*3,m*3-2:m*3)=dNX(:,n)'*s*dNX(:,m)*eye(3);
    end
    BNL(:,(n-1)*3+1:(n-1)*3+3)=...
      A*[dNX(1,n)*eye(3); dNX(2,n)*eye(3); dNX(3,n)*eye(3)];
  end
  Bt=BL+BNL;
  ke=ke+(Bt'*D(:,:,gp)*Bt+Ks)*detJ*wp(gp);
  fint=fint+Bt'*sig(:,gp)*detJ*wp(gp);
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
  xsi=[-1; 0; 1;-1; 0; 1;-1; 0; 1; -1; 0; 1;-1; 0;
           1;-1; 0; 1; -1; 0; 1;-1; 0; 1;-1; 0; 1]*g3;
  eta=[-1;-1;-1;-1;-1;-1;-1;-1;-1;  0; 0; 0; 0; 0;
           0; 0; 0; 0;  1; 1; 1; 1; 1; 1; 1; 1; 1]*g3;
  zet=[-1;-1;-1; 0; 0; 0; 1; 1; 1; -1;-1;-1; 0; 0;
           0; 1; 1; 1; -1;-1;-1; 0; 0; 0; 1; 1; 1]*g3;
  w(:,1)=[w2;w1;w2;w2;w1;w2;w2;w1;w2; w2;w1;w2;w2;w1;w2;w2;w1;w2;
                                      w2;w1;w2;w2;w1;w2;w2;w1;w2];
  w(:,2)=[w2;w2;w2;w2;w2;w2;w2;w2;w2; w1;w1;w1;w1;w1;w1;w1;w1;w1;
                                      w2;w2;w2;w2;w2;w2;w2;w2;w2];
  w(:,3)=[w2;w2;w2;w1;w1;w1;w2;w2;w2; w2;w2;w2;w1;w1;w1;w2;w2;w2;
                                      w2;w2;w2;w1;w1;w1;w2;w2;w2];
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