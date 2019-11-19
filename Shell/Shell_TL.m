function [ke,fint,epsE,Vn,sig,L]=...
         Shell_TL(nodeData,uvw,duvw,ngp,epsEn,D,dxr,sigN,oL)
nen=size(nodeData,1); lay=size(nodeData,2)-3;
coord=nodeData(:,1:3); t=nodeData(:,4:4+lay-1); T=sum(t,2);
epsE=zeros(6,ngp,lay); sig=epsE; L=zeros(ngp,9,lay);
[wp,GpLoc]=GpPos(ngp); Vn=zeros(nen,3);
ke=zeros(nen*5); fint=zeros(nen*5,1);
ey=[0 1 0].'; ez=[0 0 1].';
V1=zeros(nen,3); V2=V1; g1=V1; g2=V1;
xsi=[-1; -1; -1; 0; 1; 1; 1; 0; 0]; eta=[-1; 0; 1; 1; 1; 0; -1; -1; 0];
dNr=dershapefunc2D(xsi,eta);
dnr=dNr(1:2:end,:); dns=dNr(2:2:end,:);
disp=reshape(uvw-duvw,5,[])'; elcoord=coord+disp(:,1:3);
for n=1:nen
  Vn(n,:)=cross((dnr(n,:)*elcoord)/norm(dnr(n,:)*elcoord),...
                (dns(n,:)*elcoord)/norm(dns(n,:)*elcoord));
end
for n=1:lay
  for node=1:nen
    if abs(ey)-abs(Vn(node,:)')<1e-6
      V=cross(ez,Vn(node,:).').';
    else
      V=cross(ey,Vn(node,:).').';
    end
    V1(node,:)=V/norm(V);
    V2(node,:)=cross(Vn(node,:).',V1(node,:).').';
    g1(node,:)=-0.5*T(n)*V2(node,:);
    g2(node,:)= 0.5*T(n)*V1(node,:);
  end
  for gp=1:ngp
    xsi=GpLoc(gp,1); eta=GpLoc(gp,2); zeta=GpLoc(gp,3);
    N  =shapefunc(xsi,eta);
    tgp=N*t(:,n); Tgp=N*T;
    zet=-1+(2*sum(N*t(:,1:n))-N*(t(:,n).*(1-zeta)))/(N*T);
    dNr=dershapefunc2D(xsi,eta);
    detJ=det(dxr(:,:,gp,n));
    [B, GNL, epsEt,L(gp,:,n)]=...
      formBmatrix(N,dNr,zet,coord,dxr(:,:,gp,n),g1,g2,duvw,oL(gp,:,n));
    epsE(:,gp,n)=epsEt+epsEn(:,gp,n);
    sig(:,gp,n)=(D(:,:,gp,n)*epsEt)+sigN(:,gp,n);
    H=[sig(1,gp,n)*eye(3) sig(4,gp,n)*eye(3) sig(6,gp,n)*eye(3);
       sig(4,gp,n)*eye(3) sig(2,gp,n)*eye(3) sig(5,gp,n)*eye(3);
       sig(6,gp,n)*eye(3) sig(5,gp,n)*eye(3) sig(3,gp,n)*eye(3)];
    ke=ke+((B.'*D(:,:,gp,n)*B)+(GNL.'*H*GNL))*detJ*wp(gp)*(tgp/Tgp);
    fint=fint+B.'*sig(:,gp,n)*detJ*wp(gp)*(tgp/Tgp);
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

function [B,GNL,epsE,L]=formBmatrix(N,dNr,zet,coord,dxr,g1,g2,duvw,oL)
nen=size(coord,1);
invJ=inv(dxr);
dNx=dxr\[dNr; zeros(1,nen)];
G=zet*(invJ(:,1:2)*dNr)+invJ(:,3)*N;
B9=zeros(9,nen*5); BNL=zeros(6,nen*5); GNL=zeros(9,nen*5); L=zeros(1,9);
for n=1:3
  for m=1:nen
    L(:,(n-1)*3+1:(n-1)*3+3)=L(:,(n-1)*3+1:(n-1)*3+3)+...
     ([dNx(:,m) g1(m,n).*ones(3,1).*G(:,m) g2(m,n).*ones(3,1).*G(:,m)]*...
      [duvw(((m-1)*5)+n); duvw(((m-1)*5)+4); duvw(((m-1)*5)+5)])';
  end
end
epsE=[0.5*(L(1,1)+L(1,1)+oL(1,1)*L(1,1)+oL(1,4)*L(1,4)+oL(1,7)*L(1,7)+...
      L(1,1)*oL(1,1)+L(1,4)*oL(1,4)+L(1,7)*oL(1,7)+...
      L(1,1)*L(1,1)+L(1,4)*L(1,4)+L(1,7)*L(1,7));
      0.5*(L(1,5)+L(1,5)+oL(1,2)*L(1,2)+oL(1,5)*L(1,5)+oL(1,8)*L(1,8)+...
      L(1,2)*oL(1,2)+L(1,5)*oL(1,5)+L(1,8)*oL(1,8)+...
      L(1,2)*L(1,2)+L(1,5)*L(1,5)+L(1,8)*L(1,8));
      0.5*(L(1,9)+L(1,9)+oL(1,3)*L(1,3)+oL(1,6)*L(1,6)+oL(1,9)*L(1,9)+...
      L(1,3)*oL(1,3)+L(1,6)*oL(1,6)+L(1,9)*oL(1,9)+...
      L(1,3)*L(1,3)+L(1,6)*L(1,6)+L(1,9)*L(1,9));
      1.0*(L(1,2)+L(1,4)+oL(1,1)*L(1,2)+oL(1,4)*L(1,5)+oL(1,7)*L(1,8)+...
      L(1,1)*oL(1,2)+L(1,4)*oL(1,5)+L(1,7)*oL(1,8)+...
      L(1,1)*L(1,2)+L(1,4)*L(1,5)+L(1,7)*L(1,8));
      1.0*(L(1,8)+L(1,6)+oL(1,3)*L(1,2)+oL(1,6)*L(1,5)+oL(1,9)*L(1,8)+...
      L(1,3)*oL(1,2)+L(1,6)*oL(1,5)+L(1,9)*oL(1,8)+...
      L(1,3)*L(1,2)+L(1,6)*L(1,5)+L(1,9)*L(1,8));
      1.0*(L(1,3)+L(1,7)+oL(1,1)*L(1,3)+oL(1,4)*L(1,6)+oL(1,7)*L(1,9)+...
      L(1,1)*oL(1,3)+L(1,4)*oL(1,6)+L(1,7)*oL(1,9)+...
      L(1,1)*L(1,3)+L(1,4)*L(1,6)+L(1,7)*L(1,9))];
B9([1 4 9],1:5:end)=dNx;
B9([5 2 6],2:5:end)=dNx;
B9([8 7 3],3:5:end)=dNx;
B9(:,      4:5:end)=g1(:,[1 2 3 1 2 2 3 3 1])'.*G([1 2 3 2 1 3 2 1 3],:);
B9(:,      5:5:end)=g2(:,[1 2 3 1 2 2 3 3 1])'.*G([1 2 3 2 1 3 2 1 3],:);
BL0=B9([1:3 5 7 9],:); BL0(4:6,:)=BL0(4:6,:)+B9([4 6 8],:); tL=oL+L;
A=[tL(:,1:3:9)  zeros(1,3)  zeros(1,3);
    zeros(1,3) tL(:,2:3:9)  zeros(1,3);
    zeros(1,3)  zeros(1,3) tL(:,3:3:9);
   tL(:,2:3:9) tL(:,1:3:9)  zeros(1,3);
    zeros(1,3) tL(:,3:3:9) tL(:,2:3:9);
   tL(:,3:3:9)  zeros(1,3) tL(:,1:3:9)];
for n=1:nen
  GNL(:,(n-1)*5+1:(n-1)*5+5)=...
     [dNx(1,n)*eye(3) g1(n,:)'.*G(1,n).*ones(3,1) g2(n,:)'.*G(1,n).*ones(3,1);
      dNx(2,n)*eye(3) g1(n,:)'.*G(2,n).*ones(3,1) g2(n,:)'.*G(2,n).*ones(3,1);
      dNx(3,n)*eye(3) g1(n,:)'.*G(3,n).*ones(3,1) g2(n,:)'.*G(3,n).*ones(3,1)];
  BNL(:,(n-1)*5+1:(n-1)*5+5)=A*GNL(:,(n-1)*5+1:(n-1)*5+5);
end
B=BL0+BNL;