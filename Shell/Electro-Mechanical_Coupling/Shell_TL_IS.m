function [ke,fint,fP,epsE,Vn,sig,L,kp,mp,Rp,dRpdp]=Shell_TL_IS(layern,nodeData,uvw,duvw,ngp,epsEn,oVn,E,nu,k,sigN,oL,Tn,Told,rn,iex)
nen=size(nodeData,1); coord=nodeData(:,1:3); t=nodeData(:,4:4+layern-1);
T=sum(t,2); epsE=zeros(6,ngp); [wp,GpLoc]=GpPos(ngp);
ke=zeros(nen*5); fint=zeros(nen*5,1); dxr=zeros(3);
fP=fint; kp=zeros(nen); mp=kp; Rp=zeros(nen,1); dRpdp=zeros(nen);
ex=[1 0 0].'; ey=[0 1 0].'; ez=[0 0 1].';
V1=zeros(nen,3); V2=V1; g1=V1; g2=V1;
xsi=[-1; -1; -1; 0; 1; 1; 1; 0; 0]; eta=[-1; 0; 1; 1; 1; 0; -1; -1; 0];
dNr=dershapefunc2D(xsi,eta);
dnr=dNr(1:2:end,:); dns=dNr(2:2:end,:);
disp=reshape(uvw-duvw,5,[])'; elcoord=coord+disp(:,1:3);
for n=1:9
  Vn(n,:)=cross((dnr(n,:)*elcoord)/norm(dnr(n,:)*elcoord),...
                (dns(n,:)*elcoord)/norm(dns(n,:)*elcoord));
end
for n=1:layern
  for node=1:nen
    if norm(cross(ey,Vn(node,:)'))<1e-6          
      V=cross(ez,Vn(node,:)')';         
    else         
      V=cross(ey,Vn(node,:)')';     
    end
    V1(node,:)=V/norm(V);
    V2(node,:)=cross(Vn(node,:).',V1(node,:).').';
    g1(node,:)=-0.5*T(n)*V2(node,:);
    g2(node,:)= 0.5*T(n)*V1(node,:);
  end
  D=(E(n)/(1-nu^2))*[[1 nu; nu 1; 0 0] zeros(3,4);...
                     zeros(3,3) eye(3)*k*(1-nu)*0.5];
  D(4,4)=D(4,4)/k;
  Dp=zeros(6); Dp(1)=1; %Dp([1 8])=3.225e-4;
  for gp=1:ngp
    xsi=GpLoc(gp,1); eta=GpLoc(gp,2); zeta=GpLoc(gp,3);
    N  =shapefunc(xsi,eta);
    tgp=N*t(:,n); Tgp=N*T;
    zet=-1+(2*sum(N*t(:,1:n))-N*(t(:,n).*(1-zeta)))/(N*T);
    dNr=dershapefunc2D(xsi,eta);
    dxr(1:2,:)=(dNr*(coord+0.5*zet*(oVn.*(T*ones(1,3)))));
    dxr(3,:)=(0.5*(N.*T')*oVn);
    dNx=dxr\[dNr; zeros(1,nen)];
    detJ=det(dxr);
    VnGp=(N*oVn).';
    sdir=dxr(2,:)'/norm(dxr(2,:)');
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
    Dpsh=Q'*Dp*Q; Dpsh=Dpsh(1:3,1:3);
    Tl=N*Told;
    if Tl<0; Tl=0; elseif Tl>1; Tl=1; end
    epsP=zeros(6); epsP(1)=0; %0.15*Tl;
    epsPsh=Q'*epsP*Q;
    [B, GNL, epsEt,L(gp,:,n)]=...
      formBmatrix(xsi, eta, zet,coord,oVn,T,g1,g2,duvw,oL(gp,:,n));
    epsE(:,gp,n)=epsEt+epsEn(:,gp,n);
    sig(:,gp,n)=(D*epsEt)+sigN(:,gp,n);
    sigP=D*(-epsPsh(1:7:end)');
    fP=fP+B'*sigP*detJ*wp(gp)*(tgp/Tgp);
    sigT=sig+D*(epsPsh(1:7:end)');
    H=[sigT(1,gp,n)*eye(3) sigT(4,gp,n)*eye(3) sigT(6,gp,n)*eye(3);
       sigT(4,gp,n)*eye(3) sigT(2,gp,n)*eye(3) sigT(5,gp,n)*eye(3);
       sigT(6,gp,n)*eye(3) sigT(5,gp,n)*eye(3) sigT(3,gp,n)*eye(3)];
    ke=ke+((B.'*D*B)+(GNL.'*H*GNL))*detJ*wp(gp)*(tgp/Tgp); 
    fint=fint+B.'*sigT(:,gp,n)*detJ*wp(gp)*(tgp/Tgp);
    
    c=8; alp=0.05; gamma=0.002; miu1=0.2; miu2=0.3; beta=0.15; dt=0.05;
    r=rn(gp); Rr=1e9; Pl=N*Told; Ptol=1e-8;
    while abs(Rr)>=Ptol
      Rr=r-rn(gp)-...
         ((gamma+((miu1*r)/(miu2+Pl)))*(-r-c*Pl*(Pl-beta-1)))*dt;
      Kr=1+(gamma+(miu1/(miu2+Pl))*(2*r+c*Pl*(Pl-beta-1)))*dt;
      r=r-Kr\Rr;
    end
    rn(gp)=r;
    dRrdP=((gamma+((miu1*r)/(miu2+Pl)))*c*(2*Pl-beta-1)-...
           ((miu1*r)/(miu2+Pl)^2)*(r+c*Pl*(Pl-beta-1)))*dt;
    drdP=-Kr\dRrdP;
    fV=-c*Pl*(Pl-alp)*(Pl-1)-r*Pl   + N*iex;
    dfV=c*(-3*Pl^2+2*(1-alp)*Pl+alp)-r-Pl*drdP;
    Rp=Rp+(N'*((N*Tn-N*Told)/dt)-N'*fV)*detJ*wp(gp)*(tgp/Tgp);
    dRpdp=dRpdp+(N'*(1/dt)*N-N'*dfV*N)*detJ*wp(gp)*(tgp/Tgp);
    
    kp=kp+dNx'*Dpsh*dNx*wp(gp)*det(dxr)*(tgp/Tgp);
    mp=mp+N'*N*wp(gp)*det(dxr)*(tgp/Tgp);
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

function [B,GNL,epsE,L]=formBmatrix(xsi,eta,zet,coord,oVn,t,g1,g2,duvw,oL)
nen=size(coord,1);
N  =shapefunc(xsi,eta);
dNr=dershapefunc2D(xsi,eta);
dxr=zeros(3);
dxr(1:2,:)=(dNr*(coord+0.5*zet*(oVn.*(t*ones(1,3)))));
dxr(3,:)=(0.5*(N.*t')*oVn);
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