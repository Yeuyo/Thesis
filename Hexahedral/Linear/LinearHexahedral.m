[E,v,coord,etopol,bc,f]=cantilever_endload;
nodes=size(coord,1); nDoF=3*nodes; nels=size(etopol,1);
eDoF=zeros(24,1); DoF=reshape(1:nDoF,3,nodes)'; K=zeros(nDoF);
uvw=zeros(nDoF,1); Gpsig=[]; Iz=zeros(6); Iz(1:3,1:3)=ones(3);
I2=eye(6); I2(4:6,4:6)=0.5*eye(3); wp=ones(8,1); dN=zeros(24,8);
D=E/((1+v)*(1-2*v))*(v*Iz+((1-2*v)*I2));
Gpcoord=zeros(nels,8,3); BeGp=zeros(nels,8,6,24);
xsi=[-1 -1  1  1 -1 -1  1  1]'/sqrt(3);
eta=[-1 -1 -1 -1  1  1  1  1]'/sqrt(3);
zet=[-1  1  1 -1 -1  1  1 -1]'/sqrt(3);
for n=1:8; sxsi=sign(xsi(n)); seta=sign(eta(n)); szet=sign(zet(n));
  for Gp=1:8
    N(Gp,n)=(1+sxsi*xsi(Gp))*(1+seta*eta(Gp))*(1+szet*zet(Gp));
    dN((3*Gp)-2,n)=sxsi*(1+seta*eta(Gp))*(1+szet*zet(Gp));
    dN((3*Gp)-1,n)=seta*(1+sxsi*xsi(Gp))*(1+szet*zet(Gp));
    dN((3*Gp)-0,n)=szet*(1+sxsi*xsi(Gp))*(1+seta*eta(Gp));
  end
end; N=N/8; dN=dN/8;
for nel=1:nels
  for n=1:8
    eDoF((3*n)-2:3*n)=DoF(etopol(nel,n),:);
  end
  JT=dN*coord(etopol(nel,:),:); ke=zeros(24);
  for i=1:8; indx=[(3*i)-2:3*i]';
    Gpcoord(nel,i,1:3)=N(i,:)*coord(etopol(nel,:),1:3);
    detJ=det(JT(indx,:)); dNx=(JT(indx,:))\dN(indx,:);
    B=zeros(6,24); B(1,1:3:end)=dNx(1,:); B(2,2:3:end)=dNx(2,:);
    B(3,3:3:end)=dNx(3,:); B(4,1:3:end)=dNx(2,:);
    B(4,2:3:end)=dNx(1,:); B(5,2:3:end)=dNx(3,:);
    B(5,3:3:end)=dNx(2,:); B(6,1:3:end)=dNx(3,:);
    B(6,3:3:end)=dNx(1,:); BeGp(nel,i,:,:)=B;
    ke=ke+(B'*D*B*detJ*wp(i));
  end
  K(eDoF,eDoF)=K(eDoF,eDoF)+ke;
end
fDoF=[1:nDoF]'; pDoF=bc(:,1); d(pDoF)=bc(:,2); fDoF(pDoF)=[];
uvw(fDoF)=K(fDoF,fDoF)\(f(fDoF)-K(fDoF,pDoF)*d(pDoF));
