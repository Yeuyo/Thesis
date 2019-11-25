clear; clc;
[coord,etpl,fext0,bc,ngp,lstps,NRitmax,NRtol,ndim,E,nu,ka,oVn]=exc;
[nels,nen]=size(etpl); nodes=size(coord,1); nDoF=nodes*ndim;
neDoF=(nen*ndim)^2; krow=zeros(neDoF*nels,1); kcol=krow; kval=krow;
uvw=zeros(nDoF,1); uvwold=uvw; fint=uvw; react=uvw;
fd=(1:nDoF); fd(bc(:,1))=[]; lay=size(coord,2)-3;
epsEn=zeros(6,ngp,lay,nels); epsE=epsEn; sigN=epsEn; sig=epsEn;
Vn=zeros(nen,3,nels); oL=zeros(ngp,9,lay,nels); Told=zeros(nodes,1);
dt=0.05; Kp=zeros(nodes); Mp=Kp; Ptol=1e-8; rn=zeros(ngp,nels); T=Told; dT=Told;
c=8; alp=0.05; gamma=0.002; miu1=0.2; miu2=0.3; beta=0.15;
Rp=Told; dRpdp=Kp;
Iex=15; iex=zeros(nodes,1);
for lstp=0:lstps
  if lstp==1
    iex(1:5)=Iex;
  end
  if lstp==5
    iex=zeros(nodes,1);
  end
  fext=(lstp/lstps)*fext0; oobf=react+fext-fint; oobfnorm=2*NRtol; oobP=oobfnorm; NRit=0; 
  while oobP>Ptol %((NRit<NRitmax)&&(oobfnorm>NRtol))
    NRit=NRit+1; fint=zeros(nDoF,1); dreact=fint; dduvw=fint; fP=fint;
    if lstp>=1
      Kt=sparse(krow,kcol,kval,nDoF,nDoF);
      dduvw(bc(:,1))=(1+sign(1-NRit))*bc(:,2)/lstps;
      dduvw(fd)=Kt(fd,fd)\(oobf(fd)-Kt(fd,bc(:,1))*dduvw(bc(:,1)));
      dreact(bc(:,1))=Kt(bc(:,1),:)*dduvw-oobf(bc(:,1));
      dT=dT-dRpdp\Rp;
%       if NRit==1
%         R=Mp+0.5*dt*Kp; S=Mp-0.5*dt*Kp;
%         Tnew=R\(S*Told);
%         Told=Tnew;
%       end
    end
    uvw=uvw+dduvw; react=react+dreact; duvw=uvw-uvwold; T=T+dT;
    if NRit==1
%       for node=1:nodes
%         r=rn(node); Rr=1e9; Pl=Told(node); Pln=Told(node);
%         while abs(Rr)>=Ptol
%           Rr=r-rn(node)-...
%              ((gamma+((miu1*r)/(miu2+Pl)))*(-r-c*Pl*(Pl-beta-1)))*dt;
%           Kr=1+(gamma+(miu1/(miu2+Pl))*(2*r+c*Pl*(Pl-beta-1)))*dt;
%           r=r-Kr\Rr;
%         end
%         rn(node)=r;
%         dRrdP=((gamma+((miu1*r)/(miu2+Pl)))*c*(2*Pl-beta-1)-...
%                ((miu1*r)/(miu2+Pl)^2)*(r+c*Pl*(Pl-beta-1)))*dt;
%         drdP=-Kr\dRrdP;
%         fV=-c*Pl*(Pl-alp)*(Pl-1)-r*Pl   + iex(node);
%         dfV=c*(-3*Pl^2+2*(1-alp)*Pl+alp)-r-Pl*drdP;
%         Rp=1e9;
%         while abs(Rp)>=Ptol
%           Rp=((Pl-Pln)/dt)-fV;
%           dRpdp=(1/dt)-dfV;
%           Pl=Pl-dRpdp\Rp;
%         end
%         Told(node)=Pl;
%       end
    end
    for nel=1:nels
      ed=reshape(ones(ndim,1)*etpl(nel,:)*ndim-...
                 (ndim-1:-1:0).'*ones(1,nen),1,nen*ndim);
      [ke,felem,fPel,epsE(:,:,:,nel),Vn(:,:,nel),sig(:,:,:,nel),...
       L(:,:,:,nel),kp,mp,rp,drpdp]=Shell_TL_IS(lay,coord(etpl(nel,:),:),...
         uvw(ed),duvw(ed),ngp,epsEn(:,:,:,nel),oVn(:,:,nel),E,nu,ka,...
         sigN(:,:,:,nel),oL(:,:,:,nel),T(etpl(nel,:)),Told(etpl(nel,:)),rn(:,nel),iex(etpl(nel,:)));
      if lstp==0; krow((nel-1)*neDoF+1:nel*neDoF)=...
                    reshape(ed.'*ones(1,nen*ndim),neDoF,1);
                  kcol((nel-1)*neDoF+1:nel*neDoF)=...
                    reshape(ones(nen*ndim,1)*ed  ,neDoF,1);
      end
      kval((nel-1)*neDoF+1:nel*neDoF)=reshape(ke,neDoF,1);
      fint(ed)=fint(ed)+felem; fP(ed)=fP(ed)+fPel;
      Kp(etpl(nel,:),etpl(nel,:))=Kp(etpl(nel,:),etpl(nel,:))+kp;
      Mp(etpl(nel,:),etpl(nel,:))=Mp(etpl(nel,:),etpl(nel,:))+mp;
      Rp(etpl(nel,:))=Rp(etpl(nel,:))+rp;
      dRpdp(etpl(nel,:),etpl(nel,:))=dRpdp(etpl(nel,:),etpl(nel,:))+drpdp;
    end
    oobf=fext+react-fint; oobfnorm=norm(oobf)/norm(fext+fP+react+eps); oobP=norm(Rp);
    fprintf('%4i %4i %6.3e %6.3e\n',lstp,NRit,oobfnorm,oobP);
  end
  uvwold=uvw; epsEn=epsE; sigN=sig; oL=oL+L; Told=T;
  if lstp>=1
    thermometer(lstp,:)=Told;
    dishis(lstp,:)=uvw;
  end
end
% thermometer(thermometer<0)=0; thermometer(thermometer>1)=1; kolor='r-'; tplot=100;
% for nel=1:nels
%   ed=reshape(ones(ndim,1)*etpl(nel,:)*ndim-(ndim-1:-1:0).'*ones(1,nen),1,nen*ndim);
%   subplot(2,2,3)
%   plotShell(lay,coord(etpl(nel,:),:),zeros(45,1),zeros(45,1),zeros(1,9)','k-');
%   plotShell(lay,coord(etpl(nel,:),:),dishis(tplot,ed)',zeros(45,1),thermometer(tplot,etpl(nel,:))',kolor);
%   view(0,0)
%   subplot(2,2,4)
%   plotShell(lay,coord(etpl(nel,:),:),zeros(45,1),zeros(45,1),zeros(1,9)','k-');
%   plotShell(lay,coord(etpl(nel,:),:),dishis(tplot,ed)',zeros(45,1),thermometer(tplot,etpl(nel,:))',kolor);
%   view(0,90)
% end

plot([1:lstps], thermometer(:,1), 'k-'); hold on;
plot([1:lstps], thermometer(:,6), 'b-');
plot([1:lstps], thermometer(:,12), 'r-');
plot([1:lstps], thermometer(:,18), 'm-');
plot([1:lstps], thermometer(:,24), 'c-');
plot([1:lstps], thermometer(:,30), 'y-');
plot([1:lstps], thermometer(:,36), 'g-');

% h=gcf;
% set(h,'PaperOrientation','landscape');
% set(h,'PaperUnits','normalized');
% set(h,'PaperPosition', [0 0 1 1]);
% print(gcf, '-dpdf', 'print0.pdf');