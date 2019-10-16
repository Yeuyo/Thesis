[coord,etpl,fext0,bc,lam,miu,ngp,lstps,NRitmax,NRtol]=CookMembrane;
[nels,nen]=size(etpl); nodes=size(coord,1);
ndim=3; nDoF=nodes*ndim; neDoF=(nen*ndim)^2;
krow=zeros(neDoF*nels,1); kcol=krow; kval=krow;
uvw=zeros(nDoF,1); uvwold=uvw; fint=uvw; react=uvw;
fd=(1:nDoF); fd(bc(:,1))=[];
epsE=zeros(6,ngp,nels); sig=epsE;
Fold=zeros(3,3,ngp,nels); F=Fold;
for n=1:nels; for m=1:ngp; Fold(:,:,m,n)=eye(3); end; end
for lstp=0:lstps
  fext=(lstp/lstps)*fext0; oobf=react+fext-fint; oobfnorm=2*NRtol; NRit=0; 
  while ((NRit<NRitmax)&&(oobfnorm>NRtol))
    NRit=NRit+1; fint=zeros(nDoF,1); dreact=fint; dduvw=fint;
    if lstp>=1
      Kt=sparse(krow,kcol,kval,nDoF,nDoF);
      dduvw(bc(:,1))=(1+sign(1-NRit))*bc(:,2)/lstps;
      dduvw(fd)=Kt(fd,fd)\(oobf(fd)-Kt(fd,bc(:,1))*dduvw(bc(:,1)));
      dreact(bc(:,1))=Kt(bc(:,1),:)*dduvw-oobf(bc(:,1));
    end
    uvw=uvw+dduvw; react=react+dreact; duvw=uvw-uvwold;
    for nel=1:nels
      ed=reshape(ones(ndim,1)*etpl(nel,:)*ndim-(ndim-1:-1:0).'*ones(1,nen),1,nen*ndim);
      [ke,felem,epsE(:,:,nel),sig(:,:,nel),F(:,:,:,nel)]=ULFE(coord(etpl(nel,:),:),lam,miu,uvw(ed),duvw(ed),ngp,Fold(:,:,:,nel));
      if lstp==0; krow((nel-1)*neDoF+1:nel*neDoF)=reshape(ed.'*ones(1,nen*ndim),neDoF,1);
                  kcol((nel-1)*neDoF+1:nel*neDoF)=reshape(ones(nen*ndim,1)*ed  ,neDoF,1);
      end
      kval((nel-1)*neDoF+1:nel*neDoF)=reshape(ke,neDoF,1);
      fint(ed)=fint(ed)+felem;
    end
    oobf=fext+react-fint; oobfnorm=norm(oobf)/norm(fext+react+eps);
    fprintf('%4i %4i %6.3e\n',lstp,NRit,oobfnorm);
  end
  uvwold=uvw; Fold=F;
end
