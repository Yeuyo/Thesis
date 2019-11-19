function [coord,etopol,fext,bc,D,ngp,lstps,NRitmax,NRtol]=...
                                      GNLcantilever_endload
E=1.2e6; v=0;
ngp=8; lstps=20; NRitmax=50; NRtol=1.0e-6;
nels=300; len=10; wid=1; thk=0.1;
Iz=zeros(6); Iz(1:3,1:3)=ones(3); I2=eye(6); I2(4:6,4:6)=0.5*eye(3);
D=zeros(6,6,ngp,nels);
etopol=zeros(nels,8);
for nel=1:nels+1
  coord((nel-1)*4+1:(nel)*4,:)=[(nel-1)*(len/nels) 0 0;
                                (nel-1)*(len/nels) 0 thk;
                                (nel-1)*(len/nels) wid thk;
                                (nel-1)*(len/nels) wid 0];
  n1=(nel-1)*4+1; n2=n1+1; n3=nel*4+2; n4=n3-1;
  n5=n2+2; n6=n5-1; n7=n3+1; n8=n7+1;
  etopol(nel,:)=[n1 n2 n3 n4 n5 n6 n7 n8];
  for gp=1:ngp
    D(:,:,gp,nel)=E/((1+v)*(1-2*v))*(v*Iz+((1-2*v)*I2));
  end
end
etopol(end,:)=[];
bc=[[1:12]' zeros(12,1)];
fext=zeros(3*size(coord,1),1);
fext(3*size(coord,1)-9:3:3*size(coord,1))=1;