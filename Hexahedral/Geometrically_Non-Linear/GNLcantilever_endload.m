function [coord,etopol,fext,bc,E,v,ngp,lstps,NRitmax,NRtol]=...
                                      GNLcantilever_endload
E=1.2e6; v=0;
ngp=8; lstps=20; NRitmax=50; NRtol=1.0e-6;
nels=300; len=10; wid=1; thk=0.1;
etopol=zeros(nels,8);
for nel=1:nels+1
  coord((nel-1)*4+1:(nel)*4,:)=[(nel-1)*(len/nels) 0 0;
                                (nel-1)*(len/nels) 0 thk;
                                (nel-1)*(len/nels) wid thk;
                                (nel-1)*(len/nels) wid 0];
  n1=(nel-1)*4+1; n2=n1+1; n3=nel*4+2; n4=n3-1;
  n5=n2+2; n6=n5-1; n7=n3+1; n8=n7+1;
  etopol(nel,:)=[n1 n2 n3 n4 n5 n6 n7 n8];
end
etopol(end,:)=[];
bc=[[1:12]' zeros(12,1)];
fext=zeros(3*size(coord,1),1);
fext(3*size(coord,1)-9:3:3*size(coord,1))=1;