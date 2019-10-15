function [E,v,coord,etopol,bc,f]=cantilever_endload
E=1e9; v=0; nels=50; coord=zeros((nels-1)*4+8,3);
len=1; wid=1; thk=0.1; etopol=zeros(nels,8);
for nel=1:nels+1
  coord((nel-1)*4+1:(nel)*4,:)=[(nel-1)*(1/nels) 0 0;
                                (nel-1)*(1/nels) 0 thk;
                                (nel-1)*(1/nels) wid thk;
                                (nel-1)*(1/nels) wid 0];
  n1=(nel-1)*4+1; n2=n1+1; n3=nel*4+2; n4=n3-1;
  n5=n2+2; n6=n5-1; n7=n3+1; n8=n7+1;
  etopol(nel,:)=[n1 n2 n3 n4 n5 n6 n7 n8];
end
etopol(end,:)=[];
bc=[[1:12]' zeros(12,1)];
f=zeros(3*size(coord,1),1);
f(3*size(coord,1)-9:3:3*size(coord,1))=-0.25e+3;
