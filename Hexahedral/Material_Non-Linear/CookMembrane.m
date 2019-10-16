function [coord,etopol,fext,bc,lam,miu,ngp,lstps,NRitmax,NRtol]=CookMembrane
K=0.9154; miu=0.4225; lam=K-((2*miu)/3);
ngp=8; lstps=10; NRitmax=50; NRtol=1.0e-4;
nelx=2; nely=2; etopol=[]; bc=[];
coord=zeros((nelx+1)*(nely+1)*2,3);
for elz=1:2
  for ely=1:nely+1
    for elx=1:nelx+1
      coord(elx+(ely-1)*(nelx+1)+(elz-1)*((nelx+1)*(nely+1)),:)=[(elx-1)*48/nelx (elx-1)*(44/nelx)+((nelx+1-elx)*(((ely-1)*((44/nely)-(16/nely)))/nelx))+(ely-1)*(16/nely) elz-1];
    end
  end
end
%etopol
for ely=1:nely
  for elx=1:nelx
    etopol=[etopol; elx+(ely-1)*(nelx+1) elx+(ely-1)*(nelx+1)+((nelx+1)*(nely+1)) elx+1+(ely-1)*(nelx+1)+((nelx+1)*(nely+1)) elx+1+(ely-1)*(nelx+1)...
                    elx+(ely*(nelx+1)) elx+(ely*(nelx+1))+((nelx+1)*(nely+1)) elx+1+(ely*(nelx+1))+((nelx+1)*(nely+1)) elx+1+(ely*(nelx+1))];
  end
end
%bc
for n=1:size(coord,1)
  if coord(n,1)==0
    bc=[bc; n*3-2 0; n*3-1 0];
  end
  bc=[bc; n*3 0];
end
%fext
load=1/(4*nely);
fext=zeros(size(coord,1)*3,1);
for n=1:size(coord,1)
  if coord(n,1)==48
    if coord(n,2)==44 || coord(n,2)==60
      fext(n*3-1)=load;
    else
      fext(n*3-1)=2*load;
    end
  end
end
