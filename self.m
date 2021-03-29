%self energy of graphene
%%%%%%%%%%%%%%%%%%%%%%%%
function [sLr,sRr,wL,wR]=self(e0,H0,H1)
%-----------------------
h0_l=H0; h1_l =H1;
%--------------------------------------------------------------------------
ene=e0+complex(0,0.0000001);

[totL,tottL]=transfer(h0_l,h1_l,ene);
grL=green(totL,tottL,h0_l,h1_l,ene,-1);
sLr=h1_l'*grL*h1_l;
[totR,tottR]=transfer(h0_l,h1_l,ene);
grR=green(totR,tottR,h0_l,h1_l,ene,1);
sRr=h1_l*grR*h1_l';
gL=1i*(sLr-sLr');
gR=1i*(sRr-sRr');
%-----------------------------------
[VL,DL]=eig(full(gL));
numberL=length(DL);
icL=0;
for i1=1:numberL
    if DL(i1,i1)<0.01
        icL=icL+0;
    else
        icL=icL+1;
        lamdaL(icL)=DL(i1,i1);
        psaiL(:,icL)=VL(:,i1);
    end
end
omigaL = zeros(length(gL),icL);
wL = zeros(length(gL),icL);
for j1=1:icL
    omigaL(:,j1)=sqrt(lamdaL(j1))*psaiL(:,j1);
    wL(:,j1)=omigaL(:,j1);
end
%------------------------
[VR,DR]=eig(full(gR));
numberR=length(DR);
icR=0;
for i1=1:numberR
    if DR(i1,i1)<0.01
        icR=icR+0;
    else
        icR=icR+1;
        lamdaR(icR)=DR(i1,i1);
        psaiR(:,icR)=VR(:,i1);
    end
end
omigaR = zeros(length(gR),icR);
wR = zeros(length(gR),icR);
for j1=1:icR
    omigaR(:,j1)=sqrt(lamdaR(j1))*psaiR(:,j1);
    wR(:,j1)=omigaR(:,j1);
end