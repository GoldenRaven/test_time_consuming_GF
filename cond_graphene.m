% square lattice
% x is the bound direction
clear all
tic;
%-------------------------
wid = 90;
nlen=90;
t = 1;
%-------------------------
E = 0.2; % G = 1

numW = 1;
ww = linspace(0.1,0.1,numW);
Ndis = 1; % ensemble size
cond = zeros(Ndis,numW);

dosL = zeros(Ndis,numW);
dosR = zeros(Ndis,numW);
%-------------------------
Nsite = wid*nlen;
left = 1:wid;
right = Nsite-wid+1:Nsite;
%-------------------------main loop
h0=4*t*eye(wid)-t*diag(ones(1,wid-1),1)-t*diag(ones(1,wid-1),-1);
h1=-t*diag(ones(1,wid));
%---------------------selfenergy
[sLr,sRr,wL,wR]=self(E,h0,h1);

%----------------------------------
Ham0 = kron(speye(nlen),h0);
Ham1 = kron(diag(ones(1,nlen-1), 1),h1);
HamC = sparse(Ham0+Ham1+Ham1');
clear Ham0 Ham1 h0 h1

aa = E*ones(1,Nsite);
aa = sparse(1:Nsite,1:Nsite,aa,Nsite,Nsite);
aaa3 = aa - HamC;
aaa3(left,left)   = aaa3(left,left)   - sLr;
aaa3(right,right) = aaa3(right,right) - sRr;
aaa3 = sparse(aaa3);
clear aa HamC
%-----------------------

% xt=clock;
% xt1=xt(6);
% xt2=100*(xt1);
% iu=0;
% while iu-xt2<=0
%     iu=iu+1;xt3=rand(1);
% end

%---------------------
tic;
for nW=1:100
    [LL1,UU1,PP1,QQ1]=lu(aaa3);
end
toc;

tic;
for nW=1:100
    wwL=sparse(PP1(:,left)*wL);
    c11=LL1\wwL;
    wwR=sparse(wR'*QQ1(right,:));
    c12=wwR/UU1;
    c14=c12*c11;
    TT = c14*c14';
    cond=real(trace(TT));
end
toc;
% fprintf('T=%-8f\n',cond);
%---------------------
%         c1 = sparse(PP1(:,left)*wL);
%         c2 = LL1\c1;
%         c3 = UU1\c2;
%         c4 = QQ1*c3; % Gr!wL>
tic;
for nW=1:100
    c3 = UU1\c11;
    c4 = QQ1*c3; % Gr!wL>
    dosL = full(sum(sum(c4.*conj(c4))));
end
toc;
% fprintf('dosL=%-8f\n',dosL);
%--------------------
tic;
for nW=1:100
    c5 = sparse(PP1(:,right)*wR);
    c6 = LL1\c5;
    c7 = UU1\c6;
    c8 = QQ1*c7; % Gr!wR>
    dosR = full(sum(sum(c8.*conj(c8))));
end
toc;
% fprintf('dosR=%-8f\n',dosR);