% square lattice
% x is the bound direction
clear all
% tic;
profile on
%-------------------------
for wid = 100
    nlen = wid;
    t = 1;
    %-------------------------
    E = 0.02; % G = 1

    numW = 1;
    Ndis = 100; % ensemble size
               % cond = zeros(Ndis,numW);

    % dosL = zeros(Ndis,numW);
    % dosR = zeros(Ndis,numW);
    %-------------------------
    Nsite = wid*nlen;
    left = 1:wid;
    right = Nsite-wid+1:Nsite;
    %-------------------------main loop
    h0=4*t*eye(wid)-t*diag(ones(1,wid-1),1)-t*diag(ones(1,wid-1),-1);
    h1=-t*diag(ones(1,wid));
    %---------------------selfenergy
    [sLr,sRr,wL,wR]=self(E,h0,h1);
    % lu_T_dos(h0, h1, sLr, sRr, wL, wR, nlen, E, Nsite, left, right, Ndis, wid);
    %=========================================
    recursive_dos(E, h0, h1, wid, nlen, sLr, sRr, Ndis);
    %========================================
    % dirct_inverse_dos(h0, h1, wid, nlen, left, right, sLr, sRr, gR, E, Ndis);
    %----------------------------------
end
p = profile('info');
profsave(p,'profile_results'); % save profile

function [gintr1, gintr2, X]=maxinverse(enep,h0,h1,wid,nlen,sLr,sRr)
    %%%%%%%%%%%%%%%%%%%%%%%%% defining matrix
    h0=sparse(h0);h1=sparse(h1);
    unitm=eye(wid);h2=h1';
    for nx=1:nlen
        F{nx}=zeros(wid);X{nx}=zeros(wid);
        invF{nx}=zeros(wid);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    eminh=sparse(enep*unitm-h0);
    F{1}=eminh-sLr; % left self energy
    invF{1}=inv(F{1});
    for i=1:nlen-2
        F{i+1}=eminh-h2*invF{i}*h1;
        invF{i+1}=inv(F{i+1});
    end
    F{nlen}=eminh-sRr-h2*invF{nlen-1}*h1; %  left self energy
    invF{nlen}=inv(F{nlen});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X{nlen}=invF{nlen};
    for j=(nlen-1):-1:1
        X{j}=invF{j}*h1*X{j+1};
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    gintr1=X{1};
    gintr2=X{nlen};
end

function lu_T_dos(h0, h1, sLr, sRr, wL, wR, nlen, E, Nsite, left, right, Ndis, wid)
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
    for nW=1:Ndis
        [LL1,UU1,PP1,QQ1]=lu(aaa3);
    end
    time1 = toc;

    tic;
    for nW=1:Ndis
        wwL=sparse(PP1(:,left)*wL);
        c11=LL1\wwL;
        wwR=sparse(wR'*QQ1(right,:));
        c12=wwR/UU1;
        c14=c12*c11;
        TT = c14*c14';
        % cond=real(sum(sum(c14.*conj(c14f))));
        cond=real(sum(sum(TT)));
    end
    time2 = toc;
    fprintf('T=%-8f\n',full(cond));
    tic;
    for nW=1:Ndis
        c3 = UU1\c11;
        c4 = QQ1*c3; % Gr!wL>
        dosL = sum(sum(c4.*conj(c4)));
        % dosL = full(sum(sum(c4.*conj(c4))));
    end
    time3 = toc;
    % fprintf('dosL=%-8f\n',full(dosL));
    %--------------------
    tic;
    for nW=1:Ndis
        c5 = sparse(PP1(:,right)*wR);
        c6 = LL1\c5;
        c7 = UU1\c6;
        c8 = QQ1*c7; % Gr!wR>
        dosR = sum(sum(c8.*conj(c8)));
    end
    time4 = toc;
    fprintf('dosR=%-8f\n',full(dosR));
    fprintf('%-5g%-10.5f%-10.5f%-10.5f%-10.5f\n',wid, time1, time2, time3, time4);
end


function dirct_inverse_dos(h0, h1, wid, nlen, left, right, sLr, sRr, gR, E, Ndis)
    tic;
    for nW = 1:Ndis
        Ham0 = kron(speye(nlen),h0);
        Ham1 = kron(diag(ones(1,nlen-1), 1),h1);
        h = Ham0+Ham1+Ham1';
        h(left,left) = h(left,left) + sLr;
        h(right,right) = h(right,right) + sRr;
        gintr = inv(E*eye(wid*nlen)-h);
        gammaR = zeros(wid*nlen, wid*nlen);
        gammaR(right,right) = gR;
        DR2 = gintr*gammaR*gintr';
        dosR2 = real(trace(DR2));
        % diag(gintr*gammaR)
    end
    toc;
    dosR2
end

function dos = recursive_dos(E, h0, h1, wid, nlen, sLr, sRr, Ndis)
    tic;
    for nW = 1:Ndis
        % [gr1, gr2, X] = maxinverse(E,h0,h1,wid,nlen,sLr,sRr);
        % gL = 1i*(sLr - sLr');
        gR = 1i*(sRr - sRr');
        % T = gr1'*gL;
        % cond = real(trace(T));
        % cond = sum(sum(T.*conj(T)));
        %----------------------------------
        [gr1, gr2, X] = maxinverse(E,h0,h1,wid,nlen,sLr,sRr);
        dos = 0;
        for count = 1:wid
            gf = X{count};
            DR = gf*gR;
            dosR = sum(sum(DR.*conj(gf)));
            dos = dos + dosR;
        end
    end
    time5 = toc;
    % cond
    fprintf('recursive get dos:%-10.5f  time: %-10.5f\n',dos, time5);
end