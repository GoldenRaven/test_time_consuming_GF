%green fuction
%renwei 20mar 2002
%
function g=green(tot,tott,h0,h1,ene,igreen)
% wid=size(h0);
% wid=wid(1);
switch igreen
    
case 1    
igreen=1;
s1=h1*tot;
eh0=ene*eye(size(h0))-h0-s1;
g=inv(eh0);

otherwise %-1
igreen=-1;
s1=h1'*tott;
eh0=ene*eye(size(h0))-h0-s1;
g=inv(eh0);

end




 
