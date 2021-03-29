% transfer function
% 20 Mar 2000  Renwei
% 03 nov 2005
function [tot,tott]=transfer(h00,h01,ene)
%
t12=ene*eye(size(h00))-h00;
t11=inv(t12);                        %    check1=t11*t12,
     for my=1:2
tau(:,:,my)=zeros(size(h00)); taut(:,:,my)=zeros(size(h00));
     end 	
t12=h01;
tau(:,:,1)=t11*(t12');   taut(:,:,1)=t11*t12;   
tot=tau(:,:,1); tsum=taut(:,:,1);
tott=taut(:,:,1); tsumt=tau(:,:,1);
%
% main loop
%
for m=1:1:40,
t11=zeros;  t12=zeros;
t11=tau(:,:,1)*taut(:,:,1);  t12=taut(:,:,1)*tau(:,:,1);
s1=eye(size(h00))-(t11+t12);
s2=zeros;
s2=inv(s1);                    % check2=s2*s1,
t11=zeros;   t12=zeros;
t11=tau(:,:,1)*tau(:,:,1);    t12=taut(:,:,1)*taut(:,:,1);
tau(:,:,2)=s2*t11;   taut(:,:,2)=s2*t12;
%
t11=zeros;   s1=zeros;
t11=tsum*tau(:,:,2);  s1=tsum*taut(:,:,2);
s2=t11+tot;
tot=s2;  tsum=s1;
t11=zeros;   s1=zeros;
t11=tsumt*taut(:,:,2);   s1=tsumt*tau(:,:,2);
s2=t11+tott;
tott=s2;  tsumt=s1;
tau(:,:,1)=tau(:,:,2);   taut(:,:,1)=taut(:,:,2);
%
conv=0;   conv2=0; 
conv=conv+abs(tau(:,:,2));conv2=conv2+abs(taut(:,:,2));
if ((conv<1e-7) & (conv2<1e-7))
    break
end;

end
