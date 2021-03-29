%self energy of graphene
%%%%%%%%%%%%%%%%%%%%%%%%
function [h0_l,h1_l]=Ham(num_u,t1,vso,B,uyc)
%-----------------------
sz=[1 0; 0 -1]; id=eye(2);

iso = 1i*vso*sz;
of=(sqrt(3)/4)*linspace(5/2,5/2+(4*num_u-1)*3,4*num_u)*(2*pi*B);

beta=2;
%-----------------------
N_ucell = 8*num_u;
h0_l_hop=zeros(2*N_ucell);h1_l_hop=zeros(2*N_ucell);
h0_l_iso=zeros(2*N_ucell);h1_l_iso=zeros(2*N_ucell);
%-----------------------
%for the left leads, h0 is the hamiltonian of the monolayer
%h1 is the coupling of the left lead to the center.
%for the nearest hopping
ic=0;
for i1=1:4:8*num_u-7
    h0_l_hop(2*i1-1:2*i1,2*(i1+1)-1:2*(i1+1))=-t1*id*exp(-1i*of(1+2*ic)); % 1->2
%     h0_l_hop(2*(i1+1)-1:2*(i1+1),2*(i1+2)-1:2*(i1+2))=-t1*id; % 2->3
    h0_l_hop(2*(i1+2)-1:2*(i1+2),2*(i1+3)-1:2*(i1+3))=-t1*id*exp(+1i*of(2+2*ic)); % 3->4
%     h0_l_hop(2*(i1+3)-1:2*(i1+3),2*(i1+4)-1:2*(i1+4))=-t1*id; % 4->5

%     h0_l_hop(2*i1-1:2*i1,2*(i1+1)-1:2*(i1+1))=-t1*(1+beta*0.25*uyc*(0.25+3*ic))*id; % 1->2
    h0_l_hop(2*(i1+1)-1:2*(i1+1),2*(i1+2)-1:2*(i1+2))=-t1*(1+beta*uyc*(1+3*ic))*id; % 2->3
%     h0_l_hop(2*(i1+2)-1:2*(i1+2),2*(i1+3)-1:2*(i1+3))=-t1*(1+beta*0.25*uyc*(1.75+3*ic))*id; % 3->4
    h0_l_hop(2*(i1+3)-1:2*(i1+3),2*(i1+4)-1:2*(i1+4))=-t1*(1+beta*uyc*(2.5+3*ic))*id; % 4->5
    
%     h0_l_hop(2*i1-1:2*i1,2*(i1+1)-1:2*(i1+1))=-t1*exp(beta*0.25*uyc*(0.5+3*ic))*id; % 1->2
%     h0_l_hop(2*(i1+1)-1:2*(i1+1),2*(i1+2)-1:2*(i1+2))=-t1*exp(beta*uyc*(1.5+3*ic))*id; % 2->3
%     h0_l_hop(2*(i1+2)-1:2*(i1+2),2*(i1+3)-1:2*(i1+3))=-t1*exp(beta*0.25*uyc*(2+3*ic))*id; % 3->4
%     h0_l_hop(2*(i1+3)-1:2*(i1+3),2*(i1+4)-1:2*(i1+4))=-t1*exp(beta*uyc*(3+3*ic))*id; % 4->5
    
    ic=ic+1;
end
h0_l_hop(2*(8*num_u-3)-1:2*(8*num_u-3),2*(8*num_u-2)-1:2*(8*num_u-2))=-t1*id*exp(-1i*of(4*num_u-1));
% h0_l_hop(2*(8*num_u-2)-1:2*(8*num_u-2),2*(8*num_u-1)-1:2*(8*num_u-1))=-t1*id; % 6->7
h0_l_hop(2*(8*num_u-1)-1:2*(8*num_u-1),2*(8*num_u-0)-1:2*(8*num_u-0))=-t1*id*exp(+1i*of(4*num_u));

% h0_l_hop(2*(8*num_u-3)-1:2*(8*num_u-3),2*(8*num_u-2)-1:2*(8*num_u-2))=-t1*(1+beta*0.25*uyc*(0.25+3+6*(num_u-1)))*id; % N-3 -> N-2
h0_l_hop(2*(8*num_u-2)-1:2*(8*num_u-2),2*(8*num_u-1)-1:2*(8*num_u-1))=-t1*(1+beta*uyc*(1+3+6*(num_u-1)))*id; % N-2 -> N-1
% h0_l_hop(2*(8*num_u-1)-1:2*(8*num_u-1),2*(8*num_u-0)-1:2*(8*num_u-0))=-t1*(1+beta*0.25*uyc*(1.75+3+6*(num_u-1)))*id; % N-1 -> N

% h0_l_hop(2*(8*num_u-3)-1:2*(8*num_u-3),2*(8*num_u-2)-1:2*(8*num_u-2))=-t1*exp(beta*0.25*uyc*(0.5+3+6*(num_u-1)))*id; % N-3 -> N-2
% h0_l_hop(2*(8*num_u-2)-1:2*(8*num_u-2),2*(8*num_u-1)-1:2*(8*num_u-1))=-t1*exp(beta*uyc*(1.5+3+6*(num_u-1)))*id; % N-2 -> N-1
% h0_l_hop(2*(8*num_u-1)-1:2*(8*num_u-1),2*(8*num_u-0)-1:2*(8*num_u-0))=-t1*exp(beta*0.25*uyc*(2+3+6*(num_u-1)))*id; % N-1 -> N
%%%%%%%%%%%%
ic=0;
for i1=1:4:8*num_u-7
    h0_l_iso(2*i1-1:2*i1,2*(i1+2)-1:2*(i1+2))=-iso; % 1 -> 3
    h0_l_iso(2*(i1+1)-1:2*(i1+1),2*(i1+3)-1:2*(i1+3))=-iso; % 2 -> 4
    h0_l_iso(2*(i1+2)-1:2*(i1+2),2*(i1+4)-1:2*(i1+4))=+iso; % 3 -> 5
    h0_l_iso(2*(i1+3)-1:2*(i1+3),2*(i1+5)-1:2*(i1+5))=+iso; % 4 -> 6
    ic=ic+1;
end
h0_l_iso(2*(8*num_u-3)-1:2*(8*num_u-3),2*(8*num_u-1)-1:2*(8*num_u-1))=-iso; % 5->7
h0_l_iso(2*(8*num_u-2)-1:2*(8*num_u-2),2*(8*num_u-0)-1:2*(8*num_u-0))=-iso; % 6->8
h0_l=h0_l_iso+h0_l_iso'+h0_l_hop+h0_l_hop';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end of h0_l
%for the coupling of the h1_l
ic=0;
for i1=1:4:8*num_u
    h1_l_hop(2*i1-1:2*i1,2*(i1+1)-1:2*(i1+1))=-t1*id*exp(+1i*of(1+2*ic));
    h1_l_hop(2*(i1+3)-1:2*(i1+3),2*(i1+2)-1:2*(i1+2))=-t1*id*exp(+1i*of(2+2*ic));

%     h1_l_hop(2*i1-1:2*i1,2*(i1+1)-1:2*(i1+1))=-t1*(1+beta*0.25*uyc*(0.25+3*ic))*id; % 1 -> 2
%     h1_l_hop(2*(i1+3)-1:2*(i1+3),2*(i1+2)-1:2*(i1+2))=-t1*(1+beta*0.25*uyc*(1.75+3*ic))*id; % 4 -> 3
    
%     h1_l_hop(2*i1-1:2*i1,2*(i1+1)-1:2*(i1+1))=-t1*exp(beta*0.25*uyc*(0.5+3*ic))*id; % 1 -> 2
%     h1_l_hop(2*(i1+3)-1:2*(i1+3),2*(i1+2)-1:2*(i1+2))=-t1*exp(beta*0.25*uyc*(2+3*ic))*id; % 4 -> 3
 
    h1_l_iso(2*i1-1:2*i1,2*i1-1:2*i1)=-iso; % 1->1
    h1_l_iso(2*(i1+1)-1:2*(i1+1),2*(i1+1)-1:2*(i1+1))=+iso; % 2->2
    h1_l_iso(2*(i1+2)-1:2*(i1+2),2*(i1+2)-1:2*(i1+2))=-iso; % 3->3
    h1_l_iso(2*(i1+3)-1:2*(i1+3),2*(i1+3)-1:2*(i1+3))=+iso; % 4->4
    ic=ic+1;
end
ic=0;
for i1=1:4:8*num_u-7
    h1_l_iso(2*(i1+0)-1:2*(i1+0),2*(i1+2)-1:2*(i1+2))=+iso; % 1->3
    h1_l_iso(2*(i1+3)-1:2*(i1+3),2*(i1+1)-1:2*(i1+1))=-iso; % 4->2
    h1_l_iso(2*(i1+3)-1:2*(i1+3),2*(i1+5)-1:2*(i1+5))=-iso; % 4->6
    h1_l_iso(2*(i1+4)-1:2*(i1+4),2*(i1+2)-1:2*(i1+2))=+iso; % 5->3
    ic=ic+1;
end
h1_l_iso(2*(8*num_u-3)-1:2*(8*num_u-3),2*(8*num_u-1)-1:2*(8*num_u-1))=+iso; % 5->7
h1_l_iso(2*(8*num_u-0)-1:2*(8*num_u-0),2*(8*num_u-2)-1:2*(8*num_u-2))=-iso; % 8->6
%%%%%%%%%%%%%%%%%%%%%%%
h1_l=h1_l_iso+h1_l_hop;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end of the left lead