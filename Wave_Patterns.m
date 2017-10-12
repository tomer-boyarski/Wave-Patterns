% time step and duration 
dt=1e-6;
Nt=1e8;
% space array
dx=0.1;
x_max=5;
x=0:dx:x_max;
Nx=length(x);
 
A=zeros(4,Nx); % concentration on previous time step
B=zeros(4,Nx); % concentration on present time step
 
% initial conditions
A(1,:)=80;
A(2,:)=50;
A(3,:)=[1000*ones(1,11),zeros(1,21),1000*ones(1,19)]-10;
A(4,:)=[1000*ones(1,18),zeros(1,33)];
 
% making sure demands on initial conditions are met
total_D=trapz(x,A(1,:)+A(3,:)+A(4,:)) 
total_E=trapz(x,A(2,:)+A(4,:))
 
% various constants
v=-0.0768;
Dm=0.013;
DD = 12.5;
DE = 12.5;
Dde=Dm;
Dd=Dm;
l_dD=9.3e-4;
l_de=0.125;
l_D=0.0013;
l_E=3.8e-5;
l_e=8e-9;
 
% matrix to represent numerical derivation
M = diag(-2.*ones(1,Nx)) + diag(ones(1,Nx-1),1) + diag(ones(1,Nx-1),-1);
M(1,end)=1;
M(end,1)=1;
 
?

% changing in time 
for i=1:Nt
    B(1,:)=dt*(l_de.*A(4,:)...
        -(l_D+l_dD.*A(3,:)).*A(1,:)...
        +transpose(DD*M*transpose(A(1,:))./dx^2))...
        +A(1,:);
    B(2,:)=dt*(l_de.*A(4,:)...
        -(l_E+l_e.*A(4,:).^2).*A(3,:).*A(2,:)...
        +transpose(DE*M*transpose(A(2,:))/dx^2))...
        +A(2,:);
    B(3,:)=dt*((l_D+l_dD.*A(3,:)).*A(1,:)...
        -(l_E+l_e.*A(4,:).^2).*A(3,:).*A(2,:)...
        +transpose(Dd*M*transpose(A(3,:))/dx^2))...
        +A(3,:);
    B(4,:)=dt*((l_E+l_e.*A(4,:).^2).*A(3,:).*A(2,:)...
        -l_de.*A(4,:)...
        +transpose(Dde*M*transpose(A(4,:))/dx^2))...
        +A(4,:);
    A=B;
end
%re-arranging data
[G,I] = max(B(3,:));
B(:,:)=[B(:,I:end),B(:,1:(I-1))];
% plotting
plot(x,10*B(1,:),'k',x,10*B(2,:),'r',x,B(3,:),'y',x,B(4,:),'b','linewidth',1.5);
legend('D','E','d','de');
ylim([0,2000]);
ylabel('Concentration','fontsize',16)
xlabel('x','fontsize',16)
