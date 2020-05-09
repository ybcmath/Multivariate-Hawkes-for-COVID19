function output=seir(N,T,is,rt,r0,gamma,mu,mort)
% initial cases
I=[];
S=[];
R=[];
D=[];
E=[];
dI=[];
dD=[];
I(1)=is; % this still seems high
S(1)=N-is-rt;
R(1)=rt;
E(1)=0;
NewCase=[];

beta=(r0*gamma);

% IR simulation
for i = 2:T
    I(i)=I(i-1)+mu*E(i-1)-gamma*I(i-1);
    dI(i)=mu*E(i-1)-gamma*I(i-1);
    dD(i)=mort*gamma*I(i-1);
    R(i)=R(i-1)+gamma*I(i-1);
    E(i)=E(i-1)+beta*I(i-1)*S(i-1)/N-mu*E(i-1);
    S(i)=S(i-1)-beta*I(i-1)*S(i-1)/N;
    D(i)=mort*R(i);
    % new case each day
    NewCase(i)=mu*E(i-1)-gamma*I(i-1);
end

output.dI=dI;
output.dD=dD;
output.I=I;
output.muE=mu*E;
  

