function output=sir(N,T,is,rt,r0,gamma,mort)
  
  %initial cases
  I=[];
  S=[];
  R=[];
  D=[];
  dI=[];
  dD=[];
  I(1)=is; % this still seems high
  S(1)=N-is-rt;
  R(1)=rt;
  NewCase=[];
  beta=(r0*gamma); 
  
  
  % IR simulation
  for i = 2:T 
    I(i)=I(i-1)+beta*S(i-1)*I(i-1)/N-gamma*I(i-1);
    dI(i)=beta*S(i-1)*I(i-1)/N;
    dD(i)=mort*gamma*I(i-1);
    R(i)=R(i-1)+gamma*I(i-1);
    S(i)=S(i-1)-beta*I(i-1)*S(i-1)/N;
    D(i)=mort*R(i);
    % new case each day
    NewCase(i)=beta*I(i-1)*S(i-1)/N;
  end
  output.dI=dI;
  output.dD=dD;
  output.I=I;
  
