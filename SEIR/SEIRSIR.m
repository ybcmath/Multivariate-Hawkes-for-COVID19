function [RMSE_sir,RMSE_seir]=SEIRSIR(TS,dt,popu,is_death,D_gt,n_states)
% grid search parameters
I_sir=[];
I_seir=[];
T=50;
if is_death==1
%     i_start=[0,.001,.005,.01,.05,.1,.5,1,5,10,100,200,500,1000]; %for ca
%     i_start=[0,.001,.005,.01,.05,.1,.5,1,5,10,100,500,1000,1500]; % for simu us death
    i_start=[0,.001,.005,.01,.05,.1,.5,1,5,10,100,500,1000,1500,...
        3000,3500,5000,7500,10000]; % for eu

%     i_start=[0,.001,.005,.01,.05,.1,.5,1,5,10,100]; % for ny death
    gammas=[0,.07,.02:.02:.45];
    mus=[0,0.02:0.02:0.45];
    rs=[0,2.1,1.4:.2:6.5];
    fs=[0,.05:.05:1];
    mort=0.0001;
%     mort=[0.01,0.005,0.002,0.001,0.0005,0.0001]; % mortality

else
    i_start=[.1,.5,1,5,10,20,30,50,60,70,100,1000]; %for SEIR infection
%     i_start=[.1,.5,1,5,10,20,30,50,100]; %for SIR infection
    gammas=[.07,.02:.02:.4];
    mus=[0.02:0.02:0.4];
    rs=[2.1,1.4:.2:6];
    fs=[.05:.05:1];
    mort=.01; % mortality
end



for i=1:n_states
% y=table2array(TS(i,d_start:d_end));
% state_name=table2array(TS(i,1));

y=TS(i,:);
state_name=i;

N=popu(i);
M=length(y);

% SIR
best=-300000000;
for mt = 1:length(mort)
for r0 = 1:length(rs)
  for is = 1:length(i_start)
    for gamma = 1:length(gammas)
      
      output_sir=sir(N,T,i_start(is),0,rs(r0),gammas(gamma),mort(mt));
      dD=output_sir.dD;
      dI=output_sir.dI;
      I=output_sir.I;
      
      if (is_death==1)
        sirPred=dD;
      else
        sirPred=dI;  
      end
      
      tbest=sum(y.*log(sirPred(1:M)+10^-10))-sum(sirPred(1:M));
      
      if(tbest>best)
        r_best=rs(r0);
        g_best=gammas(gamma);
        i_best=i_start(is);
        mt_best=mort(mt);
        best=tbest;
        dDbest=dD;
        dIbest=dI;
        sirPredbest(i,:)=sirPred(M+1:M+dt);
        Ibest=I;
%         D_sir(i)=sum(dDbest(1:M+dt-1));
%         dI_sir(i,:)=dIbest(M+1:M+dt);
        [imax,ix]=max(dIbest);
        peak_best=ix;
        AIC_sir(i)=2*3-2*tbest;
      end
     end
   end
end
end
{'sir',r_best,g_best,i_best,best,peak_best,mt_best,state_name}

% SEIR
best=-300000000;
for mt = 1:length(mort)
for r0 = 1:length(rs)
  for is = 1:length(i_start)
    for gamma = 1:length(gammas)
      for mu = 1:length(mus)
      
          output_seir=seir(N,T,i_start(is),0,rs(r0),gammas(gamma),mus(mu),mort(mt));
          dD=output_seir.dD;
          dI=output_seir.dI;
          I=output_seir.I;
          dmuE=output_seir.muE;
          
          if (is_death==1)
            seirPred=dD;
          else
            seirPred=dmuE;  
          end          
          tbest=sum(y.*log(seirPred(1:M)+10^-10))-sum(seirPred(1:M));
      
          if (tbest>best)
            r_best=rs(r0);
            g_best=gammas(gamma);
            i_best=i_start(is);
            mt_best=mort(mt);
            best=tbest;
            dDbest=dD;
            dIbest=dI;
            dmuEbest=dmuE;
            Ibest=I;
            mu_best=mus(mu);
            seirPredbest(i,:)=seirPred(M+1:M+dt);
            [imax,ix]=max(dIbest);
            peak_best=ix;%dates(ix);
            AIC_seir(i)=2*4-2*tbest;
          end
      end
    end
  end
end
end
{'seir',r_best,g_best,i_best,best,peak_best,mu_best,mt_best,state_name}
end
%
% RMSE_sir=sqrt(mean((D_sir - D_gt).^2));
% RMSE_seir=sqrt(mean((D_seir - D_gt).^2));

RMSE_sir=[];
for i=1:100
    rand_idx=randsample(n_states*dt,n_states*dt,true);
    RMSE_sir(i)=sqrt(mean((sirPredbest(rand_idx) - D_gt(rand_idx)).^2));
end

RMSE_seir=[];
for i=1:100
    rand_idx=randsample(n_states*dt,n_states*dt,true);
    RMSE_seir(i)=sqrt(mean((seirPredbest(rand_idx) - D_gt(rand_idx)).^2));
end

% 
% sum(AIC_seir)
% sum(AIC_sir)