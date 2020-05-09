H=[];
simu=seirsim;
N=0;
n_states=length(simu);
d_start=1;
d_end=12;
for i=d_start:d_end
    n_case = simu(:,i);
    for j=1:n_states
        H=[H;repmat([j,i-d_start],n_case(j),1)];
    end
end

T=2;
dt=3;
n_simu=100;
D_gt=simu(:,d_end+1:d_end+dt);
% univariate
R_us=[];
D_us=[];
for i=1:n_states
    ind_s = (H(:,1)==i);
    D_us(i,:)=zeros(1,dt);
    if sum(ind_s)==0
        R_us(i)=0;
    else
        [res]=npThawkes(H(ind_s,2)',ones(1,sum(ind_s)),1,T);
        R_us(i)=res.K;
        start=[];
        start.T=d_end-d_start;
        start.n=length(H(ind_s,2));
        start.t=H(ind_s,2)';
        start.type=ones(1,sum(ind_s));
        for j=1:n_simu
            y=simu_spetas_import(dt,res.r,res.K,res.g,res.mid_t,start);
            for k=1:dt
                D_us(i,k)=D_us(i,k)+sum((y.t<=d_end-d_start+k)&(d_end-d_start+k-1<=y.t));
            end
        end
    end
end
D_us=D_us./n_simu;
RMSE=[];
for i=1:n_simu
    rand_idx=randsample(n_states*dt,n_states*dt,true);
    RMSE(i)=sqrt(mean((D_us(rand_idx) - D_gt(rand_idx)).^2));
end
RMSE_se=std(RMSE);
RMSE_avg=mean(RMSE);
% Multivariate 
H_f = H;
[res_all]=npThawkes(H_f(:,2)',H_f(:,1),n_states,0);
start=[];
start.T=H_f(end,2);
start.n=length(H_f);
start.t=H_f(:,2)';
start.type=H_f(:,1)';
D_us_m=zeros(n_states,dt);
for i=1:n_simu
    y_all=simu_spetas_import(dt,res_all.r,res_all.K,res_all.g,res_all.mid_t,start);
    for j=1:n_states
        state_ts=y_all.t(y_all.type==j);
        for k=1:dt
            D_us_m(j,k)=D_us_m(j,k)+sum((state_ts<=d_end-d_start+k)&(d_end-d_start+k-1<=state_ts));
        end
    end
end
D_us_m=D_us_m./n_simu;
RMSE_m=[];
for i=1:n_simu
    rand_idx=randsample(n_states*dt,n_states*dt,true);
    RMSE_m(i)=sqrt(mean((D_us_m(rand_idx) - D_gt(rand_idx)).^2));
end
RMSE_m_se=std(RMSE_m);
RMSE_m_avg=mean(RMSE_m);
%% simulation
% A1=importdata('A1.txt')'; %transpose due to different definition
% mu1=importdata('mu1.txt')';
% y=simu_exp(1e4,mu1,A1,1);
% [res_all]=npThawkes(y.t',y.type,10,0);

