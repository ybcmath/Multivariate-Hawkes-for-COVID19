function [RMSE,RMSE_m,D_ca_cox,D_ca_m_cox,res_all_cox]=CoxHawkesE(H,T,dt,D_gt,d_start,d_end,n_states,n_simu)
% univariate
R_ca_cox={};
D_ca_cox=[];
for i=1:n_states
    ind_s = (H(:,1)==i);
    D_ca_cox(i,:)=zeros(1,dt);
    if sum(ind_s)==0
        R_ca_cox{i}=0;
    else
        [res_cox]=npCoxhawkes(H(ind_s,2)',ones(1,sum(ind_s))',1,T);
        R_ca_cox{i}=res_cox.K;
        start=[];
        start.T=d_end-d_start;
        start.n=length(H(ind_s,2));
        start.t=H(ind_s,2)';
        start.type=ones(1,sum(ind_s));
        for j=1:n_simu
            y=simu_spetas_import(dt,res_cox.r,res_cox.K(:,:,end),res_cox.g,res_cox.mid_t,start);
            for k=1:dt
                D_ca_cox(i,k)=D_ca_cox(i,k)+sum((y.t<=d_end-d_start+k)&(d_end-d_start+k-1<=y.t));
            end
        end
    end
end
D_ca_cox=D_ca_cox./n_simu;
RMSE=[];
for i=1:n_simu
    rand_idx=randsample(n_states*dt,n_states*dt,true);
    RMSE(i)=sqrt(mean((D_ca_cox(rand_idx) - D_gt(rand_idx)).^2));
end
% D_gt=sum(D_gt,2);
% D_ca_cox=sum(D_ca_cox,2);
% for i=1:n_simu
%     rand_idx=randsample(n_states,n_states,true);
%     RMSE(i)=sqrt(mean((D_ca_cox(rand_idx) - D_gt(rand_idx)).^2));
% end

% Multivariate 
H_f = H;
[res_all_cox]=npCoxhawkes(H_f(:,2)',H_f(:,1),n_states,T);

start.T=H_f(end,2);
start.n=length(H_f);
start.t=H_f(:,2)';
start.type=H_f(:,1)';
D_ca_m_cox=zeros(n_states,dt);
for i=1:n_simu
    y_all=simu_spetas_import(dt,res_all_cox.r,res_all_cox.K(:,:,end),res_all_cox.g,res_all_cox.mid_t,start);
    for j=1:n_states
        state_ts=y_all.t(y_all.type==j);
        for k=1:dt
            D_ca_m_cox(j,k)=D_ca_m_cox(j,k)+sum((state_ts<=d_end-d_start+k)&(d_end-d_start+k-1<=state_ts));
        end
    end
end
D_ca_m_cox=D_ca_m_cox./n_simu;
RMSE_m=[];
for i=1:n_simu
    rand_idx=randsample(n_states*dt,n_states*dt,true);
    RMSE_m(i)=sqrt(mean((D_ca_m_cox(rand_idx) - D_gt(rand_idx)).^2));
end
% D_ca_m_cox=sum(D_ca_m_cox,2);
% for i=1:n_simu
%     rand_idx=randsample(n_states,n_states,true);
%     RMSE_m(i)=sqrt(mean((D_ca_m_cox(rand_idx) - D_gt(rand_idx)).^2));
% end