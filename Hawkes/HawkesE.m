function [RMSE,RMSE_m,D_us,D_us_m,res_all]=HawkesE(H,Th,dt,D_gt,d_start,d_end,n_states,n_simu)

% univariate
R_us=[];
D_us=[];
for i=1:n_states
    ind_s = (H(:,1)==i);
    D_us(i,:)=zeros(1,dt);
    if sum(ind_s)==0
        R_us(i)=0;
    else
        [res]=npThawkes(H(ind_s,2)',ones(1,sum(ind_s)),1,Th);
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

% for j=1:n_states
% for i=1:n_simu
%     rand_idx=randsample(dt,dt,true);
%     RMSE(j,i)=sqrt(mean((D_us(j,rand_idx) - D_gt(j,rand_idx)).^2));
% end
% end

% D_us=sum(D_us,2);
% D_gt=sum(D_gt,2);
% for i=1:n_simu
%     rand_idx=randsample(n_states,n_states,true);
%     RMSE(i)=sqrt(mean((D_us(rand_idx) - D_gt(rand_idx)).^2));
% end

% Multivariate 
H_f = H;
[res_all]=npThawkes(H_f(:,2)',H_f(:,1),n_states,Th);
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

% for j=1:n_states
% for i=1:n_simu
%     rand_idx=randsample(dt,dt,true);
%     RMSE_m(j,i)=sqrt(mean((D_us_m(j,rand_idx) - D_gt(j,rand_idx)).^2));
% end
% end

% D_us_m=sum(D_us_m,2);
% for i=1:n_simu
%     rand_idx=randsample(n_states,n_states,true);
%     RMSE_m(i)=sqrt(mean((D_us_m(rand_idx) - D_gt(rand_idx)).^2));
% end
