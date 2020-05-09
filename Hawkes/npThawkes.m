function [res]=npThawkes(time,topics,dim,Th)
% multivariate nonparametric temporal Hawkes process
% time: 1*n
% topic: labels of events same order as time totally dim kinds of topics

% Baichuan Yuan 5/2/2020

% Hyperparameters: 
nbins_t=10; %temporal grid nbins_t
lambda=0; %L1 regularization
min_dt=1; %no connections if dt of events <= min_dt

T=time(end);
Th=T-Th;
N=length(time);
Nh=sum(time<=Th);

%calculate temporal difference 
deltat = triu(bsxfun(@minus, time, repmat(time', [1 N])));
deltat = triu(max(deltat,1.0e-6),1);

% log scale grid
% delta_t  = 10.^linspace(-6,log10(max(max(deltat))+1),nbins_t+1);

% linear grid
delta_t = linspace(0,max(max(deltat))+1,nbins_t+1);

% calculate grid index for each point
ind_t=zeros(N);
mid_t=zeros(nbins_t,1);
for i=1:nbins_t
    ind_t=ind_t+i*(deltat>=delta_t(i) & deltat <delta_t(i+1));
    mid_t(i)=(delta_t(i)+delta_t(i+1))*0.5;
end

u_ind_t=triu(ind_t,1);
u_ind_t=u_ind_t(u_ind_t>0);

ind_bin_t=cell(1,nbins_t);
for i=1:nbins_t
    ind_bin_t{i}=(ind_t==i);
end

bin_t=diff(delta_t);

%% E-M algorithm
K=zeros(dim);
r_k=zeros(1,dim);
g=zeros(1,nbins_t);
P=triu(rand(N));
deno=sum(P);
P=bsxfun(@times,P,1./deno);

%topic index
topic_ind=cell(1,dim);
for i=1:dim
	topic_ind{i}=(topics==i);
end
topic_num=[];
for i=1:dim
    topic_num(i)=sum(topics(1:Nh)==i);
end

for kk=1:1000
    % M step
    Pnd=P-spdiags(diag(P), 0, N, N);
	%update K
    for i=1:dim
        for j=1:dim
            if (topic_num(i)<=lambda) %||(i~=j) diag-only version
                K(i,j)=0;
            else
                K(i,j)=sum(sum(Pnd(topic_ind{i}, topic_ind{j})))/(topic_num(i)-lambda);
            end
        end
    end
	%update g
    sum_p=sum(sum(triu(P,1)));
    for i=1:nbins_t
        g(i)=sum(P(ind_bin_t{i}))/(sum_p*bin_t(i));
    end
    %cal background
    diag_p=diag(P);
    for i=1:dim
       r_k(i)=sum(diag_p(topic_ind{i}))/Th;
    end
    % E step

	%update P
    Pnew=zeros(N);
    for j=1:N
        for i=1:j-1
            Pnew(i,j)=K(topics(i),topics(j))*g(ind_t(i,j)); 
        end
        if j<=Nh
            Pnew(j,j)=r_k(topics(j));
        else
            Pnew(j,j)=0;
        end
    end
    temp_diag=diag(Pnew);
	%remove near events in time
    Pnew(deltat<=min_dt)=0;
    Pnew=Pnew+diag(temp_diag);
    deno=sum(Pnew);
    deno(deno==0)=1;
    Pnew=bsxfun(@times,Pnew,1./deno);
	%determine convergence 
    err= max(max(abs(Pnew-P)));
	nb=sum(diag(Pnew));
	fprintf('iter %d: error = %g, # of background = %g\n', kk, err, nb);
	P=Pnew;
    if err<1e-03
        break
    end
end

%% output result
Pnd=P-spdiags(diag(P), 0, N, N);
for i=1:dim
	for j=1:dim
		theta_k(i,j)=sum(sum(Pnd(topics==i, topics==j)))/sum(sum(Pnd(topics==i, :)));
		var_k(i,j)=sum(sum(Pnd(topics==i, :)))*theta_k(i,j)*(1-theta_k(i,j))/(sum(topics==i))^2;
	end
end
res.delta_t=delta_t;
res.g=g;
res.K=K;
res.P=P;
res.r=r_k;
res.mid_t=mid_t;
res.var_k=var_k;
res.theta_k=theta_k;

