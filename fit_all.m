%% simu-SIR
H=[];
simu=sirsim;
N=0;
n_states=length(simu);
d_start=1;
d_end=14;
for i=d_start:d_end
    n_case = simu(:,i);
    for j=1:n_states
        H=[H;repmat([j,i-d_start],n_case(j),1)];
    end
end
T=1;
dt=14;
n_simu=100;
D_gt=simu(:,d_end+1:d_end+dt);

% Hawkes
[RMSE_simu,RMSE_simu_m]=HawkesE(H,T,dt,D_gt,d_start,d_end,n_states,n_simu);
RMSE_se=std(RMSE_simu);
RMSE_avg=mean(RMSE_simu);
RMSE_m_se=std(RMSE_simu_m);
RMSE_m_avg=mean(RMSE_simu_m);

% SIRs
TS=sirsim(:,d_start:d_end);
is_death=0;
[RMSE_simu_sir,RMSE_simu_seir]=SEIRSIR(TS,dt,popu,is_death,D_gt,n_states);
RMSE_sir_se=std(RMSE_simu_sir);
RMSE_sir_avg=mean(RMSE_simu_sir);
RMSE_seir_se=std(RMSE_simu_seir);
RMSE_seir_avg=mean(RMSE_simu_seir);

%% simu-SEIR
H=[];
simu=seirsim;
N=0;
n_states=length(simu);
d_start=1;
d_end=14;
for i=d_start:d_end
    n_case = simu(:,i);
    for j=1:n_states
        H=[H;repmat([j,i-d_start],n_case(j),1)];
    end
end
T=2;
dt=14;
n_simu=100;
D_gt=simu(:,d_end+1:d_end+dt);

% Hawkes
[RMSE_simu,RMSE_simu_m]=HawkesE(H,T,dt,D_gt,d_start,d_end,n_states,n_simu);
RMSE_se=std(RMSE_simu);
RMSE_avg=mean(RMSE_simu);
RMSE_m_se=std(RMSE_simu_m);
RMSE_m_avg=mean(RMSE_simu_m);

% SIRs
TS=seirsim(:,d_start:d_end);
is_death=0;
[RMSE_simu_sir,RMSE_simu_seir]=SEIRSIR(TS,dt,popu,is_death,D_gt);
RMSE_sir_se=std(RMSE_simu_sir);
RMSE_sir_avg=mean(RMSE_simu_sir);
RMSE_seir_se=std(RMSE_simu_seir);
RMSE_seir_avg=mean(RMSE_simu_seir);
%% EU
% ts_eu=timeseriescovid19deathseuro;
H=[];
N=0;
n_states=height(ts_eu);
d_start=30;
d_end=60;
for i=d_start:d_end
    n_case = table2array(ts_eu(:,i))-table2array(ts_eu(:,i-1));
    for j=1:n_states
        H=[H;repmat([j,i-d_start],n_case(j),1)];
    end
end

T=2; %2 for early stage
dt=14;
n_simu=100;
%hawkes
D_gt=diff(table2array(ts_eu(:,d_end:dt+d_end)),1,2);
[RMSE_ca,RMSE_ca_m,D_eu,D_eu_m,res_all]=HawkesE(H,T,dt,D_gt,d_start,d_end,n_states,n_simu);
RMSE_se=std(RMSE_ca');
RMSE_avg=mean(RMSE_ca');
RMSE_m_se=std(RMSE_ca_m');
RMSE_m_avg=mean(RMSE_ca_m');
% Italy
RMSE_m_Italy=sqrt(mean((D_eu_m(10,:) - D_gt(10,:)).^2));
RMSE_Italy=sqrt(mean((D_eu(10,:) - D_gt(10,:)).^2));
i_state=10;
plot(D_us_m(i_state,:))
hold on
plot(D_us(i_state,:))
plot(D_gt(i_state,:))

%CoxHawkes
[RMSE_ca_cox,RMSE_ca_cox_m,D_eu_cox,D_eu_m_cox,res_all_cox]=CoxHawkesE(H,T,dt,D_gt,d_start,...
d_end,n_states,n_simu);
RMSE_cox_se=std(RMSE_ca_cox);
RMSE_cox_avg=mean(RMSE_ca_cox);
RMSE_cox_m_se=std(RMSE_ca_cox_m);
RMSE_cox_m_avg=mean(RMSE_ca_cox_m);
% Italy
RMSE_mcox_Italy=sqrt(mean((D_eu_m_cox(10,:) - D_gt(10,:)).^2));
RMSE_cox_Italy=sqrt(mean((D_eu_cox(10,:) - D_gt(10,:)).^2));

% SIRs
is_death=1;
TS=table2array(ts_eu(:,d_start:d_end)); 
% popu=table2array(counties(counties.State=='CA',7));
% popu=popu(2,:);
[RMSE_ca_sir,RMSE_ca_seir]=SEIRSIR(TS,dt,popu,is_death,D_gt,n_states);
RMSE_sir_se=std(RMSE_ca_sir);
RMSE_sir_avg=mean(RMSE_ca_sir);
RMSE_seir_se=std(RMSE_ca_seir);
RMSE_seir_avg=mean(RMSE_ca_seir);
%networks
delta_k = linspace(d_start,d_end,4);
for i=1:3
    ncase=table2array(ts_eu(1:n_states,delta_k(i+1)));
    adj=res_all_cox.K(:,:,i).*repmat(ncase, [1 n_states]);
    adj(adj<1)=0;
    G=digraph(adj);
    ncolor=zeros(n_states,3);
    ncolor(ncase>0,3)=1;
    nname=ts_eu.CountryRegion;
    nname(ncase==0)='';
    figure
    LWidths = 5*G.Edges.Weight/max(G.Edges.Weight);
    plot(G,'LineWidth',LWidths,'XData',ts_eu.Long,'YData',ts_eu.Lat,...
        'NodeLabel',nname,'NodeColor',ncolor)
end
%% US
ts_us = grpstats(timeseriescovid19deathsUS,'Province_State','sum');
H=[];
N=0;
n_states=height(ts_us);
d_start=41;
d_end=61;
for i=d_start:d_end
    n_case = table2array(ts_us(:,i))-table2array(ts_us(:,i-1));
    for j=1:n_states
        H=[H;repmat([j,i-d_start],n_case(j),1)];
    end
end
T=1;
dt=14;
n_simu=1000;
% Hawkes
D_gt=diff(table2array(ts_us(:,d_end:dt+d_end)),1,2);
[RMSE_us,RMSE_us_m,D_us,D_us_m]=HawkesE(H,T,dt,D_gt,d_start,d_end,n_states,n_simu);
RMSE_se=std(RMSE_us);
RMSE_avg=mean(RMSE_us);
RMSE_m_se=std(RMSE_us_m);
RMSE_m_avg=mean(RMSE_us_m);
% CA
i_state=33;
RMSE_m_CA=sqrt(mean((D_us_m(i_state,:) - D_gt(i_state,:)).^2));
RMSE_CA=sqrt(mean((D_us(i_state,:) - D_gt(i_state,:)).^2));
plot(D_us_m(i_state,:))
hold on
plot(D_us(i_state,:))
plot(D_gt(i_state,:))

%CoxHawkes
[RMSE_us_cox,RMSE_us_cox_m,res_all_cox]=CoxHawkesE(H,T,dt,D_gt,d_start,...
    d_end,n_states,n_simu);
RMSE_cox_se=std(RMSE_us_cox);
RMSE_cox_avg=mean(RMSE_us_cox);
RMSE_cox_m_se=std(RMSE_us_cox_m);
RMSE_cox_m_avg=mean(RMSE_us_cox_m);
% SIRs
is_death=1;
TS_us=table2array(ts_us(:,d_start:d_end));
[RMSE_us_sir,RMSE_us_seir]=SEIRSIR(TS_us,dt,popu,is_death,D_gt,n_states);
RMSE_sir_se=std(RMSE_us_sir);
RMSE_sir_avg=mean(RMSE_us_sir);
RMSE_seir_se=std(RMSE_us_seir);
RMSE_seir_avg=mean(RMSE_us_seir);
% network 
% for i=1:4
%     adj=res_all_cox.K(:,:,i).*repmat(table2array(ts_us(1:n_states,d_end)), [1 n_states]);
%     adj(adj<1e-2)=0;
%     % adj(adj>10)=0;
%     % imagesc(log(adj))
%     % [u,A,w,lkh,p,para,aic] = tempestim(H(1:10000,:),33);
%     G=digraph(adj);
%     LWidths = 5*G.Edges.Weight/max(G.Edges.Weight);
%     figure
%     plot(G,'LineWidth',LWidths,'XData',long,'YData',lat,'NodeLabel',states)
% end
% network 
delta_k = linspace(d_start,d_end,5);
for i=1:4
    ncase=table2array(ts_us(1:n_states,delta_k(i+1)));
    adj=res_all_cox.K(:,:,i).*repmat(table2array(ts_us(1:n_states,d_end)), [1 n_states]);
    adj(adj<1e-2)=0;
    G=digraph(adj);
    ncolor=zeros(n_states,3);
    ncolor(ncase>0,3)=1;
    nname=states;
    nname(ncase==0)='';
    figure
    LWidths = 5*G.Edges.Weight/max(G.Edges.Weight);
    plot(G,'LineWidth',LWidths,'XData',long,'YData',lat,...
        'NodeLabel',nname,'NodeColor',ncolor)
end
%% CA
TS_ca = timeseriescovid19deathsUS(timeseriescovid19deathsUS.Province_State=="California",:);
TS_ca(59:end,:) = [];
H=[];
N=0;
n_states=height(TS_ca);
d_start=53;
d_end=82;
% d_end=84;
for i=d_start:d_end
    n_case = table2array(TS_ca(:,i))-table2array(TS_ca(:,i-1));
    for j=1:n_states
        H=[H;repmat([j,i-d_start],n_case(j),1)];
    end
end

T=2;
dt=14;
n_simu=1000;
%hawkes
D_gt=diff(table2array(TS_ca(:,d_end:dt+d_end)),1,2);
[RMSE_ca,RMSE_ca_m,D_ca,D_ca_m,res_all]=HawkesE(H,T,dt,D_gt,d_start,d_end,n_states,n_simu);
RMSE_se=std(RMSE_ca);
RMSE_avg=mean(RMSE_ca);
RMSE_m_se=std(RMSE_ca_m);
RMSE_m_avg=mean(RMSE_ca_m);
% LA
RMSE_m_la=sqrt(mean((D_ca_m(19,:) - D_gt(19,:)).^2));
RMSE_la=sqrt(mean((D_ca(19,:) - D_gt(19,:)).^2));
plot(D_ca_m(19,:))
hold on
plot(D_gt(19,:))
plot(D_ca(19,:))
%CoxHawkes
[RMSE_ca_cox,RMSE_ca_cox_m,res_all_cox]=CoxHawkesE(H,T,dt,D_gt,d_start,...
d_end,n_states,n_simu);
RMSE_cox_se=std(RMSE_ca_cox);
RMSE_cox_avg=mean(RMSE_ca_cox);
RMSE_cox_m_se=std(RMSE_ca_cox_m);
RMSE_cox_m_avg=mean(RMSE_ca_cox_m);
% SIRs
is_death=1;
TS=table2array(TS_ca(:,d_start:d_end)); 
% popu=table2array(counties(counties.State=='CA',7));
% popu=popu(2,:);
[RMSE_ca_sir,RMSE_ca_seir]=SEIRSIR(TS,dt,popu,is_death,D_gt,n_states);
RMSE_sir_se=std(RMSE_ca_sir);
RMSE_sir_avg=mean(RMSE_ca_sir);
RMSE_seir_se=std(RMSE_ca_seir);
RMSE_seir_avg=mean(RMSE_ca_seir);

% network
%     ds=ceil((d_end-d_start)/3)*i+d_start;
ncase=table2array(TS_ca(1:n_states,d_end));
adj=res_all.K.*repmat(ncase, [1 n_states]);
adj(adj<2)=0;
G=digraph(adj);
ncolor=zeros(n_states,3);
ncolor(ncase>0,3)=1;
nname=cacounties;
nname(ncase==0)='';
figure
LWidths = 5*G.Edges.Weight/max(G.Edges.Weight);
plot(G,'LineWidth',LWidths,'XData',TS_ca.Long_,'YData',TS_ca.Lat,...
    'NodeLabel',nname,'NodeColor',ncolor)
%% NY
%data
TS_ny = timeseriescovid19deathsUS(timeseriescovid19deathsUS.Province_State=="New York",:);
TS_ny(62:end,:) = [];
H=[];
N=0;
n_states=height(TS_ny);
d_start=60;
d_end=80;
% d_end=84;
for i=d_start:d_end
    n_case = table2array(TS_ny(:,i))-table2array(TS_ny(:,i-1));
    for j=1:n_states
        H=[H;repmat([j,i-d_start],n_case(j),1)];
    end
end
T=1;
dt=14;
D_gt=diff(table2array(TS_ny(:,d_end:dt+d_end)),1,2);
n_simu=1000;

%Hawkes
[RMSE_ny,RMSE_ny_m,D_ny,D_ny_m,res_all]=HawkesE(H,T,dt,D_gt,d_start,d_end,n_states,n_simu);
RMSE_se=std(RMSE_ny);
RMSE_avg=mean(RMSE_ny);
RMSE_m_se=std(RMSE_ny_m);
RMSE_m_avg=mean(RMSE_ny_m);
%CoxHawkes
[RMSE_ny_cox,RMSE_ny_cox_m]=CoxHawkesE(H,T,dt,D_gt,d_start,d_end,n_states,n_simu);
RMSE_cox_se=std(RMSE_ny_cox);
RMSE_cox_avg=mean(RMSE_ny_cox);
RMSE_cox_m_se=std(RMSE_ny_cox_m);
RMSE_cox_m_avg=mean(RMSE_ny_cox_m);
%SIRs
TS=table2array(TS_ny(:,d_start:d_end)); 
% popu=table2array(counties(counties.State=='NY',7));
% popu=popu(2:end);
is_death=1;
[RMSE_sir,RMSE_seir]=SEIRSIR(TS,dt,popu,is_death,D_gt,n_states);
RMSE_sir_se=std(RMSE_sir);
RMSE_sir_avg=mean(RMSE_sir);
RMSE_seir_se=std(RMSE_seir);
RMSE_seir_avg=mean(RMSE_seir);
%STHawkes
[RMSE_st_ny,RMSE_st_ny_m]=STHawkesE(H,T,dt,D_gt,d_start,d_end,n_states,TS_ny.Lat,TS_ny.Long_);
RMSE_st_se=std(RMSE_st_ny);
RMSE_st_avg=mean(RMSE_st_ny);
RMSE_st_m_se=std(RMSE_st_ny_m);
RMSE_st_m_avg=mean(RMSE_st_ny_m);

