function y=simu_spetas_import(T,mu,A,g,mid_t,start)
%simulate multivariate ETAS in IMA paper
%background points given
%output in chronological order
%start: historical data structure including n,t,type,T

%step1: background points
n_types=length(mu);
bp = struct('n',0,'t',[],'type',[],'father',[]);
for i=1:n_types
    n_type = poissrnd(mu(i)*T); 
    if n_type>0.5
        bp.n = bp.n + n_type;
        bp.t = [bp.t, sort(rand(1,n_type)*T)];
        bp.type = [bp.type,ones(1,n_type)*i];
    end
end
[bp.t,temp_ind] = sort(bp.t);
bp.type = bp.type(temp_ind);


bp.n = bp.n + start.n;
bp.father=zeros(1,bp.n);
bp.type=[start.type,bp.type];
bp.t=[start.t,start.T+bp.t];
flag=0;
if bp.n<0.5
    flag=2;
end
w=bp;

%step2: aftershocks
af_types=sum(A,2)';
while flag<1
    af=struct('n',[],'t',[],'type',[]);
    n2=poissrnd(af_types(w.type),1,w.n);
    af.n=sum(n2);
    af.type=[];
    af.t=[];
    af.father=[];

    for i=1:w.n
        if n2(i)>0.5  
            b1=datasample(mid_t,n2(i),'Weights',g)';
            af.t=[af.t,b1+w.t(i)];
            temp=rand(1,n2(i));
            temp_type=ones(1,n2(i));
            for j=1:n_types-1
                temp_type((temp>sum(A(w.type(i),1:j))/af_types(w.type(i)) & temp<=sum(A(w.type(i),1:j+1))/af_types(w.type(i))))=j+1;
            end
            af.type=[af.type,temp_type];
            af.father=[af.father,(i+bp.n-w.n)*ones(1,n2(i))];
        end
    end
    %combine
    if af.n>0.5
        temp_ind = af.t>start.T; %between start.T and start.T-dt?
        af.n = sum(temp_ind);
        af.t = af.t(temp_ind);
        af.type = af.type(temp_ind);
        af.father = af.father(temp_ind);
        
        bp.n=bp.n+af.n;
        bp.t=[bp.t,af.t];
        bp.type=[bp.type,af.type];
        bp.father=[bp.father,af.father];
        w=af;
        if min(af.t)>=T+start.T
            flag=2;
        end
    end
    
    if af.n<0.5 || af.n>10000
        flag=2;
    end     
end
y=struct('n',[],'t',[],'type',[]);
[y.t,ind]=sort(bp.t);
y.id=[1:1:bp.n];
y.id=y.id(ind);
y.n=bp.n;
y.t=y.t';
y.type=bp.type(ind)';
y.father=bp.father(ind);
% y.father(y.father>0)=ind(y.father(y.father>0));

end   
