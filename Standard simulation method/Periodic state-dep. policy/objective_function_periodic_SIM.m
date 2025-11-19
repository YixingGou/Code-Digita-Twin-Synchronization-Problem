%% Estimate the average cost from the sample path
function [AVG_cost,AVG_SynNUM]=objective_function_periodic_SIM(quadruple)
%% Input/decision variables:
%H is a row vector that includes h1 and h0. h1 and h0 are the synchronization decisions depend on the system state, which is 1 if DT is synchronized and 0 otherwise
%% Output/expected cost:
%%Exp_cost and Exp_NUMSyn are the expected cost and the expected number of synchronization for all the possible scenarios of system state,respectively
%% Given that when to synchronize, calculate the expected cost during the whole periods based on the sample path
global M_num p r N_c N Delta c_H c_B_1 c_B_2 SIM L y y1 TH wip y_observed REP_SIM
C_B=zeros(1,N);
C_H=zeros(1,N);
C=zeros(1,N);
cost=zeros(L-N,1);
cost_syn=zeros(L-N,1);
cost_bias=zeros(L-N,1);
NUMSyn=zeros(L-N,1);
Y=zeros(L-N,N);%Y is all the possible combinations(L-N+1 combinations) of system states over N periods in the unique sample path
Y0=zeros(L-N,1);%Y0 is the initial state vector at period 0 in the L-N+1 combinations of states 
E_N_true=zeros(L-N,N);
E_N_predict=zeros(L-N,N);
Exp_cost=zeros(SIM,1);
Exp_cost_syn=zeros(SIM,1);
Exp_cost_bias=zeros(SIM,1);
Exp_NUMSyn=zeros(SIM,1);
e_N=zeros();
%% given that when to synchronize during the whole period

%The periodic state-dependent synchronization policy with 
%the parameters (w1,m1,w0,m0) sets the synchronization decision H 
%depending on the observation of the machine state Yn being either one or zero 
%according to the following equation:
% quadruple=[0,0,0,2];%(w1,w0,m1,m0)
w=[quadruple(1);quadruple(3)];%Represent w by a column vector with two rows and one column, where the first row is w1 when the state is 1 and the second row is w0 when the state is 0. w1,w0 \in {0,1}
m=[quadruple(2);quadruple(4)];%Represent m by a column vector with two rows and one column, where the first row is m1 when the state is 1 and the second row is m0 when the state is 0. m1,m0 \in {0,1,2,…}
H=zeros(2,N);%the first row is the synchronization policy when the state is 1 and the second row is the policy when the state is 0
for x=1:2
    for n=1:N
        if m(x,1)==0
            H(x,n)=w(x,1);
        elseif m(x,1)>0 && mod(n,m(x,1))==0
            H(x,n)=w(x,1);
        elseif m(x,1)>0 && mod(n,m(x,1))~=0
            H(x,n)=1-w(x,1);
        end
    end
end
h1=H(1,:);
h0=H(2,:);
h=zeros(L-N,N);

%% For each case, we do SIM experiments
for i=1:SIM
    %% estimate the average cost of L-N combinations
    for j=1:L-N
        for n=1:N
            Y0(j)=y_observed(1,(j-1)*Delta+2);
            Y(j,n)=y_observed(1,(j+n-1)*Delta+2);%the state of system at time n*Delta in the jth combination
            %% Estimation of true performance (best estimate with the most recent observation) during the interval [n*Delta,(n+1)*Delta) from the simulation
            [~,~,~,TH_No_basic_SIM,~,~]=simulation_sample_path(M_num,p,r,N_c,1,Delta,0,0,Delta,REP_SIM,y(:,(j+n-1)*Delta+2),wip(:,(j+n-1)*Delta+1));
            E_N_true(j,n)=TH_No_basic_SIM(M_num);
            if Y(j,n)==1
                h(j,n)=h1(n);
            else
                h(j,n)=h0(n);
            end
            %% Estimation of the predicted performance during the interval [n*Delta,(n+1)*Delta) based on the sample path
            if h(j,n)==1
                x1=find(y1==Y(j,n));%select sample paths that the state of system at time n*Dleta is 1(a1 is a column vector)
                size_x1=size(x1,2);%the number of sample paths that the state of system at time n*Dleta is 1
                for x=1:size_x1
                    e_N(x)=sum(TH(1,(x1(x)-1)*Delta+2:(x1(x))*Delta+1));
                end
                E_N_predict(j,n)=sum(e_N)/size_x1;%the predicted performance estimated by using sample paths that the state at time n*Dleta is 1
                e_N=zeros();
            else
                k=find(h(j,1:n),1,'last');%k is the period of last synchronization
                if isempty(k)==1%if k is empty, it means that there is no synchronization during the interval[Delta,n*Delta).The predicted performance is estimated based on the initial state.
                    k=0;
                    x1=find(y1(1,1:L-(n-k))==Y0(j));
                    size_x1=size(x1,2);
                    for x=1:size_x1
                        e_N(x)=sum(TH(1,(x1(x)-1+n-k)*Delta+2:(x1(x)+n-k)*Delta+1));
                    end
                    E_N_predict(j,n)=sum(e_N)/size_x1;
                    e_N=zeros();                
                else%if k is not empty, the prediction performance is estimated based on the state of the last synchronization at time k*Delta
                    x1=find(y1(1,1:L-(n-k))==Y(j,k));
                    size_x1=size(x1,2);
                    for x=1:size_x1
                        e_N(x)=sum(TH(1,(x1(x)-1+n-k)*Delta+2:(x1(x)+n-k)*Delta+1));
                    end
                    E_N_predict(j,n)=sum(e_N)/size_x1;
                    e_N=zeros();                    
                end
            end
            if E_N_predict(j,n)>=E_N_true(j,n)
                c_B=c_B_1;
            else
                c_B=c_B_2;
            end
            C_B(j,n)=c_B*(E_N_true(j,n)-E_N_predict(j,n))^2;%the bias cost depends on the difference between the true performance and predicted performance during[n*Delta,(n+1)*Delta)
            C_H(j,n)=c_H*h(j,n);
            C(j,n)=C_H(j,n)+C_B(j,n);%the cost at time n*Delta
        end
        cost_syn(j)=sum(C_H(j,:));%the synchronization cost of jth scenario for the whole periods(N periods)
        cost_bias(j)=sum(C_B(j,:));%the bias cost of jth scenario for the whole periods(N periods)
        cost(j)=sum(C(j,:));%the cost of jth scenario for the whole periods(N periods)
        NUMSyn(j)=sum(h(j,:));%the number of synchronization for jth scenario
    end
    %% Output: Average cost of L-N combinations
    Exp_cost(i)=mean(cost);%expected total cost
    Exp_cost_syn(i)=mean(cost_syn);%expected synchronization cost
    Exp_cost_bias(i)=mean(cost_bias);%expected bias cost
%     Exp_cost_test=Exp_cost_syn(i)+Exp_cost_bias(i);
    Exp_NUMSyn(i)=mean(NUMSyn);
end
AVG_cost=mean(Exp_cost);%Average cost 
AVG_SynNUM=mean(Exp_NUMSyn);%Average number of synchronization
end