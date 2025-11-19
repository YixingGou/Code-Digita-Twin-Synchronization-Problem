%% Estimate the average cost by the digital twin with J sample paths 
function [AVG_cost,AVG_SynNUM]=objective_function_state_dep_new(H)
%% Input/decision variables:
%H is a row vector that includes h1, h2, h3 and h4, which are the synchronization decisions depend on the system state(observe 2 machines-4 states), which is 1 if DT is synchronized and 0 otherwise
%% Output/expected cost:
%%Exp_cost and Exp_NUMSyn are the expected cost and the expected number of synchronization for all the possible scenarios of system state,respectively
%% Given that when to synchronize, calculate the expected cost for all the possible scenarios during the whole periods
global N Delta c_H c_B_1 c_B_2 SIM L R1 R2 R3 R4 y y1 TH

C_B=zeros(1,N);
C=zeros(1,N);
cost=zeros(L-N,1);
NUMSyn=zeros(L-N,1);
Y_BN=zeros(L-N,N);%Y is all the possible combinations(L-N+1 combinations) of the states of bottleneck over N periods in the unique sample path
Y0_BN=zeros(L-N,1);%Y0 is the initial state vector at period 0 in the L-N+1 combinations of the states of bottleneck 
Y_LM=zeros(L-N,N);%Y is all the possible combinations(L-N+1 combinations) of the states of the last machine over N periods in the unique sample path
Y0_LM=zeros(L-N,1);%Y0 is the initial state vector at period 0 in the L-N+1 combinations of states of the last machine
Exp_cost=zeros(SIM,1);
Exp_NUMSyn=zeros(SIM,1);
e_N=zeros();
%% given that when to synchronize during the whole period
h1=H(1,1:N);
h2=H(1,N+1:2*N);
h3=H(1,2*N+1:3*N);
h4=H(1,3*N+1:4*N);
h=zeros(L-N,N);

%% For each case, we do SIM experiments
for i=1:SIM
    %% estimate the average cost of L-N combinations
    for j=1:L-N
        for n=1:N
            Y0_BN(j)=y(1,(j-1)*Delta+2);
            Y_BN(j,n)=y(1,(j+n-1)*Delta+2);%the state of bottleneck at time n*Delta in the jth combination
            Y0_LM(j)=y(2,(j-1)*Delta+2);
            Y_LM(j,n)=y(2,(j+n-1)*Delta+2);%the state of last machine at time n*Delta in the jth combination
            %% Estimation of the best-estimated performance of the physical system during the interval [n*Delta,(n+1)*Delta) based on the history of physical system
            if Y_BN(j,n)==1 && Y_LM(j,n)==1
                E_N_true=R1;
                h(j,n)=h1(n);
            elseif Y_BN(j,n)==1 && Y_LM(j,n)==0
                E_N_true=R2;
                h(j,n)=h2(n);               
            elseif Y_BN(j,n)==0 && Y_LM(j,n)==1
                E_N_true=R3;
                h(j,n)=h3(n);
            else
                E_N_true=R4;
                h(j,n)=h4(n);
            end
            if h(j,n)==1
                x1=find(y1(1,:)==Y_BN(j,n) & y1(2,:)==Y_LM(j,n));%select sample paths that the state of system at time n*Dleta is 1(a1 is a column vector)
                size_x1=size(x1,2);%the number of sample paths that the state of system at time n*Dleta is 1
                for x=1:size_x1
                    e_N(x)=sum(TH(1,(x1(x)-1)*Delta+2:(x1(x))*Delta+1));
                end
                E_N_predict=sum(e_N)/size_x1;%the predicted performance estimated by using sample paths that the state at time n*Dleta is 1
                e_N=zeros();
            else
                k=find(h(j,1:n),1,'last');%k is the period of last synchronization
                if isempty(k)==1%if k is empty, it means that there is no synchronization during the interval[Delta,n*Delta).The predicted performance is estimated based on the initial state.
                    k=0;
                    x1=find(y1(1,1:L-(n-k))==Y0_BN(j) & y1(2,1:L-(n-k))==Y0_LM(j));
                    size_x1=size(x1,2);
                    for x=1:size_x1
                        e_N(x)=sum(TH(1,(x1(x)-1+n-k)*Delta+2:(x1(x)+n-k)*Delta+1));
                    end
                    E_N_predict=sum(e_N)/size_x1;%the predicted performance is estimated with the average of the number of parts within the periods that the states are equal to 1 at the beginning of these periods 
                    e_N=zeros();
                else%if k is not empty, the prediction performance is estimated based on the state of the last synchronization at time k*Delta
                    x1=find(y1(1,1:L-(n-k))==Y_BN(j,k) & y1(2,1:L-(n-k))==Y_LM(j,k));
                    size_x1=size(x1,2);
                    for x=1:size_x1
                        e_N(x)=sum(TH(1,(x1(x)-1+n-k)*Delta+2:(x1(x)+n-k)*Delta+1));
                    end
                    E_N_predict=sum(e_N)/size_x1;
                    e_N=zeros();
                end
            end
            if E_N_predict>=E_N_true
                c_B=c_B_1;
            else
                c_B=c_B_2;
            end
            C_B(j,n)=c_B*(E_N_true-E_N_predict)^2;%the bias cost depends on the difference between the true performance and predicted performance during[n*Delta,(n+1)*Delta)
            C(j,n)=c_H*h(j,n)+C_B(j,n);%the cost at time n*Delta
        end
        cost(j)=sum(C(j,:));%the cost of jth scenario for the whole periods(N periods)
        NUMSyn(j)=sum(h(j,:));%the number of synchronization for jth scenario
    end
    %% Output: Average cost of L-N combinations
    Exp_cost(i)=mean(cost);
    Exp_NUMSyn(i)=mean(NUMSyn);
end
AVG_cost=mean(Exp_cost);%Average cost 
AVG_SynNUM=mean(Exp_NUMSyn);%Average number of synchronization
end