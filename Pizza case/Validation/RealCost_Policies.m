%%%state-dependent policy%%%
%% Estimation of the average cost for a given synchronization policy based only on the simulation
%%Both R* and R^DT are estimated using simulation based on the full observations
%% Input/decision variables:
%%h1 and h0 are the synchronization decisions depending not only on the period n but also on the system state（1 or 0）,which are 1 if DT is synchronized and 0 otherwise
%% Output/objective function:
%%AVG_cost and Avg_No_Syn are the average cost and the average number of synchronization during the whole prediction periods
clear
clc
overallStartTime = tic;% Record the start time of the entire program
%% parameters of the pizza line:
load('Parameters_PizzaLine_anyDist_new.mat','M_num','TBF_dists','TTR_dists','MTBF','MTTR','q','demand','assemble_all','N_c','E','BN')
M_state_0=ones(1,M_num);%the inital states of all machines
WIP_0=zeros(1,M_num-1);%the inital buffer levels of all buffers

N=7;%the number of prediction periods
Delta=8*60;%Delta is the number of cycles in each period. N*Delta is the whole prediction interval. 
c_B_1=1000;%c_B_1 is the bias cost coefficient when the digital twin overestimates, which represents a holding cost of excess resources
c_B_2=2000;%c_B_2 is the bias cost coefficient when the digital twin underestimates, which represents a penalty for production loss
L=1000;%L*Delta is the length of the unique sample path(real data)
N_val=10;%the number of experiments(the number of sample paths that used to eatimate the real cost)
REP_val=20;%the number of simulation replications for estimating R* and R_DT in each period

% %% generate N_val=100 fixed sample paths used to calculate the true cost
% for i=11:N_val
%     % [M_state,TH_all,throughput,TH_M,TH_No_basic,WIP]=simulation_PizzaLine_anyDist(M_num,TTR,TBF,q,demand,assemble_all,N_c,L,Delta,0,0,L*Delta,1,M_state_0,WIP_0);
%     [M_state,TH_all,throughput,TH_No_basic,WIP,TTR,TBF]=simulation_PizzaLine_anyDist(M_num,TTR_dists,TBF_dists,q,demand,assemble_all,N_c,L,Delta,0,0,L*Delta,1,M_state_0,WIP_0);
%     filename = ['VAL_PizzaLine_anyDist_', num2str(i), '.mat'];
%     % save(filename,'M_num','TBF','TTR','q','demand','assemble_all','E','N_c','BN','M_state','TH_all','throughput','TH_No_basic','WIP');
%     save(filename,'M_num','TBF_dists','TTR_dists','MTBF','MTTR','TBF','TTR','q','demand','assemble_all','E','N_c','BN','M_state','TH_all','throughput','TH_No_basic','WIP');
% end


y1=zeros(1,L);
r0=zeros();%r0 is used to record the performance at period n when the machine state is down(0)
C_B=zeros(L-N,N);
C_H=zeros(L-N,N);
C=zeros(L-N,N);
C_test=zeros(L-N,N);
cost=zeros(L-N,1);
cost_syn=zeros(L-N,1);
cost_bias=zeros(L-N,1);
NUMSyn=zeros(L-N,1);
Y=zeros(L-N,N);%Y is all the possible combinations(L-N+1 combinations) of system states over N periods in the unique sample path
E_N_true=zeros(L-N,N);
E_N_true_test=zeros(L-N,N);
E_N_predict=zeros(L-N,N);
E_N_predict_test=zeros(L-N,N);
h=zeros(L-N,N);%h is the synchronization decision which is 1 if DT is synchronized and 0 otherwise
%% the filename of recording the the real cost of the optimal policies
%%%%%%%%%%%%%%%%%%%% the parameters to be changed: the folder name related to ①-2 the production line:BL,UBL1,UBL2 and 2-② policy: state-dep, periodic %%%%%%%%%%%%%%%%
c_H=50;%c_H is the cost for each synchronization
str='Pizza_M1_never_c_H=50';

%% read the optimal synchronization policies
currentFolder = pwd;%the current folder0
readFilePath = fullfile(currentFolder, 'R_optimization results_PizzaLine-anyDist_new_M1.xls');%the current file path
H1 = readmatrix(readFilePath, 'Sheet','periodic','Range','E37:K39');
H0 = readmatrix(readFilePath, 'Sheet','periodic','Range','L37:R39');
Num_Policy=size(H1,1);%Number of policies

Exp_cost=zeros(N_val,Num_Policy);
Exp_cost_test=zeros(N_val,Num_Policy);
Exp_cost_syn=zeros(N_val,Num_Policy);
Exp_cost_bias=zeros(N_val,Num_Policy);
Exp_cost_bias_over=zeros(N_val,Num_Policy);
Exp_cost_bias_under=zeros(N_val,Num_Policy);
Exp_NUMSyn=zeros(N_val,Num_Policy);
elapsedTimes = zeros(Num_Policy,1);% Preallocate an array to store the elapsed time of each iteration
for w=1:Num_Policy
    w
    iterationStartTime = tic;% Record the start time of the current iteration
    h1=H1(w,:);
    h0=H0(w,:);
    %% the machine being observed: 
    Position_machine=1; %observed machine
    %% For each given policy, we do SIM experiments (using SIM sample paths to estimate the real cost for each synchronization policy)
    for i=1:N_val
        i
        %% load the 100 fixed sample paths with L*Delta time units
        filename=['VAL_PizzaLine_anyDist_', num2str(i), '.mat'];
        load(filename);
        y=M_state';%the state of all the machines in the sample path
        wip=WIP';%the buffer level of all the buffers in the sample path
        y_observed=M_state(:,Position_machine)';%the state of the machine that is selected to observe
        %% estimate the average cost of L-N combinations
        cost_bias_over=zeros(L-N,1);
        cost_bias_under=zeros(L-N,1);
        for j=1:L-N
            %%the first predicted performance of the digital twin
            [~,~,~,TH_M_SIM,TH_No_basic_SIM,~]=simulation_PizzaLine_anyDist(M_num,TTR_dists,TBF_dists,q,demand,assemble_all,N_c,N+1,Delta,0,0,(N+1)*Delta,REP_val,y(:,(j-1)*Delta+2),wip(:,(j-1)*Delta+1));
            TH_first=TH_M_SIM;%the number of parts produced by the production line during [0,(N+1)*Delta) under REP_SIM simulation replications depending on the inital state
            th_last=mean(sum(TH_M_SIM,2));%the average number fo parts produced by the production line during [0,(N+1)*Delta) over REP_SIM simulation replications
            th_last_test=TH_No_basic_SIM(M_num);
            C_B_1=zeros(L-N,N);
            C_B_2=zeros(L-N,N);
            for n=1:N
                Y(j,n)=y_observed(1,(j+n-1)*Delta+2);%the state of observed machine at time n*Delta in the jth combination
                %% Estimation of true performance (best estimate with the most recent observation) during the interval [n*Delta,(n+1)*Delta) using simulation
                [~,~,~,TH_M_SIM,TH_No_basic_SIM,~]=simulation_PizzaLine_anyDist(M_num,TTR_dists,TBF_dists,q,demand,assemble_all,N_c,N+1,Delta,n*Delta,n*Delta,(n+1)*Delta,REP_val,y(:,(j+n-1)*Delta+2),wip(:,(j+n-1)*Delta+1));
                TH_real=TH_M_SIM;%the number of parts produced by the production line during [n*Delta,(n+1)*Delta) under REP_SIM simulation replications
                E_N_true(j,n)=TH_No_basic_SIM(M_num);%the best-estimated throughput of the production line(average number of parts over REP_SIM simulation replications)
                E_N_true_test(j,n)=mean(sum(TH_real(:,n*Delta+2:(n+1)*Delta+1),2));
                if Y(j,n)==1
                    h(j,n)=h1(n);
                else
                    h(j,n)=h0(n);
                end
                %% Estimation of the predicted performance during the interval [n*Delta,(n+1)*Delta) using simulation
                if h(j,n)==1
                    [~,~,~,TH_M_SIM,TH_No_basic_SIM,~]=simulation_PizzaLine_anyDist(M_num,TTR_dists,TBF_dists,q,demand,assemble_all,N_c,N+1,Delta,n*Delta,n*Delta,(n+1)*Delta,REP_val,y(:,(j+n-1)*Delta+2),wip(:,(j+n-1)*Delta+1));
                    TH_syn=TH_M_SIM;%the number of parts produced by the production line under REP_SIM simulation replications
                    E_N_predict(j,n)=E_N_true(j,n);
                else
                    k=find(h(j,1:n),1,'last');%k is the period of last synchronization
                    if isempty(k)==1%if k is empty, it means that there is no synchronization during the interval[Delta,n*Delta).The predicted performance is estimated based on the initial state.
                        k=0;
                        No_M=sum(TH_first(:,n*Delta+2:(n+1)*Delta+1),2);
                        E_N_predict(j,n)=mean(No_M);
                    else%if k is not empty, the prediction performance is estimated based on the state of the last synchronization at time k*Delta
                        No_M=sum(TH_syn(:,n*Delta+2:(n+1)*Delta+1),2);
                        E_N_predict(j,n)=mean(No_M);
                    end
                end
                if E_N_predict(j,n)>=E_N_true(j,n)
                    c_B=c_B_1;
                    C_B_1(j,n)=c_B_1*(E_N_true(j,n)-E_N_predict(j,n))^2;
                else
                    c_B=c_B_2;
                    C_B_2(j,n)=c_B_2*(E_N_true(j,n)-E_N_predict(j,n))^2;
                end
                C_B(j,n)=c_B*(E_N_true(j,n)-E_N_predict(j,n))^2;%the bias cost depends on the difference between the true performance and predicted performance during[n*Delta,(n+1)*Delta)
                C_H(j,n)=c_H*h(j,n);
                C(j,n)=C_H(j,n)+C_B(j,n);%the total cost at time n*Delta
                C_test(j,n)=C_H(j,n)+C_B_1(j,n)+C_B_2(j,n);%the total cost at time n*Delta
            end
            cost_syn(j)=sum(C_H(j,:));%the synchronization cost of the jth state combination for the whole periods(N periods)
            cost_bias_over(j)=sum(C_B_1(j,:));%the bias cost of the jth state combination for the whole periods(N periods)
            cost_bias_under(j)=sum(C_B_2(j,:));%the bias cost of the jth state combination for the whole periods(N periods)
            cost_bias(j)=sum(C_B(j,:));%the bias cost of the jth state combination for the whole periods(N periods)           
            cost(j)=sum(C(j,:));%the cost of the jth state combination for the whole periods(N periods)
            NUMSyn(j)=sum(h(j,:));%the number of synchronization for the jth state combination
        end
        %% Output: Average cost of L-N combinations
        Exp_cost(i,w)=mean(cost);%expected total cost
        Exp_cost_syn(i,w)=mean(cost_syn);%expected synchronization cost
        Exp_cost_bias(i,w)=mean(cost_bias);%expected bias cost
        Exp_cost_bias_over(i,w)=mean(cost_bias_over);%expected bias cost for overestimation
        Exp_cost_bias_under(i,w)=mean(cost_bias_under);%expected bias cost for underestimation
        Exp_cost_test(i,w)=Exp_cost_syn(i,w)+Exp_cost_bias(i,w);
        Exp_NUMSyn(i,w)=mean(NUMSyn);
    end
    elapsedTimes(w) = toc(iterationStartTime);% Record the elapsed time of the current iteration
    save_file = ['intermediateResult_', str, '.mat'];
    if mod(w, 1) == 0
        save(save_file, 'Exp_cost', 'Exp_cost_syn', 'Exp_cost_bias', 'Exp_cost_bias_over', 'Exp_cost_bias_under', 'Exp_cost_test', 'Exp_NUMSyn', 'elapsedTimes');
    end
end
AVG_cost=mean(Exp_cost)';%Average cost
AVG_cost_syn=mean(Exp_cost_syn)';
AVG_cost_bias=mean(Exp_cost_bias)';
AVG_cost_bias_over=mean(Exp_cost_bias_over)';
AVG_cost_bias_under=mean(Exp_cost_bias_under)';
AVG_cost_test=mean(Exp_cost_test)';
AVG_No_Syn=mean(Exp_NUMSyn)';%Average number of synchronization
overallElapsedTime = toc(overallStartTime);% Record the elapsed time of the entire program
savePath = fullfile(currentFolder, 'Val_results', strcat(str, ".xls"));

writematrix(AVG_cost, savePath, 'Sheet',1,'Range','A2');
writematrix(AVG_cost_syn, savePath, 'Sheet',1,'Range','B2');
writematrix(AVG_cost_bias, savePath, 'Sheet',1,'Range','C2');
writematrix(AVG_cost_bias_over, savePath, 'Sheet',1,'Range','D2');
writematrix(AVG_cost_bias_under, savePath, 'Sheet',1,'Range','E2');
writematrix(AVG_No_Syn, savePath, 'Sheet',1,'Range','F2');
writematrix(AVG_cost_test, savePath, 'Sheet',1,'Range','H2');
writematrix(overallElapsedTime, savePath, 'Sheet',2,'Range','A2');
writematrix(elapsedTimes, savePath, 'Sheet',2,'Range','C2');
