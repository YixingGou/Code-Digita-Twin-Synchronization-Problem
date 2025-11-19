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
%% parameters of the production line:
%%%%%%%%%%%%%%%%%%the parameters shouble be changed for different lines%%%%%%%%%%%%%%%%%%%%%%%%
%Balanced lines: M_num=11, N_c=50, all-p=0.01,r=0.09
%Unbalanced lines: M_num=11, N_c=50, BN-p=0.03,r=0.09(position: Middle,Last); others-p=0.01,r=0.09
M_num=8;%the number of machines in the production line
% p=0.01*ones(1,M_num);%Failure probabilities
% r=0.09*ones(1,M_num);%Repair probabilities
% % p(ceil(M_num/2))=0.03;
% % p(M_num)=0.03;
% N_c=50*ones(1,M_num-1);%Buffer capacities: small(2), large(50)
% E=r./(p+r);%the efficiency of each machine
% M_state_0=ones(1,M_num);%the inital states of all machines
% WIP_0=zeros(1,M_num-1);%the inital buffer levels of all buffers
N=7;%the number of prediction periods
Delta=16;%Delta is the number of cycles in each period. N*Delta is the whole prediction interval. 
c_H=20;%c_H is the cost for each synchronization
c_B_1=20;%c_B_1 is the bias cost coefficient when the digital twin overestimates, which represents a holding cost of excess resources
c_B_2=20;%c_B_2 is the bias cost coefficient when the digital twin underestimates, which represents a penalty for production loss
L=1000;%L*Delta is the length of the unique sample path(real data)
SIM=100;%the number of experiments(the number of sample paths that used to eatimate the real cost)
REP_SIM=20;%the number of simulation replications for estimating R* and R_DT in each period

% %% generate SIM=100 fixed sample paths used to calculate the true cost
% for i=1:SIM
%     [M_state,TH_all,throughput,TH_M,TH_No_basic,BN,WIP]=simulation_sample_path(M_num,p,r,N_c,L,Delta,0,0,L*Delta,1,M_state_0,WIP_0);
%     filename = ['VAL_UBL2_', num2str(i), '.mat'];%the filename of SIM sample paths that are used to calculate the real cost
%     save(filename,'M_num','p','r','E','N_c','M_state','TH_all','throughput','TH_M','TH_No_basic','BN','WIP');
% end

y1=zeros(1,L);
r1=zeros();%r1 is used to record the performance at period n when the machine state is up(1)
r0=zeros();%r0 is used to record the performance at period n when the machine state is down(0)
C_B=zeros(1,N);
C_H=zeros(1,N);
C=zeros(1,N);
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
str = ['REAL_SP_dep_c_B_2=',num2str(c_B_2)];

%% read the optimal synchronization policies
%%%%%%%%%%%%%%%%%%%% the parameters to be changed: the sheet and its range related to ①-1 the production line:BL,UBL1,UBL2 and 2-① policy: state-dep, periodic %%%%%%%%%%%%%%%%
currentFolder = pwd;%the current folder
readFilePath = fullfile(currentFolder, 'R_optimization results_RealCase_fixedSP_realcost.xls');%the current file path
H1 = readmatrix(readFilePath, 'Sheet','change_c_B_SP','Range','E58:K62');
H0 = readmatrix(readFilePath, 'Sheet','change_c_B_SP','Range','L58:R62');
Num_Policy=size(H1,1);%Number of policies

Exp_cost=zeros(SIM,Num_Policy);
Exp_cost_test=zeros(SIM,Num_Policy);
Exp_cost_syn=zeros(SIM,Num_Policy);
Exp_cost_bias=zeros(SIM,Num_Policy);
Exp_NUMSyn=zeros(SIM,Num_Policy);
elapsedTimes = zeros(Num_Policy,1);% Preallocate an array to store the elapsed time of each iteration
for w=1:Num_Policy
    w
    iterationStartTime = tic;% Record the start time of the current iteration
    h1=H1(w,:);
    h0=H0(w,:);
    %% the machine being observed: Middle(ceil(M_num/2)) or Last(M_num)
%     if w<=Num_Policy/2 %the first 10 policies read from Excel were obtained by observing the middle machine,and the last 10 policies were obtained by observing the last machine
%         Position_machine=ceil(M_num/2); %observing the machine in the middle
%     else
%         Position_machine=M_num; %observing the last machine
%     end
    Position_machine=M_num; %observed machine
    %% For each given policy, we do SIM experiments (using SIM sample paths to estimate the real cost for each synchronization policy)
    E_N_predict_all = cell(SIM, 1);
    E_N_true_all = cell(SIM, 1);
    for i=1:SIM
        % i
        %% load the 100 fixed sample paths with L*Delta time units
        %%%%%%%%%%%%%%%%%%%% the parameters to be changed: ①-3 filename of the production line %%%%%%%%%%%%%%%%
        filename=['VAL_REAL_', num2str(i), '.mat'];
        load(filename);
        y=M_state';%the state of all the machines in the sample path
        wip=WIP';%the buffer level of all the buffers in the sample path
        y_observed=M_state(:,Position_machine)';%the state of the machine that is selected to observe
        %% estimate the average cost of L-N combinations
        for j=1:L-N
            %%the first predicted performance of the digital twin
            [~,~,~,TH_M_SIM,TH_No_basic_SIM,~,~]=simulation_sample_path(M_num,p,r,N_c,N+1,Delta,0,0,(N+1)*Delta,REP_SIM,y(:,(j-1)*Delta+2),wip(:,(j-1)*Delta+1));
            TH_first=TH_M_SIM;%the number of parts produced by the production line during [0,(N+1)*Delta) under REP_SIM simulation replications depending on the inital state
            th_last=mean(sum(TH_M_SIM,2));%the average number fo parts produced by the production line during [0,(N+1)*Delta) over REP_SIM simulation replications
            th_last_test=TH_No_basic_SIM(M_num);
            for n=1:N
                Y(j,n)=y_observed(1,(j+n-1)*Delta+2);%the state of observed machine at time n*Delta in the jth combination
                %% Estimation of true performance (best estimate with the most recent observation) during the interval [n*Delta,(n+1)*Delta) using simulation
                [~,~,~,TH_M_SIM,TH_No_basic_SIM,~,~]=simulation_sample_path(M_num,p,r,N_c,N+1,Delta,n*Delta,n*Delta,(n+1)*Delta,REP_SIM,y(:,(j+n-1)*Delta+2),wip(:,(j+n-1)*Delta+1));
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
                    [M_state_SIM,TH_SIM,throughput_SIM,TH_M_SIM,TH_No_basic_SIM,BN_SIM,WIP_SIM]=simulation_sample_path(M_num,p,r,N_c,N+1,Delta,n*Delta,n*Delta,(n+1)*Delta,REP_SIM,y(:,(j+n-1)*Delta+2),wip(:,(j+n-1)*Delta+1));
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
                else
                    c_B=c_B_2;
                end
                C_B(j,n)=c_B*(E_N_true(j,n)-E_N_predict(j,n))^2;%the bias cost depends on the difference between the true performance and predicted performance during[n*Delta,(n+1)*Delta)
                C_H(j,n)=c_H*h(j,n);
                C(j,n)=C_H(j,n)+C_B(j,n);%the total cost at time n*Delta
            end
            cost_syn(j)=sum(C_H(j,:));%the synchronization cost of the jth state combination for the whole periods(N periods)
            cost_bias(j)=sum(C_B(j,:));%the bias cost of the jth state combination for the whole periods(N periods)
            cost(j)=sum(C(j,:));%the cost of the jth state combination for the whole periods(N periods)
            NUMSyn(j)=sum(h(j,:));%the number of synchronization for the jth state combination
        end
        %% Output: Average cost of L-N combinations
        Exp_cost(i,w)=mean(cost);%expected total cost
        Exp_cost_syn(i,w)=mean(cost_syn);%expected synchronization cost
        Exp_cost_bias(i,w)=mean(cost_bias);%expected bias cost
        Exp_cost_test(i,w)=Exp_cost_syn(i,w)+Exp_cost_bias(i,w);
        Exp_NUMSyn(i,w)=mean(NUMSyn);
        %%save the intermediate variable as a Mat file
        E_N_predict_all{i} = E_N_predict;
        E_N_true_all{i} = E_N_true;
    end
    elapsedTimes(w) = toc(iterationStartTime);% Record the elapsed time of the current iteration
    save_file = ['intermediateResult_', str, '.mat'];
    if mod(w, 1) == 0
        save(save_file, 'Exp_cost', 'Exp_cost_syn', 'Exp_cost_bias', 'Exp_cost_test', 'Exp_NUMSyn', 'elapsedTimes');
    end
end
save(['E_N_all_', str, '.mat'], 'E_N_predict_all', 'E_N_true_all', '-v7');
AVG_cost=mean(Exp_cost)';%Average cost
AVG_cost_syn=mean(Exp_cost_syn)';
AVG_cost_bias=mean(Exp_cost_bias)';
AVG_cost_test=mean(Exp_cost_test)';
AVG_No_Syn=mean(Exp_NUMSyn)';%Average number of synchronization
overallElapsedTime = toc(overallStartTime);% Record the elapsed time of the entire program
savePath = fullfile(currentFolder, strcat(str, ".xls"));
writematrix(AVG_cost, savePath, 'Sheet',1,'Range','A2');
writematrix(AVG_cost_syn, savePath, 'Sheet',1,'Range','B2');
writematrix(AVG_cost_bias, savePath, 'Sheet',1,'Range','C2');
writematrix(AVG_No_Syn, savePath, 'Sheet',1,'Range','D2');
writematrix(AVG_cost_test, savePath, 'Sheet',1,'Range','F2');
writematrix(overallElapsedTime, savePath, 'Sheet',2,'Range','A2');
writematrix(elapsedTimes, savePath, 'Sheet',2,'Range','C2');
