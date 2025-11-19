%% State-dependent policy
%% Using enumeration approach to solve the simulation-based optimization model that minimize the expected cost
%%decision variables:h1 and h0 are the synchronization decisions at each period which is 1 if DT is synchronized and 0 otherwise
%%objective function:minimize the expected cost during the whole prediction periods
clear
clc
overallStartTime = tic;% Record the start time of the entire program

global M_num p r N_c N Delta c_H c_B_1 c_B_2 SIM L y y1 TH wip y_observed REP_SIM

%% parameters of the production line:
%%%%%%%%%%%%%%%% the parameters should be changed for different production lines %%%%%%%%%%%%%%%%
p=[0.0374, 0.087, 0.0473, 0.171, 0.121, 0.0278, 0.164, 0.655];%Failure probabilities
r=[0.145, 0.975, 0.522, 0.882, 0.88, 0.0826, 0.836, 0.256];%Repair probabilities
N_c=[26, 15, 28, 21, 21, 10, 8];%Buffer capacities
E=r./(p+r);%the efficiency of each machine
[~, BN] = min(E);%find the location of the bottleneck
M_num=size(p,2);%the number of machines in the production line
M_state_0=ones(1,M_num);%the inital states of all machines
WIP_0=zeros(1,M_num-1);%the inital buffer levels of all buffers
N=7;%the number of prediction periods
Delta=16;%Delta is the number of cycles in each period. N*Delta is the whole prediction interval. 
c_B_1=20;%c_B_1 is the bias cost coefficient when the digital twin overestimates, which represents a holding cost of excess resources
c_B_2=40;%c_B_2 is the bias cost coefficient when the digital twin underestimates, which represents a penalty for production loss
L=1000;%L*Delta is the length of the unique sample path(real data)

SIM=1;%the number of experiments for one sample path
REP_SIM=100;%the number of simulation replications when we estimate the true performance by using simulation
y1=zeros(2,L);%to record the state of machine at period l(l=0,1,…,L-1)


%% c_H is the parameter to be selected from set{10,50,150,2000}
c_H=20;%c_H is the cost for each synchronization
str='state_dep_c_H=20';

%% generate one sample paths and randomly replicate REP times
REP=5;
%%record the results of REP optimizations for REP different sample paths
min_Cost=zeros(REP,1);
AVG_No_Syn=zeros(REP,1);
Syn_Policy=zeros(REP,4*N);%record the synchronization policies for REP different sample paths
elapsed_time=zeros(REP,1);
for rep=1:REP
    rep
    iterationStartTime = tic;% Record the start time of the current iteration
    %%%%%%%%%%%%%%%%%%%% the parameters to be changed: ①-2 filename of the production line %%%%%%%%%%%%%%%%
    filename=['OPT_REAL_', num2str(rep), '.mat'];
    load(filename);
    y=M_state';%the state of all the machines
    wip=WIP';%the buffer level of all the buffers
    y_observed=M_state(:,[BN,M_num])';%the state of the bottleneck machine and last machine is selected to observe
    TH=TH_all(:,M_num)';%the throughput of the production line is equal to the the number of parts produced by the last machine
    for l=1:L
        y1(:,l)=y_observed(:,(l-1)*Delta+2);%the state of the machine being observed at period l(l=0,1,…,L-1)(time l*Delta+1)
    end

    %% Optimization: surrogateoptopt
    %the synchronization policy:sol
    %the corresponding minimum average cost:fval

    lb=zeros(1,4*N);%lower bounds for the decision variables h
    ub=ones(1,4*N);%upper bounds for the decision variables h
    fun=@objective_function_state_dep_new;%the objective function
    intcon=1:4*N;%Indices of variables that must be integer
    % options = optimoptions('surrogateopt','MaxFunctionEvaluations',120);
    options.MaxFunctionEvaluations=100;
    [sol,fval]=surrogateopt(fun,lb,ub,intcon,options);%searches for the global minimum of fun
    [~,AVG_SynNUM]=objective_function_state_dep_new(sol);
    % Syn_Policy=sol;
    Syn_Policy(rep,:)=sol;
    min_Cost(rep,:)=fval;
    AVG_No_Syn(rep,:)=AVG_SynNUM;
    
    elapsed_time(rep) = toc(iterationStartTime);% Record the elapsed time of the current iteration
    save_file = ['intermediateResult_', str, '.mat'];
    if mod(rep, 1) == 0
        save(save_file, 'min_Cost', 'AVG_No_Syn', 'Syn_Policy', 'elapsed_time');
    end
end

%% If the production line changes, the corresponding filename ("REAL") should be modified.
%%%%%%%%%%%%%%%%%%%% the parameters to be changed: the folder name related to ①-3 the production line: REAL,BL,UBL1,UBL2 %%%%%%%%%%%%%%%%
writematrix(min_Cost,".\dep_results_SIM\"+str+".xls",'Sheet',1,'Range','A2');
writematrix(AVG_No_Syn,".\dep_results_SIM\"+str+".xls",'Sheet',2,'Range','A2');
writematrix(Syn_Policy,".\dep_results_SIM\"+str+".xls",'Sheet',2,'Range','C2');
writematrix(elapsed_time,".\dep_results_SIM\"+str+".xls",'Sheet',3,'Range','A2');

overallElapsedTime = toc(overallStartTime);% Record the elapsed time of the entire program
fprintf('The entire program took %.6f seconds.\n', overallElapsedTime);