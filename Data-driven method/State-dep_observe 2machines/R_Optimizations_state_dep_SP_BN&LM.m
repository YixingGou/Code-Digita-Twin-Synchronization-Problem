%% observing the state of both the bottleneck(BN) and the last machine
%% State-dependent policy
%% Using the surrogateopt solver to solve the sample path-based optimization model that minimize the expected cost
%%generate one sample path to optimize the synchronization policy and replicate REP times randomly
%%decision variables:h1, h2, h3 and h4 are the synchronization decisions depend on the system state(observe 2 machines-4 states), which is 1 if DT is synchronized and 0 otherwise
%%objective function: minimize the average synchronization and bias costs during the whole prediction periods
clear
clc
overallStartTime = tic;% Record the start time of the entire program

global N Delta c_H c_B_1 c_B_2 SIM L R1 R2 R3 R4 y y1 TH

%% parameters of the production line:
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
c_B_1=20;%c_B_1 is the bias cost coefficient when the digital twin overestimates
c_B_2=40;%c_B_2 is the bias cost coefficient when the digital twin underestimates
L=1000;%L*Delta is the length of the unique sample path(real data)
SIM=1;%the number of experiments for one sample path

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
    y=M_state(:,[BN,M_num])';%the state of both the bottleneck machine and last machine is selected to observe
    TH=TH_all(:,M_num)';%the throughput of the production line is equal to the the number of parts produced by the last machine

    %% estimate the performance with most recent observation based on the real data
    y1=zeros(2,L);%the state of machine at period l(l=0,1,…,L)
    r1=zeros(1,L);%r1 is used to record the performance at period n when the machine state is up(1)
    r2=zeros(1,L);%r0 is used to record the performance at period n when the machine state is down(0)
    r3=zeros(1,L);
    r4=zeros(1,L);
    a=0;
    b=0;
    c=0;
    d=0;
    for l=1:L
        y1(:,l)=y(:,(l-1)*Delta+2);%the state of machine at period l(l=0,1,…,L-1)
        if y1(1,l)==1 && y1(2,l)==1
            a=a+1;
            r1(l)=sum(TH(:,(l-1)*Delta+2:l*Delta+1));%the number of parts to be produced during the interval [n*Delta+1,(n+1)*Delta) when the machine is up
        elseif y1(1,l)==1 && y1(2,l)==0
            b=b+1;
            r2(l)=sum(TH(:,(l-1)*Delta+2:l*Delta+1));%the number of parts to be produced during the interval [n*Delta+1,(n+1)*Delta) when the machine is down
        elseif y1(1,l)==0 && y1(2,l)==1
            c=c+1;
            r3(l)=sum(TH(:,(l-1)*Delta+2:l*Delta+1));
        else
            d=d+1;
            r4(l)=sum(TH(:,(l-1)*Delta+2:l*Delta+1));        
        end
    end
    R1=sum(r1)/a;%the expected performance of physical system when machine is up(1)
    R2=sum(r2)/b;%the expected performance of physical system when machine is down(0)
    R3=sum(r3)/c;
    R4=sum(r4)/d;

    %% Optimization: surrogateoptopt
    %the synchronization policy:sol
    %the corresponding minimum average cost:fval
    lb=zeros(1,4*N);%lower bounds for the decision variables h
    ub=ones(1,4*N);%upper bounds for the decision variables h
    fun=@objective_function_state_dep_new;%the objective function
    intcon=1:4*N;%Indices of variables that must be integer
    options.MaxFunctionEvaluations=100;
    [sol,fval]=surrogateopt(fun,lb,ub,intcon,options);%searches for the global minimum of fun
    [~,AVG_SynNUM]=objective_function_state_dep_new(sol);
    Syn_Policy(rep,:)=sol;
    min_Cost(rep,:)=fval;
    AVG_No_Syn(rep,:)=AVG_SynNUM;
    
    elapsed_time(rep) = toc(iterationStartTime);% Record the elapsed time of the current iteration
    save_file = ['intermediateResult_', str, '.mat'];
    if mod(rep, 1) == 0
        save(save_file, 'min_Cost', 'AVG_No_Syn', 'Syn_Policy', 'elapsed_time');
    end
end

%% save the results in the Excel files:
writematrix(min_Cost,".\dep_results_SP\"+str+".xls",'Sheet',1,'Range','A2');
writematrix(AVG_No_Syn,".\dep_results_SP\"+str+".xls",'Sheet',2,'Range','A2');
writematrix(Syn_Policy,".\dep_results_SP\"+str+".xls",'Sheet',2,'Range','C2');
writematrix(elapsed_time,".\dep_results_SP\"+str+".xls",'Sheet',3,'Range','A2');

overallElapsedTime = toc(overallStartTime);% Record the elapsed time of the entire program
fprintf('The entire program took %.6f seconds.\n', overallElapsedTime);