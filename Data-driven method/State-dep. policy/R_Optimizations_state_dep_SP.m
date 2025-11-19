%% State-dependent policy
%% observing the state of the bottleneck(BN)
%% Using the surrogateopt solver to solve the sample path-based optimization model that minimize the expected cost
%%generate one sample path to optimize the synchronization policy and replicate REP times randomly
%%decision variables:h1 and h0 are the synchronization decisions at each period which is 1 if DT is synchronized and 0 otherwise
%%objective function: minimize the average synchronization and bias costs during the whole prediction periods
clear
clc
overallStartTime = tic;% Record the start time of the entire program

global N Delta c_H c_B_1 c_B_2 SIM L R1 R0 y_observed y1 TH

%% parameters of the production line:
%%%%%%%%%%%%%%%% the parameters should be changed for different production lines %%%%%%%%%%%%%%%%
%Balanced lines: M_num=11, N_c=50, all-p=0.01,r=0.09
%Unbalanced lines: M_num=11, N_c=50, BN-p=0.03,r=0.09(position: Middle,Last); others-p=0.01,r=0.09
M_num=8;%the number of machines in the production line
% p=0.001*ones(1,M_num);%Failure probabilities
% r=0.009*ones(1,M_num);%Repair probabilities
% N_c=50*ones(1,M_num-1);%Buffer capacities: small(2), large(50)
% E=r./(p+r);%the efficiency of each machine
% M_state_0=ones(1,M_num);%the inital states of all machines
% WIP_0=zeros(1,M_num-1);%the inital buffer levels of all buffers
N=7;%the number of prediction periods
Delta=16;%Delta is the number of cycles in each period. N*Delta is the whole prediction interval. 
%c_H=20;%c_H is the cost of each synchronization
c_B_1=20;%c_B_1 is the bias cost coefficient when the digital twin overestimates
c_B_2=40;%c_B_2 is the bias cost coefficient when the digital twin underestimates
L=1000;%L*Delta is the length of the unique sample path(real data)

%% For optimization, the experiments are repeated REP times at each value of c_H (REP sample paths for optimization are fixed for different c_H)
REP=5;
% REP_SIM=100;%the number of simulation replications when we estimate R* using simulation(SP-based method does not use simulation)
SIM=1;%the repetitions for one sample path

% %generate REP fixed sample paths for optimization
% for i=1:REP
% %     [M_state,TH_all,throughput,TH_M,TH_No_basic,BN,WIP]=simulation_sample_path(M_num,p,r,N_c,L,Delta,0,0,L*Delta,1,M_state_0,WIP_0);
%     [M_state,TH_all,throughput,TH_No_basic,BN,WIP]=simulation_sample_path(M_num,p,r,N_c,L,Delta,0,0,L*Delta,1,M_state_0,WIP_0);
%     filename = ['OPT_UBL5M1-3_', num2str(i), '.mat'];
%     save(filename,'M_num','p','r','E','N_c','BN','M_state','TH_all','throughput','TH_No_basic','WIP');
% end


%% the machine being observed:
 %%%%%%%%%%%%%%%%%%%% the parameters to be changed for different observed machine %%%%%%%%%%%%%%%
OM=ceil(M_num);%OM /in {1,2,…,M}
for x2=1:size(OM,2)
    Position_machine=OM(x2);
    %% c_H is the parameter to be selected from set C_H=[1,50,250]
    C_H=20;% C_H=[1,10,50,250,350,450];
    for x1=1:size(C_H,2)
%         iterationStartTime = tic;% Record the start time of the current iteration
        c_H=C_H(x1);%c_H is the cost for each synchronization
    %     str=['state_dep_c_H=', num2str(c_H)]; 
        str = ['state_dep_M', num2str(Position_machine), '_c_H=',num2str(c_H)];
        %% REP fixed sample paths are used to optimize for different c_H
        %%Each sample path consists of L periods that includes the state of the machine at each cycle (L*Delta cycles)
        %%record the results of REP optimizations for REP different sample paths
        min_Cost=zeros(REP,1);
        AVG_No_Syn=zeros(REP,1);
        Syn_Policy=zeros(REP,2*N);%record the synchronization policies for REP different sample paths
        elapsed_time=zeros(REP,1);
        for rep=1:REP
            rep
            iterationStartTime = tic;% Record the start time of the current iteration
            %%%%%%%%%%%%%%%%%%%% the parameters to be changed: ①-2 filename of the production line %%%%%%%%%%%%%%%%
            filename=['OPT_REAL_', num2str(rep), '.mat'];
            load(filename);
            y=M_state';%the state of all the machines in the sample path
            wip=WIP';%the buffer level of all the buffers in the sample path
            y_observed=M_state(:,Position_machine)';%the state of the machine that is selected to observe
            TH=TH_all(:,M_num)';%the throughput of the production line is equal to the the number of parts produced by the last machine

            %% estimate the performance with most recent observation based on the real data
            y1=zeros(1,L);
            r1=zeros(1,L);%r1 is used to record the performance at period n when the machine state is up(1)
            r0=zeros(1,L);%r0 is used to record the performance at period n when the machine state is down(0)
            a=0;
            b=0;
            for l=1:L
                y1(l)=y_observed((l-1)*Delta+2);%the state of machine at period l(l=0,1,…,L-1)
                if y1(l)==1
                    a=a+1;
                    r1(l)=sum(TH(:,(l-1)*Delta+2:l*Delta+1));%the number of parts to be produced during the interval [n*Delta+1,(n+1)*Delta) when the machine is up
                else
                    b=b+1;
                    r0(l)=sum(TH(:,(l-1)*Delta+2:l*Delta+1));%the number of parts to be produced during the interval [n*Delta+1,(n+1)*Delta) when the machine is down
                end
            end
            R1=sum(r1)/a;%the expected performance of physical system when machine is up(1)
            R0=sum(r0)/b;%the expected performance of physical system when machine is down(0)

            %% Optimization: surrogateoptopt
            %the synchronization policy:sol
            %the corresponding minimum average cost:fval
            lb=zeros(1,2*N);%lower bounds for the decision variables h
            ub=ones(1,2*N);%upper bounds for the decision variables h
            fun=@objective_function_state_dep_SP;%the objective function
            intcon=1:2*N;%Indices of variables that must be integer
            options.MaxFunctionEvaluations=100;
            [sol,fval]=surrogateopt(fun,lb,ub,intcon,options);%searches for the global minimum of fun
            [~,AVG_SynNUM]=objective_function_state_dep_SP(sol);
            Syn_Policy(rep,:)=sol;
            min_Cost(rep,:)=fval;
            AVG_No_Syn(rep,:)=AVG_SynNUM;
            
            elapsed_time(rep) = toc(iterationStartTime);% Record the elapsed time of the current iteration
            save_file = ['intermediateResult_', str, '.mat'];% Record the intermediate results of each rep
            if mod(rep, 1) == 0
                save(save_file, 'min_Cost', 'AVG_No_Syn', 'Syn_Policy', 'elapsed_time');
            end
        end
        %% If the production line changes, the corresponding filename ("REAL") should be modified.
        %%%%%%%%%%%%%%%%%%%% the parameters to be changed: the folder name related to ①-3 the production line: REAL,BL,UBL1,UBL2 %%%%%%%%%%%%%%%%
        writematrix(min_Cost,".\dep_results_REAL_SP\"+str+".xls",'Sheet',1,'Range','A2');
        writematrix(AVG_No_Syn,".\dep_results_REAL_SP\"+str+".xls",'Sheet',2,'Range','A2');
        writematrix(Syn_Policy,".\dep_results_REAL_SP\"+str+".xls",'Sheet',2,'Range','C2');
        writematrix(elapsed_time,".\dep_results_REAL_SP\"+str+".xls",'Sheet',3,'Range','A2');
    end
end
overallElapsedTime = toc(overallStartTime);% Record the elapsed time of the entire program
fprintf('The entire program took %.6f seconds.\n', overallElapsedTime);