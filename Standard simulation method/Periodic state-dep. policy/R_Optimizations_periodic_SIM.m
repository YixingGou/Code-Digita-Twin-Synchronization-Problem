%% Periodic state-dependent policy
%% observing the state of the bottleneck(BN)
%% Using the surrogateopt solver to solve the sample path-based optimization model that minimize the expected total cost
%%generate one sample path to optimize the synchronization policy and replicate REP times randomly
%%decision variables: h1 and h0 are the synchronization decisions at each period which is 1 if DT is synchronized and 0 otherwise
%%objective function: minimize the average synchronization and bias costs during the whole prediction periods
clear
clc
overallStartTime = tic;% Record the start time of the entire program

global M_num p r N_c N Delta c_H c_B_1 c_B_2 SIM L y y1 TH wip y_observed REP_SIM 

%% parameters of the production line:
%%%%%%%%%%%%%%%% the parameters should be changed for different production lines %%%%%%%%%%%%%%%%
%Balanced lines: M_num=11, N_c=50, all-p=0.01,r=0.09
%Unbalanced lines: M_num=11, N_c=50, BN-p=0.03,r=0.09(position: Middle,Last); others-p=0.01,r=0.09
M_num=11;%the number of machines in the production line
p=0.01*ones(1,M_num);%Failure probabilities
r=0.09*ones(1,M_num);%Repair probabilities
p(ceil(M_num/2))=0.03;
N_c=50*ones(1,M_num-1);%Buffer capacities: small(2), large(50)
E=r./(p+r);%the efficiency of each machine
% M_state_0=ones(1,M_num);%the inital states of all machines
% WIP_0=zeros(1,M_num-1);%the inital buffer levels of all buffers
N=7;%the number of prediction periods
Delta=16;%Delta is the number of cycles in each period. N*Delta is the whole prediction interval. 
c_B_1=10;%c_B_1 is the bias cost coefficient when the digital twin overestimates, which represents a holding cost of excess resources
c_B_2=20;%c_B_2 is the bias cost coefficient when the digital twin underestimates, which represents a penalty for production loss
L=1000;%L*Delta is the length of the sample path(real data)

y1=zeros(1,L);%to record the state of machine at period l(l=0,1,…,L-1)

%% For optimization, the experiments are repeated REP times at each value of c_H (REP sample paths for optimization are fixed for different c_H)
REP=5;
REP_SIM=100;%the number of simulation replications when we estimate R* using simulation(SP-based method does not use simulation)
SIM=1;%the repetitions for one sample path

% %generate REP fixed sample paths for optimization
% for i=1:REP
% %     [M_state,TH_all,throughput,TH_M,TH_No_basic,BN,WIP]=simulation_sample_path(M_num,p,r,N_c,L,Delta,0,0,L*Delta,1,M_state_0,WIP_0);
%     [M_state,TH_all,throughput,TH_No_basic,BN,WIP]=simulation_sample_path(M_num,p,r,N_c,L,Delta,0,0,L*Delta,1,M_state_0,WIP_0);
%     filename = ['OPT_BL4_', num2str(i), '.mat'];
%     save(filename,'M_num','p','r','E','N_c','BN','M_state','TH_all','throughput','TH_No_basic','WIP');
% end

%% the machine being observed:
 %%%%%%%%%%%%%%%%%%%% the parameters to be changed for different observed machine %%%%%%%%%%%%%%%
OM=[ceil(M_num/2),M_num];
for x2=1:size(OM,2)
    Position_machine=OM(x2);
    %% c_H is the parameter to be selected from set C_H=[1,50,250]
    C_H=50;% C_H=[1,10,50,250,350,450];
    for x1=1:size(C_H,2)
        c_H=C_H(x1);%c_H is the cost for each synchronization
        str = ['periodic_M', num2str(Position_machine), '_c_H=',num2str(c_H)];
        %% REP fixed sample paths are used to optimize for different c_H
        %%Each sample path consists of L periods that includes the state of the machine at each cycle (L*Delta cycles)
        %%record the results of REP optimizations for REP different sample paths
        quadruple=zeros(REP,4);
        min_Cost=zeros(REP,1);
        AVG_No_Syn=zeros(REP,1);
        Syn_Policy=zeros(REP,2*N);%record the synchronization policies for REP different sample paths
        elapsed_time=zeros(REP,1);
        for rep=1:REP
            rep
            iterationStartTime = tic;% Record the start time of the current iteration
            %%%%%%%%%%%%%%%%%%%% the parameters to be changed: ①-2 filename of the production line %%%%%%%%%%%%%%%%
            filename=['OPT_UBL1_', num2str(rep), '.mat'];
            load(filename);
            y=M_state';%the state of all the machines in the sample path
            wip=WIP';%the buffer level of all the buffers in the sample path 
            y_observed=M_state(:,Position_machine)';%the state of the machine that is selected to observe
            TH=TH_all(:,M_num)';%TH is the throughput of the production line that is equal to the the number of parts produced by the last machine
            for l=1:L
                y1(l)=y_observed(1,(l-1)*Delta+2);%the state of the observed machine at period l(l=0,1,…,L-1)(time l*Delta+1)
            end
            %% Optimization: surrogateopt
            %the quadruple(w1,m1,w0,m0): sol
            %the corresponding minimum average cost: fval
            nvars=4;%the number of decision variables
            lb=[0,0,0,0];%lower bounds for the decision variables quadruple
            ub=[1,floor(N/2),1,floor(N/2)];%upper bounds for the decision variables quadruple
            intcon=1:nvars;
            fun=@objective_function_periodic_SIM;%the objective function
            options.MaxFunctionEvaluations=100;
            [sol,fval]=surrogateopt(fun,lb,ub,intcon,options);%searches for the global minimum of fun
            quadruple(rep,:)=sol;%(w1,m1,w0,m0)
            min_Cost(rep,:)=fval;
            [~,AVG_SynNUM]=objective_function_periodic_SIM(sol);
            AVG_No_Syn(rep,:)=AVG_SynNUM;
            %% Restore the periodic synchronization policy according to the quadruple(w1,m1,w0,m0)
            w=[quadruple(rep,1);quadruple(rep,3)];%Represent w by a column vector with two rows and one column, where the first row is w1 when the state is 1 and the second row is w0 when the state is 0. w1,w0 \in {0,1}
            m=[quadruple(rep,2);quadruple(rep,4)];%Represent m by a column vector with two rows and one column, where the first row is m1 when the state is 1 and the second row is m0 when the state is 0. m1,m0 \in {0,1,2,…}
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
            Syn_Policy(rep,:)=[H(1,:),H(2,:)];
            
            elapsed_time(rep) = toc(iterationStartTime);% Record the elapsed time of the current iteration
            save_file = ['intermediateResult_', str, '.mat'];
            if mod(rep, 1) == 0
                save(save_file, 'min_Cost', 'AVG_No_Syn', 'Syn_Policy', 'elapsed_time');
            end
        end

        %% If the production line changes, the corresponding filename ("UBL1") should be modified.
        %%%%%%%%%%%%%%%%%%%% the parameters to be changed: the folder name related to ①-3 the production line: REAL,BL,UBL1,UBL2 %%%%%%%%%%%%%%%%      
        writematrix(min_Cost,".\per_results_UBL1_SIM\"+str+".xls",'Sheet',1,'Range','A2');
        writematrix(quadruple,".\per_results_UBL1_SIM\"+str+".xls",'Sheet',1,'Range','C2');
        writematrix(AVG_No_Syn,".\per_results_UBL1_SIM\"+str+".xls",'Sheet',2,'Range','A2');
        writematrix(Syn_Policy,".\per_results_UBL1_SIM\"+str+".xls",'Sheet',2,'Range','C2');
        writematrix(elapsed_time,".\per_results_UBL1_SIM\"+str+".xls",'Sheet',3,'Range','A2');
    end
end
overallElapsedTime = toc(overallStartTime);% Record the elapsed time of the entire program
fprintf('The entire program took %.6f seconds.\n', overallElapsedTime);