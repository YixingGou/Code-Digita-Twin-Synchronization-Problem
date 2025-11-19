function [M_state,TH,throughput,TH_No_basic,WIP,TTR,TBF]=simulation_PizzaLine_anyDist(M_num,TTR_dists,TBF_dists,q,demand,assemble_all,N_c,N,Delta,t0,t1,t2,REP_No,M_state_0,WIP_0)
%仿真实验所需参数
%%需要根据实验要求初始化的参数

% function [M_state,TH,throughput,TH_M,TH_No_basic,BN,WIP]=simulation_sample_path(M_num,p,r,N_c,N,Delta,t0,t1,t2,REP_No,M_state_0,WIP_0)

%% Input Parameters
% M_num = 6; % Number of machines
% N = 10; % Number of prediction periods
% Delta = 10; % Cycles in each prediction period
% REP_No = 1; % Number of repetitions of the simulation
% N_c = 20 * ones(1, M_num - 1); % Buffer capacities
% demand = 0.9; % Market demand rate
% q = [0.05, 0, 0.05, 0.05, 0.05, 0]; % Scrap rate of each machine
% r_pack=6;%Assembly ratio of pack to pizzas: Assemble r_pack pizzas into one pack
% r_box=3;%Assembly ratio of box to packs: Assemble r_box packe into one box
% assemble_all=[1,1,r_pack,r_box,1,1];%Assembly ratio in each machine
% t0 = 0; % Start time of the simulation
% t1 = 0; % Beginning of prediction interval
% t2 = N * Delta; % End of prediction interval
% M_state_0 = ones(1, M_num); % Initial states of all machines
% WIP_0 = zeros(1, M_num - 1); % Initial buffer levels of all buffers


% %% Define arbitrary distributions for TTR and TBF (machine-specific)
% TTR_dists = {
%     makedist('Exponential', 'mu', oven_mttr), ... % Machine 1
%     makedist('Exponential', 'mu', refrigerator_mttr), ...              % Machine 2
%     makedist('Uniform', 'Lower', packer_mttr_LB, 'Upper', packer_mttr_UB), ...     % Machine 3
%     makedist('Uniform', 'Lower', boxer_mttr_LB, 'Upper', boxer_mttr_UB), ... % Machine 4
%     makedist('Exponential', 'mu', robot_mttr), ...              % Machine 5
%     % makedist('Normal', 'mu', 4, 'sigma', 1)         % Machine 6
% };
% 
% TBF_dists = {
%     makedist('Exponential', 'mu', oven_mtbf), ...             % Machine 1
%     makedist('Exponential', 'mu', refrigerator_mtbf), ... % Machine 2
%     makedist('Uniform', 'Lower', packer_mtbf_LB, 'Upper', packer_mtbf_UB), ...    % Machine 3
%     makedist('Uniform', 'Lower', boxer_mtbf_LB, 'Upper', boxer_mtbf_UB) , ... % Machine 4
%     makedist('Exponential', 'mu', robot_mtbf), ...             % Machine 5
%     % makedist('Normal', 'mu', 10, 'sigma', 1)         % Machine 6
% };

%% Generate TTR and TBF for each machine
TTR = cell(1, M_num);
TBF = cell(1, M_num);

for m = 1:M_num-1
    TTR{m} = random(TTR_dists{m}, N * Delta, 1); % Generate TTR sequence
    TBF{m} = random(TBF_dists{m}, N * Delta, 1); % Generate TBF sequence
end
TBF{M_num}=ones(N * Delta, 1);
TTR{M_num}=zeros(N * Delta, 1);


%% Simulation Variables Initialization
TH_No = zeros(REP_No, M_num); % Throughput of each machine during the whole simulation period over REP_No replications
TH_M = zeros(REP_No, N * Delta + 1); % Throughput of the last machine during N*Delta + 1 time units
No_scrap = zeros(1, M_num); % Count of scrapped products

%% Simulation Loop
for k = 1:REP_No
    % State tracking
    M_state = ones(N * Delta + 1, M_num); % Machine states
    M_state(1, :) = M_state_0;
    WIP = zeros(N * Delta + 1, M_num - 1); % Buffer levels
    WIP(1, :) = WIP_0;
    TH = zeros(N * Delta + 1, M_num); % Throughput tracking
    
    % Time tracking for TTR and TBF
    time_in_current_state = zeros(1, M_num); % Time in the current state
    current_ttr = zeros(1, M_num); % Current TTR for each machine
    current_tbf = zeros(1, M_num); % Current TBF for each machine

    % Initialize TTR and TBF
    for m = 1:M_num
        current_tbf(m) = TBF{m}(1);
        current_ttr(m) = TTR{m}(1);
    end
    
    for t = t0 + 2:N * Delta + 1
        for m = 1:M_num
            % Update machine states based on TTR and TBF
            if M_state(t - 1, m) == 1 % Machine is operational
                if time_in_current_state(m) >= current_tbf(m) && current_ttr(m)>0
                    M_state(t, m) = 0; % Transition to failure
                    time_in_current_state(m) = 0;
                    current_ttr(m) = TTR{m}(mod(t, length(TTR{m})) + 1); % Update TTR
                else
                    time_in_current_state(m) = time_in_current_state(m) + 1;
                    M_state(t, m) = 1; % Continue operational
                end
            else % Machine is in failure
                if time_in_current_state(m) >= current_ttr(m)
                    M_state(t, m) = 1; % Transition to operational
                    time_in_current_state(m) = 0;
                    current_tbf(m) = TBF{m}(mod(t, length(TBF{m})) + 1); % Update TBF
                else
                    time_in_current_state(m) = time_in_current_state(m) + 1;
                    M_state(t, m) = 0; % Continue in failure
                end
            end
        end

        % Production logic for each machine
        for m=M_num:-1:1
            if m==M_num
                w1=(rand(1,1) < demand);%determine whether there is a market demand: 1-yes; 0-no
                if M_state(t,m)==1 && WIP(t-1,m-1)~=0 && w1==1
                    TH(t,m)=1;
                    WIP(t,m-1)=WIP(t-1,m-1)-1;
                else
                    TH(t,m)=0;
                    WIP(t,m-1)=WIP(t-1,m-1);
                end
            elseif m==1
                w2=(rand(1,1)>q(m));%determine whether the part produced by machine m is scrapped:0-scrap defective part,1-good quality
                if w2==0
                    No_scrap(m)=1;%count the number of scrapped pizzas
                end
                if M_state(t,m)==1 && WIP(t-1,m)~=N_c(1,m) && w2==1
                    TH(t,m)=1;
                    WIP(t,m)=WIP(t,m)+1;
                else
                    TH(t,m)=0;
                    WIP(t,m)=WIP(t,m);
                end
            else
                w3=(rand(1,1)>q(m));%determine whether the part produced by machine m is scrapped:0-scrap defective part,1-good quality
                if w3==0
                    No_scrap(m)=1;%count the number of scrapped pizzas
                end
                %%Machine m will prduce a pizza/pack/box only if it is up, the product is not scrapped, 
                % the upstream buffer level is more than assemble_all(m) and the downstream buffer is not full.
                if M_state(t,m)==1 && w3==1 && WIP(t-1,m)~=N_c(1,m) && WIP(t-1,m-1)>=assemble_all(m)
                    TH(t,m)=1;
                    WIP(t,m)=WIP(t,m)+1;
                    WIP(t,m-1)=WIP(t-1,m-1)-assemble_all(m);
                else
                    TH(t,m)=0;
                    WIP(t,m)=WIP(t,m);
                    WIP(t,m-1)=WIP(t-1,m-1);
                end
            end
        end

    end
    % Record throughput
    throughput = sum(TH(t1 + 2:t2 + 1, :) , 1);
    TH_No(k, :) = throughput;%the throughput of each machine during time interval [t1,t2) at kth simulation (each row: 1-row M_num-column vector)
    TH_M(k, :) = TH(:, M_num)';%the throughput of the last machine during N*Delta + 1 time units at kth simulation (each row: 1-row (N*Delta + 1)-column vector)
end
%% Output Results
TH_No_basic = sum(TH_No, 1) / REP_No; % Average throughput of each machine over REP_No replications

