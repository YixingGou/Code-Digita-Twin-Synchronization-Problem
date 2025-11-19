function [M_state,TH,throughput,TH_No_basic,BN,WIP]=simulation_sample_path(M_num,p,r,N_c,N,Delta,t0,t1,t2,REP_No,M_state_0,WIP_0)
%% Input:
% M_num is the number of machines in the geometric production line
% p is a vector of the failure probabilities of all machines
% r is a vector of the repair probabilities of all machines
% N_c is a vector of the buffer capacity
% N is the number of prediction periods
% Delta is the number of cycles in each prediction period
% t0 is the start time of the simulation
% t1 is the beginning of prediction interval[t1,t2)
% t2 is the end of prediction interval[t1,t2)
% REP_No is the number of repetitions of the simulation
% M_state_0 is the inital states of all machines at time t0
% WIP_0 is the inital buffer levels of all buffers
%% Output:
% rng(1) %For reproducibility
E=r./(p+r);%Machine efficiencies
[~, BN] = min(E);%find the location of the bottleneck
% Fs_No=zeros(REP_No,M_num);%to record the probability of starvation for each machine under Exp_No simulation repetitions
TH_No=zeros(REP_No,M_num);%to record the production rate of each machine under Exp_No simulation repetitions
for k=1:REP_No
    M_state=ones(N*Delta+1,M_num);%to record the states of all machines at each time unit
    M_state(1,:)=M_state_0;%the 1st row in array M_state is used to record the inital states
    M_state(2,:)=M_state(1,:);%the 2nd row in array M_state is the states at time 1, which is setted to be the same as the inital states
    WIP=zeros(N*Delta+1,M_num-1);%to record the buffer levels of all buffers at each time unit
    WIP(1,:)=WIP_0;%the 1st row in array M_state is used to record the inital buffer levels
    TH=zeros(N*Delta+1,M_num);%to record the number of parts produced by each machine at each time unit in a simulation(one of REP simulation repetitions)
    Fs=zeros(N*Delta+1,M_num-1);%to record whether each machine is starved at each time unit in a simulation(one of REP simulation repetitions)
    %%%%%%%%%%%%%%%%%%
    N_empty=zeros(1,M_num-1);
    M_oppo=ones(1,M_num);
    %%%%%%%%%%%%%%%%%%%%%Main of the simulation%%%%%%%%%%%%%%%%%
    %for t=2:T
    for t=t0+2:N*Delta+1
%         t
        for m=M_num:-1:1
            if m==M_num
                if M_state(t,m)==1 && WIP(t-1,m-1)~=0
                    TH(t,m)=1;
                    WIP(t,m-1)=WIP(t-1,m-1)-1;
                else
                    TH(t,m)=0;
                    WIP(t,m-1)=WIP(t-1,m-1);
                end
            elseif m==1
                if M_state(t,m)==1 && WIP(t-1,m)~=N_c(1,m)
                    TH(t,m)=1;
                    WIP(t,m)=WIP(t,m)+1;
                else
                    TH(t,m)=0;
                    WIP(t,m)=WIP(t,m);
                end
            else
                if M_state(t,m)==1 && WIP(t-1,m)~=N_c(1,m) && WIP(t-1,m-1)~=0
                    TH(t,m)=1;
                    WIP(t,m)=WIP(t,m)+1;
                    WIP(t,m-1)=WIP(t-1,m-1)-1;
                else
                    TH(t,m)=0;
                    WIP(t,m)=WIP(t,m);
                    WIP(t,m-1)=WIP(t-1,m-1);
                end
            end

        end
        a=rand(1,M_num);
        b1=(a>p);
        b2=(a<=r);
        b3=b1.*M_state(t,:);
        b4=b2.*(M_oppo-M_state(t,:));
        M_state(t+1,:)=b3+b4;
        c=(WIP(t,:)==N_empty);
        [~,d2]=find(c==1);
        e=size(d2,2);
        for i=1:e
            M_state(t+1,d2(1,i)+1)=1;
        end
        f=(WIP(t,:)==N_c);
        [~,g2]=find(f==1);
        h=size(g2,2);
        for j=1:h
            M_state(t+1,g2(1,j))=1;
        end
        for v=2:M_num
            if WIP(t,v-1)==1 && M_state(t+1,v-1)==0 && M_state(t+1,v)==1
                Fs(t,v)=1;
            else
                Fs(t,v)=0;
            end
        end
    end
    throughput=sum(TH(t1+2:t2+1,:));
    TH_No(k,:)=throughput;
    % E_No(k,:)=throughput/(t2-t1);
    % Fs_No(k,:)=sum(Fs(t1+2:t2+1,:))/(t2-t1); 
end
TH_No_basic=sum(TH_No,1)/REP_No;
% Fs_No_basic=sum(Fs_No,1)/REP_No;
