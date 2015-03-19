% Vignesh Ramchandran
% Ion Channel -- Final Project
% Worked with Rajiv Deshpande and Kush Gupta

clear all
clc


%PART I (a)

%Voltage is starting at -20mV

%Define State Variables

V = -20; %mV
alpha = 0.835*(exp(0.027*(V - 35))); %rate exiting state 1
beta = 0.033*(exp(-0.093*(V - 35))); %rate exiting state 2


t = 0;                      %starting time

check_1 = 1;                %Checks if channel is in state 1
time = [];                  %array of transition times
state = [];                 %array of what states it transitions to (either 1 or 2)

number_in_1 = 0;            %Number of times in 1
number_in_2 = 0;            %Number of  times in 2

time_in_1 = 0;              %Time spent in 1
time_in_2 = 0;              %Time spent in 2


%Simulation


while t <= 20000
    
    if check_1 == 1                        %Channel is in state 1
        x = exprnd(1/alpha);               %Create a random exponential variate with parameter 1/alpha (1/lambda)
        
        time = [time; t];                  %note down time at transition
        state = [state; 1];                %note down state at transition
        
        number_in_1 = number_in_1 + 1;     %increase number of times in state 1 by 1
        time_in_1 = time_in_1 + x;         %increase time spent in state 1 by x
        t = t + x;                         %increase total time by x
        check_1 = 0;                       %switch states
        
    elseif check_1 == 0                    %Channel is in state 2
        x2 = exprnd(1/beta);               %Create a random exponential variate with parameter 1/alpha (1/lambda)
        
        time = [time; t];                  %note down time at transition
        state = [state; 2];                %note down time at transition
                
        number_in_2 = number_in_2 + 1;     %increase number of times in state 2 by 1
        time_in_2 = time_in_2 + x2;        %increase time in state 2 by x2
        t = t + x2;                        %increase total time by x2
        check_1 = 1;                       %switch states
        
    end
end



figure(1)                                               %Create Figure 1
stairs(time, state)                                     %Stairs aligns the two vectors and plots them as bars
xlim([0 300])                                           %Setting an x range
ylim([0 3])                                             %Setting a y range
xlabel('Time (ms)')                                     %Axis label for x
ylabel('Current State')                                 %Axis label for y
title('Single Channel Simulation at V = -20 mV')        %Title


fprintf('Average Dwell Time in 2 for -20mV')            
average_time_in_2 = time_in_2 / number_in_2             %Computes average time spent in state 2


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


%Now we change the voltage to 5 mV

%Define State Variables

V = 5; %mV
alpha = 0.835*(exp(0.027*(V - 35))); %rate exiting state 1
beta = 0.033*(exp(-0.093*(V - 35))); %rate exiting state 2


t = 0;                      %starting time

check_1 = 1;                %Checks if channel is in state 1
time = [];                  %array of transition times
state = [];                 %array of what states it transitions to (either 1 or 2)

number_in_1 = 0;            %Number of times in 1
number_in_2 = 0;            %Number of  times in 2

time_in_1 = 0;              %Time spent in 1
time_in_2 = 0;              %Time spent in 2


%Simulation


while t <= 20000
    
    if check_1 == 1                        %Channel is in state 1
        x = exprnd(1/alpha);               %Create a random exponential variate with parameter 1/alpha (1/lambda)
        
        time = [time; t];                  %note down time at transition
        state = [state; 1];                %note down state at transition
        
        number_in_1 = number_in_1 + 1;     %increase number of times in state 1 by 1
        time_in_1 = time_in_1 + x;         %increase time spent in state 1 by x
        t = t + x;                         %increase total time by x
        check_1 = 0;                       %switch states
        
    elseif check_1 == 0                    %Channel is in state 2
        x2 = exprnd(1/beta);               %Create a random exponential variate with parameter 1/alpha (1/lambda)
        
        time = [time; t];                  %note down time at transition
        state = [state; 2];                %note down time at transition
                
        number_in_2 = number_in_2 + 1;     %increase number of times in state 2 by 1
        time_in_2 = time_in_2 + x2;        %increase time in state 2 by x2
        t = t + x2;                        %increase total time by x2
        check_1 = 1;                       %switch states
        
    end
end



figure(2)                                               %Create Figure 1
stairs(time, state)                                     %Stairs aligns the two vectors and plots them as bars
xlim([0 300])                                           %Setting an x range
ylim([0 3])                                             %Setting a y range
xlabel('Time (ms)')                                     %Axis label for x
ylabel('Current State')                                 %Axis label for y
title('Single Channel Simulation at V = 5 mV')          %Title


fprintf('Average Dwell Time in 2 for 5 mV')            
average_time_in_2 = time_in_2 / number_in_2             %Computes average time spent in state 2


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


%Now we change the voltage to 40 mV

%Define State Variables

V = 40; %mV
alpha = 0.835*(exp(0.027*(V - 35))); %rate exiting state 1
beta = 0.033*(exp(-0.093*(V - 35))); %rate exiting state 2


t = 0;                      %starting time

check_1 = 1;                %Checks if channel is in state 1
time = [];                  %array of transition times
state = [];                 %array of what states it transitions to (either 1 or 2)

number_in_1 = 0;            %Number of times in 1
number_in_2 = 0;            %Number of  times in 2

time_in_1 = 0;              %Time spent in 1
time_in_2 = 0;              %Time spent in 2


%Simulation


while t <= 20000
    
    if check_1 == 1                        %Channel is in state 1
        x = exprnd(1/alpha);               %Create a random exponential variate with parameter 1/alpha (1/lambda)
        
        time = [time; t];                  %note down time at transition
        state = [state; 1];                %note down state at transition
        
        number_in_1 = number_in_1 + 1;     %increase number of times in state 1 by 1
        time_in_1 = time_in_1 + x;         %increase time spent in state 1 by x
        t = t + x;                         %increase total time by x
        check_1 = 0;                       %switch states
        
    elseif check_1 == 0                    %Channel is in state 2
        x2 = exprnd(1/beta);               %Create a random exponential variate with parameter 1/alpha (1/lambda)
        
        time = [time; t];                  %note down time at transition
        state = [state; 2];                %note down time at transition
                
        number_in_2 = number_in_2 + 1;     %increase number of times in state 2 by 1
        time_in_2 = time_in_2 + x2;        %increase time in state 2 by x2
        t = t + x2;                        %increase total time by x2
        check_1 = 1;                       %switch states
        
    end
end



figure(3)                                               %Create Figure 1
stairs(time, state)                                     %Stairs aligns the two vectors and plots them as bars
xlim([0 300])                                           %Setting an x range
ylim([0 3])                                             %Setting a y range
xlabel('Time (ms)')                                     %Axis label for x
ylabel('Current State')                                 %Axis label for y
title('Single Channel Simulation at V = 40 mV')         %Title


fprintf('Average Dwell Time in 2 for 40mV')            
average_time_in_2 = time_in_2 / number_in_2             %Computes average time spent in state 2



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Part I (b)


%SWITCH FROM -100mV to 5mV

%Define State Variables

V = -100; %mV
alpha = 0.835*(exp(0.027*(V - 35))); %rate exiting state 1
beta = 0.033*(exp(-0.093*(V - 35))); %rate exiting state 2


t = 0;                      %starting time

check_1 = 1;                %Checks if channel is in state 1
time = [];                  %array of transition times
state = [];                 %array of what states it transitions to (either 1 or 2)

number_in_1 = 0;            %Number of times in 1
number_in_2 = 0;            %Number of  times in 2

time_in_1 = 0;              %Time spent in 1
time_in_2 = 0;              %Time spent in 2


%Simulation


while t <= 200
    
    if t >= 100                               %at 100 ms -- switch voltage
        V = 5;
        alpha = 0.835*(exp(0.027*(V - 35)));  %re-define rate exiting state 1
        beta = 0.033*(exp(-0.093*(V - 35)));  %re-define rate exiting state 2
    end
    
    if check_1 == 1                        %Channel is in state 1
        x = exprnd(1/alpha);               %Create a random exponential variate with parameter 1/alpha (1/lambda)
        
        time = [time; t];                  %note down time at transition
        state = [state; 1];                %note down state at transition
        
        number_in_1 = number_in_1 + 1;     %increase number of times in state 1 by 1
        time_in_1 = time_in_1 + x;         %increase time spent in state 1 by x
        t = t + x;                         %increase total time by x
        check_1 = 0;                       %switch states
        
    elseif check_1 == 0                    %Channel is in state 2
        x2 = exprnd(1/beta);               %Create a random exponential variate with parameter 1/alpha (1/lambda)
        
        time = [time; t];                  %note down time at transition
        state = [state; 2];                %note down time at transition
                
        number_in_2 = number_in_2 + 1;     %increase number of times in state 2 by 1
        time_in_2 = time_in_2 + x2;        %increase time in state 2 by x2
        t = t + x2;                        %increase total time by x2
        check_1 = 1;                       %switch states
        
    end
end

figure(4)
stairs(time, state)                                     %Stairs aligns the two vectors and plots them as bars
xlim([0 200])                                           %Setting an x range
ylim([0 3])                                             %Setting a y range
xlabel('Time (ms)')                                     %Axis label for x
ylabel('Current State')                                 %Axis label for y
title('Single Channel Simulation: Switch to 5mV')       %Title


fprintf('Average Dwell Time in 2 for Switch to 5 mV')            
average_time_in_2 = time_in_2 / number_in_2             %Computes average time spent in state 2


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


%SWTICH FROM -100mV to 40mV

%Define State Variables

V = -100; %mV
alpha = 0.835*(exp(0.027*(V - 35))); %rate exiting state 1
beta = 0.033*(exp(-0.093*(V - 35))); %rate exiting state 2


t = 0;                      %starting time

check_1 = 1;                %Checks if channel is in state 1
time = [];                  %array of transition times
state = [];                 %array of what states it transitions to (either 1 or 2)

number_in_1 = 0;            %Number of times in 1
number_in_2 = 0;            %Number of  times in 2

time_in_1 = 0;              %Time spent in 1
time_in_2 = 0;              %Time spent in 2


%Simulation


while t <= 200
    
    if t >= 100                               %at 100 ms -- switch voltage
        V = 40;
        alpha = 0.835*(exp(0.027*(V - 35)));  %re-define rate exiting state 1
        beta = 0.033*(exp(-0.093*(V - 35)));  %re-define rate exiting state 2
    end
    
    if check_1 == 1                        %Channel is in state 1
        x = exprnd(1/alpha);               %Create a random exponential variate with parameter 1/alpha (1/lambda)
        
        time = [time; t];                  %note down time at transition
        state = [state; 1];                %note down state at transition
        
        number_in_1 = number_in_1 + 1;     %increase number of times in state 1 by 1
        time_in_1 = time_in_1 + x;         %increase time spent in state 1 by x
        t = t + x;                         %increase total time by x
        check_1 = 0;                       %switch states
        
    elseif check_1 == 0                    %Channel is in state 2
        x2 = exprnd(1/beta);               %Create a random exponential variate with parameter 1/alpha (1/lambda)
        
        time = [time; t];                  %note down time at transition
        state = [state; 2];                %note down time at transition
                
        number_in_2 = number_in_2 + 1;     %increase number of times in state 2 by 1
        time_in_2 = time_in_2 + x2;        %increase time in state 2 by x2
        t = t + x2;                        %increase total time by x2
        check_1 = 1;                       %switch states
        
    end
end

figure(5)
stairs(time, state)                                     %Stairs aligns the two vectors and plots them as bars
xlim([0 200])                                           %Setting an x range
ylim([0 3])                                             %Setting a y range
xlabel('Time (ms)')                                     %Axis label for x
ylabel('Current State')                                 %Axis label for y
title('Single Channel Simulation: Switch to 40mV')      %Title


fprintf('Average Dwell Time in 2 for Switch to 40 mV')            
average_time_in_2 = time_in_2 / number_in_2             %Computes average time spent in state 2


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%PART III (CHALLENGE PROBLEM)

[t, p] = ode45('ode100', [0 100], 0);
[t2, p_200] = ode45('ode200', [100 200], p(end));

p1 = p(end) * 100;
p2 = p_200(end) * 100;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%PART II (a)


%SWITCH FROM -100mV to 5mV
%Population of 50 channels



%Define State Variables

V = -100; %mV
alpha = 0.835*(exp(0.027*(V - 35))); %rate exiting state 1
beta = 0.033*(exp(-0.093*(V - 35))); %rate exiting state 2


t = 0;                      %starting time

check_1 = 1;                %Checks if channel is in state 1
time = [];                  %array of transition times
state = [];                 %array of what states it transitions to (either 1 or 2)

number_in_1 = 0;            %Number of times in 1
number_in_2 = 0;            %Number of  times in 2

time_in_1 = 0;              %Time spent in 1
time_in_2 = 0;              %Time spent in 2

Total_time = [];            %Collects synchronized time series data
Big_Matrix = [];
%Simulation



for i = 1:5
    while t <= 200


        if t >= 100                               %at 100 ms -- switch voltage
            V = 5;
            alpha = 0.835*(exp(0.027*(V - 35)));  %re-define rate exiting state 1
            beta = 0.033*(exp(-0.093*(V - 35)));  %re-define rate exiting state 2
        end

        if check_1 == 1                        %Channel is in state 1
            ax = 1/alpha;
            x = exprnd(ax);                    %Create a random exponential variate with parameter 1/alpha (1/lambda)

            time = [time; t];                  %note down time at transition
            state = [state; 1];                %note down state at transition

            number_in_1 = number_in_1 + 1;     %increase number of times in state 1 by 1
            time_in_1 = time_in_1 + x;         %increase time spent in state 1 by x
            t = t + x;                         %increase total time by x
            check_1 = 0;                       %switch states

        
        elseif check_1 == 0                    %Channel is in state 2
            bx = 1/beta;
            x2 = exprnd(bx);                   %Create a random exponential variate with parameter 1/alpha (1/lambda)

            time = [time; t];                  %note down time at transition
            state = [state; 2];                %note down time at transition

            number_in_2 = number_in_2 + 1;     %increase number of times in state 2 by 1
            time_in_2 = time_in_2 + x2;        %increase time in state 2 by x2
            t = t + x2;                        %increase total time by x2
            check_1 = 1;                       %switch states
        end

    end

Total_time = time_series(time)

time = [];
state = [];
t = 0;
check_1 = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







