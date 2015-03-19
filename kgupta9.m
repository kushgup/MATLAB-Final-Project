% Kush Gupta
% Final Project - Ion Channels
% Worked with Rajiv Deshpande and Vignesh Ramchandaran


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART Ia: Run a Simulation for 20000 milliseconds (20s)

i = 1;   %counter to index the figures containing graphs

for v = [-20 5 40]  %in mV
    
    alpha = .835 * exp (.027 * (v - 35) );	%equivalent to lambda for state 1
    beta = .033 * exp (-.093 * (v - 35) );	%equivalent to lambda for state 2

    % Initializing Variables

    t = 0;  %time
    state = 1;  %state
    transTimes = [];  %list of transition times
    listStates = []; %list of states
    numTime1 = 0;   %number of times in state 1
    totTime1 = 0;   %total time in state 1
    numTime2 = 0;   %number of times in state 2
    totTime2 = 0;   %total time in state 2

    while t < 20000
        if (state == 1)
            u = rand;
            dwell = -(1/alpha)*log(u);

            transTimes = [transTimes; t];
            listStates = [listStates; state];

            totTime1 = totTime1 + dwell;
            numTime1 = numTime1 + 1;

            t = t + dwell;
            state = 2;
        end

        if (state == 2)
            u = rand;
            dwell = -(1/beta)*log(u);

            transTimes = [transTimes; t];
            listStates = [listStates; state];

            totTime2 = totTime2 + dwell;
            numTime2 = numTime2 + 1;

            t = t + dwell;
            state = 1;
        end


    end

    figure(i)
    stairs(transTimes,listStates)   %stairs keeps the function at the current y-value until the next is called, creating 'stairs'.
    axis([0 200 0 3])
    %xlim([0 200])   %reduced time (x-axis) to 1/100th of the time for readability of graph
    %ylim([0 3])
    xlabel('Time in ms')
    ylabel('Markov State')

    title(sprintf('Single 2-State Model Simulation at %d mV', v))   %sprintf method helps create formatted text strings for the title
    
    fprintf('At %d volts:\n', v)
    fprintf('Theoretical mean occupancy time for State 2: %f\n', 1/beta)
    fprintf('Actual mean occupancy time for State 2: %f\n\n', totTime2/numTime2)
 
    i = i + 1;
    
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART Ib: Run a Simulation for 200 milliseconds (.2s) with voltage switch
% at 100 milliseconds


%SIMULATION #1: -100mV to 5mV

% Initializing Variables
t = 0;
state = 1;  %state
transTimes = [];  %list of transition times
listStates = []; %list of states
numTime1 = 0;   %number of times in state 1
totTime1 = 0;   %total time in state 1
numTime2 = 0;   %number of times in state 2
totTime2 = 0;   %total time in state 2

while t <= 200
    
    if t < 100
        v = -100;
    else
        v = 5;
    end
    
    alpha = .835 * exp (.027 * (v - 35) );	%equivalent to lambda for state 1
    beta = .033 * exp (-.093 * (v - 35) );	%equivalent to lambda for state 2

    if (state == 1)
        u = rand;
        dwell = -(1/alpha)*log(u);

        transTimes = [transTimes; t];
        listStates = [listStates; state];

        totTime1 = totTime1 + dwell;
        numTime1 = numTime1 + 1;

        t = t + dwell;
        state = 2;
    end

    if (state == 2)
        u = rand;
        dwell = -(1/beta)*log(u);

        transTimes = [transTimes; t];
        listStates = [listStates; state];

        totTime2 = totTime2 + dwell;
        numTime2 = numTime2 + 1;

        t = t + dwell;
        state = 1;
    end
    
end

figure(4)
stairs(transTimes,listStates)   %stairs keeps the function at the current y-value until the next is called, creating 'stairs'
axis([0 200 0 3])
%xlim([0 200])   %reduced time (x-axis) to 1/100th of the time for readability of graph
%ylim([0 3])
xlabel('Time in ms')
ylabel('Markov State')
title('Single 2-State Model Simulation with switch from -100 mV to 5 mV at 100ms')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SIMULATION #2: -100mV to 40mV

% Initializing Variables
t = 0;
state = 1;  %state
transTimes = [];  %list of transition times
listStates = []; %list of states
numTime1 = 0;   %number of times in state 1
totTime1 = 0;   %total time in state 1
numTime2 = 0;   %number of times in state 2
totTime2 = 0;   %total time in state 2

while t <= 200
    
    if t < 100
        v = -100;
    else
        v = 40;
    end
    
    alpha = .835 * exp (.027 * (v - 35) );	%equivalent to lambda for state 1
    beta = .033 * exp (-.093 * (v - 35) );	%equivalent to lambda for state 2

    if (state == 1)
        u = rand;
        dwell = -(1/alpha)*log(u);

        transTimes = [transTimes; t];
        listStates = [listStates; state];

        totTime1 = totTime1 + dwell;
        numTime1 = numTime1 + 1;

        t = t + dwell;
        state = 2;
    end

    if (state == 2)
        u = rand;
        dwell = -(1/beta)*log(u);

        transTimes = [transTimes; t];
        listStates = [listStates; state];

        totTime2 = totTime2 + dwell;
        numTime2 = numTime2 + 1;

        t = t + dwell;
        state = 1;
    end
    
end

figure(5)
stairs(transTimes,listStates)   %stairs keeps the function at the current y-value until the next is called, creating 'stairs'
axis([0 200 0 3])
%xlim([0 200])   %reduced time (x-axis) to 1/100th of the time for readability of graph
%ylim([0 3])
xlabel('Time in ms')
ylabel('Markov State')
title('Single 2-State Model Simulation with switch from -100 mV to 40 mV at 100ms')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHALLENGE PROBLEM

[t0, p2_0] = ode45('initialProb', [0 200] , 0); %function initialProb should initially have voltage at -100mV. At t=0, state is 1, so p2 = 0


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART IIa -- simulate 50 channels

% construct master matrix
a = 0:.5:200;
b = zeros(401,50);
master = [a' b];

for i = 2:51
    %run simulation:
    t = 0;
    state = 1;  %state
    transTimes = [];  %list of transition times
    listStates = []; %list of states
    while t <= 200

        if t < 100
            v = -100;
        else
            v = 5;
        end

        alpha = .835 * exp (.027 * (v - 35) );	%equivalent to lambda for state 1
        beta = .033 * exp (-.093 * (v - 35) );	%equivalent to lambda for state 2

        if (state == 1)
            u = rand;
            dwell = -(1/alpha)*log(u);

            transTimes = [transTimes; t];
            listStates = [listStates; state];

            t = t + dwell;
            state = 2;
        end

        if (state == 2)
            u = rand;
            dwell = -(1/beta)*log(u);

            transTimes = [transTimes; t];
            listStates = [listStates; state];

            t = t + dwell;
            state = 1;
        end

    end

    %get size of vectors
    len = length(transTimes);

    %fill in appropriate column in master with the state of the markov
    %model
    
    j = 1;
    for k = 1:401
        if j>len
            master(k,i) = listStates(j-1);
        elseif master(k,1) < transTimes(j)
            master(k,i) = listStates(j-1);
        elseif master(k,1) >= transTimes(j)
            master(k,i) = listStates(j);
            j = j+1;
        end
            
    end

end

% create state-2 occupancy probability vector
occProb2 = [];
for i = 1:401
    count = 0;
    for j = 2:51
        if master(i,j) == 2
            count = count + 1;
        end
    end
    occProb2 = [occProb2 (count/50)];
end

time = master(:,1);

figure(6)
plot(t0,p2_0,'r')
hold on
plot(time,occProb2,'b');
axis([0 200 0 .6])
xlabel('Time in ms')
ylabel('State 2 Occupancy Probability')
title('Theoretical vs. 50 Channel Estimate of Occupancy Probability')
legend('Theoretical Probability','Calculated Probability')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART IIb -- simulate 500 channels

% construct master matrix
a = 0:.5:200;
b = zeros(401,500);
master = [a' b];

for i = 2:501
    %run simulation:
    t = 0;
    state = 1;  %state
    transTimes = [];  %list of transition times
    listStates = []; %list of states
    while t <= 200

        if t < 100
            v = -100;
        else
            v = 5;
        end

        alpha = .835 * exp (.027 * (v - 35) );	%equivalent to lambda for state 1
        beta = .033 * exp (-.093 * (v - 35) );	%equivalent to lambda for state 2

        if (state == 1)
            u = rand;
            dwell = -(1/alpha)*log(u);

            transTimes = [transTimes; t];
            listStates = [listStates; state];

            t = t + dwell;
            state = 2;
        end

        if (state == 2)
            u = rand;
            dwell = -(1/beta)*log(u);

            transTimes = [transTimes; t];
            listStates = [listStates; state];

            t = t + dwell;
            state = 1;
        end

    end

    %get size of vectors
    len = length(transTimes);

    %fill in appropriate column in master with the state of the markov
    %model

    j = 1;
    for k = 1:401
        if j>len
            master(k,i) = listStates(j-1);
        elseif master(k,1) < transTimes(j)
            master(k,i) = listStates(j-1);
        elseif master(k,1) >= transTimes(j)
            master(k,i) = listStates(j);
            j = j+1;
        end
            
    end

end

% create state-2 occupancy probability vector
occProb2 = [];
for i = 1:401
    count = 0;
    for j = 2:501
        if master(i,j) == 2
            count = count + 1;
        end
    end
    occProb2 = [occProb2 (count/500)];
end

time = master(:,1);

figure(7)
plot(t0,p2_0,'r')
hold on
plot(time,occProb2,'b');
axis([0 200 0 .6])
xlabel('Time in ms')
ylabel('State 2 Occupancy Probability')
title('Theoretical vs. 500 Channel Estimate of Occupancy Probability')
legend('Theoretical Probability','Calculated Probability')
hold off

