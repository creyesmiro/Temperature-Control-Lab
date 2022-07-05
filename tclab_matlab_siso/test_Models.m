%% test heaters and temp read out
clc; clear all; close all; format shortg

% include tclab
tclab

%% adjustable constants
Ta = 20.4 + 273.15;       % K
U = 10.0;                 % W/m^2 K
alpha = 0.01;             % W / % heater

%% run and simulate heater
% simulation time
time = 10;           % time in min
loops = time*60;        % simulation time in seconds

% Predictions
% preallocate a t x t matrix with T1C() reading
% this matrix will be updated from the Arduino
% and used for supplying the initial T value to
% the RK function
T1 = ones(loops) * T1C();       % measured T (sensor 1)
Tp_fnd = ones(loops) * T1C();       % pred T (fundamental model)
Tp_fod = ones(loops) * T1C();

% Error & Power allocation
error_fnd = ones(loops);          % error matrix for fundamental model
error_fod = zeros(loops);         % error matrix for FOPDT model
Qlog = ones(loops); 

% get Temp from Arduino
start_time = clock;
prev_time = start_time;

% dynamic plot (note: subplots needs to be declared here first)
figure(1)
subplot(2,1,1)
hold on, grid on
anexp = animatedline('LineStyle','-', 'Color', 'g', 'LineWidth', 2);
anpredfnd = animatedline('LineStyle','--','Color', 'b','LineWidth', 2);
anpredfod = animatedline('LineStyle','--','Color', 'r', 'LineWidth', 2);
xlabel('time (sec)')
ylabel('Temperature \circC')
legend('T_1 Measured', 'Energy Bal Pred', 'FOPDT Pred', 'Location', 'northwest')
title('Temperature Simulation')
subplot(2,1,2)
hold on, grid on
yyaxis left
anerrorfnd = animatedline('LineStyle','-', 'Color', 'b', 'LineWidth', 2);
anerrorfod = animatedline('LineStyle','-', 'Color', 'r' , 'LineWidth', 2);
xlabel('Time (sec)')
ylabel('Cumulative Error')
yyaxis right
title('Step and Error Simulation')
anQ = animatedline('LineStyle','-', 'Color', 'k', 'LineWidth', 3);
ylabel('Power Level Q (%)')
xlabel('time (sec)')
legend('Energy Balance Error', 'FOPDT Error', 'Location', 'northwest')
title('Fundemental vs. FOPDT Error')
% set levels !!!!update FOPDT function manually !!!!!
level = ones(loops,1);
level(1:9) = 0;
level(10:300) = 90;
level(241:end) = 0;

for ii = 1:loops
    % adjust power level
    Q = level(ii);
    h1(Q);
    Qlog(ii) = level(ii);
    
    % Pause Sleep time
    pause_max = 1.0;
    pause_time = pause_max - etime(clock,prev_time);
    if pause_time >= 0.0
        pause(pause_time - 0.01)
    else
        pause(0.01)
    end
    
    % Record time and change in time
    t = clock;
    dt = etime(t,prev_time);
    prev_time = t;
    
    % non-linear energy balance
    jj = ii+1;
    [tsim, Tnext_fnd] = ode45(@(tsim,x)energy_bal(tsim,x,Q,...
        alpha,Ta,U), [0 dt], Tp_fnd(jj-1)+273.15);
    Tp_fnd(jj) = Tnext_fnd(end) - 273.15;
    
    % FOPDT model
    Tnext_fod = fopdt(Ta,dt,Tp_fod(jj-1),ii,loops,level);
    Tp_fod(jj) = Tnext_fod;
    
    % read and record from temperature controller
    T1(ii) = T1C();
    
    % calculate error
    error_fnd(jj) = error_fnd(ii) + abs(T1(ii) - Tp_fnd(ii));
    error_fod(jj) = error_fod(ii) + abs(T1(ii)-Tp_fod(ii));
    fnd_fod = error_fnd(end,1) - error_fod(end,1);
    
    % plot
    addpoints(anexp,ii,T1(ii))
    addpoints(anpredfnd,ii,Tp_fnd(ii))
    addpoints(anpredfod,ii,Tp_fod(ii))
    addpoints(anerrorfnd,ii,error_fnd(ii))
    addpoints(anerrorfod,ii,error_fod(ii))
    addpoints(anQ,ii,Qlog(ii))
    drawnow
    
end
disp(['Energy Balance Cumulative Error =', num2str(error_fnd(end,1))]);
disp(['FOPDT Cumulative Error =', num2str(error_fod(end,1))]);
disp(['Heat Transfer Coefficient (U) value =', num2str(U)]);
level = 0;
h1(level);
disp(['Heater 1 off: ', num2str(level)  '%'])
% turn off heater but keep led on if T > 50
if T1C() > 50
    led(1)
    disp(['Warning, heater temperature =', num2str(T1C())])
else
    led(0)
end