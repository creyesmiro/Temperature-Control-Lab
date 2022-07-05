function Tp_fod = fopdt(Tss,dt,Tp_prev,ii,loops,level)
%   FOPDT model
Kp = 0.7;               % degC/%  // proportional (gain)
tauP = 120.0;           % seconds // 
thetaP = 10;            % seconds (integer) // dead time
Qss = 0;                % heater initial power
TssC = Tss - 273.15;    % convert from K to deg C

% Simulate one time step with FOPDT model
z = exp(-dt/tauP);
Tp_fod = (Tp_prev - TssC) * z + level(max(1,ii-thetaP-1)-Qss)...
    *(1-z)*Kp + TssC;
end

