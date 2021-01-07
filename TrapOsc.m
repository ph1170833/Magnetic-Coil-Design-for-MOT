function dydt=TrapOsc(t, y, I, Current)
    M = 2.207e-25; %Mass of Cesium atom
    %alpha = 3.684e-22;
    %kappa = 3.293e-19;
    %I = 1e3; %Intensity of trap laser
    %Current= 1; %Current through coils
    [alpha, kappa] = TrapPara(I, Current);
    % We have differential equation as:
    %    Mx'' = -alpha*x'-kappa*x
    % This can be broken down as two couple first order DEs:
    %   1. x' = v
    %   2. Mv' = -alpha*v-kappa*x
    % We will solve using these 2 DEs
    dydt = [y(2); -(alpha/M)*y(2)-(kappa/M)*y(1)];
end