function dCdx = steadyReactor(x,C,v0,k0_forward,k0_backward,Ea_forward, ...
    Ea_backward,Acs,rho_bulk,T)
% This function defines the material balance equations 
% the derivative terms are initialized

dCdx = zeros(4,1);
% the input concentrations for each component are defined 
CA = C(1); % mol/m^3
CB = C(2); % mol/m^3
CC = C(3); % mol/m^3
CD = C(4); % mol/m^3

Rm = 8.314; % universal gas constant (J/mol/K)

% rate constants are calculated using the Arrhenius equation
kin1 = k0_forward * exp(-Ea_forward/Rm/(T+273.15)); % mol g^-1 h^-1 atm^-2
kin2 = k0_backward * exp(-Ea_backward/Rm/(T+273.15));

% the rates are calculated using mass-action kinetics. Since rate constant 
% values are given for partial pressures, unit conversion is carried out.

rate(1) = kin1 * CA * CB * (8.205e-5 * (T+273.15))^2; % mol g^-1 h^-1
rate(2) = kin2 * CC * CD * (8.205e-5 * (T+273.15))^2;

% the derivative of concentrations calculated using mole balance equation

dCdx(1) = Acs / v0 * rho_bulk * 1000 / 3600 * (-rate(1) + rate(2)); 
dCdx(2) = Acs / v0 * rho_bulk * 1000 / 3600 * (-rate(1) + rate(2));
dCdx(3) = Acs / v0 * rho_bulk * 1000 / 3600 * (rate(1) - rate(2));
dCdx(4) = Acs / v0 * rho_bulk * 1000 / 3600 * (rate(1) - rate(2));

% hr to second unit by /3600
% kg to g conversion by *1000

end