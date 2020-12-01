function [cA, cB, cC, cD] = runModel(T, ghsv, ratio, k0_forward, ...
    k0_backward, Ea_forward, Ea_backward, mw, rhos, Ts)
% This function calculates the exit concentrations and the conversion based
% on given input conditions and rate parameters

%% Reactor geometry, catalyst properties

Rm = 8.314; % universal gas constant (J/mol/K)
dr_h = 0.0127; % diameter of the reactor (m)
Acs = pi / 4 * dr_h^2; % cross sectional area of the reactor (m2)
Lr = 0.3048; % length of the reactor (m)
Vr = Acs * Lr; % volume of the reactor (m)
rho_bulk = 800; % bulk density of the catalyst (kg/m3)

%% Physical property
% the index of the densities are determined using the first column of the
% density array and the input temperature. 
idx = find(Ts == T);
rho = rhos(idx,:); % the density array of component (kg / m3)

%% Inlet concentration calculations

% mole fraction calculation for inlet components 
y_A = 1 / (1 + ratio); 
y_B = 1 - y_A; 

% using volume of reactor and the gas-hour-space-velocity, volumetric flow
% rate is calculated. 
v0 = Vr * ghsv / 3600; % m3 / s

% molar volume of each component
mv = mw ./ rho' / 1000; % m3 / mol

% total molar volume is calculated using mole fractions of inlet parameters
mv_total = y_A * mv(1) + y_B * mv(2); % m3 / mol

% inital molar flow estimated using the molar volume and the volumetric
% flow rate
n_In_total = v0 /  mv_total; % total molar flow (mol/s)

% total molar inlet flow is defined as a vector
n_in = [n_In_total*y_A; n_In_total*y_B; 0 ; 0]; % molar flow of components 
%(mol/s)

% concentrations are calculated using the inlet molar flows and the
% volumetric flow rate
c1 = n_in(1) / v0; % mol/m3
c2 = n_in(2) / v0; % mol/m3

%% Reactor simulation

% the inlet concentration vector defined for ode15s function

init(1) = c1; % concentration of A at entrance
init(2) = c2; % concentration of B at entrance
init(3) = 1e-30; % concentration of C at entrance (1e-30 used instead of 0)
init(4) = 1e-30; % concentration of D at entrance

xspan = [0 Lr]; % ODE is integrated through the length of reactor

% The ODE integration 
[x,C] = ode15s(@(x,C) steadyReactor(x,C,v0,k0_forward,k0_backward, ...
    Ea_forward,Ea_backward,Acs,rho_bulk,T),xspan,init); 

% After integration, output concentrations at the exit are returned 
cA = C(end,1); % mol/m3
cB = C(end,2); % mol/m3
cC = C(end,3); % mol/m3
cD = C(end,4); % mol/m3

end

