function [rel_sens_matrix, percent_rel_sens] = ...
    calculateRelativeSensitivity(T, ghsv, ratio, parameters, delta_h, ...
    mw, Ts, rhos, dr_h, Lr, rho_bulk)
% This function calculates the relative sensitivity of output parameters
% with respect to changes in the kinetic parameter values

%% Kinetic parameters
k0_forward = parameters(1); % specified pre-exponential factor for forward reaction
k0_backward = parameters(2); % specified pre-exponential factor for backward reaction
Ea_forward = parameters(3); % specified activation energy of forward reaction (J/mol)
Ea_backward = Ea_forward + delta_h; % The calculated activation energy of
% backward reaction

%% Inlet calculations

Acs = pi / 4 * dr_h^2; % cross sectional area of the tube (m2)
Vr = Acs * Lr; % volume of the tube calculated from cross-sectional area 
% and length (m3)

% the index of the densities are determined using the first column of the
% density array and the input temperature. 
idx = find(Ts == T);
rho = rhos(idx,:); % the density array of component (kg / m3)

% mole fraction calculation for inlet components 
y_A = 1 / (1 + ratio); 
y_B = 1 - y_A;

% using volume of reactor and the gas-hour-space-velocity, volumetric flow
% rate is calculated. 
v0 = Vr * ghsv / 3600; % m3 / s

% molar volume of each component at given temperature.
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

% the inlet concentration vector defined for ode15s function
init(1) = c1; % concentration of A at entrance
init(2) = c2; % concentration of B at entrance
init(3) = 1e-30; % concentration of C at entrance (1e-30 used instead of 0)
init(4) = 1e-30; % concentration of D at entrance

%% Sensitivity evaluation 

xspan = linspace(0,Lr,50); % ODE is integrated through the length of reactor
% the output variables are calculated at every discrete point specified.
% 50 nodes is selected for the demo.

% The concentration profile of the reactor is calculated at parameter
% values 1% higher than the original parameter values. 
[~,C0] = ode15s(@(x,C) steadyReactor(x,C,v0,k0_forward,k0_backward,Ea_forward,Ea_backward,Acs,rho_bulk,T),xspan,init); 
[~,Cdelta1] = ode15s(@(x,C) steadyReactor(x,C,v0,k0_forward * 1.01,k0_backward,Ea_forward,Ea_backward,Acs,rho_bulk,T),xspan,init); 
[~,Cdelta2] = ode15s(@(x,C) steadyReactor(x,C,v0,k0_forward,k0_backward * 1.01,Ea_forward,Ea_backward,Acs,rho_bulk,T),xspan,init); 
[~,Cdelta3] = ode15s(@(x,C) steadyReactor(x,C,v0,k0_forward,k0_backward,Ea_forward * 1.01,Ea_forward * 1.01 + delta_h,Acs,rho_bulk,T),xspan,init); 


% The derivative of the concentration values with respect to small
% parameter changes are calculated using Taylor's approximation.
dCdk0f = (Cdelta1 - C0) / (k0_forward * 0.01);
dCdk0b = (Cdelta2 - C0) / (k0_backward * 0.01);
dCdEaf = (Cdelta3 - C0) / (Ea_forward * 0.01);

% The derivatives normalized using the inlet concentration and the original
% parameter values. 
norm_dCdk0f = dCdk0f * k0_forward ./ C0;
norm_dCdk0b = dCdk0b * k0_backward ./ C0;
norm_dCdEaf = dCdEaf * Ea_forward ./ C0;

% Relative sensitivity calculation
for i=1:size(norm_dCdk0f,2)% for each  species (for matrices dimension 2 is columns)
    RS_k0f_C(i) = sqrt(sum(norm_dCdk0f(:,i).^2)) / size(norm_dCdk0f(:,i),1);
    RS_k0b_C(i) = sqrt(sum(norm_dCdk0b(:,i).^2)) / size(norm_dCdk0b(:,i),1);
    RS_Eaf_C(i) = sqrt(sum(norm_dCdEaf(:,i).^2)) / size(norm_dCdEaf(:,i),1);
end

% The output parameters are prepared 
rel_sens_matrix = [RS_k0f_C; RS_k0b_C; RS_Eaf_C];
% The percent relative sensitivity is calculated on the basis of the total
% relative sensitivity of all variables must be add up to unity.
percent_rel_sens = rel_sens_matrix ./ sum(rel_sens_matrix) * 100;

end