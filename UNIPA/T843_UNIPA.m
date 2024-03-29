%% MULTI-OBJECTIVE OPTIMIZATION FOR THE PLANNING OF A SINGLE BUILDING REC
% This script allows to optimize the synthesis, design and operation stages
% through a multi-carrier and multi-component Energy Hub model.
% The model allows managing six final demands (electricity, 
% high-temperature and low-temperature heating, cooling, water, and 
% hydrogen).
% The building is linked to six main grids (electricity, district heating, 
% district cooling, water, natural gas, and hydrogen).
% The model can be used for the simplified modeling of a generic energy 
% system (building, cluster of buildings, microgrids, etc.) and perfectly
% fits for the planning problem of a single-building with many apartments 
% making up a Renewable Energy Community (REC).
% The building can provide ancillary services to all the grids exploiting
% the flexibility of the demands.
% The final solution is considered robust since it is based on the 
% uncertainty related to the availability of renewable energy sources and
% to the final demands.
% The optimization is performed over a standard year according to five 
% objective functions (economic, primary energy, carbon emissions,
% grid interactions, critical raw materials) and is based on a 
% mixed-integer linear programming mathematical model.

clear 
clc


%% DEFINE AND EVALUATE SIMULATION PARAMETERS
% General Parameters
numScenarios = 10;                                                        % Number of Monte Carlo simulations
K_hour = 24;                                                               % Simulation hours in a day [hours/day]; can be 1, 2, 4, 6, 12 or 24
K_day = 1;                                                                 % Simulation days in a year [days/year]
T = K_hour*K_day;                                                          % Number of timesteps

% Multi-Objective Optimization Parameters
% Order of objective functions: 1) costs; 2) primary energy consumption;
% 3) carbon emissions; 4) grid interactions; 5) critical raw materials

% Multi-Objective Method
Method = 's';                                                                % 's' for scalarization technique, 'g' for genetic algorithm

% Multi-Objective Weighting factors (for scalarization technique)
w = [0.25 0.25 0.25 0.25 0];

% Multi-Objective Normalization factors (for scalarization technique)
n = [10^6 10^6 10^5 10^7 1];

% Components and grids inclusion (1 to include, 0 to exclude)
E_grid_in = 1;                                                             % Electrical grid
NG_grid_in = 1;                                                            % Natural gas network
H_grid_in = 1;                                                             % District heating network
F_grid_in = 1;                                                             % District cooling network
W_grid_in = 1;                                                             % District water network
H2_grid_in = 1;                                                            % District hydrogen network
E_grid_out = 1;                                                            % Electrical grid
H_grid_out = 1;                                                            % District heating network
F_grid_out = 1;                                                            % District cooling network
W_grid_out = 1;                                                            % District water network
H2_grid_out = 1;                                                           % District hydrogen network
GB = 1;                                                                    % Natural Gas-fired Boiler
CHP = 0;                                                                   % Combined Heat and Power
HP = 1;                                                                    % Electrical Heat Pump
EC = 1;                                                                    % Electrical Chiller (modeled as part of the electrical heat pump)
AC = 1;                                                                    % Absorption Chiller
EL = 1;                                                                    % Water Electrolyzer
FC = 1;                                                                    % Cogenerative Fuel Cell
WT = 1;                                                                    % Wind turbine
PV = 1;                                                                    % Photovoltaic plant
CSP = 1;                                                                   % Concentrating Solar Power
STC = 1;                                                                   % Solar thermal collector
HSS = 1;                                                                   % Hydrogen Storage System
WSS = 1;                                                                   % Water Storage System
ESS = 1;                                                                   % Electricity Storage System
TSS_H = 1;                                                                 % Heating Storage System
TSS_F = 1;                                                                 % Cooling Storage System

%% DEFINE THE REC DEMANDS
% Final loads uncertainty is assessed assuming a Gaussian distribution
% function
Trends = readtable('Data_Input.xlsx','sheet','Loads');

mu_E_dem =      Trends.mu_E_dem;
mu_H_HT_dem =   Trends.mu_H_HT_dem;
mu_H_LT_dem =   Trends.mu_H_LT_dem;
mu_F_dem =      Trends.mu_F_dem;
mu_H2_dem =     Trends.mu_H2_dem;
mu_W_dem =      Trends.mu_W_dem;

sigma_E_dem =       Trends.sigma_E_dem;
sigma_H_HT_dem =    Trends.sigma_H_HT_dem;
sigma_H_LT_dem =    Trends.sigma_H_LT_dem;
sigma_F_dem =       Trends.sigma_F_dem;
sigma_H2_dem =      Trends.sigma_H2_dem;
sigma_W_dem =       Trends.sigma_W_dem;

E_out_flex_TOT = 18;                                                       % Flexible Demand of Electricity [kWh_e]
H_out_HT_flex_TOT = 0;                                                     % Flexible Demand of High-temperature Heating [kWh_h]
H_out_LT_flex_TOT = (1.067*150+36.67)*1.162*(45-15)/1000*33;               % Flexible Demand of Low-temperature Heating [kWh_h]
F_out_flex_TOT = 0;                                                        % Flexible Demand of Cooling [kWh_f]
H2_out_flex_TOT = 5*15;                                                    % Flexible Demand of Hydrogen [kg]
W_out_flex_TOT = 0;                                                        % Flexible Demand of Water [kg]


%% DEFINE THE RENEWABLE ENERGY SOURCES AVAILABILITY
% Air temperature uncertainty is assessed assuming a Gaussian distribution
% function. Solar radiation uncertainty is assessed assuming a Beta
% distribution function. Wind speed uncertainty is assessed assuming a
% Weibull distribution function. 
Sol = readtable('Data_Input.xlsx','sheet','Solar');
G_0 = 1366;                                                                % Solar constant [W/m^2]

mu_T_air =      Sol.mu_T_air;
mu_G_b =        Sol.mu_G_b/G_0;
mu_G_d =        Sol.mu_G_d/G_0;
mu_G_r =        Sol.mu_G_r/G_0;

sigma_T_air =   Sol.sigma_T_air;
sigma_G_b =     Sol.sigma_G_b/G_0;
sigma_G_d =     Sol.sigma_G_d/G_0;
sigma_G_r =     Sol.sigma_G_r/G_0;

k_w = ncread("microclimate-Point_Torino.nc","weib_k_combined");                   % Weibull function shape factor at 50, 100, and 200 m
A_w = ncread("microclimate-Point_Torino.nc","weib_A_combined");                   % Weibull function scale factor at 50, 100, and 200 m


%% DEFINE OBJECTIVE FUNCTIONS PARAMETERS

% Order of energy vectors: E, DH, DC, H2, W
R_flex = -[0.001 0.001 0.001 0.001 0.001];                                 % Flexibility remuneration [euro/kWh or euro/kg]
R_REC_ele = -0.11;                                                         % REC electrical energy sharing remuneration
R_REC_th = -0.09;                                                          % REC thermal energy sharing remuneration

% Order of equipment: NGB, CHP, HP, AC, EL, FC, WT, PV, CSP, STC, HSS, WSS, ESS, TSS_H, TSS_F
Equip = readtable('Data_Input.xlsx','sheet','Equipment');
Equip = table2array(Equip(:,2:end));
Equip(isnan(Equip)) = 0;

% Objective functions parameters for equipment
i_def = 0.052;                                                             % Real interest rate
UL = Equip(1,:)';                               	                       % Useful life of equipment in years
C_capex = Equip(2,:)';                           	                       % First order investment costs for equipment [euro/size]
C_capex_0 = Equip(3,:)';                          	                       % Zero-th order investment costs for equipment [euro]
C_opex = Equip(4,:)';                             	                       % Operating costs for equipment [euro/size/year]
PE_embodied = Equip(5,:)';                             	                   % Cumulative Energy Demand for equipment production [MJ/unit]
GWP_embodied = Equip(6,:)';                             	               % Global Warming Potential for equipment production [kg CO2-eq/unit]
CRM_embodied = Equip(7,:)';                                                % Critical Raw Materials content for equipment [kg/unit]
UCRF = i_def*((1+i_def).^UL)./(((1+i_def).^UL)-1);	                       % Uniform Series Capital Recovery Factors

% Order of grids: E, NG, DH, DC, H2, W
Grids = readtable('Data_Input.xlsx','sheet','Grids');
Grids = table2array(Grids(:,2:end));
Grids(isnan(Grids)) = 0;

% Objective functions parameters for grids
C_grids =  ones(T,6).*Grids(1,:);                                          % Operating costs for purchases from grids [euro/kWh or euro/kg]
R_grids = -ones(T,6).*Grids(2,:);                                          % Operating incomes for purchases to grids [euro/kWh or euro/kg]
PE_grids = Grids(3,:)';                                                    % Cumulative Energy Demand [MJ/kWh or MJ/kg]
GWP_grids = Grids(4,:)';                                                   % Global Warming Potential [kg CO2-eq/kWh or kg CO2-eq/kg]
GII_grids_in = Grids(5,:)';                                                % Penalty factor [-]
GII_grids_out = Grids(6,:)';                                               % Penalty factor [-]


%% DEFINE TECHNICAL SCENARIO PARAMETERS

% Auxiliary parameters
Q_E = 10*(1+max(mu_E_dem)+E_out_flex_TOT);                                 % Upper bound for electrical energy flows
Q_H = 10*(1+max(mu_H_HT_dem)+max(mu_H_LT_dem)+H_out_HT_flex_TOT+H_out_LT_flex_TOT); % Upper bound for heating energy flows
Q_F = 10*(1+max(mu_F_dem)+F_out_flex_TOT);                                 % Upper bound for cooling energy flows
Q_H2 = 10*(1+max(mu_H2_dem)+H2_out_flex_TOT);                              % Upper bound for hydrogen mass flows
Q_W = 10*(1+max(mu_W_dem)+W_out_flex_TOT);                                 % Upper bound for water mass flows

% Grids and network parameters
K_E_grid = 1;                                                              % Electrical grid efficiency [-]
K_H_grid = 1;	                                                           % District heating network efficiency [-]
K_F_grid = 1;	                                                           % District cooling network efficiency [-]
K_W_grid = 1;                  	                                           % Water network efficiency [-]
K_H2_grid = 1;	                                                           % Hydrogen network efficiency [-]

% Converters parameters
K_GB_gh = 0.90;                                                            % Boiler gas-to-heating efficiency [-]

K_CHP_ge = 0.34;                                                           % Combined Heat and Power gas-to-electricity efficiency [-]
K_CHP_gh = 0.49;                                                           % Combined Heat and Power gas-to-heating efficiency [-]

K_HP_eh = 5.7;                                                             % Heat Pump electricity-to-heating coefficient of performance [-]
K_EC_ef = 4.1;                                                             % Chiller electricity-to-cooling coefficient of performance [-]

K_AC_hf = 0.9;                                                             % Absorption Chiller heating-to-cooling coefficient of performance [-]

K_EL_eH2 = 0.016;                                                          % Water Electrolyzer electricity-to-hydrogen conversion factor [kg/kWh_e]
K_EL_H2w = 8.55;                                                           % Water Electrolyzer hydrogen-to-water conversion factor [l/kg]

K_FC_H2h = 0.4;                                                            % Fuel Cell hydrogen-to-heating conversion factor [kg/kWh_th]
K_FC_H2e = 12.23;                                                          % Fuel Cell hydrogen-to-electricity conversion factor [kWh_e/kg]
K_FC_H2w = 9.47;                                                           % Fuel Cell hydrogen-to-water conversion factor [kg/kg]

H_GB_min = 0;                                                              % Boiler minimum partialization rate (with respect to the rated size) [-]
E_CHP_min = 0;                                                             % Combined Heat and Power minimum partialization rate (with respect to the rated size) [-]
F_HP_min = 0;                                                              % Heat Pump (cooling mode) minimum partialization rate (with respect to the rated size) [-]
H_HP_min= 0;                                                               % Heat Pump (heating mode) minimum partialization rate (with respect to the rated size) [-]
F_AC_min = 0;                                                              % Absorption Chiller minimum partialization rate (with respect to the rated size) [-]
E_EL_min = 0;                                                              % Water Electrolyzer minimum partialization rate (with respect to the rated size) [-]
E_FC_min = 0;                                                              % Fuel Cell  minimum partialization rate (with respect to the rated size) [-]
RU_GB = 1;                                                                 % Boiler ramp-up rate (with respect to the rated size) [-]
RU_CHP = 1;                                                                % Combined Heat and Power ramp-up rate (with respect to the rated size) [-]
RU_EC = 1;                                                                 % Heat Pump ramp-up rate (with respect to the rated size) [-]
RU_AC = 1;                                                                 % Absorption Chiller ramp-up rate (with respect to the rated size) [-]
RU_EL = 1;                                                                 % Water Electrolyzer ramp-up rate (with respect to the rated size) [-]
RU_FC = 1;                                                                 % Fuel Cell ramp-up rate (with respect to the rated size) [-]
RD_GB = 1;                                                                 % Boiler ramp-down rate (with respect to the rated size) [-]
RD_CHP = 1;                                                                % Combined Heat and Power ramp-down rate (with respect to the rated size) [-]
RD_EC = 1;                                                                 % Heat Pump ramp-down rate (with respect to the rated size) [-]
RD_AC = 1;                                                                 % Absorption Chiller ramp-down rate (with respect to the rated size) [-]
RD_EL = 1;                                                                 % Water Electrolyzer ramp-down rate (with respect to the rated size) [-]
RD_FC = 1;                                                                 % Fuel Cell ramp-down rate (with respect to the rated size) [-]

% Generators parameters
K_T_w = 6;                                                                 % number of wind measures per hour
h_WT = 80;                                                                 % wind turbine hub height [m]
h_ref = 100;                                                               % anemometer height [m] - can be 50, 100, or 200 m
v_ci = 4;                                                                  % wind turbine cut-in wind speed [m/s]
v_co = 22;                                                                 % wind turbine cut-off wind speed [m/s]
v_r = 10;                                                                  % wind turbine rated wind speed [m/s]
P_r_WT = 400;                                                              % wind turbine rated power [kW]

NOCT = 45;                                                                 % PV nominal operating cell temperature [°C]
G_NOCT = 800;                                                              % reference solar radiance [W/m^2]
T_a_NOCT = 20;                                                             % reference air temperature [°C]
A_PV = 1.722 * 1.134;                                                      % One PV module surface [m^2]
K_PV_BOP = 0.95;                                                           % PV system Balance Of Plant efficiency [-]
K_PV_se = 0.21;                                                            % PV system solar-to-electricity efficiency  [-]
T_ref = 25;                                                                % reference cell temperature [°C]
beta_PV = -2.75*10^-3 ;                                                    % PV cell temperature coefficient [°C^-1]
P_r_PV = A_PV*K_PV_se;                                                     % PV system rated power [kW]

K_CSP_se = 0.1394;                                                         % Concentrating Solar Power solar-to-electricity efficiency  [-]
K_CSP_sh = 0.3964;                                                         % Concentrating Solar Power solar-to-heating efficiency  [-]
A_CSP = 400;                                                               % Concentrating Solar Power surface [m^2]
P_r_CSP = 1000;                                                            % CSP system rated power [kW]

eta_0 = 0.734;                                                             % Solar Thermal Collector zero-loss effificency [-]
a_1 = 1.529;                                                               % Solar Thermal Collector first order heat loss coefficient [W/m^2 K]
a_2 = 0.0166;                                                              % Solar Thermal Collector second order heat loss coefficient [W/m^2 K^2]
A_STC = 1.867;                                                             % Solar Thermal Collector aperture area [m^2]
T_m = 40;                                                                  % Mean temperature of heat transfer fluid [°C]

% Storage systems parameters
K_HSS_in = 1;                                                              % Hydrogen storage charge efficiency [-]
K_HSS_out = 1;                                                             % Hydrogen storage discharge efficiency [-]
phi_HSS = 0.02;                                                            % Hydrogen storage self-discharge rate [-] (articolo Optimizing Renewable Power Management in Transmission Congestion. An Energy Hub Model Using Hydrogen Storage)
DoD_HSS = 0;                                                               % Hydrogen storage Depth Of Discharge [-]
DoC_HSS = 1;                                                               % Hydrogen storage Depth Of Charge [-]

K_WSS_in = 0.99;                                                           % Water storage charge efficiency [-]
K_WSS_out = 0.99;                                                          % Water storage discharge efficiency [-]
phi_WSS = 0.01;                                                            % Water storage self-discharge rate [-]
DoD_WSS = 0;                                                               % Water storage Depth Of Discharge [-]
DoC_WSS = 1;                                                               % Water storage Depth Of Charge [-]

K_ESS_in = 0.97;                                                           % Electrical storage charge efficiency [-] (articolo Francesco Urban energy hubs economic optimization and environmental comparison in Italy and Vietnam)
K_ESS_out = 0.97;                                                          % Electrical storage discharge efficiency [-]
phi_ESS = 0.01;                                                            % Electrical storage self-discharge rate [-]
DoD_EES = 0.2;                                                             % Electrical storage Depth Of Discharge [-] (articolo Francesco Urban energy hubs economic optimization and environmental comparison in Italy and Vietnam)
DoC_EES = 1;                                                               % Electrical storage Depth Of Charge [-]

K_TSS_H_in = 1;                                                            % Thermal storage (heating) charge efficiency [-] (articolo Multi-Objective Optimization of Urban Microgrid Energy Supply According to Economic and Environmental Criteria) 
K_TSS_H_out = 1;                                                           % Thermal storage (heating) discharge efficiency [-]
K_TSS_F_in = 1;                                                            % Thermal storage (cooling) charge efficiency [-]
K_TSS_F_out = 1;                                                           % Thermal storage (cooling) discharge efficiency [-]
phi_TSS_H = 0.01;                                                          % Thermal storage (heating) self-discharge rate [-] (articolo Francesco Multi-Objective Optimization of Urban Microgrid Energy Supply According to Economic and Environmental Criteria)
phi_TSS_F = 0.01;                                                          % Thermal storage (cooling) self-discharge rate [-]
DoD_TSS_H = 0;                                                             % Thermal storage (heating) Depth Of Discharge [-] 
DoD_TSS_F = 0;                                                             % Thermal storage (cooling) Depth Of Discharge [-] 
DoC_TSS_H = 1;                                                             % Thermal storage (heating) Depth Of Charge [-]
DoC_TSS_F = 1;                                                             % Thermal storage (cooling) Depth Of Charge [-]


%% EVALUATE UNCERTAIN LOADS
E_out_fix = nan(24/K_hour*T,numScenarios);                                 % Electricity requirement [kWh]
H_out_HT_fix = nan(24/K_hour*T,numScenarios);                              % High temperature heating requirement [kWh]
H_out_LT_fix = nan(24/K_hour*T,numScenarios);                              % Low temperature heating requirement [kWh]
F_out_fix = nan(24/K_hour*T,numScenarios);                                 % Cooling requirement [kWh]
H2_out_fix = nan(24/K_hour*T,numScenarios);                                % Hydrogen requirement [kWh]
W_out_fix = nan(24/K_hour*T,numScenarios);                                 % Freshwater requirement [kWh]

for j = 1:numScenarios 
    for i = 1:24/K_hour*T
        E_out_fix(i,j) = random('Normal',mu_E_dem(i),sigma_E_dem(i));
        H_out_HT_fix(i,j) = random('Normal',mu_H_HT_dem(i),sigma_H_HT_dem(i));
        H_out_LT_fix(i,j) = random('Normal',mu_H_LT_dem(i),sigma_H_LT_dem(i));
        F_out_fix(i,j) = random('Normal',mu_F_dem(i),sigma_F_dem(i));
        H2_out_fix(i,j) = random('Normal',mu_H2_dem(i),sigma_H2_dem(i));
        W_out_fix(i,j) = random('Normal',mu_W_dem(i),sigma_W_dem(i));
    end
end

clear Trends i j mu_E_dem mu_H_HT_dem mu_H_LT_dem mu_F_dem mu_H2_dem mu_W_dem sigma_E_dem sigma_H_HT_dem sigma_H_LT_dem sigma_F_dem sigma_H2_dem sigma_W_dem

%% EVALUATE UNCERTAIN RENEWABLE ENERGY AVAILABILITY

% Wind Turbine
v_aver_50 = A_w(1)*gamma(1+1/k_w(1));
v_aver_100 = A_w(2)*gamma(1+1/k_w(2));
alpha_w = log(v_aver_100/v_aver_50)/log(100/50);                           % Wind shear coefficient [-]

if h_ref == 50
    k_w = k_w(1);
    A_w = A_w(1);
elseif h_ref == 100
    k_w = k_w(2);
    A_w = A_w(2);
else
    k_w = k_w(3);
    A_w = A_w(3);
end

P_WT = nan(K_T_w*24/K_hour*T,numScenarios);                                % Wind turbine power output [kW]
e_WT = nan(T,numScenarios);                                                % Wind turbine energy output per timestep [kWh]

v_w_ref = random('Weibull',A_w,k_w,K_T_w*24/K_hour*T,numScenarios);        % Wind speed [m/s]
v_w = v_w_ref*(h_WT/h_ref)^alpha_w;                                        % Wind speed at the turbine rotor height [m/s]

for i = 1:K_T_w*24/K_hour*T*numScenarios
    if v_w(i) < v_ci
        P_WT(i) = 0;
    elseif v_w(i) >= v_co
        P_WT(i) = 0;
    elseif v_w(i) >= v_ci
        if v_w(i) < v_r
        P_WT(i) = P_r_WT*((v_w(i)-v_ci)/(v_r-v_ci));
        else
            P_WT(i) = P_r_WT;
        end
    end
end

for i = 1:24/K_hour*T*numScenarios
    e_WT(i) = sum(P_WT(i*K_T_w-K_T_w+1:i*K_T_w))/K_T_w;
end

clear k_w A_w h_ref K_T_w v_w v_ci v_co v_r P_WT i v_w_ref v_aver_50 v_aver_100 h_WT alpha_w


% Solar systems
beta_G_b = (1-mu_G_b).*((mu_G_b.*(1-mu_G_b))./sigma_G_b.^2-1);
beta_G_d = (1-mu_G_d).*((mu_G_d.*(1-mu_G_d))./sigma_G_d.^2-1);
beta_G_r = (1-mu_G_r).*((mu_G_r.*(1-mu_G_r))./sigma_G_r.^2-1);
beta_G_b(isnan(beta_G_b)) = 0;
beta_G_d(isnan(beta_G_d)) = 0;
beta_G_r(isnan(beta_G_r)) = 0;

alpha_G_b = mu_G_b.*beta_G_b./(1-mu_G_b);
alpha_G_d = mu_G_d.*beta_G_d./(1-mu_G_d);
alpha_G_r = mu_G_r.*beta_G_r./(1-mu_G_r);
alpha_G_b(isnan(alpha_G_b)) = 0;
alpha_G_d(isnan(alpha_G_d)) = 0;
alpha_G_r(isnan(alpha_G_r)) = 0;

T_air = nan(24/K_hour*T,numScenarios);                                     % Outdoor temperature [°C]
G_b = nan(24/K_hour*T,numScenarios);                                       % Direct irradiance on the inclined plane [W/m^2]
G_d = nan(24/K_hour*T,numScenarios);                                       % Diffuse irradiance on the inclined plane [W/m^2]
G_r = nan(24/K_hour*T,numScenarios);                                       % Reflected irradiance on the inclined plane [W/m^2]

for j = 1:numScenarios
    for i = 1:24/K_hour*T
        T_air(i,j) = random('Normal',mu_T_air(i),sigma_T_air(i));
        G_b(i,j) = G_0*random('Beta',beta_G_b(i),alpha_G_b(i));
        G_d(i,j) = G_0*random('Beta',beta_G_d(i),alpha_G_d(i));
        G_r(i,j) = G_0*random('Beta',beta_G_r(i),alpha_G_r(i));
    end
end

G_b(isnan(G_b)) = 0;
G_d(isnan(G_d)) = 0;
G_r(isnan(G_r)) = 0;

G_sol = G_b + G_d + G_r;                                                   % Global solar irradiance on the inclined plane [W/m^2]

T_cell = T_air + G_sol*(NOCT-T_a_NOCT)/G_NOCT;                             % PV cell temperature [°C]
eta_PV = K_PV_se*(1+beta_PV*(T_cell-T_ref));                               % PV cell efficiency [-]
e_PV = A_PV*G_sol*K_PV_BOP.*eta_PV/1000;                                   % PV electricity output [kWh]

e_CSP = A_CSP*G_b*K_CSP_se;                                                % CSP electricity output [kWh]
h_CSP = A_CSP*G_b*K_CSP_sh;                                                % CSP heating output [kWh]

h_STC = A_STC*(G_sol*eta_0 - a_1*(T_m-T_air) - a_2*(T_m-T_air).^2);        % STC heating output [kWh]
h_STC(h_STC<0) = 0;

clear Sol G_0 NOCT T_a_NOCT G_NOCT A_PV K_PV_BOP eta_PV_ref T_ref beta mu_T_air mu_G_b mu_G_d
clear mu_G_r sigma_T_air sigma_G_b sigma_G_d sigma_G_r beta_G_b beta_G_d beta_G_r beta_G_b beta_G_d beta_G_r alpha_G_b
clear alpha_G_d alpha_G_r alpha_G_b alpha_G_d alpha_G_r eta_PV T_air T_cell i j
clear a_1 a_2 eta_0 beta_PV T_m A_CSP K_CSP_se K_CSP_sh K_PV_se

%% Data checks and calculation of other simulation parameters
if sum(w) ~= 1
    error('Check weighting factors');
end

if max(w) == 1
    n = ones(1,6);
end

Ngen = 4;                                                                  % number of generators (components converting a renewable energy source in an energy vector)
Nsto = 5;                                                                  % number of storage systems 
Ncomp = 15;                                                                % number of components

nvars = 63*T + 2*Ncomp;                                                    % number of variables
eqconst = 11*T + Nsto + 6;                                                 % number of equality constraints  
ineqconst = 63*T + Ncomp;                                                  % number of inequality constraints 


%% BUILD OBJECTIVE FUNCTION VECTORS 
% Order of variables
% x = [ E_grid,in (0), E_grid,out (1), E_CHP (2), E_FC (3), E_ESS,in (4),
% E_ESS,out (5), E_HP (6), E_EC (7), E_EL (8), E_out,fl (9), H_grid,in (10),
% H_HT,grid,out (11), H_LT,grid,out (12), H_GB,HT (13), H_GB,LT (14),
% H_CHP,HT (15), H_FC,HT (16), H_CSP,HT (17), H_STC,HT (18), H_TSS,HT,in (19),
% H_TSS,HT,out (20), H_TSS,LT,out (21),  H_out,HT,fl (22), H_out,LT,fl (23),
% F_grid,in (24), F_grid,out (25), F_AC (26), F_TSS,in (27), F_TSS,out (28),
% F_out,fl (29), H_2,grid,in (30), H_2,grid,out (31), H_2,HSS,in (32),
% H_2,HSS,out (33), H_2,out,fl (34), W_grid,in (35), W_grid,out (36),
% W_WSS,in (37), W_WSS,out (38), W_out,fl (39), SOC_ESS (40), SOC_TSS,H (41),
% SOC_TSS,F (42), SOC_HSS (43), SOC_WSS (44), lambda_GB (45), lambda_CHP (46), 
% lambda_HP (47), lambda_AC (48), lambda_EL (49), lambda_FC (50),
% delta_ESS,in (51), delta_ESS,out (52), delta_TSS,H,in (53),
% delta_TSS,H,out (54), delta_TSS,F,in (55), delta_TSS,F,out (56),
% delta_HSS,in (57), delta_HSS,out (58), delta_WSS,in (59), delta_WSS,out (60),
% delta_HP (61), delta_EC (62), S_GB, SCHP, SHP, SAC, SEL, SFC, NWT, NPV,
% NCSP, NSTC, SHSS, SWSS, SESS, STSS,H, STSS,F, theta_GB, theta_CHP,
% theta_HP, theta_AC, theta_EL, theta_FC, theta_WT, theta_PV, theta_CSP,
% theta_STC, theta_HSS, theta_WSS, theta_ESS, theta_TSS,H, theta_TSS,F ]

% Cost function
f_cost = zeros(nvars,1);
f_cost(1+ 0*T: 1*T) = 365/K_day*C_grids(:,1);
f_cost(1+ 1*T: 2*T) = 365/K_day*(R_grids(:,1)-R_REC_ele);                  
f_cost(1+ 2*T: 3*T) = 365/K_day*(C_grids(:,2)/K_CHP_ge);
f_cost(1+ 3*T: 4*T) = 365/K_day*(R_REC_ele+K_FC_H2h/K_FC_H2e*R_REC_th);    
f_cost(1+ 6*T: 7*T) = 365/K_day*K_HP_eh*R_REC_th;                          
f_cost(1+ 7*T: 8*T) = 365/K_day*K_EC_ef*R_REC_th;                         
f_cost(1+ 9*T:10*T) = 365/K_day*R_flex(1);
f_cost(1+10*T:11*T) = 365/K_day*C_grids(:,3);
f_cost(1+11*T:12*T) = 365/K_day*(R_grids(:,3)-R_REC_th);                  
f_cost(1+12*T:13*T) = 365/K_day*(R_grids(:,3)-R_REC_th);                  
f_cost(1+13*T:14*T) = 365/K_day*(C_grids(:,2)/K_GB_gh);
f_cost(1+14*T:15*T) = 365/K_day*(C_grids(:,2)/K_GB_gh);
f_cost(1+22*T:24*T) = 365/K_day*R_flex(2);
f_cost(1+24*T:25*T) = 365/K_day*C_grids(:,4);
f_cost(1+25*T:26*T) = 365/K_day*R_grids(:,4);
f_cost(1+29*T:30*T) = 365/K_day*R_flex(3);
f_cost(1+30*T:31*T) = 365/K_day*C_grids(:,5);
f_cost(1+31*T:32*T) = 365/K_day*R_grids(:,5);
f_cost(1+34*T:35*T) = 365/K_day*R_flex(4);
f_cost(1+35*T:36*T) = 365/K_day*C_grids(:,6);
f_cost(1+36*T:37*T) = 365/K_day*R_grids(:,6);
f_cost(1+39*T:40*T) = 365/K_day*R_flex(4);
f_cost(end-29:end-15) = C_capex.*UCRF + C_opex;
f_cost(end-14:end) = C_capex_0.*UCRF;

% Primary Energy function
f_energy = zeros(nvars,1);
f_energy(1+ 0*T: 1*T) = 365/K_day*PE_grids(1);
f_energy(1+ 2*T: 3*T) = 365/K_day*(PE_grids(2)/K_CHP_ge);
f_energy(1+10*T:11*T) = 365/K_day*PE_grids(3);
f_energy(1+13*T:15*T) = 365/K_day*(PE_grids(2)/K_GB_gh);
f_energy(1+24*T:25*T) = 365/K_day*PE_grids(4);
f_energy(1+30*T:31*T) = 365/K_day*PE_grids(5);
f_energy(1+35*T:36*T) = 365/K_day*PE_grids(6);
f_energy(end-29:end-15) = PE_embodied./UL;

% Carbon function
f_carbon = zeros(nvars,1);
f_carbon(1+ 0*T: 1*T) = 365/K_day*GWP_grids(1);
f_carbon(1+ 2*T: 3*T) = 365/K_day*(GWP_grids(2)/K_CHP_ge);
f_carbon(1+10*T:11*T) = 365/K_day*GWP_grids(3);
f_carbon(1+13*T:15*T) = 365/K_day*(GWP_grids(2)/K_GB_gh);
f_carbon(1+24*T:25*T) = 365/K_day*GWP_grids(4);
f_carbon(1+30*T:31*T) = 365/K_day*GWP_grids(5);
f_carbon(1+35*T:36*T) = 365/K_day*GWP_grids(6);
f_carbon(end-29:end-15) = GWP_embodied./UL;

% Grid Interaction function
f_grid = zeros(nvars,1);
f_grid(1+ 0*T: 1*T) = GII_grids_in(1);                                     % Egrid,in
f_grid(1+ 1*T: 2*T) = GII_grids_out(1);                                    % Egrid,out
f_grid(1+ 2*T: 3*T) = GII_grids_in(2);                                     % CHP,in
f_grid(1+10*T:11*T) = GII_grids_in(3);                                     % Hgrid,in
f_grid(1+11*T:12*T) = GII_grids_out(3);                                    % HHT,grid,out
f_grid(1+12*T:13*T) = GII_grids_out(3);                                    % HLT,grid,out
f_grid(1+13*T:14*T) = GII_grids_in(2);                                     % GB,in
f_grid(1+14*T:15*T) = GII_grids_in(2);                                     % GB,in
f_grid(1+24*T:25*T) = GII_grids_in(4);                                     % Fgrid,in
f_grid(1+25*T:26*T) = GII_grids_out(4);                                    % Fgrid,out
f_grid(1+30*T:31*T) = GII_grids_in(5);                                     % H2 grid,in
f_grid(1+31*T:32*T) = GII_grids_out(5);                                    % H2 grid,out
f_grid(1+35*T:36*T) = GII_grids_in(6);                                     % Wgrid,in 
f_grid(1+36*T:37*T) = GII_grids_out(6);                                    % Wgrid,out

% CRM function
f_CRM = zeros(nvars,1);
f_CRM(end-29:end-15) = CRM_embodied./UL;

f = [f_cost f_energy f_carbon f_grid f_CRM];
F = sum(w.*f./n,2);
fitnessfcn = @(sol) sum(sol*f,1);

%% OPTIMIZATION ALGORITHM INPUTS
A_eq = zeros(eqconst,nvars);
b_eq = zeros(eqconst,1);
A = zeros(ineqconst,nvars);       
b = zeros(ineqconst,1);           
lb = zeros(nvars,1);
ub = inf(nvars,1);
intcon = [45*T+1:63*T, 63*T+7: 63*T+10, 63*T+Ncomp+1:nvars];

% Upper bounds
ub(1+ 0*T: 1*T) = E_grid_in*Q_E;                                           % E_grid_in
ub(1+ 1*T: 2*T) = E_grid_out*Q_E;                                          % E_grid_out
ub(1+ 2*T: 3*T) = NG_grid_in*CHP*Q_E;                                      % E_CHP
ub(1+ 3*T: 4*T) = FC*Q_H2*K_FC_H2e;                                        % E_FC
ub(1+ 4*T: 5*T) = ESS*Q_E;                                                 % E_ESS_in
ub(1+ 5*T: 6*T) = ESS*Q_E;                                                 % E_ESS_out
ub(1+ 6*T: 7*T) = HP*Q_F/K_EC_ef;                                          % E_HP
ub(1+ 7*T: 8*T) = EC*Q_E;                                                  % E_EC
ub(1+ 8*T: 9*T) = EL*Q_H2/K_EL_eH2;                                        % E_EL
ub(1+ 9*T:10*T) = Q_E;                                                     % E_out_fl
ub(1+10*T:11*T) = H_grid_in*Q_H;                                           % H_grid_in
ub(1+11*T:12*T) = H_grid_out*Q_H;                                          % H_HT_grid_out
ub(1+12*T:13*T) = H_grid_out*Q_H;                                          % H_LT_grid_out
ub(1+13*T:14*T) = NG_grid_in*GB*Q_H;                                       % H_GB_HT
ub(1+14*T:15*T) = NG_grid_in*GB*Q_H;                                       % H_GB_LT
ub(1+15*T:16*T) = CHP*Q_H;                                                 % H_CHP_HT
ub(1+16*T:17*T) = FC*Q_H;                                                  % H_FC_HT
ub(1+17*T:18*T) = CSP*Q_H;                                                 % H_CSP_HT
ub(1+18*T:19*T) = STC*Q_H;                                                 % H_STC_HT
ub(1+19*T:20*T) = TSS_H*Q_H;                                               % H_TSS_HT_in
ub(1+20*T:21*T) = TSS_H*Q_H;                                               % H_TSS_HT_out
ub(1+21*T:22*T) = TSS_H*Q_H;                                               % H_TSS_LT_out
ub(1+22*T:23*T) = Q_H;                                                     % H_out_HT_fl
ub(1+23*T:24*T) = Q_H;                                                     % H_out_LT_fl
ub(1+24*T:25*T) = F_grid_in*Q_F;                                           % F_grid_in
ub(1+25*T:26*T) = F_grid_out*Q_F;                                          % F_grid_out
ub(1+26*T:27*T) = AC*Q_F;                                                  % F_AC
ub(1+27*T:28*T) = TSS_F*Q_F;                                               % F_TSS_in
ub(1+28*T:29*T) = TSS_F*Q_F;                                               % F_TSS_out
ub(1+29*T:30*T) = Q_F;                                                     % F_out_fl
ub(1+30*T:31*T) = H2_grid_in*Q_H2;                                         % H2_grid_in
ub(1+31*T:32*T) = H2_grid_out*Q_H2;                                        % H2_grid_out
ub(1+32*T:33*T) = HSS*Q_H2;                                                % H2_HSS_in
ub(1+33*T:34*T) = HSS*Q_H2;                                                % H2_HSS_out
ub(1+34*T:35*T) = Q_H2;                                                    % H2_out_fl
ub(1+35*T:36*T) = W_grid_in*Q_W;                                           % W_grid_in
ub(1+36*T:37*T) = W_grid_out*Q_W;                                          % W_grid_out
ub(1+37*T:38*T) = WSS*Q_W;                                                 % W_WSS_in
ub(1+38*T:39*T) = WSS*Q_W;                                                 % W_WSS_out
ub(1+39*T:40*T) = Q_W;                                                     % W_out_fl
ub(1+40*T:41*T) = ESS*Q_E;                                                 % SOC_ESS
ub(1+41*T:42*T) = TSS_H*Q_H;                                               % SOC_TSS_H
ub(1+42*T:43*T) = TSS_F*Q_F;                                               % SOC_TSS_F
ub(1+43*T:44*T) = HSS*Q_H2;                                                % SOC_HSS
ub(1+44*T:45*T) = WSS*Q_W;                                                 % SOC_WSS
ub(1+45*T:63*T) = 1;                                                       % lambdas and deltas
ub(end-29) = GB*Q_H;                                                       % S_GB
ub(end-28) = CHP*Q_E;                                                      % S_CHP
ub(end-27) = HP*Q_F;                                                       % S_HP
ub(end-26) = AC*Q_F;                                                       % S_AC
ub(end-25) = EL*Q_E;                                                       % S_EL
ub(end-24) = FC*Q_E;                                                       % S_FC
ub(end-23) = WT*Q_E;                                                       % N_WT
ub(end-22) = PV*Q_E;                                                       % N_PV
ub(end-21) = CSP*Q_E;                                                      % N_CSP
ub(end-20) = STC*Q_H;                                                      % N_STC
ub(end-19) = HSS*Q_H2;                                                     % S_HSS
ub(end-18) = WSS*Q_W;                                                      % S_WSS
ub(end-17) = ESS*Q_E;                                                      % S_ESS
ub(end-16) = TSS_H*Q_H;                                                    % S_TSS_H
ub(end-15) = TSS_F*Q_F;                                                    % S_TSS_F
ub(end-14:end) = 1;                                                        % thetas

% Inequality constraints vector
b(2*T+1:3*T) = 1;                                                          % electrical storage status
b(7*T+1:8*T) = 1;                                                          % hot thermal storage status
b(12*T+1:13*T) = 1;                                                        % cold thermal storage status
b(17*T+1:18*T) = 1;                                                        % hydrogen storage status
b(22*T+1:23*T) = 1;                                                        % water storage status
b(32*T+1:33*T) = 1;                                                        % HP 

% Inequality constraints coefficients matrix
for i = 1:T

% ESS
A(i+0*T,i+4*T) = 1;
A(i+0*T,i+51*T) = -Q_E;
A(i+1*T,i+5*T) = 1;
A(i+1*T,i+52*T) = -Q_E;
A(i+2*T,i+51*T) = 1;
A(i+2*T,i+52*T) = 1;
A(i+3*T,63*T+13) = 24/K_hour*DoD_EES;
A(i+3*T,i+40*T) = -1;
A(i+4*T,i+40*T) = 1;
A(i+4*T,63*T+13) = -24/K_hour*DoC_EES;

% TSS,H
A(i+5*T,i+19*T) = 1;
A(i+5*T,i+53*T) = -Q_H;
A(i+6*T,i+20*T) = 1;
A(i+6*T,i+21*T) = 1;
A(i+6*T,i+54*T) = -Q_H;
A(i+7*T,i+53*T) = 1;
A(i+7*T,i+54*T) = 1;
A(i+8*T,63*T+14) = 24/K_hour*DoD_TSS_H;
A(i+8*T,i+41*T) = -1;
A(i+9*T,i+41*T) = 1;
A(i+9*T,63*T+14) = -24/K_hour*DoC_TSS_H;

% TSS,F
A(i+10*T,i+27*T) = 1;
A(i+10*T,i+55*T) = -Q_F;
A(i+11*T,i+28*T) = 1;
A(i+11*T,i+56*T) = -Q_F;
A(i+12*T,i+55*T) = 1;
A(i+12*T,i+56*T) = 1;
A(i+13*T,63*T+15) = 24/K_hour*DoD_TSS_F;
A(i+13*T,i+42*T) = -1;
A(i+14*T,i+42*T) = 1;
A(i+14*T,63*T+15) = -24/K_hour*DoC_TSS_F;

% HSS
A(i+15*T,i+32*T) = 1;
A(i+15*T,i+57*T) = -Q_H2;
A(i+16*T,i+33*T) = 1;
A(i+16*T,i+58*T) = -Q_H2;
A(i+17*T,i+57*T) = 1;
A(i+17*T,i+58*T) = 1;
A(i+18*T,63*T+11) = 24/K_hour*DoD_HSS;
A(i+18*T,i+43*T) = -1;
A(i+19*T,i+43*T) = 1;
A(i+19*T,63*T+11) = -24/K_hour*DoC_HSS;

% WSS
A(i+20*T,i+37*T) = 1;
A(i+20*T,i+59*T) = -Q_W;
A(i+21*T,i+38*T) = 1;
A(i+21*T,i+60*T) = -Q_W;
A(i+22*T,i+59*T) = 1;
A(i+22*T,i+60*T) = 1;
A(i+23*T,63*T+12) = 24/K_hour*DoD_WSS;
A(i+23*T,i+44*T) = -1;
A(i+24*T,i+44*T) = 1;
A(i+24*T,63*T+12) = -24/K_hour*DoC_WSS;

% CHP
A(i+25*T,i+46*T) = E_CHP_min;
A(i+25*T,i+2*T)=-1;
A(i+26*T,i+2*T)=1;
A(i+26*T,i+46*T)=-Q_E;
A(i+27*T,i+2*T)=1;
A(i+27*T,63*T+2)=-1;
A(i+28*T,i+2*T)=1;
A(i+28*T,i+2*T+1)=-1;
A(i+28*T,63*T+2)=-RD_CHP;
A(i+29*T,i+2*T+1)=1;
A(i+29*T,i+2*T)=-1;
A(i+29*T,63*T+2)=-RU_CHP;

% HP
A(i+30*T,i+6*T)=1;
A(i+30*T,i+61*T)=-Q_F;
A(i+31*T,i+7*T)=1;
A(i+31*T,i+62*T)=-Q_F;
A(i+32*T,i+61*T)=1;
A(i+32*T,i+62*T)=1;
A(i+33*T,i+47*T)=F_HP_min;
A(i+33*T,i+7*T)=-K_EC_ef;
A(i+34*T,i+7*T)=1;
A(i+34*T,i+47*T)=-Q_F;
A(i+35*T,i+7*T)=K_EC_ef;
A(i+35*T,63*T+3)=-1;
A(i+36*T,i+47*T)=H_HP_min;
A(i+36*T,i+6*T)=-K_HP_eh;
A(i+37*T,i+6*T)=1;
A(i+37*T,i+47*T)=-Q_F;
A(i+38*T,i+6*T)=K_EC_ef;
A(i+38*T,63*T+3)=-1;
A(i+39*T,i+7*T)=1;
A(i+39*T,i+7*T+1)=-1;
A(i+39*T,63*T+3)=-RD_EC;
A(i+40*T,i+7*T+1)=1;
A(i+40*T,i+7*T)=-1;
A(i+40*T,63*T+3)=-RU_EC;
A(i+41*T,i+6*T)=1;
A(i+41*T,i+6*T+1)=-1;
A(i+41*T,63*T+3)=-RD_EC*K_HP_eh/K_EC_ef;
A(i+42*T,i+6*T+1)=1;
A(i+42*T,i+6*T)=-1;
A(i+42*T,63*T+3)=-RU_EC*K_HP_eh/K_EC_ef;

% AC
A(i+43*T,i+48*T)=F_AC_min;
A(i+43*T,i+26*T)=-1;
A(i+44*T,i+26*T)=1;
A(i+44*T,i+48*T)=-Q_F;
A(i+45*T,i+26*T)=1;
A(i+45*T,63*T+4)=-1;
A(i+46*T,i+26*T)=1;
A(i+46*T,i+26*T+1)=-1;
A(i+46*T,63*T+4)=-RD_AC;
A(i+47*T,i+26*T+1)=1;
A(i+47*T,i+26*T)=-1;
A(i+47*T,63*T+4)=-RU_AC;

% EL
A(i+48*T,i+49*T)=E_EL_min;
A(i+48*T,i+8*T)=-1;
A(i+49*T,i+8*T)=1;
A(i+49*T,i+49*T)=-Q_E;
A(i+50*T,i+8*T)=1;
A(i+50*T,63*T+5)=-1;
A(i+51*T,i+8*T)=1;
A(i+51*T,i+8*T+1)=-1;
A(i+51*T,63*T+5)=-RD_EL;
A(i+52*T,i+8*T+1)=1;
A(i+52*T,i+8*T)=-1;
A(i+52*T,63*T+5)=-RU_EL;

% FC
A(i+53*T,i+50*T)=E_FC_min;
A(i+53*T,i+3*T)=-1;
A(i+54*T,i+3*T)=1;
A(i+54*T,i+50*T)=-Q_E;
A(i+55*T,i+3*T)=1;
A(i+55*T,63*T+6)=-1;
A(i+56*T,i+3*T)=1;
A(i+56*T,i+3*T+1)=-1;
A(i+56*T,63*T+6)=-RD_FC;
A(i+57*T,i+3*T+1)=1;
A(i+57*T,i+3*T)=-1;
A(i+57*T,63*T+6)=-RU_FC;

% GB
A(i+58*T,i+45*T)=H_GB_min;
A(i+58*T,i+13*T)=-1;
A(i+58*T,i+14*T)=-1;
A(i+59*T,i+13*T)=1;
A(i+59*T,i+14*T)=1;
A(i+59*T,i+45*T)=-Q_H;
A(i+60*T,i+13*T)=1;
A(i+60*T,i+14*T)=1;
A(i+60*T,63*T+1)=-1;
A(i+61*T,i+13*T)=1;
A(i+61*T,i+14*T)=1;
A(i+61*T,i+13*T+1)=-1;
A(i+61*T,i+14*T+1)=-1;
A(i+61*T,63*T+1)=-RD_GB;
A(i+62*T,i+13*T+1)=1;
A(i+62*T,i+14*T+1)=1;
A(i+62*T,i+13*T)=-1;
A(i+62*T,i+14*T)=-1;
A(i+62*T,63*T+1)=-RU_GB;

end

% Link between design and synthesis variables
for i=1:Ncomp
    A(i+63*T,63*T+i)=1;
end
A(1 +63*T,63*T+16) = -Q_H;                                                 % GB
A(2 +63*T,63*T+17) = -Q_E;                                                 % CHP
A(3 +63*T,63*T+18) = -Q_F;                                                 % HP
A(4 +63*T,63*T+19) = -Q_F;                                                 % AC
A(5 +63*T,63*T+20) = -Q_E;                                                 % EL
A(6 +63*T,63*T+21) = -Q_E;                                                 % FC
A(7 +63*T,63*T+22) = -Q_E;                                                 % WT
A(8 +63*T,63*T+23) = -Q_E;                                                 % PV
A(9 +63*T,63*T+24) = -Q_E;                                                 % CSP
A(10+63*T,63*T+25) = -Q_H;                                                 % STC
A(11+63*T,63*T+26) = -Q_H2;                                                % HSS
A(12+63*T,63*T+27) = -Q_W;                                                 % WSS
A(13+63*T,63*T+28) = -Q_E;                                                 % ESS
A(14+63*T,63*T+29) = -Q_H;                                                 % TSS,H
A(15+63*T,63*T+30) = -Q_F;                                                 % TSS,F

% Corrections to the A matrix
A(29*T,2*T+1)=-1;
A(29*T,3*T+1)=0;
A(30*T,3*T+1)=0;
A(30*T,2*T+1)=1;

A(40*T,7*T+1)=-1;
A(40*T,8*T+1)=0;
A(41*T,7*T+1)=1;
A(41*T,8*T+1)=0;

A(42*T,6*T+1)=-1;
A(42*T,7*T+1)=0;
A(43*T,6*T+1)=1;
A(43*T,7*T+1)=0;

A(47*T,26*T+1)=-1;
A(47*T,27*T+1)=0;
A(48*T,26*T+1)=1;
A(48*T,27*T+1)=0;

A(52*T,8*T+1)=-1;
A(52*T,9*T+1)=0;
A(53*T,8*T+1)=1;
A(53*T,9*T+1)=0;

A(57*T,3*T+1)=-1;
A(57*T,4*T+1)=0;
A(58*T,3*T+1)=1;
A(58*T,4*T+1)=0;

A(62*T,13*T+1)=-1;
A(62*T,15*T+1)=0;
A(63*T,13*T+1)=1;
A(63*T,15*T+1)=0;

if Method == 's'
    x_sol = nan(numScenarios,nvars);
    x_fval = nan(numScenarios,1);
elseif Method == 'g'
    x_sol = nan(1,nvars);
end
x_exitflag = nan(numScenarios,1);

for p = 1:numScenarios
p                                                                          % no semicolon to check the current iteration during the optimisation process

% Equality constraints vector
for i = 1:T
b_eq(i+0*T) = sum(E_out_fix(i*24/K_hour-24/K_hour+1:i*24/K_hour,p));
b_eq(i+1*T) = sum(H_out_HT_fix(i*24/K_hour-24/K_hour+1:i*24/K_hour,p));
b_eq(i+2*T) = sum(H_out_LT_fix(i*24/K_hour-24/K_hour+1:i*24/K_hour,p));
b_eq(i+3*T) = sum(F_out_fix(i*24/K_hour-24/K_hour+1:i*24/K_hour,p));
b_eq(i+4*T) = sum(H2_out_fix(i*24/K_hour-24/K_hour+1:i*24/K_hour,p));
b_eq(i+5*T) = sum(W_out_fix(i*24/K_hour-24/K_hour+1:i*24/K_hour,p));
end

b_eq(end-5) = E_out_flex_TOT;
b_eq(end-4) = H_out_HT_flex_TOT;
b_eq(end-3) = H_out_LT_flex_TOT;
b_eq(end-2) = F_out_flex_TOT;
b_eq(end-1) = H2_out_flex_TOT;
b_eq(end) = W_out_flex_TOT;


% Equality constraints coefficients matrix
for i = 1:T
   
% Electrical energy balance
A_eq(i+0*T,i+0*T) = K_E_grid;
A_eq(i+0*T,i+1*T) = - 1/K_E_grid;
A_eq(i+0*T,i+2*T) = 1;
A_eq(i+0*T,i+3*T) = 1;
A_eq(i+0*T,i+4*T) = -1;
A_eq(i+0*T,i+5*T) = 1;
A_eq(i+0*T,i+6*T) = -1;
A_eq(i+0*T,i+7*T) = -1;
A_eq(i+0*T,i+8*T) = -1;
A_eq(i+0*T,i+9*T) = -1;
A_eq(i+0*T,end-23) = e_WT(i,p);
A_eq(i+0*T,end-22) = e_PV(i,p);
A_eq(i+0*T,end-21) = e_CSP(i,p);

% High temperature heating energy balance
A_eq(i+1*T,i+10*T) = K_H_grid;
A_eq(i+1*T,i+11*T) = - 1/K_H_grid;
A_eq(i+1*T,i+13*T) = 1;
A_eq(i+1*T,i+15*T) = 1;
A_eq(i+1*T,i+16*T) = 1;
A_eq(i+1*T,i+17*T) = 1;
A_eq(i+1*T,i+18*T) = 1;
A_eq(i+1*T,i+19*T) = -1;
A_eq(i+1*T,i+20*T) = 1;
A_eq(i+1*T,i+6*T) = K_HP_eh;
A_eq(i+1*T,i+26*T) = -1/ K_AC_hf;
A_eq(i+1*T,i+22*T) = -1;

% Low temperature heating energy balance
A_eq(i+2*T,i+12*T) = - 1/K_H_grid;
A_eq(i+2*T,i+14*T) = 1;
A_eq(i+2*T,i+2*T) = K_CHP_gh/K_CHP_ge;
A_eq(i+2*T,i+15*T) = -1;
A_eq(i+2*T,i+3*T) = K_FC_H2h/K_FC_H2e;
A_eq(i+2*T,i+16*T) = -1;
A_eq(i+2*T,end-21) = h_CSP(i,p);
A_eq(i+2*T,i+17*T) = -1;
A_eq(i+2*T,end-20) = h_STC(i,p);
A_eq(i+2*T,i+18*T) = -1;
A_eq(i+2*T,i+21*T) = 1;
A_eq(i+2*T,i+23*T) = -1;

% Cooling energy balance
A_eq(i+3*T,i+24*T) = K_F_grid;
A_eq(i+3*T,i+25*T) = - 1/K_F_grid;
A_eq(i+3*T,i+7*T) = K_EC_ef;
A_eq(i+3*T,i+26*T) = 1;
A_eq(i+3*T,i+27*T) = -1;
A_eq(i+3*T,i+28*T) = 1;
A_eq(i+3*T,i+29*T) = -1;

% Hydrogen mass balance
A_eq(i+4*T,i+30*T) = K_H2_grid;
A_eq(i+4*T,i+31*T) = -1/K_H2_grid;
A_eq(i+4*T,i+8*T) = K_EL_eH2;
A_eq(i+4*T,i+3*T) = -1/K_FC_H2e;
A_eq(i+4*T,i+32*T) = -1;
A_eq(i+4*T,i+33*T) = 1;
A_eq(i+4*T,i+34*T) = -1;

% Water mass balance
A_eq(i+5*T,i+35*T) = K_W_grid;
A_eq(i+5*T,i+36*T) = -1/K_W_grid;
A_eq(i+5*T,i+8*T) = -K_EL_eH2*K_EL_H2w;
A_eq(i+5*T,i+3*T) = K_FC_H2w/K_FC_H2e;
A_eq(i+5*T,i+37*T) = -1;
A_eq(i+5*T,i+38*T) = 1;
A_eq(i+5*T,i+39*T) = -1;

% Electrical storage hourly balance
A_eq(i+6*T,i+40*T-1) = 1-phi_ESS;
A_eq(i+6*T,i+40*T) = -1;
A_eq(i+6*T,i+4*T) = K_ESS_in;
A_eq(i+6*T,i+5*T) = -1/K_ESS_out;

% Hot thermal storage hourly balance
A_eq(i+7*T,i+41*T-1) = 1-phi_TSS_H;
A_eq(i+7*T,i+41*T) = -1;
A_eq(i+7*T,i+19*T) = K_TSS_H_in;
A_eq(i+7*T,i+20*T) = -1/K_TSS_H_out;
A_eq(i+7*T,i+21*T) = -1/K_TSS_H_out;

% Cold thermal storage hourly balance
A_eq(i+8*T,i+42*T-1) = 1-phi_TSS_F;
A_eq(i+8*T,i+42*T) = -1;
A_eq(i+8*T,i+27*T) = K_TSS_F_in;
A_eq(i+8*T,i+28*T) = -1/K_TSS_F_out;

% Hydrogen storage balance:
A_eq(i+9*T,i+43*T-1) = 1-phi_HSS;
A_eq(i+9*T,i+43*T) = -1;
A_eq(i+9*T,i+32*T) = K_HSS_in;
A_eq(i+9*T,i+33*T) = -1/K_HSS_out;

% Water storage balance:
A_eq(i+10*T,i+44*T-1) = 1-phi_WSS;
A_eq(i+10*T,i+44*T) = -1;
A_eq(i+10*T,i+37*T) = K_WSS_in;
A_eq(i+10*T,i+38*T) = -1/K_WSS_out;

end

% Electrical storage period balance for first and last hours
A_eq(1+6*T,41*T) = 1 - phi_ESS;
A_eq(1+6*T,40*T) = 0;
A_eq(1+11*T,1+40*T) = 1;
A_eq(1+11*T,41*T) = -1;

% Hot thermal storage period balance for first and last hours
A_eq(1+7*T,42*T) = 1 - phi_TSS_H;
A_eq(1+7*T,41*T) = 0;
A_eq(2+11*T, 1+41*T) = 1;
A_eq(2+11*T,42*T) = -1;

% Cold thermal storage period balance for first and last hours
A_eq(1+8*T,43*T) = 1 - phi_TSS_F;
A_eq(1+8*T,42*T) = 0;
A_eq(3+11*T, 1+42*T) = 1;
A_eq(3+11*T,43*T) = -1;

% Hydrogen storage period balance for first and last hours
A_eq(1+9*T,44*T) = 1 - phi_HSS;
A_eq(1+9*T,43*T) = 0;
A_eq(4+11*T, 1+43*T) = 1;
A_eq(4+11*T,44*T) = -1;

% Water storage period balance for first and last hours
A_eq(1+10*T,45*T) = 1 - phi_WSS;
A_eq(1+10*T,44*T) = 0;
A_eq(5+11*T,1+44*T) = 1;
A_eq(5+11*T,45*T) = -1;

% Flexible loads balance
A_eq(6 +11*T,1+ 9*T:10*T) = 1;
A_eq(7 +11*T,1+22*T:23*T) = 1;
A_eq(8 +11*T,1+23*T:24*T) = 1;
A_eq(9 +11*T,1+29*T:30*T) = 1;
A_eq(10+11*T,1+34*T:35*T) = 1;
A_eq(11+11*T,1+39*T:40*T) = 1;

% Cost objective function
f_cost(end-23) = f_cost(end-23)+365/K_day*sum(e_WT(:,p))*R_REC_ele;        
f_cost(end-22) = f_cost(end-22)+365/K_day*sum(e_PV(:,p))*R_REC_ele;        
f_cost(end-21) = f_cost(end-21)+365/K_day*(sum(e_CSP(:,p))*R_REC_ele+sum(h_CSP(:,p))*R_REC_th);
f_cost(end-20) = f_cost(end-20)+365/K_day*sum(h_STC(:,p))*R_REC_th;      

% Solutions
if Method == 's'
    A = sparse(A);
    A_eq = sparse(A_eq);
    [sol,fval,exitflag,output] = intlinprog(F,intcon,A,b,A_eq,b_eq,lb,ub);
    if numScenarios ~= 1
        x_sol(p,:) = sol';
        x_fval(p) = fval;
        x_exitflag(p) = exitflag;
        x_output(p) = output;
    else
        FVAL = sol'*f;
    end
elseif Method == 'g'
    [sol,fval,exitflag,output,population,scores] = gamultiobj(fitnessfcn,nvars,A,b,A_eq,b_eq,lb,ub,[],intcon);
    if numScenarios ~= 1
        [m,~] = size(sol);
        x_sol(p*m-m+1:p*m,:) = sol;
        x_exitflag(p) = exitflag;
        x_output(p) = output;
    end
end

end

%% OPTIMIZATION ALGORITHM OUTPUTS

if numScenarios ~= 1
        x_FVAL = x_sol*f;
end

if Method == 's'
    if numScenarios == 1
        synthesis_solution = sol(end-14:end);
        design_solution = sol(end-29:end-15);

        figure
        plot(FVAL(:,1),FVAL(:,2),'ro')
        title('Objective Functions Results 1')
        xlabel('Annualized Cost [Euro]')
        ylabel('Annualized Life Cycle Primary Energy Use [MJ]')

        figure
        plot(FVAL(:,1),FVAL(:,3),'go')
        title('Objective Functions Results 2')
        xlabel('Annualized Cost [Euro]')
        ylabel('Annualized Life Cycle Carbon Emissions [kg CO_2-eq]')

        figure
        plot(FVAL(:,2),FVAL(:,3),'bo')
        title('Objective Functions Results 3')
        xlabel('Annualized Life Cycle Primary Energy Use [MJ]')
        ylabel('Annualized Life Cycle Carbon Emissions [kg CO_2-eq]')

        figure
        plot(sol(end-22),sol(end-19),'c*')
        title('PV against ESS Sizes')
        xlabel('Photovoltaic System Rated Power [kWp]')
        ylabel('Electrical Energy Storage Rated Capacity [kWh]')

                
        figure
        scatter3(sol(end-22),sol(end-16),sol(end-25),'mx')
        title('PV against HSS against EL Sizes')
        xlabel('Photovoltaic System Rated Power [kWp]')
        ylabel('Hydrogen Storage Rated Capacity [kg]')
        zlabel('Water Electrolyzer Rated Capacity [kW]')
    else
        synthesis_solution = zeros(4,15);
        design_solution = zeros(4,15);
        for i = 1:15
            synthesis_solution(1,i) = max(x_sol(:,end-15+i));
            synthesis_solution(2,i) = mean(x_sol(:,end-15+i));
            synthesis_solution(3,i) = mode(x_sol(:,end-15+i));
            synthesis_solution(4,i) = median(x_sol(:,end-15+i));
            design_solution(1,i) = max(x_sol(:,end-30+i));
            design_solution(2,i) = mean(x_sol(:,end-30+i));
            design_solution(3,i) = mode(x_sol(:,end-30+i));
            design_solution(4,i) = median(x_sol(:,end-30+i));
        end

        figure
        plot(x_FVAL(:,1),x_FVAL(:,2),'ro')
        title('Objective Functions Results 1')
        xlabel('Annualized Cost [Euro]')
        ylabel('Annualized Life Cycle Primary Energy Use [MJ]')

        figure
        plot(x_FVAL(:,1),x_FVAL(:,3),'go')
        title('Objective Functions Results 2')
        xlabel('Annualized Cost [Euro]')
        ylabel('Annualized Life Cycle Carbon Emissions [kg CO_2-eq]')

        figure
        plot(x_FVAL(:,2),x_FVAL(:,3),'bo')
        title('Objective Functions Results 3')
        xlabel('Annualized Life Cycle Primary Energy Use [MJ]')
        ylabel('Annualized Life Cycle Carbon Emissions [kg CO_2-eq]')

        figure
        scatter3(x_FVAL(:,1),x_FVAL(:,2),x_FVAL(:,3),'m+')
        title('Objective Functions Results 4')
        xlabel('Annualized Cost [Euro]')
        ylabel('Annualized Life Cycle Primary Energy Use [MJ]')
        zlabel('Annualized Life Cycle Carbon Emissions [kg CO_2-eq]')

        figure
        plot(x_sol(:,end-22),x_sol(:,end-19),'c*')
        title('PV against ESS Sizes')
        xlabel('Photovoltaic System Rated Power [kWp]')
        ylabel('Electrical Energy Storage Rated Capacity [kWh]')
        
        figure
        scatter3(x_sol(:,end-22),x_sol(:,end-16),x_sol(:,end-25),'mx')
        title('PV against HSS against EL Sizes')
        xlabel('Photovoltaic System Rated Power [kWp]')
        ylabel('Hydrogen Storage Rated Capacity [kg]')
        zlabel('Water Electrolyzer Rated Capacity [kW]')
    end
elseif Method == 'g'
    if numScenarios == 1
        figure
        plot(fval(:,1),fval(:,2),'ro')
        title('Objective Functions Results 1')
        xlabel('Annualized Cost [Euro]')
        ylabel('Annualized Life Cycle Primary Energy Use [MJ]')

        figure
        plot(fval(:,1),fval(:,3),'go')
        title('Objective Functions Results 2')
        xlabel('Annualized Cost [Euro]')
        ylabel('Annualized Life Cycle Carbon Emissions [kg CO_2-eq]')

        figure
        plot(fval(:,2),fval(:,3),'bo')
        title('Objective Functions Results 3')
        xlabel('Annualized Life Cycle Primary Energy Use [MJ]')
        ylabel('Annualized Life Cycle Carbon Emissions [kg CO_2-eq]')

        figure
        plot(sol(:,end-22),sol(:,end-19),'c*')
        title('PV against ESS Sizes')
        xlabel('Photovoltaic System Rated Power [kWp]')
        ylabel('Electrical Energy Storage Rated Capacity [kWh]')
        
        figure
        scatter3(x_sol(:,end-22),x_sol(:,end-16),x_sol(:,end-25),'mx')
        title('PV against HSS against EL Sizes')
        xlabel('Photovoltaic System Rated Power [kWp]')
        ylabel('Hydrogen Storage Rated Capacity [kg]')
        zlabel('Water Electrolyzer Rated Capacity [kW]')
    else
        figure
        plot(x_FVAL(:,1),x_FVAL(:,2),'ro')
        title('Objective Functions Results 1')
        xlabel('Annualized Cost [Euro]')
        ylabel('Annualized Life Cycle Primary Energy Use [MJ]')

        figure
        plot(x_FVAL(:,1),x_FVAL(:,3),'go')
        title('Objective Functions Results 2')
        xlabel('Annualized Cost [Euro]')
        ylabel('Annualized Life Cycle Carbon Emissions [kg CO_2-eq]')

        figure
        plot(x_FVAL(:,2),x_FVAL(:,3),'bo')
        title('Objective Functions Results 3')
        xlabel('Annualized Life Cycle Primary Energy Use [MJ]')
        ylabel('Annualized Life Cycle Carbon Emissions [kg CO_2-eq]')

        figure
        plot(x_sol(:,end-22),x_sol(:,end-19),'c*')
        title('PV against ESS against HSS Sizes')
        xlabel('Photovoltaic System Rated Power [kWp]')
        ylabel('Electrical Energy Storage Rated Capacity [kWh]')
    end
end

clear E_grid_in NG_grid_in H_grid_in F_grid_in W_grid_in H2_grid_in GB CHP HP EC AC EL FC WT PV CSP STC HSS WSS ESS TSS_H TSS_F
clear Q_E Q_H Q_F Q_H2 Q_W LHV_NG K_E_grid K_H_grid K_F_grid K_W_grid K_H2_grid
clear K_GB_gh K_CHP_ge K_CHP_gh K_HP_eh K_EC_ef K_AC_hf K_EL_eH2 K_EL_H2w K_FC_H2h K_FC_H2e K_FC_H2w
clear H_GB_min E_CHP_min F_HP_min H_HP_min F_AC_min E_EL_min E_FC_min RU_GB RU_CHP RU_EC RU_AC
clear RU_EL RU_FC RD_GB RD_CHP RD_EC RD_AC RD_EL RD_FC K_HSS_in K_HSS_out phi_HSS DoD_HSS DoC_HSS Q_HSS
clear K_WSS_in K_WSS_out phi_WSS DoD_WSS DoC_WSS K_ESS_in K_ESS_out phi_ESS DoD_EES DoC_EES
clear K_TSS_H_in K_TSS_H_out K_TSS_F_in K_TSS_F_out phi_TSS_H phi_TSS_F DoD_TSS_H DoD_TSS_F DoC_TSS_H DoC_TSS_F
clear Ngen Nsto Ncomp nvars eqconst ineqconst p p_CSP UL C_capex C_capex_0 C_opex C_grids 
clear R_grids R_flex R_REC_ele R_REC_th i_def UCRF PE_embodied PE_grids GWP_embodied GWP_grids
clear GII_grids_in GII_grids_out CRM_embodied CRM_grids Equip Grids intcon P_r_WT A_STC P_r_PV P_r_CSP 
clear A A_eq b b_eq E_out_fix F_out_fix H2_out_fix H_out_HT_fix H_out_LT_fix W_out_fix
clear lb ub i E_out_flex_TOT H_out_HT_flex_TOT H_out_LT_flex_TOT F_out_flex_TOT H2_out_flex_TOT W_out_flex_TOT G_b G_d G_r G_sol
