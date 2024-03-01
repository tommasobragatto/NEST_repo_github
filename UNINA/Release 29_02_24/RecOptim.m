%% SOFTWARE FOR PROBABILISTIC OPTIMAL CONTROL OF RESOURCES IN A RENEWABLE ENERGY COMMUNITIY
%  ***********************************************************************************************************************************************
%  RELEASE 29/02/2024 
%  FABIO MOTTOLA AND DANIELA PROTO UNIVERSITY OF NAPLES FEDERICO II, NAPLES (ITALY)
%  fabio.mottola@unina.it, daniela.proto@unina.it
%
%  THE SOFTWARE IS PART OF Project funded under the National Recovery and Resilience Plan (NRRP), 
%  Mission 4 Component 2 Investment 1.3 - Call for tender No. 1561 of 11.10.2022 of Ministero dell Università e della Ricerca (MUR). 
%  European Union – NextGenerationEU. Award Number: Project code PE0000021, Concession Decree No. 1561 of 11.10.2022 
%  adopted by Ministero dell Università e della Ricerca (MUR), CUP E63C22002160007 - Project title “Network 4 Energy Sustainable Transition – NEST.
%  ************************************************************************************************************************************************


clear all
clc
% function [Output_DR_Type_1,Output_DR_Type_2,Output_DR_Type_3,P_BESS_mean,SoC_mean,P_BESS_set_point] = RecOptim(P_dr_scheduled)

%%% Desrizione 



%% IMPORT INPUT DATA (see documentation on how to provide the excel file 'INPUT.xlsx')

[P_load_input,xl] = xlsread('Input.xlsx','Load');       % Input load power forecasting
[P_gen_input,xg] = xlsread('Input.xlsx','RES');         % Input generation power forecasting
Pr_sh = xlsread('Input.xlsx','Shared_En_Price');        % Purchasing Price
Pr_s = xlsread('Input.xlsx','Selling_Price');           % Selling Price
Pr_p = xlsread('Input.xlsx','Purchasing_Price');        % Purchasing Price
Sim_Data = xlsread('Input.xlsx','Sim_Info');            % Simulation parameters
P_dr_type_I = xlsread('Input.xlsx','DR_Load_type_I');	% Demand response Load type I data
P_dr_type_II = xlsread('Input.xlsx','DR_Load_type_II');	% Demand response Load type II data
P_dr_type_III = xlsread('Input.xlsx','DR_Load_type_III');	% Demand response Load type III data
range_opt = size(xlsread('Output.xlsx','BESS_SoC'));   % length of excel file previous output
load Inputs_MAT % load the values of P_L_dr_II_fixed_next and SoC next: Note that they can be []

% Time interval and iteration details
dt = Sim_Data(1,1);     % Length of the time interval
N = Sim_Data(1,2);      % Number of samples of the MonteCarlo procedure

% Identification of the number of forward control intervals
[nt_l,ct_l] = size(P_load_input);
[nt_g,ct_g] = size(P_gen_input);
[nt_sh,ct_sh] = size(Pr_sh);
[nt_p,ct_p] = size(Pr_p);
[nt_s,ct_s] = size(Pr_s);
nt = min([nt_l,nt_g,nt_p,nt_s,nt_sh]);          % Number of forward control intervals 

if nt<nt_p
    Pr_p=Pr_p(1:nt,1);
end
if nt<nt_s
    Pr_s=Pr_s(1:nt,1);
end
if nt<nt_sh
    Pr_sh=Pr_sh(1:nt,1);
end

P_L_dr_II_fixed=[P_L_dr_II_fixed_next,zeros(1,nt-length(P_L_dr_II_fixed_next))];    % note that length P_L_dr_II_fixed_next cannot be higher than nt 

% IMPORT DATA ON THE BATTERY ENERGY STORAGE SYSTEM

BESS_data = xlsread('Input.xlsx','BESS');
E_b = BESS_data(1,1);                   % battery capacity [kWh]
P_b = BESS_data(2,1);                   % battery power rating [kW]
eta_b = sqrt(BESS_data(3,1)/100);       % battery round trip efficiency
max_dod = BESS_data(4,1)/100;           % battery max Deopth of Discharge
SoC_band = BESS_data(5,1)/100;          % battery SoC available for other services
eta_i = BESS_data(6,1)/100;             % battery converter DC/AC efficiency
P_conv = BESS_data(7,1);                % Converter power rating
SoC_0 = BESS_data(8,1)/100;             % Battery initial SoC
SoH = BESS_data(9,1)/100;               % Battery State oh Health

E_b = E_b*SoH;
Pb = P_b/P_b;
E_max = (1-SoC_band);                 % Maximum energy can be stored
E_min = (1-max_dod)+SoC_band;         % Minimum energy can be stored

if SoC_next == -1                % Initial value of the State of Charge
    E_0 = SoC_0;                      % Value provided in the excel file
else
    E_0 = SoC_0;                   % Value derived from the previous software running
end

%% Scheduling of the DR_Load_Type_I 
if size(P_dr_type_I)==0
    X_dr_I = [];
    P_L_dr_I = zeros(1,nt);
    Output_DR_Type_1 = [];
else
    [r_dr,c_dr] = size(P_dr_type_I);
    n_dr_I = c_dr/2;            % number of DR Load type I % IMPORTANT: INPUT must contain data such that n_dr_I<=nt
    n_t_dr = min([r_dr,nt]);    % number of forward control intervals for DR load_type I % IMPORTANT: INPUT must contain data such that n_dr_I<=nt
    UB = [];
    P_dr = [];
    A_eq = [];
    e_dr = [];
    k=0;
    for h=1:2:n_dr_I*2
        k=k+1;
        P_dr = [P_dr,P_dr_type_I(1:n_t_dr,h)'];     % constant DR Load power    [kW]
        e_dr = [e_dr,P_dr_type_I(1,h+1)];    % required DR Load energy   [kWh]
        UB = [UB,(P_dr_type_I(1:n_t_dr,h)>0)'];
        A_eq = [A_eq;zeros(1,n_t_dr*(k-1)),ones(1,n_t_dr).*P_dr_type_I(1:n_t_dr,h)',zeros(1,n_t_dr*(n_dr_I-k))];
    end
    B_eq = e_dr';
    LB = zeros(1,n_t_dr*n_dr_I);
    
    intcon = [1:n_t_dr*n_dr_I];
    f = repmat(Pr_p(1:n_t_dr,1),n_dr_I,1);
    [X_dr_I,fval1,extflg1] = intlinprog(f,intcon,[],[],A_eq,B_eq,LB,UB);
    if extflg1<=0
        X_dr_I = zeros(size(f));
        disp('WARNING: Constraints imposed on DR Loads Type I are not feasible')
    end
    P_L_dr_I = zeros(1,nt);
    Output_DR_Type_1 = zeros(n_dr_I,nt);
    for h=1:n_dr_I
        P_L_dr_I(1,1:n_t_dr) = P_L_dr_I(1,1:n_t_dr)+(P_dr_type_I(1:n_t_dr,(2*h-1)).*X_dr_I((h-1)*n_t_dr+1:h*n_t_dr,1))';
        Output_DR_Type_1(h,1:n_t_dr) = (P_dr_type_I(1:n_t_dr,(2*h-1)).*X_dr_I((h-1)*n_t_dr+1:h*n_t_dr,1))';
    end
end

%% Scheduling of the RD_Load_Type_II
if size(P_dr_type_II)==0
    X_dr_II = [];
    P_L_dr_II = zeros(1,nt);
    Output_DR_Type_2 = [];
else
    [r_dr2,c_dr2] = size(P_dr_type_II);
    n_dr_II = c_dr2/3;            % number of DR Load type I
    n_t_dr2 = min([r_dr2,nt]);    % number of forward control intervals for DR load_type II % IMPORTANT: INPUT must contain data such that n_dr_II<=nt
    UB2 = [];
    P_dr2 = [];
    A = [];
    B=[];
    q_dr = [];
    p_dr = [];
    k=0;
    for h=1:3:n_dr_II*3
        k=k+1;
        P_dr2 = [P_dr2,P_dr_type_II(1:n_t_dr2,h)'];     % DR time intervals
        q_dr = [q_dr,P_dr_type_II(1,h+1)];              % time limit
        p_dr = [p_dr,P_dr_type_II(1,h+2)];              % constant power [kW]
        UB2 = [UB2,(P_dr_type_II(1:n_t_dr2,h)>0)'];
        for t=1:n_t_dr2-q_dr(1,k)+1
            A = [A;zeros(1,n_t_dr2*(k-1)),zeros(1,t-1),-ones(1,q_dr(1,k)),zeros(1,n_t_dr2-(t-1)-q_dr(1,k)),zeros(1,n_t_dr2*(n_dr_II-k))];
            B = [B;-p_dr(1,k)];
        end
    end
    LB2 = zeros(1,n_t_dr2*n_dr_II);
    
    intcon2 = [1:n_t_dr2*n_dr_II];
    f2 = repmat(Pr_p(1:n_t_dr2,1),n_dr_II,1);
    [X_dr_II,fval2,extflg2] = intlinprog(f2,intcon2,A,B,[],[],LB2,UB2);
    if extflg2<=0
        X_dr_II = zeros(size(f2));
        disp('WARNING: Constraints imposed on DR Loads Type II are not feasible')
    end
    P_L_dr_II = zeros(1,nt);
    Output_DR_Type_2 = zeros(n_dr_II,nt);
    for h=1:n_dr_II
        P_L_dr_II(1,1:n_t_dr2) = P_L_dr_II(1,1:n_t_dr2)+(P_dr_type_II(1:n_t_dr2,(3*h-2)).*X_dr_II((h-1)*n_t_dr2+1:h*n_t_dr2,1))';
        Output_DR_Type_2(h,1:n_t_dr2) = (P_dr_type_II(1:n_t_dr2,(3*h-2)).*X_dr_II((h-1)*n_t_dr2+1:h*n_t_dr2,1))';
    end
end
P_L_dr_II_fixed_next = P_L_dr_II(1,2:nt)+P_L_dr_II_fixed(1,2:nt);

%% DR Load Type III
if size(P_dr_type_III)==0
    X_dr_III = [];
    P_L_dr_III = zeros(1,nt);
    Output_DR_Type_3 = [];
    X03 = [];
    LB3 = [];
    UB3 = [];
    n_dr_III = 0;
    X_dr_III = zeros(size(X03'));
    e_rIII = 0;
    n_t_dr3 = nt;
else
    [r_dr3,c_dr3] = size(P_dr_type_III);
    n_dr_III = c_dr3/3;            % number of DR Load type I % IMPORTANT: INPUT must contain data such that n_dr_I<=nt
    n_t_dr3 = min([r_dr3,nt]);     % number of forward control intervals for DR load_type III % IMPORTANT: INPUT must contain data such that n_dr_III<=nt
    LB3 = [];
    UB3 = [];
    P_dr3_min = [];
    P_dr3_max = [];
    P_dr3_mean = [];
    A3 = [];
    k=0;
    for h=1:3:n_dr_III*3
        k=k+1;
        P_dr3_min = [P_dr3_min,P_dr_type_III(1:n_t_dr3,h)'];     % DR time intervals
        P_dr3_max = [P_dr3_max,P_dr_type_III(1:n_t_dr3,h+1)'];     % DR time intervals
        P_dr3_mean = [P_dr3_mean,(P_dr_type_III(1:n_t_dr3,h)'+P_dr_type_III(1:n_t_dr3,h+1)')/2];
        LB3 = [LB3,P_dr_type_III(1:n_t_dr3,h)'];
        UB3 = [UB3,P_dr_type_III(1:n_t_dr3,h+1)'];
        check_3(k,1) = sum(P_dr_type_III(1:n_t_dr3,h+1))*dt;
        e_rIII(k,1) = P_dr_type_III(1,h+2);
        A3 = [A3;zeros(1,n_t_dr3*(k-1)),-ones(1,n_t_dr3).*P_dr3_mean(1,n_t_dr3*(k-1)+1:n_t_dr3*k),zeros(1,n_t_dr3*(n_dr_III-k))];
    end
    if check_3>=e_rIII
        B3 = -e_rIII;
        P_dr3_mean = (P_dr3_max+P_dr3_min)/2;
        X03 = (LB3+UB3)/2;  % initial solution
        intcon3 = [1:n_t_dr3*n_dr_III];
        f3 = repmat(Pr_p(1:n_t_dr3,1),n_dr_III,1);
        f3 = f3.*P_dr3_mean';
        [X_dr_III,fval3,extflg3] = intlinprog(f3,intcon3,A3,B3,[],[],zeros(size(LB3)),(UB3>0)*1);
        if extflg3 <=0
            X_dr_III = ones(size(X03'));
        end
    else
        disp('WARNING: Constraints imposed on DR Loads Type III are not feasible')
        X03 = [];
        LB3 = [];
        UB3 = [];
        n_dr_III = 0;
        X_dr_III = zeros(size(X03'));
    end
end

%% Evaluation of RES and Total Load profiles

%  Load

if xl(2,1)=="Gaussian"
    for t=1:nt
        P_load(t,:) = normrnd(P_load_input(t,1),P_load_input(t,2),[1,N]);
    end
end
P_load = P_load+P_L_dr_I'+P_L_dr_II'+P_L_dr_II_fixed';
 
%  Generation

if xg(2,1)=="Gaussian"
    for t=1:nt
        P_gen(t,:) = normrnd(P_gen_input(t,1),P_gen_input(t,2),[1,N]);
    end
end

%% REC PROBABILISTIC OPTIMAL OPERATION
e_rIII = e_rIII/P_b;
X0 = [ones(1,round(nt/2)),-ones(1,nt-round(nt/2)),X03];   % Initial solution
lb = [-ones(1,nt),LB3/P_b.*X_dr_III'];
ub = [ones(1,nt),UB3/P_b.*X_dr_III'];
X(1,:) = X0*0;      % In case no feasible solution is found, null-vector is provided as solution
n_iter_aux = 0;
for n_iter=1:N
    for t=1:nt
        P_l(t,1)=P_load(t,n_iter)/P_b;
        P_g(t,1)=P_gen(t,n_iter)/P_b;
    end
    if n_iter>1
        x0=x;
    end
      
    options = optimset('Display','off','Diagnostic','off','TolX',1e-9,'TolCon',1e-9,'TolFun',1e-6,'algorithm','active-set');
    [x,fval,exit_flag(n_iter)] = fmincon('RecObj',X0,[],[],[],[],lb,ub,'RecCon',options,P_l,P_g,nt,eta_b,eta_i,dt,E_max,E_min,E_b,P_b,Pr_sh,Pr_s,Pr_p,E_0,e_rIII,n_dr_III,n_t_dr3);
    if exit_flag(n_iter)>0
        n_iter_aux = n_iter_aux+1;
        X(n_iter_aux,:) = x; % pu
        E_1 =-(x(1:nt)/eta_i.*eta_b.*(x(1:nt)<0)+x(1:nt)/eta_i/eta_b.*(x(1:nt)>0))*dt*P_b/E_b;  % SoC evaluation for all forward time intervals
        E_11 = [E_0,E_1];
        SoC(n_iter_aux,:) = cumsum(E_11);
    end
end

%% OUTPUT Elaboration
P_BESS_mean = mean(X(:,1:nt))*P_b;
P_BESS_std = std(X(:,1:nt))*P_b;
SoC_mean = mean(SoC);
SoC_std = std(SoC);
SoC_next = SoC_mean(1,2);   % SoC initial value for the next running of the software

% Inputs_MAT
save Inputs_MAT P_L_dr_II_fixed_next SoC_next

Output_DR_Type_3 = zeros(n_dr_III*2,nt);

k=0;
for h = 1:n_dr_III
    k=k+1;
    Output_DR_Type_3(k,1:n_t_dr3) = round(mean(X(:,nt+(h-1)*n_t_dr3+1:nt+h*n_t_dr3))*P_b,3);    % mean values
    k=k+1;
    Output_DR_Type_3(k,1:n_t_dr3) = round(std(X(:,nt+(h-1)*n_t_dr3+1:nt+h*n_t_dr3))*P_b,3);      % standard devation
end

Excel = actxserver('Excel.Application');
d=dir('Output.xlsx');
Workbook = Excel.Workbooks.Open(fullfile(d.folder,d.name));
WorkSheets = Excel.sheets;

TargetSheet = get(WorkSheets,'item','DR_Load_Type I');   
Activate(TargetSheet);
Excel.Range('A2:AD1001').Select;     % It is recomended to increase the number of columns in case of more than 30 loads and the number of rows in case of more than 1000 time intervals
Excel.Selection.Clear;

TargetSheet = get(WorkSheets,'item','DR_Load_Type II');
Activate(TargetSheet);
Excel.Range('A2:AD1001').Select;     % It is recomended to increase the number of columns in case of more than 30 loads and the number of rows in case of more than 1000 time intervals
Excel.Selection.Clear;

TargetSheet = get(WorkSheets,'item','DR_Load_Type III');
Activate(TargetSheet);
Excel.Range('A2:BH1001').Select;     % It is recomended to increase the number of columns in case of more than 30 loads and the number of rows in case of more than 1000 time intervals
Excel.Selection.Clear;

TargetSheet = get(WorkSheets,'item','BESS_SoC');
Activate(TargetSheet);
Excel.Range('A3:D1002').Select;     % It is recomended to increase the number of rows in case of more than 1000 time intervals
Excel.Selection.Clear;

Excel.ActiveWorkbook.Save;
Excel.ActiveWorkbook.Close(false);
delete(Excel); clear Excel

xlswrite('Output.xlsx',P_BESS_mean(1),'BESS_power','A2')
xlswrite('Output.xlsx',P_BESS_mean','BESS_SoC','A3')
xlswrite('Output.xlsx',P_BESS_std','BESS_SoC','B3')
xlswrite('Output.xlsx',SoC_mean','BESS_SoC','C3')
xlswrite('Output.xlsx',SoC_std','BESS_SoC','D3')

xlswrite('Output.xlsx',Output_DR_Type_1','DR_Load_Type I','A2')
xlswrite('Output.xlsx',Output_DR_Type_2','DR_Load_Type II','A2')
xlswrite('Output.xlsx',Output_DR_Type_3','DR_Load_Type III','A2')