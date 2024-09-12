%% Doyle-Fuller-Newman Model
%   Published June 14, 2016 by Professor Scott Moura
%   Energy, Controls, and Applications Lab (eCAL)
%   University of California, Berkeley
%   http://ecal.berkeley.edu/
clc;
clear;
tic;
delete('results*.mat');
disp('Fast DFN')
disp('%%%%%%%%')

%% Electrochemical Model Parameters
% Load Lithium Cobolt Oxide Params, adopted from DUALFOIL
run param/params_LCO
%run params_NMC_Samsung_new_iteration

%% Input charge/discharge Current Data %%
% % Current | Positive <=> Discharge, Negative <=> Charge

% Calculate C-rate in terms of [A/m^2] using low/high voltage cutoffs
[cn_low,cp_low] = init_cs(p,p.volt_min);
[cn_high,cp_high] = init_cs(p,p.volt_max);
Delta_cn = cn_high-cn_low;
Delta_cp = cp_low-cp_high;
p.OneC = min(p.epsilon_s_n*p.L_n*Delta_cn*p.Faraday/3600, p.epsilon_s_p*p.L_p*Delta_cp*p.Faraday/3600);


% %exponential  discharge, constant discharge, exp. charge
% 
% start_time = 1;
% end_time = 40;
% time_decay = linspace(start_time, end_time, end_time - start_time + 1);
% decay_rate = -log(1/5) / (end_time - start_time); % To make it reach 1C after 40 seconds
% I_decay = 5 * p.OneC * exp(-decay_rate * (time_decay - start_time));
% % Constant current at 1C
% I_const = ones(1, 80) * p.OneC;
% 
% % Exponential increase at the end
% start_time = 121;
% end_time = 181;
% time_rise = linspace(start_time, end_time, end_time - start_time + 1);
% rise_rate = log(5) / (end_time - start_time); % To make it reach 5C after 60 seconds
% I_rise = p.OneC * exp(rise_rate * (time_rise - start_time));
% 
% % Combine the three parts of the current profile
% I = [I_decay, I_const, I_rise];

%I = -1*p.OneC*ones(size(t));

% I(11:40) = 5*p.OneC;
% 
% I((40+91):(40+90+30)) = -5*p.OneC;

%I = 5*p.OneC*ones(size(t));

%%%%%%%%%%%%%%% DYNAMIC CHARGE/DISCHARGE CYCLES FROM EXPERIMENTS %%%%%%%%%%%%%%%
% load('data/DC1_batt_ObsData');
% i_custom=I*10;
% 
% 
% clear I;

%%%%%%%%%%%%%%% MANUAL INPUT WITH C-RATE %%%%%%%%%%%%%%%%%%%%%%%%%
p.delta_t = 1;
tmax=5000; %this sets the size of .mat files regularly flushed to disk. also preallocates arrays to tmax
t = 0:p.delta_t:tmax;
NT = length(t);
CTRL=true; %Controller == ON

% amplitude = p.OneC;     
% frequency = 0.01;     
% sampleRate = 1;  
% duration = 300;      
% t = 1:(1/sampleRate):(duration);
% i = -amplitude * sin(2 * pi * frequency * t);

%% Script for Determining OCV

chargemode='charge'; %Input: Charge or Discharge?
V0 = p.volt_max-0.4;%4.3;%4.35; %4.42 for 99.5% SOC volt_exp(1); Modify this to change initial SOC
%% Initial Conditions & Preallocation
% Solid concentration

[csn0,csp0] = init_cs(p,V0);

% Electrolyte concentration
ce0 = p.c_e;

% Temperature
T0 = p.T_amb;

% Vector lengths
Ncsn = p.PadeOrder * (p.Nxn-1);
Ncsp = p.PadeOrder * (p.Nxp-1);
Nce = p.Nx - 3;
Nc = Ncsn+Ncsp+Nce;
Nn = p.Nxn - 1;
Ns = p.Nxs - 1;
Np = p.Nxp - 1;
Nnp = Nn+Np;
Nx = p.Nx - 3;
Nz = 3*Nnp + Nx;

% Output Discretization params
disp('Discretization Params');
fprintf(1,'No. of FDM nodes in Anode | Separator | Cathode : %1.0f | %1.0f | %1.0f\n',p.Nxn,p.Nxs,p.Nxp);
fprintf(1,'Order of Pade Approx for Solid Concentration : %1.0f\n',p.PadeOrder);
fprintf(1,'Time Step : %2.2f sec\n',p.delta_t);
disp(' ');

c_s_n0 = zeros(p.PadeOrder,1);
c_s_p0 = zeros(p.PadeOrder,1);

%%%%% Initial condition based on Jordan form
c_s_n0(3) = csn0;
c_s_p0(3) = csp0;

c_s_n = zeros(Ncsn,NT);
c_s_p = zeros(Ncsp,NT);

c_s_n(:,1) = repmat(c_s_n0, [Nn 1]);
c_s_p(:,1) = repmat(c_s_p0, [Nn 1]);

% Electrolyte concentration
c_e = zeros(Nx,NT);
c_e(:,1) = ce0 * ones(Nx,1);

c_ex = zeros(Nx+4,NT);
c_ex(:,1) = c_e(1,1) * ones(Nx+4,1);

% Temperature
T = zeros(NT,1);
T(1) = T0;

% Current
I = zeros(NT,1);
I(1) = -0.1*p.OneC;
I(2) = -0.1*p.OneC;

% Solid Potential
Uref_n0 = refPotentialAnode(p, csn0(1)*ones(Nn,1) / p.c_s_n_max);
Uref_p0 = refPotentialCathode(p, csp0(1)*ones(Np,1) / p.c_s_p_max);

phi_s_n = zeros(Nn,NT);
phi_s_p = zeros(Np,NT);
phi_s_n(:,1) = Uref_n0;
phi_s_p(:,1) = Uref_p0;

% Electrolyte Current
i_en = zeros(Nn,NT);
i_ep = zeros(Np,NT);

% Electrolyte Potential
phi_e = zeros(Nx+2,NT);

% Molar Ionic Flux
jn = zeros(Nn,NT);
jp = zeros(Np,NT);

% Surface concentration
c_ss_n = zeros(Nn,NT);
c_ss_p = zeros(Np,NT);
c_ss_n(:,1) = repmat(csn0, [Nn 1]);
c_ss_p(:,1) = repmat(csp0, [Np 1]);

% Volume average concentration
c_avg_n = zeros(Nn,NT);
c_avg_p = zeros(Np,NT);
c_avg_n(:,1) = repmat(csn0, [Nn 1]);
c_avg_p(:,1) = repmat(csp0, [Np 1]);

% SOC (Bulk Anode SOC)
SOC = zeros(NT,1);
SOC(1) = (mean(c_avg_n(:,1)) - cn_low) / (cn_high - cn_low);

% Overpotential
eta_n = zeros(Nn,NT);
eta_p = zeros(Np,NT);

% Constraint Outputs
c_e_0n = zeros(NT,1);
c_e_0n(1) = c_ex(1,1);

c_e_0p = zeros(NT,1);
c_e_0p(1) = c_ex(end,1);

eta_s_Ln = zeros(NT,1);
eta_s_Ln(1) = phi_s_p(1,1) - phi_e(1,1);

% Voltage
Volt = zeros(NT,1);
Volt(1) = phi_s_p(end,1) - phi_s_n(1,1) - p.R_c*I(1);

% Conservation of Li-ion matter
n_Li_s = zeros(NT,1);
n_Li_e = zeros(NT,1);

n_Li_e(1) = sum(c_e(1:Nn,1)) * p.L_n*p.delta_x_n * p.epsilon_e_n * p.Area ...
     + sum(c_e(Nn+1:end-Np,1)) * p.L_s*p.delta_x_s * p.epsilon_e_s * p.Area ...
     + sum(c_e(end-Np+1:end,1)) * p.L_p*p.delta_x_p * p.epsilon_e_p * p.Area;

% Stats
newtonStats.iters = zeros(NT,1);
newtonStats.relres = cell(NT,1);
newtonStats.condJac = zeros(NT,1);

% Initial Conditions
x0 = [c_s_n(:,1); c_s_p(:,1); c_e(:,1); T(1)];

z0 = [phi_s_n(:,1); phi_s_p(:,1); i_en(:,1); i_ep(:,1);...
      phi_e(:,1); jn(:,1); jp(:,1)];

%% Preallocate states
x = zeros(length(x0), NT);
z = zeros(length(z0), NT);

x(:,1) = x0;
z(:,1) = z0;

%% Precompute data
% % Solid concentration matrices
% [A_csn,B_csn,A_csp,B_csp,C_csn,C_csp,A_csn_normalized, A_csp_normalized] = c_s_mats(p);
% p.A_csn = A_csn;
% p.A_csn_normalized= A_csn_normalized;
% p.B_csn = B_csn;
% p.A_csp = A_csp;
% p.A_csp_normalized=A_csp_normalized;
% p.B_csp = B_csp;
% p.C_csn = C_csn;
% p.C_csp = C_csp;
% 
% clear A_csn B_csn A_csp B_csp C_csn C_csp A_csn_normalized A_csp_normalized;

% Adjust Temperature Dependent Parameters, based on present temperaure
% Solid concentration matrices
p.D_s_n = p.D_s_n0 * exp(p.E.Dsn/p.R*(1/p.T_ref - 1/T(1))); % ATT:
%commented out for NMC samsung
p.D_s_p = p.D_s_n0 * exp(p.E.Dsp/p.R*(1/p.T_ref - 1/T(1)));

[A_csn,B_csn,A_csp,B_csp,C_csn,C_csp] = c_s_mats(p);
p.A_csn = A_csn;
p.B_csn = B_csn;
p.A_csp = A_csp;
p.B_csp = B_csp;
p.C_csn = C_csn;
p.C_csp = C_csp;

clear A_csn B_csn A_csp B_csp C_csn C_csp;

% Electrolyte concentration matrices
[M1n,M2n,M3n,M4n,M5n, M1s,M2s,M3s,M4s, M1p,M2p,M3p,M4p,M5p, C_ce] = c_e_mats(p);

p.ce.M1n = M1n;
p.ce.M2n = M2n;
p.ce.M3n = M3n;
p.ce.M4n = M4n;
p.ce.M5n = M5n;

p.ce.M1s = M1s;
p.ce.M2s = M2s;
p.ce.M3s = M3s;
p.ce.M4s = M4s;

p.ce.M1p = M1p;
p.ce.M2p = M2p;
p.ce.M3p = M3p;
p.ce.M4p = M4p;
p.ce.M5p = M5p;

p.ce.C = C_ce;

rM3 = [Nn; Ns; Np];
cM3 = rM3';
p.ce.M3 = sparse(blkdiagFast(rM3, cM3, p.ce.M3n, p.ce.M3s, p.ce.M3p));

clear M1n M2n M3n M4n M5n M1s M2s M3s M4s M1p M2p M3p M4p M5p C_ce rM1 cM1;

% Solid Potential
[F1_psn,F1_psp,F2_psn,F2_psp,G_psn,G_psp,...
    C_psn,C_psp,D_psn,D_psp] = phi_s_mats(p);
p.F1_psn = F1_psn;
p.F1_psp = F1_psp;
p.F2_psn = F2_psn;
p.F2_psp = F2_psp;
p.G_psn = G_psn;
p.G_psp = G_psp;
p.C_psn = C_psn;
p.C_psp = C_psp;
p.D_psn = D_psn;
p.D_psp = D_psp;

clear F1_psn F1_psp F2_psn F2_psp G_psn G_psp C_psn C_psp D_psn D_psp;

% Electrolyte Current
[F1_ien,F1_iep,F2_ien,F2_iep,F3_ien,F3_iep] = i_e_mats(p);
p.F1_ien = F1_ien;
p.F1_iep = F1_iep;
p.F2_ien = F2_ien;
p.F2_iep = F2_iep;
p.F3_ien = F3_ien;
p.F3_iep = F3_iep;

clear F1_ien F1_iep F2_ien F2_iep F3_ien F3_iep;

% Electrolyte Potential
p.M1_pen_skel = sparse(diag(ones(p.Nxn-2,1),1) + diag(-ones(p.Nxn-2,1),-1));
p.M1_pes_skel = sparse(diag(ones(p.Nxs-2,1),1) + diag(-ones(p.Nxs-2,1),-1));
p.M1_pep_skel = sparse(diag(ones(p.Nxp-2,1),1) + diag(-ones(p.Nxp-2,1),-1));

[M1_pe,M2_pe,M3_pe,M4_pe,C_pe] = phi_e_mats(p);
p.M1_pe = M1_pe;
p.M2_pe = M2_pe;
p.M3_pe = M3_pe;
p.M4_pe = M4_pe;
p.C_pe = C_pe;

clear M1_pe M2_pe M3_pe M4_pe C_pe

% Jacobian
[f_x, f_z, g_x, g_z] = jac_dfn_pre(p);
p.f_x = f_x;
p.f_z = f_z;
p.g_x = g_x;
p.g_z = g_z;
clear f_x f_z g_x g_z

%%CC-CV Controller, PID Initialisation 

Kp_voltage = 0.3*125.4;  %0.3 to avoid oscillations buidling up at the end of CV phase
Ki_voltage = 0.1*125.4;  %change integrator (too agressive)
Kd_voltage = 0.3*31.35;

integral_term = 0;
previous_error = 0;
accumulated_error = 0;

if strcmp(chargemode, 'charge')
    constant_current = -1*p.OneC;
end
if strcmp(chargemode, 'discharge')
    constant_current = 1*p.OneC;
end
desired_voltage_c = p.volt_max-0.22; %4.5; T=25°C:0.22 T=0°C:0.18 T=-20°C: +0.07 T=10°C: -0.1
desired_voltage_d = p.volt_min+0.2;%3.5;

mode= "CC"; %ATT: before: CC
modename='CC1'; %ATT: before: CC1
%desiredSOC=0.15;

% Full Model Validation Test
C_rates_val=[0.1,0.2,0.5,1,2]*p.OneC;
C_rates_val=flip(C_rates_val);
C_rates_val_start=-1*p.OneC;
k_val=0;
DDP=[zeros(1,10)+1.5, zeros(1,20)+0.4,zeros(1,5)-0.8,zeros(1,30)]*p.OneC;
k_DDP=1;

% UDDS Test

data = readtable('UDDS.csv', 'VariableNamingRule', 'preserve');
UDDS_t = data.("# Time [s]"); 
UDDS_i = data.("Current [A]")*12.706*2; %scaled by 12.706 so that I_max = 2C 
UDDS_indices = abs(UDDS_i) > 2*p.OneC;
UDDS_i(UDDS_indices) = UDDS_i(UDDS_indices) / 2;
k_UDDS=1;

%Relaxation Behaviour investigation
C_rates_relax=[0.5,1,2,4];
k_relax=1;

%define HPPC cycle
%C_rates=[0.5,1,1.5,2];
C_rates=[0.05,0.1,0.15,0.2,0.5,1,1.5,2];
HPPC_voltage_increment=0.05; %reduce SOC in 5% steps in between HPPC cycles
HPPC_cycle_count=0; %keeps track of current HPPC cycle count for later storage of data
HPPC_current=[];
HPPC_cycle=[ones(1,10)*-p.OneC,zeros(1,40),ones(1,10)*p.OneC,zeros(1,40)];
for C_rate = C_rates
    HPPC_current=[HPPC_current,HPPC_cycle*C_rate];
end
T_lim=p.T_amb+0.5;
%T_lim=26+273.15; %Cooldown Temp
HPPC_SOC_levels=[99,98,97,96,95,90,85,80,75,70,65,60,55,50,45,40,35,30,25,20,15,10,5,4,3,2,1]/100;
%HPPC_SOC_levels=[90:-10:10]/100;
OCV_levels=[99.5:-0.5:95.5,95:-1:5,4.5:-0.5:0.5]/100;%,2:1:4,5:5:95,96:1:99]/100;
%OCV_levels=[99:-1:97,96:1:99]/100;
k_OCV=1;
voltage_exceed_count=1;
current_reduction_increment=0.25*p.OneC;
v_min=false;
v_min_firstoccurence=false;
%%
C_rates_extreme={[0.05,0.1,0.15,0.2,0.5,1,1.4,1.5,2],[0.05,0.1,0.15,0.2,0.5,1,1.5,1.6,2],[0.05,0.1,0.15,0.2,0.5,1,1.5,1.8,2],[0.05,0.1,0.15,0.2,0.5,1,1.5,1.8,2],[0.05,0.1,0.15,0.2,0.5,1,1.1,1.5,2],[0.05,0.1,0.15,0.2,0.5,1,1.5,2]};
C_ex=[1.4,1.6,1.8,1.8,1.1,0.5]; %Max C for 99,98,97,3,2,1 SOC
%% just for CHEN at 0 deg C
% C_rates_extreme={[0.05,0.1,0.15,0.2,0.5,1,1.4,1.5,2],[0.05,0.1,0.15,0.2,0.5,1,1.5,1.6,2],[0.05,0.1,0.15,0.2,0.5,1,1.5,1.8,2],[0.05,0.1,0.15,0.2,0.5,1,1.5,2],[0.05,0.1,0.15,0.2,0.5,1,1.5,2],[0.05,0.1,0.15,0.2,0.5,1,1.5,2]};
% C_ex=[1.4,1.6,1.8,1,0.5,0.2]; %Max C for 99,98,97,3,2,1 SOC
%%
HPPC_current_extreme={};
for i = 1:length(C_rates_extreme)
    HPPC_current_extreme{i}=[];
    for C_rate = C_rates_extreme{i}
        if C_rate <= C_ex(i)
            HPPC_current_extreme{i}=[HPPC_current_extreme{i},HPPC_cycle*C_rate];
        elseif i <=  3 %High SOC
            HPPC_current_extreme{i}=[HPPC_current_extreme{i},ones(1,ceil((C_rate/C_ex(i))*10*2))*-p.OneC/2*C_ex(i),zeros(1,40),ones(1,10)*p.OneC*C_rate,zeros(1,40)];
        elseif i > 3 %Low SOC
            HPPC_current_extreme{i}=[HPPC_current_extreme{i},ones(1,10)*-p.OneC*C_rate,zeros(1,40),ones(1,ceil((C_rate/C_ex(i))*10*2))*p.OneC/2*C_ex(i),zeros(1,40)];    
        end
    end
end
k_HPPC_ex=1; %count at which of the extreme SOCs we are at
k_HPPC_cycle=1;%count at which of the HPPC cycles we are at
% figure()
% %hold on
% for i = 1:length(HPPC_current_extreme)
%     HPPC_current_ex = HPPC_current_extreme{i};
%     plot(1:length(HPPC_current_ex), HPPC_current_ex/p.OneC)
%     grid on
% end



%%
k_length_HPPC=length(HPPC_current);
HPPC_voltage=zeros(floor(1/HPPC_voltage_increment),length(HPPC_current)); %create a "number of HPPC cycles" x "length of HPPC cycles" matrix to store the data for later identification
HPPC_startingOCV=zeros(floor(1/HPPC_voltage_increment),1);
%plot(1:length(HPPC_current),HPPC_current/p.OneC);
%% Integrate!
disp('Simulating DFN Model...');
k_full=1;
k=1;
last_saved_batch=0;
k_end=0; %placeholder to initiate variable
constant_current=C_rates_val_start;
while CTRL == true
    
    % Current
    if(k == 1)
        Cur_vec = [I(k), I(k), I(k+1)];
    else
        Cur_vec = [I(k-1), I(k), I(k+1)];
    end
    
    % Adjust Temperature Dependent Parameters, based on present temperaure
    %ATT: commented out: NMC 
    p.k_n = p.k_n0 * exp(p.E.kn/p.R*(1/p.T_ref - 1/T(k)));
    p.k_p = p.k_p0 * exp(p.E.kp/p.R*(1/p.T_ref - 1/T(k)));
    %ATT: commented out: NMC
    % Solid concentration matrices
    p.D_s_n = p.D_s_n0 * exp(p.E.Dsn/p.R*(1/p.T_ref - 1/T(k)));
    p.D_s_p = p.D_s_p0 * exp(p.E.Dsp/p.R*(1/p.T_ref - 1/T(k)));
    
    [A_csn,B_csn,A_csp,B_csp,C_csn,C_csp] = c_s_mats(p);
    p.A_csn = A_csn;
    p.B_csn = B_csn;
    p.A_csp = A_csp;
    p.B_csp = B_csp;
    p.C_csn = C_csn;
    p.C_csp = C_csp;
    
    clear A_csn B_csn A_csp B_csp C_csn C_csp;
    
    % Step-forward in time
    [x(:,k+1), z(:,k+1), stats] = cn_dfn(x(:,k),z(:,k),Cur_vec,p);

    % Parse out States
    c_s_n(:,k+1) = x(1:Ncsn, k+1);
    c_s_p(:,k+1) = x(Ncsn+1:Ncsn+Ncsp, k+1);
    c_e(:,k+1) = x(Ncsn+Ncsp+1:Nc, k+1);
    T(k+1) = x(end, k+1);
    
    phi_s_n(:,k+1) = z(1:Nn, k+1);
    phi_s_p(:,k+1) = z(Nn+1:Nnp, k+1);
    i_en(:,k+1) = z(Nnp+1:Nnp+Nn, k+1);
    i_ep(:,k+1) = z(Nnp+Nn+1:2*Nnp, k+1);
    phi_e(:,k+1) = z(2*Nnp+1:2*Nnp+Nx+2, k+1);
    jn(:,k+1) = z(2*Nnp+Nx+3:2*Nnp+Nx+Nn+2, k+1);
    jp(:,k+1) = z(2*Nnp+Nx+Nn+3:end, k+1);
    
    newtonStats.iters(k+1) = stats.iters;
    newtonStats.relres{k+1} = stats.relres;
%     newtonStats.condJac(k+1) = stats.condJac;
    
    % Output data
    [~, ~, y] = dae_dfn(x(:,k+1),z(:,k+1),I(k+1),p);
    
    c_ss_n(:,k+1) = y(1:Nn);
    c_ss_p(:,k+1) = y(Nn+1:Nnp);
    
    c_avg_n(:,k+1) = y(Nnp+1:Nnp+Nn);
    c_avg_p(:,k+1) = y(Nnp+Nn+1 : 2*Nnp);
    SOC(k+1) = (mean(c_avg_n(:,k+1)) - cn_low) / (cn_high - cn_low);
    
    c_ex(:,k+1) = y(2*Nnp+1:2*Nnp+Nx+4);
    
    eta_n(:,k+1) = y(2*Nnp+Nx+4+1 : 2*Nnp+Nx+4+Nn);
    eta_p(:,k+1) = y(2*Nnp+Nx+4+Nn+1 : 2*Nnp+Nx+4+Nn+Np);
    
    c_e_0n(k+1) = y(end-5);
    c_e_0p(k+1) = y(end-4);
    eta_s_Ln(k+1) = y(end-3);
    
    Volt(k+1) = y(end-2);
    if strcmp(mode, 'HPPC')  % this is here so that during current can first be updated, before Volt is calculated. otherwise the whole matrix is shifted one element to the right
        %ATT
        HPPC_voltage(HPPC_cycle_count,k_HPPC)=Volt(k+1);
        if k_HPPC == 1
            HPPC_startingOCV(HPPC_cycle_count)=Volt(k); %store the OCV just before HPPC testing for later averaging
            HPPC_startingSOC(HPPC_cycle_count)=SOC(k);
        end
        k_HPPC=k_HPPC + 1;
    end
    %CC-CV charge controller: start  

    measured_voltage = Volt(k+1);

    if  strcmp(modename, 'CC1') && measured_voltage >=  desired_voltage_c
        mode = 'CV';
        modename= 'CV1';
        
    end
 
    if strcmp(modename, 'CV1') && (abs(I(k)) <= (1/20)*p.OneC || SOC(k+1)>1)  %Rest Phase 
        mode = 'Rest';
        modename='Rest1'; 
        
        k_rest_end1=k_full+500; %ATT 
        %k_end=k_full+7200; %ATT
        k_val=k_val+1;
    end
    
    % Cycle 2: discharge

    if strcmp(modename, 'Rest1')  && k_full >= k_rest_end1 && T(k+1)<T_lim %ATT : Rest1
         integral_term = 0;
         mode = 'CC';
         modename='CC2';
         constant_current = C_rates_val(k_val);

    end

    if strcmp(modename, 'CC2')  && measured_voltage  <=  desired_voltage_d 
        mode = 'CV';
        modename ='CV2';
    end 

    if strcmp(modename, 'CV2')  && (abs(I(k)) <= (1/20)*p.OneC || SOC(k+1)<0)
        mode = 'Rest';
        modename='Rest2';
        k_rest_end2=k_full+500; 
    end

    %Cycle 3: charge

    if strcmp(modename, 'Rest2')  && k_full >= k_rest_end2 && T(k+1)<T_lim
        integral_term = 0;
        mode = 'CC';
        modename='CC1'; %ATT goes back to CC1
        if k_val==length(C_rates_val)
            modename='CC3';
        end
        constant_current = -C_rates_val(k_val);

    end

    if strcmp(modename, 'CC3')  && measured_voltage >=  desired_voltage_c
        mode = 'CV';
        modename ='CV3';
    end 

    if strcmp(modename, 'CV3')  && (abs(I(k)) <= (1/20)*p.OneC || SOC(k+1)>1)
        mode = 'Rest';
        modename='Rest3';
        k_rest_end3=k_full+500;%+7200; 
        k_val=k_val+1;
    end

    % Max C-Rates Testing

    if strcmp(modename, 'Rest1ATT') && T(k+1)< T_lim
        desired_SOC=0.04;
        mode = 'CC';
        modename= 'charge C-Test';
        constant_current=-0.1*p.OneC;
    end
    
    if strcmp(modename, 'charge C-Test') && SOC(k+1) >= desired_SOC
        mode = 'Rest';
        modename ='Rest C-Test';
    end

    if strcmp(modename, 'Rest C-Test') && T(k+1) < T_lim
        mode = 'HPPC';
        modename = 'C-maxing';
        k_HPPC=1;
    end
    
    % Relaxation behaviour investigation
    
    if strcmp(modename, 'RestATT') && k_full == k_rest_end1 %ATT
         integral_term = 0;
         mode = 'CC';
         modename='CC';
         constant_current = p.OneC*C_rates_relax(k_relax);

    end
    
    if strcmp(modename, 'CC') && Volt(k+1)<desired_voltage_d
        mode = 'Rest';
        modename = 'Rest_relax';
        k_rest_relax = k_full +3600;
        if k_relax == length(C_rates_relax)
            k_end= k_full+3600;
            modename= 'End';
        end
        k_relax=k_relax+1;
    end
    
    if strcmp(modename, 'Rest_relax')  && k_full == k_rest_relax 
        mode = 'CC';
        modename='CC_after_relax';
        constant_current = -1*p.OneC;

    end

    if strcmp(modename, 'CC_after_relax')  && measured_voltage >=  desired_voltage_c
        mode = 'CV';
        modename ='CV_after_relax';
    end 

    if strcmp(modename, 'CV_after_relax')  && abs(I(k)) <= (1/20)*p.OneC
        mode = 'Rest';
        modename='Rest1';
        k_rest_end1=k_full+720;%+7200; 
    end
    

    % OCV testing

    if strcmp(modename, 'Rest1ATT') && k_full >= k_rest_end1 && T(k+1)<T_lim
        mode ='OCV';
        modename ='OCV';
        desiredSOC=OCV_levels(k_OCV);
    end
    
    tol=0.0004;
    if strcmp(modename, 'OCV') && abs(SOC(k+1)-desiredSOC) <= tol
        mode ='Rest';
        modename ='RestOCV';
        if k_OCV<length(OCV_levels)
            k_rest_OCV=k_full+100;%3600;
        end
        if k_OCV==length(OCV_levels)
            k_end=k_full+1000;%+3600;
            k_rest_OCV=k_end+10;
        end
        k_OCV=k_OCV+1;
        k_OCV_Volt=k_full+101; %make sure the if condition if OCV_Volt isnt checked  before 100 elements are here
        OCV_Volt=[]; %special Volt vector, not cleared every 5000 steps, to allow for Volt(k-100) check
        k_OCV_count=1;
    end

    if strcmp(modename,'RestOCV') && k_full >= k_rest_OCV && T(k+1)<T_lim && k_full >=k_OCV_Volt %k_OCV_count is basically überflüssig right now.
        if abs(OCV_Volt(k_OCV_count-100)-Volt(k))< 0.0005
            mode = 'OCV';
            modename= 'OCV';
            desiredSOC=OCV_levels(k_OCV);
        end
    end

    % DDP Testing

    if strcmp(modename, 'Rest3ATT') && T(k+1)< T_lim && k_full>=k_rest_end3 %%k_full == k_rest_end1 %ATT : 'OCV') && k_full == k_rest_end1
        mode='DDP';
        modename='DPP';
        desiredSOC=HPPC_SOC_levels(k_HPPC_cycle);
        %if Volt(k+1) < Volt_min 
    end

    % UDDS Testing

    if strcmp(modename, 'Rest3') && T(k+1)< T_lim && k_full>=k_rest_end3 %%k_full == k_rest_end1 %ATT : 'OCV') && k_full == k_rest_end1
        mode='UDDS';
        modename='UDDS';
    end

    % HPPC testing every 5% SOC + 1% at edges

    if strcmp(modename, 'Rest1ATT') && T(k+1)< T_lim && k_full>=k_rest_end1 %%k_full == k_rest_end1 %ATT : 'OCV') && k_full == k_rest_end1
        mode='SOC1';
        modename='SOC1';
        desiredSOC=HPPC_SOC_levels(k_HPPC_cycle);
        %desiredSOC=SOC(k+1)-HPPC_voltage_increment; %reduce SOC in 5% increments
        %desiredSOC = round(desiredSOC * 100) / 100; %takes care of slippage
        % if SOC(k+1)>0.99 
        %     desiredSOC=0.99;
        % end
        %k_end=k_rest_end1+500;         
    end

    if strcmp(modename, 'SOC1') && SOC(k+1) <= desiredSOC
        mode = 'Rest';
        modename ='RestSOC1';
        k_restSOC_end=k_full+500; %Rest for 1h
    end

    if strcmp(modename, 'RestSOC1') && T(k+1) < T_lim && k_full >= k_restSOC_end %k_full == k_restSOC_end
        mode='HPPC';
        modename='HPPC';
        k_HPPC=1;
        if SOC(k+1)>0.965 || SOC(k+1)<0.035
            if k_HPPC_ex > length(HPPC_current_extreme)
                k_HPPC_ex = length(HPPC_current_extreme);
            end
            k_HPPC_end=k_full+length(HPPC_current_extreme{k_HPPC_ex});
        else
            k_HPPC_end=k_full+k_length_HPPC;
        end
        HPPC_cycle_count = HPPC_cycle_count + 1;
    end

    if strcmp(modename, 'HPPC') && k_full==k_HPPC_end
        mode = 'SOC1';
        modename = 'SOC1';
        k_HPPC_cycle=k_HPPC_cycle+1;
        if SOC(k+1)>0.965 || SOC(k+1)<0.035
                 k_HPPC_ex=k_HPPC_ex+1;
        end
        if k_HPPC_cycle <= length(HPPC_SOC_levels)
            desiredSOC=HPPC_SOC_levels(k_HPPC_cycle);
        else
            k_end=k_full;
        end
    % if desiredSOC <= 0.01 % stop HPPC testing under 5% SOC
    %     k_end=k_full;
    end
    
    
    %%Stop condition

    
    if ~strcmp(chargemode, 'none')
        error_c = measured_voltage - desired_voltage_c;
        error_d = measured_voltage - desired_voltage_d;
        if abs(error_c) < abs(error_d)
            error = error_c; % error_c is closer to zero
        else
            error = error_d; % error_d is closer or equal to zero
        end
        derivative_term = Kd_voltage * (error - previous_error) / p.delta_t;
        integral_term = integral_term + Ki_voltage * error * p.delta_t;    
        previous_error = error;
    end

    if strcmp(mode,'CC')
        I(k+1:k+2)= constant_current; 
        integral_term=constant_current-Kp_voltage*error-derivative_term;
    elseif strcmp(mode,'CV')
        control_signal_current = Kp_voltage * error + integral_term + derivative_term;
        I(k+1:k+2) = control_signal_current;  
    elseif strcmp(mode, 'Rest')
        I(k+1:k+2) = 0;
         if strcmp(modename, 'RestOCV')
            OCV_Volt=[OCV_Volt,Volt(k)]; %special case for OCV testing. extra Volt vector to allow backcheck (Volt(k-100)). could be cleared otherwise (every 5000 steps)
            k_OCV_count=k_OCV_count+1;
         end
    elseif strcmp(mode, 'Test')
        I(k+1:k+2) = I_Test(k_Test);
        k_Test=k-k_Rest_Test+2;
    elseif strcmp(mode, 'SOC1') %discharge at some C-rate
        I(k+1:k+2) = 2*p.OneC;
        if SOC(k+1) < 0.15
            I(k+1:k+2)=1/10*I(k+1:k+2);
        end
        %if SOC(k+1)<0.05 %|| SOC(k+1)>0.92 %ATT
        %    I(k+1:k+2)=1/10*I(k+1:k+2);
        %end
    elseif strcmp(mode, 'DDP')
        I(k+1:k+2) = DDP(k_DDP);
        
    elseif strcmp(mode, 'UDDS')
        I(k+1:k+2) = UDDS_i(k_UDDS);
        if v_min_firstoccurence
            I(k+1:k+2)=1/10*I(k+1:k+2);
        end
        if v_min && v_min_firstoccurence == false %SOC(k+1) < 0.05
            I(k+1:k+2)=1/10*I(k+1:k+2);
            v_min_firstoccurence=true; % only reduce current once
            k_v_min=k_full+10; %10 second grace period for v to relax above v_min
            v_min=false;
        end
        k_UDDS=k_UDDS+1;
        if v_min && v_min_firstoccurence && k_full>k_v_min
            k_end=k_full;
        end
        if k_UDDS==length(UDDS_i)
            k_UDDS=1; %start another UDDS loop
        end
    elseif strcmp(mode, 'HPPC')
        if SOC(k+1)>0.968 || SOC(k+1)<0.038
            I(k+1:k+2) = HPPC_current_extreme{k_HPPC_ex}(k_HPPC);
        else
            I(k+1:k+2) = HPPC_current(k_HPPC);
        end
        %if SOC(k+1)>0.915 || SOC(k+1)<0.075
        %    I(k+1:k+2)=I(k+1:k+2)*1/10;
        %end

    elseif strcmp(mode, 'OCV')
        if k_OCV==1 %first cycle
            I(k+1:k+2) = 0.5*p.OneC;
        elseif OCV_levels(k_OCV)<OCV_levels(k_OCV-1)
            I(k+1:k+2) = (0.5)*p.OneC;         
        else
            I(k+1:k+2) = (-0.5)*p.OneC;
        end
        if SOC(k+1)<0.1 %|| SOC(k+1)>0.92 %ATT
            I(k+1:k+2)=1/10*I(k+1:k+2);
        end
        if v_min 
            I(k+1:k+2) = 1/2*I(k+1:k+2);%/voltage_exceed_count;
        end



        
    elseif strcmp(mode,'custom')
        if k_full+2==length(i_custom)
            k_end=k_full;
            CTRL=false; %stop controller
            t=0:p.delta_t:k_full;
        end
        I(k+1:k+2) = i_custom(k);
    end  

    if k_full == k_end
        CTRL=false; %stop controller
        t=0:p.delta_t:k_full-1; %-1 not sure why actually.
    end

    n_Li_s(k+1) = y(end-1);
    n_Li_e(k+1) = y(end);
    
    eta_s_n = phi_s_n - phi_e(1:Nn,:);
    eta_s_p = phi_s_p - phi_e(end-Np+1:end, :);
  
%   Commented by SHP.  
   fprintf(1,'Time : %3.2f sec | Mode : %s |C-rate : %2.2f | Temp : %2.1fdegC | SOC : %1.3f | Voltage : %2.3fV | Newton Iters : %2.0f\n',...
       k_full,modename,I(k+1)/p.OneC,T(k+1)-273.15,SOC(k+1),Volt(k+1),stats.iters);
    
   if v_min
       voltage_exceed_count=voltage_exceed_count+1;
   else
       voltage_exceed_count=2;
   end
   v_min=false;
    if(Volt(k+1) < p.volt_min)
        fprintf(1,'Min Voltage of %1.1fV exceeded\n',p.volt_min);
        v_min=true;
        beep;
        %break;
    elseif(Volt(k+1) > p.volt_max)
        fprintf(1,'Max Voltage of %1.1fV exceeded\n',p.volt_max);
        beep;
        %break;
    elseif(any(c_ex(:,k) < 1))
        fprintf(1,'c_e depleted below 1 mol/m^3\n');
        v_min=true;
        beep;
        %break;
    end
    

    if mod(k,tmax) == 0 || k_full == k_end
        filename = sprintf('results_%08d_%08d.mat', last_saved_batch+1, k_full);
        T_save = T(2:end);
        SOC_save = SOC(2:end);
        I_save = I(2:end-1);
        Volt_save = Volt(2:end);
        c_ss_n_save = [c_ss_n(1,2:end); c_ss_n(end,2:end)];
        c_ss_p_save = [c_ss_p(1,2:end); c_ss_p(end,2:end)];
        c_e_0n_save = c_e_0n(2:end);
        c_e_0p_save = c_e_0p(2:end);
        save(filename, 'T_save', 'SOC_save', 'I_save', 'Volt_save', 'c_ss_n_save', 'c_ss_p_save', 'c_e_0n_save', 'c_e_0p_save');
        clear T_save SOC_save I_save Volt_save c_ss_n_save c_ss_p_save c_e_0n_save c_e_0p_save; 
        
        T = T(end);
        SOC = SOC(end);
        I = I(end-1:end);
        Volt = Volt(end);
        c_e_0n = c_e_0n(end);
        c_e_0p = c_e_0p(end);
        eta_s_Ln = eta_s_Ln(end);
        n_Li_e = n_Li_e(end);
        n_Li_s = n_Li_s(end);

        x = x(:,end);
        z = z(:,end);
        c_ss_n = c_ss_n(:,end);
        c_avg_n = c_avg_n(:,end);
        c_avg_p = c_avg_p(:,end);
        c_e = c_e(:,end);
        c_ex = c_ex(:,end);
        c_s_n = c_s_n(:,end);
        c_s_p = c_s_p(:,end);
        c_ss_p = c_ss_p(:,end);
        eta_n = eta_n(:,end);
        eta_p = eta_p(:,end);
        eta_s_n = eta_s_n(:,end);
        eta_s_p = eta_s_p(:,end);
        i_en = i_en(:,end);
        i_ep = i_ep(:,end);
        jp = jp(:,end);
        jn = jn(:,end);
        phi_e = phi_e(:,end);
        phi_s_n = phi_s_n(:,end);
        phi_s_p = phi_s_p(:,end);

        k=0; %reset the within-batch k, (set to k=1 in the next lines)
        last_saved_batch = k_full;
    end
    
    k=k+1;
    k_full = k_full+1;
    

end
I(end)=[];


%% Outputs
disp('Simulating Output Vars...');
simTime = toc;
fprintf(1,'Simulation Time : %3.2f min\n',simTime/60);
disp('To plots results, run...');
disp(' plot_dfn')
disp(' animate_dfn')

HPPC_startingOCV=nonzeros(HPPC_startingOCV); %sometimes has a zero at the end, not sure why
HPPC_voltage_zeroindex = sum(abs(HPPC_voltage),2)>0;
HPPC_voltage=HPPC_voltage(HPPC_voltage_zeroindex,:);
save('HPPC_data.mat', 'HPPC_voltage', 'HPPC_startingOCV','HPPC_startingSOC','HPPC_current');

%% Save Output Data for Plotting
% out.date=date;
% out.time=t;
% out.cur=I;
% out.volt=Volt;
% out.soc=SOC;
% out.c_ss_n=c_ss_n;
% out.c_ss_p=c_ss_p;
% out.temp = T;
% % out.eta_s_Ln=eta_s_Ln;
% out.ce0n=c_e_0n;
% out.ce0p=c_e_0p;
% out.simtime=simTime;

% save('data/Int_Obs/dfn_UDDS_NCM20Q.mat', '-struct', 'out');
% save('data/new/dfn_ce_new.mat', '-struct', 'out'); %10C Discharge LiCoO2

%% Save Output Data for Sensitivity Analysis
out.date=date;
out.p = p;
out.time=t;
out.cur=I;
out.volt=Volt;
out.soc=SOC;
out.temp=T;
out.x = x;
out.z = z;
out.simtime=simTime;

% save('data/sensitivity/0C_dfn.mat', 'out');



