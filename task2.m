clear all; clc; close all;

% A line of sight free space optical communication system has a link
% with a range (L) of 1.2 km

L_in_km = 1.2; % in km
L_in_m = L_in_km * 1000; % in m

% The system is used in both clear weather and fog condition.
% In clear weather the attenuation coefficient of the atmospheric channel
% \gamma_{t}(\lambda) is 0.4 km^{-1}

gamma_t = 0.4; % in km^{-1}

% Figure 2 shows the international visibility specification under different fog condition
% The system includes a transmitter using a laser having a very
% narrow full-angle beam divergence \alpha of 0.2 degree

angle = 0.2; % in degree
angle_rad = deg2rad(angle); % in radian

% and emitted power P_{t} of 40 mW

P_t = 40; % in mW
P_t_dbm = 10 * log10(P_t); % in dBm

% At the receiver end a concentration lens with a radius r = 30 cm

r = 30; % in cm
r_in_meter = r / 100; % in meter

% is used to focus the light to the photodetector.
% The receiver has a sensitivity of -30 dBm

rec_sensitivity = -30; % in dBm

% We need to investigate the link gain and loss under both clear and foggy weather conditions.
% Please work out the following questions

% 1. In clear weather condition (no fog), determine the link budget (20%)

% Link Budget = Gain - Loss
% Link Budget = emitted power - recieved sensitivity - atmospheric loss - beam divergence loss

% STEP 1: Calculate the Atmospheric Loss
atm_loss = exp(-gamma_t * L_in_km); % in mW
atm_loss_db = 10 * log10(atm_loss); % in dB

% STEP 2: Calculate the Beam Divergence Loss / Geometry Loss

area_coverage = pi * (L_in_m * tan(angle_rad / 2))^2; % in m^2
area_collection = pi * (r_in_meter)^2;
geometry_path_loss = area_collection / area_coverage;
geometry_path_loss_dbm = 10 * log10(geometry_path_loss); % in dB

% STEP 3: Calculate the Link Budget
total_gain = P_t_dbm - rec_sensitivity;
total_loss = atm_loss_db + geometry_path_loss_dbm;

link_budget = total_gain + total_loss; % Total Loss is negative

% print the detail of every component of link budget
fprintf('Power Transmitted (dBm) = %.2f dBm\n', P_t_dbm);
fprintf('Recieved Sensitivity = %.2f dBm\n', rec_sensitivity);
fprintf('Atmospheric Loss = %.2f dB\n', atm_loss_db);
fprintf('Geometry Path Loss = %.2f dB\n', geometry_path_loss_dbm);
fprintf('Total Gain = %.2f dB\n', total_gain);
fprintf('Total Loss = %.2f dB\n', total_loss);
fprintf('Link Budget in Clear Weather Condition = %.2f dB\n', link_budget);

% NEXT QUESTION
% Using Figure 2, calculate the atmospheric power loss under moderate fog.
% Comment on whether the available link power margin is sufficient when we have moderate fog in the channel?

% Moderate Fog
% Visibility: 500 m
% Attenuation Coefficient: 28.9 dB/Km


att_coeff = 28.9; % in dB/km
fog_loss = att_coeff * L_in_km; % in dB

fprintf('\n\n');
fprintf('Atmospheric Loss in Moderate Fog = %.2f dB\n', fog_loss);
fprintf('Link Budget in Moderate Fog = %.2f dB\n', link_budget - fog_loss);

fprintf('Is it sufficient? ');
if (link_budget - fog_loss) > 0
    fprintf('Yes\n');
else
    fprintf('No\n');
end