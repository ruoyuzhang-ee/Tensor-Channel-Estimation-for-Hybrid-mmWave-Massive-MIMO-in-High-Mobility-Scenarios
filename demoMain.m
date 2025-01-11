% The demo file for high-mobility channel estimation 

% Setting
c_light = 3*1e8;
fc = 60*1e9;
K_bar = 128; 
N_BS = 64;
N_MS = 32;
L = 3;
fs = 100 * 1e6; %Sampling rate
tau_max = 320 * 1e-9; 
N_CP = 32; %CP length
% Training parameters
N = 20;
N_RF = 6;
K = 10;

SNR_dB = 20;

% Doppler generate
T_sampling = 1/(fs); 
v_max = 120;
w_max = 2*pi * fc/c_light * (v_max)/3.6 * T_sampling * (N_CP+K_bar);  


%% 2. Generate channel H
sin_AoD_true = sin(random('unif',-pi/2,pi/2,1,L));
sin_AoA_true = sin(random('unif',-pi/2,pi/2,1,L));
tau_true = tau_max * random('unif',0,1,1,L);
alpha_true = 1/sqrt(2) * (randn(1,L)+ 1i*randn(1,L));
velocity_ll = random('unif',0,v_max,1,L);
Doppler_true = 2*pi * fc/c_light * velocity_ll/3.6 * T_sampling * (N_CP+K_bar);

a_MS_matrix = 1/sqrt(N_MS)*exp(1j*pi*(0:N_MS-1).'*sin_AoA_true);
a_BS_matrix = 1/sqrt(N_BS)*exp(1j*pi*(0:N_BS-1).'*sin_AoD_true);

% true channel
H = zeros(N_MS, N_BS, K, N);
for tt = 1:N
    for kk = 1:K
        H_kk = zeros(N_MS,N_BS);
        for ll = 1:L
            H_kk = H_kk + alpha_true(ll) * exp(1i*Doppler_true(ll)*tt) * exp(-1i*2*pi*tau_true(ll)*fs*kk/K_bar) * a_MS_matrix(:,ll) * a_BS_matrix(:,ll)';
        end
        H(:,:,kk,tt) = H_kk;
    end
end



% Training
phase = 2*pi * rand(N_BS, N);
F_temp = exp(1j*phase);
F = sqrt(N)*F_temp/norm(F_temp, 'fro'); %precoder
phase = 2*pi * rand(N_MS, N_RF);
W_temp = exp(1j*phase);
W = sqrt(N_RF)*W_temp/norm(W_temp, 'fro'); %combiner

% Factor matrix
AA = W' * a_MS_matrix;
BB = zeros(N, L);
for l = 1:L
    w_list_true = Doppler_true(l) * [1:1:N];
    BB(:,l) = diag(exp(1i*w_list_true)) * F' * a_BS_matrix(:,l);
end
CC = zeros(K, L);
g_tau_true = zeros(K,L);
for ll = 1:L
    for kk = 1:K
        g_tau_true(kk,ll) = exp(-1i*2*pi*tau_true(ll)*fs*kk/K_bar);
    end
end
for ll = 1:L
    CC(:,ll) = alpha_true(ll) * g_tau_true(:,ll);
end

% Generate the tensor
U = cell(1,3);
U{1,1} = AA;
U{1,2} = BB;
U{1,3} = CC;
Y_noiseless = cpdgen(U);
Y = noisy(Y_noiseless, SNR_dB);

% CPD
options.Display = false; % Show progress on the command line.
options.Initialization = @cpd_gevd; % Select pseudorandom initialization.
options.Algorithm = @cpd_als; % Select ALS as the main algorithm.
options.AlgorithmOptions.TolFun = 1e-12; % Set function tolerance stop criterion
options.AlgorithmOptions.TolX = 1e-12; % Set step size tolerance stop criterion
[Uhat,~] = cpd(Y,L,options);


%% Parameter estimate
search_dg = 1e-5;
CEOptions.F = F;
CEOptions.N_BS = N_BS;
CEOptions.W = W;
CEOptions.N_MS = N_MS;
CEOptions.search_dg = search_dg;
CEOptions.L = L;
CEOptions.fs = fs;
CEOptions.K = K;
CEOptions.K_bar = K_bar;
CEOptions.tau_max = tau_max;
CEOptions.T = N;
optionsAltOpt.T = N;
optionsAltOpt.w_max = w_max;
optionsAltOpt.search_dg = search_dg;
optionsAltOpt.Num_Iter = 30;
[H_est] = func_CP_TVFS( Uhat, CEOptions, optionsAltOpt );

H_est_N = H_est(:,:,:,N);
H_N = H(:,:,:,N);
mse_H = norm(H_est_N(:) - H_N(:))^2 / norm(H_N(:))^2

