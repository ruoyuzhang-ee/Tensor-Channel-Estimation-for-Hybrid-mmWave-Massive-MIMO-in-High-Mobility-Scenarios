function [H_est] = func_CP_TVFS( Uhat, CEOptions, optionsAltOpt )
% Channel Estimation
F = CEOptions.F;
N_BS = CEOptions.N_BS;
W = CEOptions.W;
N_MS = CEOptions.N_MS;
search_dg = CEOptions.search_dg;
L = CEOptions.L;
fs = CEOptions.fs;
K = CEOptions.K;
K_bar = CEOptions.K_bar;
tau_max = CEOptions.tau_max;
T = CEOptions.T;

AAhat = Uhat{1,1};
BBhat = Uhat{1,2};
CChat = Uhat{1,3};


% est AoA
[ sin_AoA_est ] = func_1D_searchAngle( AAhat, W, N_MS, search_dg);

% est AoD and Doppler
[sin_AoD_est, Doppler_est] = func_JADE( BBhat, F, N_BS, optionsAltOpt );

% est tau
[ tau_est ] = func_1D_searchTau( CChat, tau_max, K, K_bar, fs, search_dg);

% compute the nonsingular diagonal matrices
a_MS_matrix_est = 1/sqrt(N_MS)*exp(1j*pi*(0:N_MS-1).'*sin_AoA_est);
a_BS_matrix_est = 1/sqrt(N_BS)*exp(1j*pi*(0:N_BS-1).'*sin_AoD_est);
AAest = W' * a_MS_matrix_est;
BBest = zeros(T, L);
for l = 1:L
    w_list_est = Doppler_est(l) * [1:1:T];
    BBest(:,l) = diag(exp(1i*w_list_est)) * F' * a_BS_matrix_est(:,l);
end
[ Singular_Diag_M1, ~ ] = func_scaling_ambiguities( AAest, AAhat );
[ Singular_Diag_M2, ~ ] = func_scaling_ambiguities( BBest, BBhat );
Singular_Diag_M3 = inv(Singular_Diag_M2) * inv(Singular_Diag_M1) * eye(L);

% est alpha
g_tau_est = exp(-1i*2*pi* [1:K]' * tau_est *fs/K_bar);
[ Singular_Diag_alpha, ~ ] = func_scaling_ambiguities( g_tau_est, CChat );
alpha_est = (diag(Singular_Diag_alpha) ./ diag(Singular_Diag_M3)).';

% est H
H_est = zeros(N_MS, N_BS, K, T);
for tt = 1:T
    for kk = 1:K
        H_kk_est = zeros(N_MS,N_BS);
        for ll = 1:L
            H_kk_est = H_kk_est + alpha_est(ll) * exp(1i*Doppler_est(ll)*tt) * exp(-1i*2*pi*tau_est(ll)*fs*kk/K_bar) * a_MS_matrix_est(:,ll) * a_BS_matrix_est(:,ll)';
        end
        H_est(:,:,kk,tt) = H_kk_est;
    end
end


end
