function [ tau_est ] = func_1D_searchTau( CChat, tau_max, K, K_bar, fs, search_dg)
% one-dimensional search for delay

[~,L] = size(CChat);
tau_est = zeros(1,L);
resol = tau_max/K_bar/4;
candiGroupMax = 1;
% 1 Coarse grid search
tau_Dict = [0+resol: resol: tau_max-resol];
g_tau_coarse = exp(-1i*2*pi* [1:K]' * tau_Dict *fs/K_bar);

tau_coarse = zeros(candiGroupMax,L);
for l = 1:L
    f_corr1 = zeros(length(tau_Dict),1);
    cc_hat = CChat(:,l);
    for ii = 1:length(tau_Dict)
        f_corr1(ii) = abs(cc_hat' * g_tau_coarse(:,ii))^2/(norm(cc_hat)^2 * norm(g_tau_coarse(:,ii))^2);
    end
    [~, pos1] = sort(f_corr1, 'descend');
    tau_coarse(:,l) = tau_Dict(pos1(1:candiGroupMax));
end

% 2 Refining the search in the vicinity of possible grid points.
f_corr_candiGroup = zeros(candiGroupMax,L);
tauest_candiGroup = zeros(candiGroupMax,L);
for l = 1:L
    cc_hat = CChat(:,l);
    for jj = 1:candiGroupMax
        search_tau_temp = [tau_coarse(jj,l)-resol; tau_coarse(jj,l)+resol]; % The initial range
        tau_grid = [search_tau_temp(1,:) : search_dg : search_tau_temp(2,:)];
        Num_Taucandi = length(tau_grid);
        f_corr2 = zeros(Num_Taucandi,1);
        g_tau_refine = exp(-1i*2*pi* [1:K]' * tau_grid *fs/K_bar);
        for ii = 1:Num_Taucandi
            f_corr2(ii) = abs(cc_hat' * g_tau_refine(:,ii))^2/(norm(cc_hat)^2 * norm(g_tau_refine(:,ii))^2);
        end
        [~, pos] = sort(f_corr2,'descend'); % Find the maxmum one
        
        f_corr_candiGroup(jj,l) = f_corr2(pos(1));
        tauest_candiGroup(jj,l) = tau_grid(pos(1));
    end
end
for l = 1:L
    f_corrFinal = f_corr_candiGroup(:,l);
    [~, posFinal] = sort(f_corrFinal,'descend');
    posFinalmax = posFinal(1);
    tau_est(l) = tauest_candiGroup(posFinalmax,l);
end
    
end