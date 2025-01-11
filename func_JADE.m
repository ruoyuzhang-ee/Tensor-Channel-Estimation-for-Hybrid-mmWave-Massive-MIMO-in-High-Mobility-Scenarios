function [sin_AoD_est, Doppler_est] = func_JADE( BBhat, F, N_BS, optionsAltOpt )

[~,L] = size(BBhat);

T = optionsAltOpt.T;

sin_AoD_est = zeros(1,L);
Doppler_est = zeros(1,L);

search_dg = optionsAltOpt.search_dg;

Num_Iter = optionsAltOpt.Num_Iter;
w_max = optionsAltOpt.w_max;
for ll = 1:L
    bb_hat = BBhat(:,ll);
    
    w_init = zeros(1,T); 
    w_prev = w_init;
    w_l1_hat_list = zeros(Num_Iter,1);
    pos_max_list = zeros(Num_Iter,1);
    for ii = 1:Num_Iter
        % Estimate of c
        [ sin_AoD ] = func_1D_searchAngle( diag(exp(-1i*w_prev)) * bb_hat, F, N_BS, search_dg);
        A_Dict = 1/sqrt(N_BS)*exp(1j*pi*(0:N_BS-1).'*sin_AoD);
        Phi = F' * A_Dict;
        h_hat = pinv(Phi) * diag(exp(-1i*w_prev)) * bb_hat;
        pos_max_list(ii) = sin_AoD;
        
        w_range = [-w_max:0.0001:w_max];
        fval = zeros(length(w_range),1);
        for w_ii = 1:length(w_range)
            ww = w_range(w_ii);
            temp_fval = real(bb_hat' * diag(exp(1i*ww*(1:1:T))) * Phi * h_hat );
            fval(w_ii) = temp_fval;
        end
        [~, pos] = max(fval);
        w_l1_hat = w_range(pos);
        
        % Update of w
        w_new = w_l1_hat * (1:1:T);
        w_prev = w_new;
        w_l1_hat_list(ii) = w_l1_hat;
    end
    sin_AoD_est(ll) = pos_max_list(end);
    Doppler_est(ll) = w_l1_hat_list(end);
end
end
