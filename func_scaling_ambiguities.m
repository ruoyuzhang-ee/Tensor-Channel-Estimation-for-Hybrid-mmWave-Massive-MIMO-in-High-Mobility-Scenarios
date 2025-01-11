function [ Singular_Diag_M1, norm_err ] = func_scaling_ambiguities( AAest, AAhat )
% Estimate the ambiguity

[~,L] = size(AAhat);

Singular_Diag_M1 = zeros(L,L);
for ll = 1:L
    Singular_Diag_M1(ll,ll) = pinv(AAest(:,ll)) * AAhat(:,ll);
end

errAA = AAhat - AAest * Singular_Diag_M1;
norm_err = norm(errAA);

end

