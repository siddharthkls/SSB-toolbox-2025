function [gamma_k,stdNonSecAsym] = computing_gamma_k(num_eig,lambda,M,N)

% Computing the sparsity parameters \gamma_k for CorrITRoCA algorithm
% 
% Last updated: 2017-03-25
% 
% INPUT
%     num_eig       number of eigenvectors to be computed
%     lambda        the sample eigenvalues of size 1-by-(num_eig)
%     N             the number of samples or sequences
%     M             the number of variables or sites
%     
% OUTPUT
%     gamma_k       1-by-(num_eig) vector with sparsity parameters
%     stdNonSecAsym standard deviation of the asymptotic distribution of
%                   nonsec coordinates
c = M / N;
gamma_k = zeros(1, num_eig);

lambda_pop = ((lambda(1:num_eig) + 1 - c) + sqrt((lambda(1:num_eig) + 1 - c).^2 - 4 * lambda(1:num_eig))) / 2;

% Asymptotic distribution of arbitrary non-sector coordinates
c_v_sample = sqrt((1 - c ./ (lambda_pop - 1).^2) ./ (1 + c ./ (lambda_pop - 1)));
stdNonSecAsym = sqrt((1 - c_v_sample.^2) ./ (M - num_eig));

% Use a cell array for Octave compatibility
progress_msgs = cell(1, num_eig);

fprintf('Starting computation of gamma_k thresholds for %d PCs...\n', num_eig);

parfor kk = 1:num_eig
    x = 0:1e-8:6 * stdNonSecAsym(kk);

    % CDF
    F = round2dp((erf(x / (stdNonSecAsym(kk) * sqrt(2)))).^M, 5);

    [F2, mask] = unique(F); % removing non-unique elements from CDF
    x2 = x(mask);

    a = cumsum(F2);
    a2 = a / max(a);

    b = find(a2 <= 0.95);
    gamma_k(kk) = x2(b(end) + 1);

    % Save progress messages as character vectors
    progress_msgs{kk} = sprintf('  Finished gamma_k computation for PC %d', kk);
end

% Display progress after the parallel loop
for kk = 1:num_eig
    fprintf('%s\n', progress_msgs{kk});
end

fprintf('All gamma_k thresholds computed.\n');

