diary('roca_output1.txt');  % Start recording output
diary on;

B ='Binary data.csv';      % CSV file name
row_start = 2;               % Starting row number
row_end   = 3001;               % Ending row number
col_start = 2;                % Starting column number
col_end   = 2001;               % Ending column number

% === READER ===
fid = fopen(B, 'r');
if fid == -1
    error('Cannot open file: %s', B);
end

% Skip header if it exists (optional)
header_line = fgetl(fid);  % comment out if no header

% Skip rows before row_start
for i = 2:row_start-1
    fgetl(fid);
end

% Initialize data matrix
num_rows = row_end - row_start + 1;
num_cols = col_end - col_start + 1;
B = zeros(num_rows, num_cols);

% Read desired rows and extract required columns
for i = 1:num_rows
    line = fgetl(fid);
    if ischar(line)
        values = str2double(strsplit(line, ','));
        B(i, :) = values(col_start:col_end);
    else
        break;
    end
end

fclose(fid);

[N, M] = size(B);  % Number of genomes x Pfam domains
disp(['Input matrix B size: ', num2str(N), ' x ', num2str(M)]);

[Bcap,lambda] = remove_phylogeny(B,1);
disp(['Bcap size: ', num2str(size(Bcap,1)), ' x ', num2str(size(Bcap,2))]);
disp(['Lambda size: ', num2str(length(lambda))]);

%Threshold to only capture the large Eigen values and remove noise
N_shuffles = 50; %Number of shuffles
thresh = 1; %percentile to define threshold for significant positive and negative correlations
rand("seed", 42); % Octave
[lambda_max_rnd,pos_thresh,neg_thresh] = ...
    computing_lambda_max_rnd(Bcap,N_shuffles,thresh);
disp(['lambda_max_rnd: ', num2str(lambda_max_rnd)]);
disp(['pos_thresh: ', num2str(pos_thresh)]);
disp(['neg_thresh: ', num2str(neg_thresh)]);

alpha = sum(lambda>lambda_max_rnd);
disp(['Alpha (number of significant eigenvalues): ', num2str(alpha)]);

gamma_k = computing_gamma_k(alpha,lambda,M,N);
disp(['Gamma_k size: ', num2str(length(gamma_k))]);

PC_roca = CorrITSPCA(Bcap, gamma_k, alpha);
disp(['PC_roca size: ', num2str(size(PC_roca,1)), ' x ', num2str(size(PC_roca,2))]);

% Load binary matrix from CSV
raw_data = dlmread('Actual data.csv', ',');  % Make sure this is the correct filename

% Subset the matrix
B1 = raw_data(row_start:row_end, col_start:col_end);

% Identify columns where all values are 0 (fully absent domains)
site_freq_no_domains = find(all(B1==0, 1));
disp(['Number of fully absent domains (all zeros columns): ', num2str(length(site_freq_no_domains))]);

ls = size(B,2);       % Get number of columns in table T
disp(['Total number of columns in B: ', num2str(ls)]);

true_indices = setdiff([1:ls],site_freq_no_domains);
disp(['Number of columns after removing absent domains: ', num2str(length(true_indices))]);

[sec_eig_roca, sec_eig_roca_true, sec_eig_roca_incl_cs, length_sec_roca] = ...
    form_sectors_roca(PC_roca,alpha,site_freq_no_domains,true_indices,ls);
disp(['Number of sectors formed by RoCA: ', num2str(length(sec_eig_roca))]);
disp(sec_eig_roca);

n_secs_roca = length(sec_eig_roca);

C_hat = compute_clean_C(Bcap,alpha,lambda); %cleaned standardized correlation matrix
disp(['Size of cleaned correlation matrix C_hat: ', num2str(size(C_hat,1)), ' x ', num2str(size(C_hat,2))]);

% B1 is your binary matrix (or any numeric matrix)
freq_bin = mean(B1, 1);  % This returns a 1 x num_columns vector of column means
disp(['Size of frequency vector freq_bin: ', num2str(length(freq_bin))]);

[mc,mean_abs_corr,per_neg_corr,per_pos_corr] = ...
    stats_sectors(C_hat,sec_eig_roca,n_secs_roca,freq_bin,pos_thresh,neg_thresh);
disp('Mean conservation (mc) per sector:');
disp(mc);

disp('Mean absolute correlation per sector:');
disp(mean_abs_corr);

disp('Percentage negative correlations per sector:');
disp(per_neg_corr);

disp('Percentage positive correlations per sector:');
disp(per_pos_corr);
diary off;  % Stop recording output
