function [params,fluctuation_regime,selection,output_foldername] = input_params(kappa)

%% input parameters

% temporal autocorrelation
% kappa = [0.6]; % IN FUNCTION ARGUMENT temporal autocorrelation parameter positive = red noise, 0 = white noise, negative = blue noise

fs_w1 = 1; % fecundity selection for residents
shyb = [0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99]; % fecundity selection for hybrids
fs1 = 1-shyb; % 1-fecundity selection coefficient for hybrids (mean)
sigma_fs1 = [0.05,0.15,0.25,0.35,0.45]; %[0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5]; % standard deviation of fecundity coefficient for hybrids
s_immi = 0.95; % fecundity selection coefficient for immigrants
fs_d = 1-s_immi; % fs_d = 1-fecundity coefficient for immigrants
tol_fs1 = sigma_fs1/100;
tol_sigma_fs1 = sigma_fs1/100;

L = 50; % size of the arena LxL
K1 = 20; % K1 carrying cap (# indivs) for the resident
q1 = 1;  % quality of the habitat for immigrants = K1*frac_qual

fluctuation_regime = 'temporal' % 'spatial_random' or 'spatial_lineargrad', 'temporal' or 'both_random' 'both_lineargrad' 'static'
selection = 'fecundity'; % 'viability' or 'fecundity' or 'both' % currently 'fecundity' only

N1_0 = L*L*K1; % initial pop size of species 1
N1_immig = round(N1_0/100);
N2_immig = 10; %[2,5,round(N1_0/5000),round(N1_0/2000),round(N1_0/1000),round(N1_0/500)]; % number of immigrants every time unit
freq_joining = 1; % every X time steps N2_0 number of 00's will join the arena. 

% dispersal
alpha = [0.5,1,2]; % inverse of mean dispersal distance
Pr_disp = 0.2; % probability of dispersal for each individual
max_breed_disp = 50; % the max number of mate searching dispersal per individual; = 1 if no redispersal

% demographic parameters
mean_fecund = 20; % mean fecundity

% viability TBA

% simulation paramters
reps = 20; % number of replicates for each param combination
reps_start = 1; reps_end = reps_start+reps-1; % reps_start = which number you want to start labeling outputs with? 
param_combos = length(K1)*length(fs1)*length(sigma_fs1)*length(fs_d)*length(N2_immig)*length(alpha)*length(Pr_disp)*length(mean_fecund); % total number of combinations of parameter values
total_reps = reps*param_combos; % grand total of simulation runs
num_gen = 100; % number of generations to simulate
output_foldername = '/uufs/chpc.utah.edu/common/home/u6033116/matlab/output_alpha/'; % adjust to your folder setting

% construct the parameter vectors for parallel computing
params = zeros(20,total_reps);

m = 1;
for i = 1:length(fs1)
    for j = 1:length(sigma_fs1)
        for k = 1:length(kappa)
            for l = 1:length(N2_immig)
                for p = 1:length(alpha)
                    for q = 1:length(Pr_disp)
                        for s = 1:length(fs_d)
            params(15,m) = fs1(i);
            params(16,m) = sigma_fs1(j);
            params(17,m) = kappa(k);
            params(7,m) = N2_immig(l);
            params(9,m) = alpha(p);
            params(10,m) = Pr_disp(q);
            params(14,m) = fs_d(s);
            params(18,m) = tol_fs1(j);
            params(19,m) = tol_sigma_fs1(j);
            m = m + 1;
                        end
                    end
                end
            end
        end
    end
end
params(7,m:total_reps) = repmat(params(7,1:m-1),1,reps-1);
params(9,m:total_reps) = repmat(params(9,1:m-1),1,reps-1);
params(10,m:total_reps) = repmat(params(10,1:m-1),1,reps-1);
params(14,m:total_reps) = repmat(params(14,1:m-1),1,reps-1);
params(15,m:total_reps) = repmat(params(15,1:m-1),1,reps-1);
params(16,m:total_reps) = repmat(params(16,1:m-1),1,reps-1);
params(17,m:total_reps) = repmat(params(17,1:m-1),1,reps-1);
params(18,m:total_reps) = repmat(params(18,1:m-1),1,reps-1);
params(19,m:total_reps) = repmat(params(19,1:m-1),1,reps-1);

params(1,:) = repmat(num_gen,1,total_reps);
params(2,:) = repmat(L,1,total_reps);
params(3,:) = repmat(K1,1,total_reps);
params(4,:) = repmat(q1,1,total_reps);
params(5,:) = repmat(N1_0,1,total_reps);
params(6,:) = repmat(N1_immig,1,total_reps);
params(8,:) = repmat(freq_joining,1,total_reps);
params(11,:) = repmat(max_breed_disp,1,total_reps);
params(12,:) = repmat(mean_fecund,1,total_reps);
params(13,:) = repmat(fs_w1,1,total_reps);
params(20,:) = reshape(repmat(reps_start:reps_end,param_combos,1),total_reps,1);

end
