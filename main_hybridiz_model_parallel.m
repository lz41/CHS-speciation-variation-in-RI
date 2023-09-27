% Linyi project
function main_hybridiz_model_parallel(seed_rng0,kappa) % input parameters: seed_rng0 = a random seed, kappa = autocorrelation parameter

%% the model

%% input parameters
[params_in,fluctuation_regime,selection,output_foldername] = input_params(kappa);

%% set up
p = parpool(23); % number of cores = number of available cores -1 
parfor runs = 1:size(params_in,2) % to use the parallel computation toolbox (PCT)
%for runs = 1:size(params_in,2) % if no PCT

    

    params = params_in(:,runs);


rng(seed_rng0+params(20));

      
    num_gen = params(1);% number of generations to simulate 

    L = params(2); % size of the arena LxL
    K1 = params(3); % K1 carrying cap (# indivs) for resident; constant or variable
    q1 = params(4);  % quality of foreign habitat = Qual_resid*frac_qual

    N1_0 = params(5); % initial pop size of species 1
    N1_immig = params(6);
    N2_immig = params(7); % number of immigrants every time unit
    freq_joining = params(8); % every X time steps N2_0 number of 00's will join the arena

    % dispersal
    alpha = params(9); % inverse of mean dispersal distance
    Pr_disp = params(10); % probability of dispersal for each individual
    max_breed_disp = params(11); % = 1 if no redispersal

    % demographic parameters
    mean_fecund = params(12); % mean fecundity

    % viability

    % fecundity selection fs
    fs_w1 = params(13);
    fs1 = params(15); % fecundity coefficient for hybrids (mean)
    sigma_fs1 = params(16); % std of fecundity coefficient for hybrids
    fs_d = params(14); % fecundity coefficient for immigrants
    tol_fs1 = params(18); % tolerance for fs1
    tol_sigma_fs1 = params(19); % tolerance for sigma_fs1

    % temporal autocorrelation
    kappa = params(17); % temporal autocorrelation parameter positive = red noise, 0 = white noise, negative = blue noise

    X = ones(L,L); % all habitat 1

    % assigning carrying capacity
    X1 = (X==1)*K1; % home habitat
    X1 = X1 + (X==0)*K1*q1; % foreign habitat

    % landscape of selectivity (fs1 or s1)
    s_fluct = correlated_noise2(kappa, num_gen, fluctuation_regime, selection, X, fs1, sigma_fs1, tol_fs1, tol_sigma_fs1);

    % individuals
    Wasps = zeros(N1_0,6);
    Wasps(:,1:2) = randi(L,N1_0,2); % location (tree, patch) x,y
    Wasps(:,4:5) = 1; % 1 = resident alleles, 0 = foreign alleles
    Wasps(:,6) = rand(size(Wasps,1),1)<0.5; %1 = females

    % distance matrix
    [r,c]=find(X1);
    dist_matrix = squareform(pdist([c,r]));
    [~, Locb] = ismember(Wasps(:,1:2),[c,r],'rows'); % which tree?
    Wasps(:,3) = Locb; % tree (patch) id

    Pr_d = exp(-alpha*(dist_matrix)); % dispersal kernel
    Pr_d = Pr_d.*(1-eye(size(Pr_d)));
    Pr_d = Pr_d./sum(Pr_d,2); % normalize = dispersal within the arena (~reflecting boundaries)

    % data output table
    OP_data = zeros(num_gen*L*L,6);
    k = 1;

    %% iteration
    t0=tic;
    for t = 1:num_gen

        % immigrant inviability (before mating)

        % natal dispersal
        [Wasps] = wasp_dispersal(Wasps, c, r, L, Pr_disp, Pr_d);

        % immigrants and more residents from source come and join -- immigrants are
        % adults, so they can survive to reproduce and mostly die (had passed critical
        % developmental period)
        [Wasps] = immigr_join(Wasps, L, N2_immig, N1_immig, c, r, t, freq_joining);

        % mating -- mating and look for mates (breeding dispersal)
        [Wasps,female_mated,num_breed_disp] = wasp_mating(Wasps, c, r, L, Pr_disp, Pr_d, max_breed_disp);

        switch fluctuation_regime
            case {'both_random','both_lineargrad'} % not implemented in this version
                s_fluct_new = correlated_noise_both(kappa, num_gen, fs1, sigma_fs1, tol_fs1, s_fluct);
                % reproduction + fecundity selection
                [Wasps] = wasp_sexualrepro(Wasps, female_mated, mean_fecund, X1, s_fluct, fs_w1, fs_d, t, fluctuation_regime, selection, s_fluct_new);
            otherwise
                [Wasps] = wasp_sexualrepro(Wasps, female_mated, mean_fecund, X1, s_fluct, fs_w1, fs_d, t, fluctuation_regime, selection);
        end
       
        % viability selection -- not implemented yet

        % output data

        trees = unique(Wasps(:,3));

        for i = trees'
            OP_data(k,1) = t; % time
            OP_data(k,2) = i; % tree
            OP_data(k,3) = sum(ismember(Wasps(Wasps(:,3)==i,3:5),[i,1,1],'rows')); % genotype 11
            OP_data(k,4) = sum(ismember(Wasps(Wasps(:,3)==i,3:5),[i,1,0],'rows')); % genotype 10
            OP_data(k,5) = sum(ismember(Wasps(Wasps(:,3)==i,3:5),[i,0,1],'rows')); % genotype 01
            OP_data(k,6) = sum(ismember(Wasps(Wasps(:,3)==i,3:5),[i,0,0],'rows')); % genotype 00

            k = k+1;
        end

        % to print out during test run
        % fprintf('time step = %1i out of %1i\n', t, num_gen)
        % fprintf('number of breeding dispersal %1i\n', num_breed_disp)
        % [sum(OP_data(OP_data(:,1)==t,3),1),sum(sum(OP_data(OP_data(:,1)==t,4:5),1)),sum(OP_data(OP_data(:,1)==t,6),1)]

    end

    OP_data(find(OP_data(:,1)==0,1):size(OP_data,1),:) = []; % remove the empty rows
    Ttoc = toc(t0)

    %% data visualization

    % pop_size = zeros(num_gen,5);
    % for tt = 1:max(OP_data(:,1))
    %     pop_size(tt,1) = tt;
    %     pop_size(tt,2:5) = sum(OP_data(OP_data(:,1)==tt,3:6),1);
    % end

    %outfilename = sprintf('/Users/etsuko/Dropbox/Etsukos stuff/Linyi project/model_outputs/rep%i_seed%i_%s_kappa%i_fs1%i_sigma_fs1%i_N2immig%i',params(20),seed_rng0,fluctuation_regime,round(kappa*10),round(fs1*100),round(sigma_fs1*100),N2_immig);
    outfilename = sprintf([output_foldername,'rep%i_seed%i_kappa%i_%s_fs1%i_sigma_fs1%i_fs_d%i_alpha%i_Pr_d%i'], ...
        params(20),seed_rng0,round(kappa*10),fluctuation_regime,round(fs1*100),round(sigma_fs1*100),round(fs_d*100),round(alpha*10),round(Pr_disp*100));
    parsave(outfilename,OP_data,Ttoc,params);

end

delete(p) % shut down the parallel
