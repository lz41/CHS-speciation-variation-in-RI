clear; % clear variables in the workspace

seed_rng0 = 3; 
fs1 = [0.05,0.2,0.4,0.6,0.8];
sigma_fs1 = [0.05,0.15,0.25,0.35,0.45];
fs_d = [0.1,0.2,0.5];
N2_immig = 10; 
num_gen = 100;
rep_num = 10;

%% static
fluctuation_regime = 'static';
folder_name = sprintf('/uufs/chpc.utah.edu/common/home/u6033116/matlab/output_fs_d/',fluctuation_regime); % folder name where the output files are

%seed_rng0 = 88; % if you used different ones for different scenarios etc. 
sigma_fs1_static = 0.2; % does not matter for static
kappa = 0.6; % does not matter for static

immig_allele_static = zeros(length(sigma_fs1_static),length(fs1));
allele0 = zeros(num_gen, rep_num);
num_indiv = zeros(num_gen, rep_num);

for i = 1:length(fs1)
    for q = 1:length(sigma_fs1_static) 
        for rep = 1:rep_num
            infilename = sprintf([folder_name,'rep%i_seed%i_%s_fs1%i_sigma_fs1%i_fs_d%i'], ...
                rep,seed_rng0,fluctuation_regime,round(fs1(i)*100),round(sigma_fs1*100),round(fs_d*100));

            load(infilename);

            allele0(:,rep) = accumarray(OP_data(:,1),OP_data(:,4))+accumarray(OP_data(:,1),OP_data(:,5))...
                +2*accumarray(OP_data(:,1),OP_data(:,6)); % adding the number of allele0
            num_indiv(:,rep) = accumarray(OP_data(:,1),sum(OP_data(:,3:6),2));

        end
        immig_allele_static(q,i) = mean(mean(allele0(51:num_gen,:),1)); % calculating the mean for the last 50 generations
    end
end

% to export allele0 and num_indiv
writematrix(allele0,[folder_name,'allele0_static.csv'])
writematrix(num_indiv,[folder_name,'num_indiv_static.csv'])

%% temporal kappa = 0 white noise
fluctuation_regime = 'temporal';
folder_name = sprintf('/uufs/chpc.utah.edu/common/home/u6033116/matlab/output_fs_d/',fluctuation_regime); % folder name where the output files are

seed_rng0 = 4; % if you used a different seed for these runs
kappa = 0;

immig_allele_temporal_0 = zeros(length(sigma_fs1),length(fs1));
allele0 = zeros(num_gen, rep_num);
num_indiv = zeros(num_gen, rep_num);

for i = 1:length(fs1) 
    for q = 1:length(sigma_fs1) 
        for rep = 1:rep_num

            infilename = sprintf([folder_name,'rep%i_seed%i_%s_fs1%i_sigma_fs1%i_fs_d%i'],...
                rep,seed_rng0,fluctuation_regime,round(kappa*10),round(fs1(i)*100),round(sigma_fs1(q)*100),N2_immig);

            load(infilename);

            allele0(:,rep) = accumarray(OP_data(:,1),OP_data(:,4))+accumarray(OP_data(:,1),OP_data(:,5))...
                +2*accumarray(OP_data(:,1),OP_data(:,6)); % adding the number of allele0
            num_indiv(:,rep) = accumarray(OP_data(:,1),sum(OP_data(:,3:6),2));

        end
        immig_allele_temporal_0(q,i) = mean(mean(allele0(51:num_gen,:),1)); % calculating the mean for the last 50 generations
    end
end

% to export allele0 and num_indiv
writematrix(allele0,[folder_name,'allele0_temporal0.csv'])
writematrix(num_indiv,[folder_name,'num_indiv_temporal0.csv'])

%% temporal kappa = 0.6 red noise
seed_rng0 = 5;
kappa = 0.6;

immig_allele_temporal_6 = zeros(length(sigma_fs1),length(fs1));
allele0 = zeros(num_gen, rep_num);
num_indiv = zeros(num_gen, rep_num);

for i = 1:length(fs1) 
    for q = 1:length(sigma_fs1) 
        for rep = 1:rep_num

            infilename = sprintf([folder_name,'rep%i_seed%i_%s_fs1%i_sigma_fs1%i_fs_d%i'],...
                rep,seed_rng0,fluctuation_regime,round(kappa*10),round(fs1(i)*100),round(sigma_fs1(q)*100),N2_immig);
            load(infilename);

            allele0(:,rep) = accumarray(OP_data(:,1),OP_data(:,4))+accumarray(OP_data(:,1),OP_data(:,5))...
                +2*accumarray(OP_data(:,1),OP_data(:,6)); % adding the number of allele0
            num_indiv(:,rep) = accumarray(OP_data(:,1),sum(OP_data(:,3:6),2));

        end
        immig_allele_temporal_6(q,i) = mean(mean(allele0(51:num_gen,:),1)); % calculating the mean for the last 50 generations
    end
end

% to export allele0 and num_indiv
writematrix(allele0,[folder_name,'allele0_temporal6.csv'])
writematrix(num_indiv,[folder_name,'num_indiv_temporal6.csv'])

%% spatial_random
fluctuation_regime = 'spatial_random';
folder_name = sprintf('/uufs/chpc.utah.edu/common/home/u6033116/matlab/output_fs_d/',fluctuation_regime); % folder name where the output files are

seed_rng0 = 1; % if you used a different seed for these runs
kappa = 0; % anything you used (kappa is not used in the simulations)

immig_allele_spatial_random = zeros(length(sigma_fs1),length(fs1));
allele0 = zeros(num_gen, rep_num);
num_indiv = zeros(num_gen, rep_num);

for i = 1:length(fs1) 
    for q = 1:length(sigma_fs1) 
        for rep = 1:rep_num

            infilename = sprintf([folder_name,'rep%i_seed%i_%s_fs1%i_sigma_fs1%i_fs_d%i'],...
                rep,seed_rng0,fluctuation_regime,round(kappa*10),round(fs1(i)*100),round(sigma_fs1(q)*100),N2_immig);
            load(infilename);

            allele0(:,rep) = accumarray(OP_data(:,1),OP_data(:,4))+accumarray(OP_data(:,1),OP_data(:,5))...
                +2*accumarray(OP_data(:,1),OP_data(:,6)); % adding the number of allele0
            num_indiv(:,rep) = accumarray(OP_data(:,1),sum(OP_data(:,3:6),2));

        end
        immig_allele_spatial_random(q,i) = mean(mean(allele0(51:num_gen,:),1)); % calculating the mean for the last 50 generations
    end
end

% to export allele0 and num_indiv
writematrix(allele0,[folder_name,'allele0_spatial_random.csv'])
writematrix(num_indiv,[folder_name,'num_indiv_spatial_random.csv'])

%% spatial linear gradient
fluctuation_regime = 'spatial_lineargrad'; 
folder_name = sprintf('/uufs/chpc.utah.edu/common/home/u6033116/matlab/output_fs_d/',fluctuation_regime); % folder name where the output files are

seed_rng0 = 2; % if you used a different seed for these runs
kappa = 0; % anything you used (kappa is not used in the simulations)

immig_allele_spatial_lineargrad = zeros(length(sigma_fs1),length(fs1));
allele0 = zeros(num_gen, rep_num);
num_indiv = zeros(num_gen, rep_num);

for i = 1:length(fs1) 
    for q = 1:length(sigma_fs1) 
        for rep = 1:rep_num

            infilename = sprintf([folder_name,'rep%i_seed%i_%s_fs1%i_sigma_fs1%i_fs_d%i'],...
                rep,seed_rng0,fluctuation_regime,round(kappa*10),round(fs1(i)*100),round(sigma_fs1(q)*100),N2_immig);
            load(infilename);

            allele0(:,rep) = accumarray(OP_data(:,1),OP_data(:,4))+accumarray(OP_data(:,1),OP_data(:,5))...
                +2*accumarray(OP_data(:,1),OP_data(:,6)); % adding the number of allele0
            num_indiv(:,rep) = accumarray(OP_data(:,1),sum(OP_data(:,3:6),2));

        end
        immig_allele_spatial_lineargrad(q,i) = mean(mean(allele0(51:num_gen,:),1)); % calculating the mean for the last 50 generations
    end
end

% to export allele0 and num_indiv
writematrix(allele0,[folder_name,'allele0_spatial_lineargrad.csv'])
writematrix(num_indiv,[folder_name,'num_indiv_spatial_lineargrad.csv'])

%% heat maps -- temporal variation
load('/Users/etsuko/Dropbox/Etsukos stuff/Linyi project/model/newmap3.mat') % to color the heat maps in the way we want

relative_temporal_0 = (immig_allele_temporal_0-repmat(immig_allele_static,11,1))./immig_allele_plot; %change in gene flow
figure(1); clf(1); subplot(1,2,1); 
imagesc((fliplr(flipud(relative_temporal_0)))); colormap(newmap4); colorbar; axis square
xlabel('fecundity selection coefficient','FontSize',16); ylabel('variance of fecundity selection coefficient','FontSize',16); 
xticks([2,4,6,8,10]); yticks([2,4,6,8,10])
xticklabels({'0.3','0.5','0.7','0.9','0.99'}) % in terms of fs coeff (in a normal way) %fs1 = [0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8];
yticklabels({'0.45','0.35','0.25','0.15','0.05'})

relative_temporal_6 = (immig_allele_temporal_6-repmat(immig_allele_static,11,1))./immig_allele_plot; 
figure(1); subplot(1,2,2);  
imagesc((fliplr(flipud(relative_temporal_6)))); colormap(newmap4); colorbar; axis square
xlabel('fecundity selection coefficient'); ylabel('variance of fecundity selection coefficient'); 
xticks([2,4,6,8,10]); yticks([2,4,6,8,10])
xticklabels({'0.3','0.5','0.7','0.9','0.99'}) % in terms of fs coeff (in a normal way) %fs1 = [0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8];
yticklabels({'0.45','0.35','0.25','0.15','0.05'})

subplot(1,2,1)
clim([min([min(relative_temporal_0(:)),min(relative_temporal_6(:))]),max([max(relative_temporal_0(:)),max(relative_temporal_6(:))])]);
subplot(1,2,2)
clim([min([min(relative_temporal_0(:)),min(relative_temporal_6(:))]),max([max(relative_temporal_0(:)),max(relative_temporal_6(:))])]);

%% heat maps -- spatial variation

relative_spatial_random = (immig_allele_spatial_random-repmat(immig_allele_static,11,1))./immig_allele_plot; %change in gene flow
figure(2); clf(2); subplot(1,2,1); 
imagesc((fliplr(flipud(relative_spatial_random)))); colormap(newmap4); colorbar; axis square
xlabel('fecundity selection coefficient','FontSize',16); ylabel('variance of fecundity selection coefficient','FontSize',16); 
xticks([2,4,6,8,10]); yticks([2,4,6,8,10])
xticklabels({'0.3','0.5','0.7','0.9','0.99'}) % in terms of fs coeff (in a normal way) %fs1 = [0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8];
yticklabels({'0.45','0.35','0.25','0.15','0.05'})

relative_spatial_lineargrad = (immig_allele_spatial_lineargrad-repmat(immig_allele_static,11,1))./immig_allele_plot; 
figure(2); subplot(1,2,2);  
imagesc((fliplr(flipud(relative_spatial_lineargrad)))); colormap(newmap4); colorbar; axis square
xlabel('fecundity selection coefficient'); ylabel('variance of fecundity selection coefficient'); 
xticks([2,4,6,8,10]); yticks([2,4,6,8,10])
xticklabels({'0.3','0.5','0.7','0.9','0.99'}) % in terms of fs coeff (in a normal way) %fs1 = [0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8];
yticklabels({'0.45','0.35','0.25','0.15','0.05'})

subplot(1,2,1)
clim([min([min(relative_temporal_0(:)),min(relative_temporal_6(:))]),max([max(relative_temporal_0(:)),max(relative_temporal_6(:))])]);
subplot(1,2,2)
clim([min([min(relative_temporal_0(:)),min(relative_temporal_6(:))]),max([max(relative_temporal_0(:)),max(relative_temporal_6(:))])]);