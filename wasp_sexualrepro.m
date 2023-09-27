function [Wasps] = wasp_sexualrepro(Wasps, female_mated, mean_fecund, X1, s_fluct, fs_w1, fs_d, t, fluctuation_regime, selection, s_fluct_new)

X1 = X1(:);

num_offsp = poissrnd(mean_fecund,1,size(Wasps,1)); % draw a number of offspring for each femal from the poisson distn with mean = mean_fecund

num_offsp(sum(Wasps(:,4:5),2)==2) = round(num_offsp(sum(Wasps(:,4:5),2)==2)*fs_w1); % residents fecundity*coefficient of fecundity fitness = fecundity fitness

if strcmp(selection,'fecundity') % hybrid's fluctuating fecundity selection
    switch fluctuation_regime
        case {"temporal","static"}
            num_offsp(sum(Wasps(:,4:5),2)==1) = round(num_offsp(sum(Wasps(:,4:5),2)==1)*(s_fluct(t))); % reduce the number of offspring for hybrid
%         case {"both_random","both_lineargrad"} % currently not used
%             s_fluct = s_fluct + s_fluct_new(t); s_fluct(s_fluct<0)=0;
%             num_offsp(sum(Wasps(:,4:5),2)==1) = round(num_offsp(sum(Wasps(:,4:5),2)==1).*(s_fluct(Wasps(sum(Wasps(:,4:5),2)==1,3)))'); % hybrids (fluctuating fecundity selection)
        case {"spatial_random","spatial_lineargrad"}
            num_offsp(sum(Wasps(:,4:5),2)==1) = round(num_offsp(sum(Wasps(:,4:5),2)==1).*(s_fluct(Wasps(sum(Wasps(:,4:5),2)==1,3)))'); % selection coeff is location specific
    end
    num_offsp(sum(Wasps(:,4:5),2)==0) = round(num_offsp(sum(Wasps(:,4:5),2)==0)*fs_d); % immigrants (strong fecundity selection no matter what)
end

num_offsp(female_mated==0) = 0; % unmated females have 0 offspring
trees = unique(Wasps(Wasps(:,6)==1,3)); % patches females occupy
for i = trees' % for each female in patch i, choose a surviving number of offspring probabilistically with prob carry cap/total number of offspring born 
    num_offsp(Wasps(:,3)==i & Wasps(:,6)==1) = binornd(num_offsp(Wasps(:,3)==i & Wasps(:,6)==1),min(1,X1(i)/sum(num_offsp(Wasps(:,3)==i & Wasps(:,6)==1))));
end

newWasps = zeros(sum(num_offsp),6); % create a matrix with the number of rows = the total number of surviving offspring

which_females = find(female_mated>0); % find mothers

n_kids = 1;
for i = which_females'
    newWasps(n_kids:(n_kids+num_offsp(i)-1),1:3) = repmat(Wasps(i,1:3),num_offsp(i),1); % location is where the mother is
    inherited_from_mom = (rand(num_offsp(i),1)<0.5); % one allele will come from father but which one is being decided here
    newWasps(n_kids:(n_kids+num_offsp(i)-1),4:5) = repmat(Wasps(i,4:5),num_offsp(i),1).*[inherited_from_mom,(1-inherited_from_mom)]...
        + repmat(Wasps(female_mated(i),4:5),num_offsp(i),1).*[(1-inherited_from_mom),inherited_from_mom];% then others allele
    newWasps(n_kids:(n_kids+num_offsp(i)-1),6) = rand(num_offsp(i),1)<0.5; % 50:50 sex ratio of offspring
    
    n_kids = n_kids + num_offsp(i); % counting the total number of offspring processed
end

Wasps = newWasps; % overwrite the Wasp population with this new generation
end