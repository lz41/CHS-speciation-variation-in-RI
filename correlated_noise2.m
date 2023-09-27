function s_fluct = correlated_noise2(kappa, num_gen, fluctuation_regime, selection, X, fs1, sigma_fs1, tol_fs1, tol_sigma_fs1)

% fecundity selection
if strcmp(selection, 'fecundity')
switch fluctuation_regime
    case 'temporal'

            s_fluct = zeros(1,num_gen);
            s_fluct(1) = normrnd(0,sigma_fs1); % draw values from normal distribution with mean 0 and std sigma_fs1

            for i = 2:num_gen
                s_fluct(i) = kappa*s_fluct(i-1)+normrnd(0,sigma_fs1)*sqrt(1-kappa^2); % AR(1), Ruokolainen et al. 2009 TREE
            end

        s_fluct2 = s_fluct+(fs1-sum(s_fluct(s_fluct>=0))/length(s_fluct)); % adjust the mean to be = fs1
        
        kk = 0; % to iteratively attain all positive values and the desired mean at the same time
        while abs(sum(s_fluct2(s_fluct2>=0))/length(s_fluct2)-fs1)>tol_fs1 % is the mean of s_fluct2 beyond the tolerance level for fs1?
            s_fluct2 = s_fluct2+(fs1-sum(s_fluct2(s_fluct2>=0))/length(s_fluct2)); % if yes, add the deviation and re-evaluate (ignore negatives)
            kk = kk + 1; % to check how many times it have to iterate
        end
        s_fluct2(s_fluct2<0)=0; % set negative values = 0
        s_fluct2 = s_fluct2.*(sigma_fs1/std(s_fluct2)); % adjust the std
        s_fluct = s_fluct2; 

    case {'spatial_random','both_random'} % both_random to be included in the future

        %tobefixed = ones(size(X)); % not needed?
        s_fluct = zeros(size(X));
        criteria_B = abs(mean(s_fluct(X==1))-fs1) > tol_fs1; % mark if beyond tolerance, recalculate
        criteria_C = abs(std(s_fluct(X==1))-sigma_fs1) > tol_sigma_fs1; % mark if beyond tolerance, recalculate

        mu = log((fs1^2)/sqrt(sigma_fs1^2+fs1^2)); sigma = sqrt(log(sigma_fs1^2/(fs1^2)+1)); % transform lognormal params to normal params

        kk = 0;
        while criteria_B || criteria_C % if the desirable conditions are not met
            
%             s_fluct = s_fluct.*(tobefixed==0) + (tobefixed==1 & X==1).*abs(lognrnd(mu,sigma,size(X))); % draw values from the lognormal distribution
% 
%             criteria_B = abs(mean(s_fluct(X==1))-fs1) > tol_fs1; % mark if beyond tolerance
%             criteria_C = abs(std(s_fluct(X==1))-sigma_fs1) > tol_sigma_fs1; % mark if beyond tolerance

            s_fluct = (X==1).*abs(lognrnd(mu,sigma,size(X))); % draw values from the lognormal distribution

            criteria_B = abs(mean(s_fluct(X==1))-fs1) > tol_fs1; % mark if beyond tolerance
            criteria_C = abs(std(s_fluct(X==1))-sigma_fs1) > tol_sigma_fs1; % mark if beyond tolerance

            kk = kk + 1; % to check how many times it had to try to meet the tolerance levels
        end
        s_fluct = s_fluct(:); % this s_fluct was 2D but linearized for indexing

    case {'spatial_lineargrad','both_lineargrad'} % both_lineargrad to be included in the future
        s_fluct = zeros(size(X));
        criteria_B = abs(mean(s_fluct(X==1))-fs1) > tol_fs1;
        criteria_C = abs(std(s_fluct(X==1))-sigma_fs1) > tol_sigma_fs1;

        mu = log((fs1^2)/sqrt(sigma_fs1^2+fs1^2)); sigma = sqrt(log(sigma_fs1^2/(fs1^2)+1));

        kk = 0;
        while criteria_B || criteria_C % if the desirable conditions are not met
         
            s_fluct = (X==1).*abs(lognrnd(mu,sigma,size(X)));

            criteria_B = abs(mean(s_fluct(X==1))-fs1) > tol_fs1;
            criteria_C = abs(std(s_fluct(X==1))-sigma_fs1) > tol_sigma_fs1;

            kk = kk + 1;
        end

        % order the values to create a linear gradient
        s_fluct = reshape(sort(s_fluct(:)),size(X));
        s_fluct = s_fluct(:); % this s_fluct was 2D but linearized for indexing

    case 'static'

        s_fluct = ones(1,num_gen).*fs1;

end

elseif strcmp(selection,'viability') 

    % to be filled in the future

end

end