function [Wasps] = wasp_dispersal(Wasps, c, r, L, Pr_disp, Pr_d, female_mated)

if nargin == 6 %natal dispersal
    who_disp = (rand(size(Wasps,1),1) < Pr_disp); % dispersers = 1
elseif nargin == 7 % breeding dispersal
    who_disp = (female_mated == 0); % dispersers = non-mated females and all males
end

trees_w_disp_wasps = unique(Wasps(who_disp ==1,3)); % find trees containing dispersing indivs
new_home = nan(size(Wasps,1),1);
for i = trees_w_disp_wasps'
    new_Trees = randsample(L^2,sum((who_disp==1)&Wasps(:,3)==i),'true',Pr_d(i,:)); % probabilistically choose a destination patch
    new_home((who_disp == 1)&Wasps(:,3)==i) = new_Trees;
end

% natal dispersal
Wasps(who_disp==1,1:2) = [c(new_home(who_disp==1)),r(new_home(who_disp==1))]; % impute new coordinates
Wasps(who_disp==1,3) = new_home(who_disp==1); % impute new patch ids

end