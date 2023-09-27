function [Wasps,female_mated,num_breed_disp] = wasp_mating(Wasps, c, r, L, Pr_disp, Pr_d, max_breed_disp)

% females mate only once but males can mate indefinitely many times. all females mate, not all males can mate
% Select a mate at the local site. If cannot find a mate, disperse again.
% Once mated, females do not move until all females are mated. Males can
% keep moving until all the females are mated. 

num_breed_disp = 0;
female_mated = zeros(size(Wasps,1),1);

% random mating
while sum(female_mated(Wasps(:,6)==1)==0)>0 && num_breed_disp <= max_breed_disp
    trees_w_wasps = unique(Wasps(Wasps(:,6)==1 & female_mated==0,3)); % trees with unmated females
    for i = trees_w_wasps'
        f = find(((Wasps(:,3)==i) & (Wasps(:,6)==1))==1 & female_mated==0); % 1 = on this tree, not mated yet
        m = find((Wasps(:,3)==i) & (Wasps(:,6)==1)==0);
        if ~isempty(f) && ~isempty(m)
            if length(m)==1
                f_mates(1:length(f)) = m;
            else
                f_mates = randsample(m,length(f),'true');
            end
            female_mated(f) = f_mates; % enter male id in female_mated
        end
        f_mates = [];
    end

    % males disperse after mating until all females are mated, mated
    % females do not move anymore after mating
    % males and unmated females disperse
    if sum(female_mated(Wasps(:,6)==1)==0)>0
    num_breed_disp = num_breed_disp + 1;
    [Wasps] = wasp_dispersal(Wasps, c, r, L, Pr_disp, Pr_d, female_mated);
    end
end

end
