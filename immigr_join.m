function [Wasps] = immigr_join(Wasps, L, N2_immig, N1_immig, c, r, t, freq_joining)

if mod(t,freq_joining) == 0 % is this the time step when immigrants join?
    immig_wasps = zeros(N2_immig,6);
    immig_wasps(:,1:2) = randi(L,N2_immig,2); % x,y, randomly distribute
    immig_wasps(:,4:5) = 0; % all alleles are foreign (0)
    immig_wasps(:,6) = rand(size(immig_wasps,1),1)<0.5; % assign sex (1 = females)
    [~, Locb] = ismember(immig_wasps(:,1:2),[c,r],'rows'); % which tree?
    immig_wasps(:,3) = Locb; % tree id
    Wasps = [Wasps;immig_wasps]; % append immigrants to existing indivs
end

if t > 1 % new residents flowing into this area from the core habitat
    new_resident_wasps = zeros(N1_immig,6);
    new_resident_wasps(:,1:2) = randi(L,N1_immig,2); % x,y, randomly distribute
    new_resident_wasps(:,4:5) = 1; % all alleles are resident (1)
    new_resident_wasps(:,6) = rand(size(new_resident_wasps,1),1)<0.5; % assign sex (1 = females)
    [~, Locb] = ismember(new_resident_wasps(:,1:2),[c,r],'rows'); % which tree?
    new_resident_wasps(:,3) = Locb; % tree id
    Wasps = [Wasps;new_resident_wasps]; % append immigrants
end

end