%initial setup
clc, clearvars, close all;
file_path3 = 'C:\Users\USER\Downloads\PSCompPars_2024.01.11_08.36.44.csv';
Data = readtable(file_path3);


%creating key parameters for filtering the data
st_temp = Data.st_teff;
disc_facility = Data.disc_facility;
pl_axis = Data.pl_orbsmax;
disc_method = Data.discoverymethod;

%using the key parameters to adjust the filtered data
filter_condition = (st_temp>4500 & st_temp<6500 & pl_axis<2 & strcmp(disc_method, 'Transit') & strcmp(disc_facility, 'Kepler'));
st_temp = st_temp(filter_condition);
disc_facility = disc_facility(filter_condition);
pl_axis = pl_axis(filter_condition);
disc_method = disc_method(filter_condition);

%using the key parameters to adjust all other filtered parameters
st_met = Data.st_met(filter_condition);
st_mass = Data.st_mass(filter_condition);
st_age = Data.st_age(filter_condition);
pl_orb_per = Data.pl_orbper(filter_condition);
pl_mass = Data.pl_bmassj(filter_condition);
pl_orb_ecc = Data.pl_orbeccen(filter_condition);
pl_radius = Data.pl_radj(filter_condition);
num_planets = Data.sy_pnum(filter_condition);
st_name = Data.hostname(filter_condition);
pl_name = Data.pl_name(filter_condition);









% Creating the system parameters (size and orbital period ratios):
SystemHosts = unique(st_name(num_planets>1));

%ratios for 2 planets systems
size_ratios_2pl = zeros(1, length(unique(st_name(num_planets==2)))-10);
orb_per_ratios_2pl = zeros(1, length(unique(st_name(num_planets==2)))-10);

%ratios for 3 planets systems
size_ratios_3pl_MinMid = zeros(1, length(unique(st_name(num_planets==3)))-10);
orb_per_ratios_3pl_MinMid = zeros(1, length(unique(st_name(num_planets==3)))-10);

size_ratios_3pl_MidMax = zeros(1, length(unique(st_name(num_planets==3)))-10);
orb_per_ratios_3pl_MidMax = zeros(1, length(unique(st_name(num_planets==3)))-10);

size_ratios_3pl_MinMax = zeros(1, length(unique(st_name(num_planets==3)))-10);
orb_per_ratios_3pl_MinMax = zeros(1, length(unique(st_name(num_planets==3)))-10);

%creating the X parameters for the systems:
SystemMet_2Pl = zeros(1, length(unique(st_name(num_planets==2)))-10);
SystemMet_3Pl = zeros(1, length(unique(st_name(num_planets==3)))-10);

SystemAge_2Pl = zeros(1, length(unique(st_name(num_planets==2)))-10);
SystemAge_3Pl = zeros(1, length(unique(st_name(num_planets==3)))-10);

SystemMass_2Pl = zeros(1, length(unique(st_name(num_planets==2)))-10);
SystemMass_3Pl = zeros(1, length(unique(st_name(num_planets==3)))-10);

k1 = 0;
k2 = 0;

for i = 1:height(SystemHosts)
    SystemHostsData = Data(strcmp(SystemHosts(i), Data.hostname), :);

    SystemPlRad = SystemHostsData.pl_radj;
    SystemOrbPer = SystemHostsData.pl_orbper;

    if height(SystemHostsData) == 2
        k1 = k1+1;

        size_ratios_2pl(k1) = SystemPlRad(find(SystemOrbPer == min(SystemOrbPer))) / SystemPlRad(find(SystemOrbPer == max(SystemOrbPer)));
        orb_per_ratios_2pl(k1) = min(SystemOrbPer) / max(SystemOrbPer);
      
        SystemMet_2Pl(k1) = SystemHostsData.st_met(1);
        SystemAge_2Pl(k1) = SystemHostsData.st_age(1);
        SystemMass_2Pl(k1) = SystemHostsData.st_mass(1);
    end

    if height(SystemHostsData) == 3
        sort(SystemOrbPer);
        k2 = k2+1;

        size_ratios_3pl_MinMid(k2) = SystemPlRad(find(SystemOrbPer == SystemOrbPer(1))) / SystemPlRad(find(SystemOrbPer == SystemOrbPer(2)));
        orb_per_ratios_3pl_MinMid(k2) = SystemOrbPer(1) / SystemOrbPer(2);
     
        size_ratios_3pl_MidMax(k2) = SystemPlRad(find(SystemOrbPer == SystemOrbPer(2))) / SystemPlRad(find(SystemOrbPer == SystemOrbPer(3)));
        orb_per_ratios_3pl_MidMax(k2) = SystemOrbPer(2) / SystemOrbPer(3);

        size_ratios_3pl_MinMax(k2) = SystemPlRad(find(SystemOrbPer == SystemOrbPer(1))) / SystemPlRad(find(SystemOrbPer == SystemOrbPer(3)));
        orb_per_ratios_3pl_MinMax(k2) = SystemOrbPer(1) / SystemOrbPer(3);

        
        SystemMet_3Pl(k2) = SystemHostsData.st_met(1);
        SystemAge_3Pl(k2) = SystemHostsData.st_age(1);
        SystemMass_3Pl(k2) = SystemHostsData.st_mass(1);
    end
end


%controling all of the figures:
set(0,'DefaultFigureVisible','off');







%advanced tests with a moving threshold with a dependence on METALLICITY:
pl_radius_distributions_met = zeros(1,1000);
pl_num_distribution_met = zeros(1,1000);

size_ratios_2pl_distribution_met = zeros(1,1000);
size_ratios_3pl_MinMid_distribution_met = zeros(1,1000);
size_ratios_3pl_MidMax_distribution_met = zeros(1,1000);
size_ratios_3pl_MinMax_distribution_met = zeros(1,1000);

orb_per_ratios_2pl_distribution_met = zeros(1,1000);
orb_per_ratios_3pl_MinMid_distribution_met = zeros(1,1000);
orb_per_ratios_3pl_MidMax_distribution_met = zeros(1,1000);
orb_per_ratios_3pl_MinMax_distribution_met = zeros(1,1000);

% ks tests on the planetery parameters
k1 = 0;
for i = linspace(min(st_met), max(st_met), 1000)
    k1 = k1+1;

    if length(pl_radius(st_met>=i))>=100 && length(pl_radius(st_met<=i))>=100
        [~,b] = kstest2(pl_radius(st_met>=i), (pl_radius(st_met<=i)));
        pl_radius_distributions_met(k1) = b;

        [~,b] = kstest2(num_planets(st_met>=i), (num_planets(st_met<=i)));
        pl_num_distribution_met(k1) = b;
    end
end

% ks tests on the 2 planet systems
k1 = 0;
for i = linspace(min(SystemMet_2Pl), max(SystemMet_2Pl), 1000)
    k1 = k1+1;

    if length(size_ratios_2pl(SystemMet_2Pl>=i))>=50 && length(size_ratios_2pl(SystemMet_2Pl<=i))>=50
        [~,b] = kstest2(size_ratios_2pl(SystemMet_2Pl>=i), (size_ratios_2pl(SystemMet_2Pl<=i)));
        size_ratios_2pl_distribution_met(k1) = b;

        [~,b] = kstest2(orb_per_ratios_2pl(SystemMet_2Pl>=i), (orb_per_ratios_2pl(SystemMet_2Pl<=i)));
        orb_per_ratios_2pl_distribution_met(k1) = b;
    end
end

% ks tests on the 3 planet systems
k1 = 0;
for i = linspace(min(SystemMet_3Pl), max(SystemMet_3Pl), 1000)
    k1 = k1+1;
    
    if true%length(size_ratios_3pl_MinMid(SystemMet_3Pl>=i))>=10 && length(size_ratios_3pl_MinMid(SystemMet_3Pl<=i))>=10

    % MinMid
    [~,b] = kstest2(size_ratios_3pl_MinMid(SystemMet_3Pl>=i), (size_ratios_3pl_MinMid(SystemMet_3Pl<=i)));
    size_ratios_3pl_MinMid_distribution_met(k1) = b;

    [~,b] = kstest2(orb_per_ratios_3pl_MinMid(SystemMet_3Pl>=i), (orb_per_ratios_3pl_MinMid(SystemMet_3Pl<=i)));
    orb_per_ratios_3pl_MinMid_distribution_met(k1) = b;

    % MidMax
    [~,b] = kstest2(size_ratios_3pl_MidMax(SystemMet_3Pl>=i), (size_ratios_3pl_MidMax(SystemMet_3Pl<=i)));
    size_ratios_3pl_MidMax_distribution_met(k1) = b;

    [~,b] = kstest2(orb_per_ratios_3pl_MidMax(SystemMet_3Pl>=i), (orb_per_ratios_3pl_MidMax(SystemMet_3Pl<=i)));
    orb_per_ratios_3pl_MidMax_distribution_met(k1) = b;

    % MinMax
    [~,b] = kstest2(size_ratios_3pl_MinMax(SystemMet_3Pl>=i), (size_ratios_3pl_MinMax(SystemMet_3Pl<=i)));
    size_ratios_3pl_MinMax_distribution_met(k1) = b;

    [~,b] = kstest2(orb_per_ratios_3pl_MinMax(SystemMet_3Pl>=i), (orb_per_ratios_3pl_MinMax(SystemMet_3Pl<=i)));
    orb_per_ratios_3pl_MinMax_distribution_met(k1) = b;

    end
end

% plotting the tests results:
ksPlot1var(pl_radius_distributions_met, st_met, 'Stellar metallicity', 'Planet radius KS')
ksPlot1var(pl_num_distribution_met, st_met, 'Stellar metallicity', 'Planet number KS')
ksPlot1var(size_ratios_2pl_distribution_met, SystemMet_2Pl, 'Stellar metallicity', '2 planet size ratios KS')
ksPlot1var(orb_per_ratios_2pl_distribution_met, SystemMet_2Pl, 'Stellar metallicity', '2 planet orbital period ratios KS')

ksPlot3var(size_ratios_3pl_MinMid_distribution_met, size_ratios_3pl_MidMax_distribution_met, size_ratios_3pl_MinMax_distribution_met, SystemMet_3Pl, 'Stellar metallicity', '3 planet size ratios KS')
ksPlot3var(orb_per_ratios_3pl_MinMid_distribution_met, orb_per_ratios_3pl_MidMax_distribution_met, orb_per_ratios_3pl_MinMax_distribution_met, SystemMet_3Pl, 'Stellar metallicity', '3 planet orbital period ratios KS')







%advanced tests with a moving threshold with a dependence on Age:
pl_radius_distributions_age = zeros(1,1000);
pl_num_distribution_age = zeros(1,1000);

size_ratios_2pl_distribution_age = zeros(1,1000);
size_ratios_3pl_MinMid_distribution_age = zeros(1,1000);
size_ratios_3pl_MidMax_distribution_age = zeros(1,1000);
size_ratios_3pl_MinMax_distribution_age = zeros(1,1000);

orb_per_ratios_2pl_distribution_age = zeros(1,1000);
orb_per_ratios_3pl_MinMid_distribution_age = zeros(1,1000);
orb_per_ratios_3pl_MidMax_distribution_age = zeros(1,1000);
orb_per_ratios_3pl_MinMax_distribution_age = zeros(1,1000);

high_met_num_dis = zeros(1,1000);
low_met_num_dis = zeros(1,1000);


crit_met = -0.1;
low_met_age=(st_age(logical(st_met<crit_met)));
high_met_age = (st_age(logical(st_met>crit_met)));
low_met_num=(num_planets(st_met<crit_met));
high_met_num = (num_planets(st_met>crit_met));

% ks tests on the planetery parameters
k1 = 0;
for i = linspace(min(st_age), max(st_age), 1000)
    k1 = k1+1;

    if length(pl_radius(st_age>=i))>=100 && length(pl_radius(st_age<=i))>=100
        [~,b] = kstest2(pl_radius(st_age>=i), (pl_radius(st_age<=i)));
        pl_radius_distributions_age(k1) = b;

        [~,b] = kstest2(num_planets(st_age>=i), (num_planets(st_age<=i)));
        pl_num_distribution_age(k1) = b;

        [~,b] = kstest2(high_met_num(high_met_age>=i), (high_met_num(high_met_age<=i)));
        high_met_num_dis(k1) = b;

        [~,b] = kstest2(low_met_num(low_met_age>=i), (low_met_num(low_met_age<=i)));
        low_met_num_dis(k1) = b;
    end
end

% ks tests on the 2 planet systems
k1 = 0;
for i = linspace(min(SystemAge_2Pl), max(SystemAge_2Pl), 1000)
    k1 = k1+1;

    if length(size_ratios_2pl(SystemAge_2Pl>=i))>=50 && length(size_ratios_2pl(SystemAge_2Pl<=i))>=50

        [~,b] = kstest2(size_ratios_2pl(SystemAge_2Pl>=i), (size_ratios_2pl(SystemAge_2Pl<=i)));
        size_ratios_2pl_distribution_age(k1) = b;

        [~,b] = kstest2(orb_per_ratios_2pl(SystemAge_2Pl>=i), (orb_per_ratios_2pl(SystemAge_2Pl<=i)));
        orb_per_ratios_2pl_distribution_age(k1) = b;

    end
end

% ks tests on the 3 planet systems
k1 = 0;
for i = linspace(min(SystemAge_3Pl), max(SystemAge_3Pl), 1000)
    k1 = k1+1;
    
    % MinMid
    [~,b] = kstest2(size_ratios_3pl_MinMid(SystemAge_3Pl>=i), (size_ratios_3pl_MinMid(SystemAge_3Pl<=i)));
    size_ratios_3pl_MinMid_distribution_age(k1) = b;

    [~,b] = kstest2(orb_per_ratios_3pl_MinMid(SystemAge_3Pl>=i), (orb_per_ratios_3pl_MinMid(SystemAge_3Pl<=i)));
    orb_per_ratios_3pl_MinMid_distribution_age(k1) = b;

    % MidMax
    [~,b] = kstest2(size_ratios_3pl_MidMax(SystemAge_3Pl>=i), (size_ratios_3pl_MidMax(SystemAge_3Pl<=i)));
    size_ratios_3pl_MidMax_distribution_age(k1) = b;

    [~,b] = kstest2(orb_per_ratios_3pl_MidMax(SystemAge_3Pl>=i), (orb_per_ratios_3pl_MidMax(SystemAge_3Pl<=i)));
    orb_per_ratios_3pl_MidMax_distribution_age(k1) = b;

    % MinMax
    [~,b] = kstest2(size_ratios_3pl_MinMax(SystemAge_3Pl>=i), (size_ratios_3pl_MinMax(SystemAge_3Pl<=i)));
    size_ratios_3pl_MinMax_distribution_age(k1) = b;

    [~,b] = kstest2(orb_per_ratios_3pl_MinMax(SystemAge_3Pl>=i), (orb_per_ratios_3pl_MinMax(SystemAge_3Pl<=i)));
    orb_per_ratios_3pl_MinMax_distribution_age(k1) = b;
end



ksPlot1var(pl_radius_distributions_age, st_age, 'Stellar age', 'Planet radius KS')
ksPlot1var(pl_num_distribution_age, st_age, 'Stellar age', 'Planet number KS')
ksPlot1var(size_ratios_2pl_distribution_age, SystemAge_2Pl, 'Stellar age', '2 planet size ratios KS')
ksPlot1var(orb_per_ratios_2pl_distribution_age, SystemAge_2Pl, 'Stellar age', '2 planet orbital period ratios KS')

ksPlot3var(size_ratios_3pl_MinMid_distribution_age, size_ratios_3pl_MidMax_distribution_age, size_ratios_3pl_MinMax_distribution_age, SystemAge_3Pl, 'Stellar age', '3 planet size ratios KS')
ksPlot3var(orb_per_ratios_3pl_MinMid_distribution_age, orb_per_ratios_3pl_MidMax_distribution_age, orb_per_ratios_3pl_MinMax_distribution_age, SystemAge_3Pl, 'Stellar age', '3 planet orbital period ratios KS')




%4 variable tests to check for metallicity bias influencing age tests:
figure();
plot(linspace(min(high_met_age), max(high_met_age), 1000), (high_met_num_dis),LineWidth=8)
set(gca, "YScale", "log","FontSize",25);
xlabel('Stellar age','FontSize',30);
ylabel('Number of planets KS',FontSize=30);
title('Number of planets KS for high metallicities',FontSize=30);


figure();
plot(linspace(min(low_met_age), max(low_met_age), 1000), (low_met_num_dis),LineWidth=8)
set(gca, "YScale", "log","FontSize",25);
xlabel('Stellar age','FontSize',30);
ylabel('Number of planets KS',FontSize=30);
title('Number of planets KS for low metallicities',FontSize=30);





%advanced tests with a moving threshold with a dependence on MASS:
pl_radius_distributions_mass = zeros(1,1000);
pl_num_distribution_mass = zeros(1,1000);

size_ratios_2pl_distribution_mass = zeros(1,1000);
size_ratios_3pl_MinMid_distribution_mass = zeros(1,1000);
size_ratios_3pl_MidMax_distribution_mass = zeros(1,1000);
size_ratios_3pl_MinMax_distribution_mass = zeros(1,1000);

orb_per_ratios_2pl_distribution_mass = zeros(1,1000);
orb_per_ratios_3pl_MinMid_distribution_mass = zeros(1,1000);
orb_per_ratios_3pl_MidMax_distribution_mass = zeros(1,1000);
orb_per_ratios_3pl_MinMax_distribution_mass = zeros(1,1000);

average_3pl_orb_per_ratios_distribution_mass = zeros(1,1000);

% ks tests on the planetery parameters
k1 = 0;
for i = linspace(min(st_mass), max(st_mass), 1000)
    k1 = k1+1;

    if length(pl_radius(st_mass>=i))>=100 && length(pl_radius(st_mass<=i))>=100
        [~,b] = kstest2(pl_radius(st_mass>=i), (pl_radius(st_mass<=i)));
        pl_radius_distributions_mass(k1) = b;

        [~,b] = kstest2(num_planets(st_mass>=i), (num_planets(st_mass<=i)));
        pl_num_distribution_mass(k1) = b;
    end
end

% ks tests on the 2 planet systems
k1 = 0;
for i = linspace(min(SystemMass_2Pl), max(SystemMass_2Pl), 1000)
    k1 = k1+1;

    if length(size_ratios_2pl(SystemMass_2Pl>=i))>=50 && length(size_ratios_2pl(SystemMass_2Pl<=i))>=50

        [~,b] = kstest2(size_ratios_2pl(SystemMass_2Pl>=i), (size_ratios_2pl(SystemMass_2Pl<=i)));
        size_ratios_2pl_distribution_mass(k1) = b;

        [~,b] = kstest2(orb_per_ratios_2pl(SystemMass_2Pl>=i), (orb_per_ratios_2pl(SystemMass_2Pl<=i)));
        orb_per_ratios_2pl_distribution_mass(k1) = b;
    end
end

% ks tests on the 3 planet systems
k1 = 0;
for i = linspace(min(SystemMass_3Pl), max(SystemMass_3Pl), 1000)
    k1 = k1+1;
    
    % MinMid
    [~,b] = kstest2(size_ratios_3pl_MinMid(SystemMass_3Pl>=i), (size_ratios_3pl_MinMid(SystemMass_3Pl<=i)));
    size_ratios_3pl_MinMid_distribution_mass(k1) = b;

    [~,b] = kstest2(orb_per_ratios_3pl_MinMid(SystemMass_3Pl>=i), (orb_per_ratios_3pl_MinMid(SystemMass_3Pl<=i)));
    orb_per_ratios_3pl_MinMid_distribution_mass(k1) = b;

    % MidMax
    [~,b] = kstest2(size_ratios_3pl_MidMax(SystemMass_3Pl>=i), (size_ratios_3pl_MidMax(SystemMass_3Pl<=i)));
    size_ratios_3pl_MidMax_distribution_mass(k1) = b;

    [~,b] = kstest2(orb_per_ratios_3pl_MidMax(SystemMass_3Pl>=i), (orb_per_ratios_3pl_MidMax(SystemMass_3Pl<=i)));
    orb_per_ratios_3pl_MidMax_distribution_mass(k1) = b;

    % MinMax
    [~,b] = kstest2(size_ratios_3pl_MinMax(SystemMass_3Pl>=i), (size_ratios_3pl_MinMax(SystemMass_3Pl<=i)));
    size_ratios_3pl_MinMax_distribution_mass(k1) = b;

    [~,b] = kstest2(orb_per_ratios_3pl_MinMax(SystemMass_3Pl>=i), (orb_per_ratios_3pl_MinMax(SystemMass_3Pl<=i)));
    orb_per_ratios_3pl_MinMax_distribution_mass(k1) = b;

    %just expiramenting
    average_3pl_orb_per_ratios_distribution_mass(k1) = mean([orb_per_ratios_3pl_MinMid_distribution_mass(k1),orb_per_ratios_3pl_MidMax_distribution_mass(k1),orb_per_ratios_3pl_MinMax_distribution_mass(k1)]);
end

% plotting the tests results:
ksPlot1var(pl_radius_distributions_mass, st_mass, 'Stellar mass', 'Planet radius KS')
ksPlot1var(pl_num_distribution_mass, st_mass, 'Stellar mass', 'Planet number KS')
ksPlot1var(size_ratios_2pl_distribution_mass, SystemMass_2Pl, 'Stellar mass', '2 planet size ratios KS')
ksPlot1var(orb_per_ratios_2pl_distribution_mass, SystemMass_2Pl, 'Stellar mass', '2 planet orbital period ratios KS')

ksPlot3var(size_ratios_3pl_MinMid_distribution_mass, size_ratios_3pl_MidMax_distribution_mass, size_ratios_3pl_MinMax_distribution_mass, SystemMass_3Pl, 'Stellar mass', '3 planet size ratios KS')
ksPlot3var(orb_per_ratios_3pl_MinMid_distribution_mass, orb_per_ratios_3pl_MidMax_distribution_mass, orb_per_ratios_3pl_MinMax_distribution_mass, SystemMass_3Pl, 'Stellar mass', '3 planet orbital period ratios KS')







%Basic tests, scatter plots and some hists:

%controling all of the figures:
set(0,'DefaultFigureVisible','off');



%all threshold conditions:
metal_rich_condition = st_met > 0;
metal_poor_condition = st_met < 0;
median_st_age = median(st_age(~isnan(st_age)));
median_st_mass = median(st_mass(~isnan(st_mass)));



%basic tests when the dependence is on stellar metallicity:
figure;
scatter(st_met, pl_radius, 50, 'b', 'filled');
set(gca, "YScale", "log", "FontSize", 25);
xlabel('Stellar metallicity', FontSize=30);
ylabel('Planet radius', FontSize=30);
%title('Planet radius vs Stellar Metallicity', FontSize=30);



figure;
[n,x] = hist(num_planets(metal_rich_condition));
bar(n/sum(n));
hold on;
[n,x] = hist(num_planets(metal_poor_condition));
bar(n/sum(n));
hold off


figure;
[n,x] = hist(num_planets(metal_rich_condition));
bar(x,n/sum(n));
hold on
[n,x] = hist(num_planets(metal_poor_condition));
bar(x,n/sum(n));
xlabel('Number of Planets in System');
ylabel('Fraction');
title('Distribution of Number of Planets');
legend('Metal-Rich Systems', 'Metal-Poor Systems');
hold off;




figure
hold on
[n,x]= hist(pl_radius(metal_rich_condition),30);
bar(x,n/sum(n))
[n2,~] = hist(pl_radius(metal_poor_condition),30);
bar(x,n2/sum(n2))




figure
hold on
[n,x]= hist(pl_radius(st_met>-0.1),30);
[n2,~] = hist(pl_radius(st_met<-0.1),30);
bar(x,n/sum(n+n2))
bar(x,n2/sum(n+n2), 'BarWidth',0.6)




figure
hold on
[n,x]= hist(num_planets(metal_rich_condition),8);
[n2,~] = hist(num_planets(metal_poor_condition),8);
bar(n/sum(n+n2))
bar(n2/sum(n+n2), 'BarWidth',0.6)





figure();
scatter(SystemMet_2Pl, size_ratios_2pl, 50, 'b', 'filled')
set(gca, "YScale", "log","FontSize",25);
xlabel('Stellar metallicity','FontSize',30);
ylabel('Planet radius ratio','FontSize',30);
%title('Planet Radius Ratio of 2 Planet Systems vs Stellar Metallicity',FontSize=30);

figure();
scatter(SystemMet_2Pl, orb_per_ratios_2pl, 50, 'b', 'filled')
set(gca, "YScale", "log","FontSize",25);
xlabel('Stellar metallicity','FontSize',30);
ylabel('Orbital period ratio','FontSize',30);
%title('Orbital Period Ratio of 2 Planet Systems vs Stellar Metallicity',FontSize=30);


figure();
scatter(SystemMet_3Pl, size_ratios_3pl_MinMid, 50, 'blue', 'filled')
hold on;
scatter(SystemMet_3Pl, size_ratios_3pl_MidMax,50, 'red', 'filled')
hold on;
scatter(SystemMet_3Pl, size_ratios_3pl_MinMax, 50, 'green', 'filled')
set(gca, "YScale", "log","FontSize",25);
xlabel('Stellar metallicity','FontSize',30);
ylabel('Planet radius ratio',FontSize=30);
legend('Closest planet/Mid planet', 'Mid planet/Furthest planet', 'Closest planet/Furthest planet',fontsize=15);
%title('Planet Radius Ratio of 3 Planet Systems vs Stellar Metallicity',FontSize=30);
hold off;

figure();
scatter(SystemMet_3Pl, orb_per_ratios_3pl_MinMid, 50, 'blue', 'filled')
hold on;
scatter(SystemMet_3Pl, orb_per_ratios_3pl_MidMax, 50, 'red', 'filled')
hold on;
scatter(SystemMet_3Pl, orb_per_ratios_3pl_MinMax, 50, 'green', 'filled')
set(gca, "YScale", "log","FontSize",25);
xlabel('Stellar metallicity','FontSize',30);
ylabel('Orbital period ratio','FontSize',30);
legend('Closest planet/Mid planet', 'Mid planet/Furthest planet', 'Closest planet/Furthest planet',fontsize=15);
%title('Orbital Period Ratio of 3 Planet Systems vs Stellar Metallicity',FontSize=30);
hold off;



%basic tests when the dependence is on stellar AGE
figure;
scatter(st_age, pl_radius, 50, 'b', 'filled');
set(gca, "YScale", "log","FontSize",25);
xlabel('Stellar age','FontSize',30);
ylabel('Planet radius','FontSize',30);
%title('Planet Radius vs Stellar Age',FontSize=30);

figure;
[n,x] = hist(num_planets(st_age>median_st_age));
bar(x, n/sum(n));
hold on;
[n2,x] = hist(num_planets(st_age<median_st_age));
bar(x, n2/sum(n2));
xlabel('Number of Planets in System');
ylabel('Amount');
title('Distribution of Number of Planets');
legend('Older Stars', 'Younger Stars');
hold off;

figure();
scatter(SystemAge_2Pl, size_ratios_2pl, 50, 'b', 'filled')
set(gca, "YScale", "log","FontSize",30);
xlabel('Stellar age','FontSize',30);
ylabel('Planet radius ratio','FontSize',30);
%title('Planet Radius Ratio of 2 Planet Systems vs Stellar Age',FontSize=30);

figure();
scatter(SystemAge_2Pl, orb_per_ratios_2pl, 50, 'b', 'filled')
set(gca, "YScale", "log","FontSize",25);
xlabel('Stellar age','FontSize',30);
ylabel('Orbital period ratio','FontSize',30);
%title('Orbital Period Ratio of 2 Planet Systems vs Stellar Age',FontSize=30);


figure();
scatter(SystemAge_3Pl, size_ratios_3pl_MinMid, 50, 'blue', 'filled')
hold on;
scatter(SystemAge_3Pl, size_ratios_3pl_MidMax,50, 'red', 'filled')
hold on;
scatter(SystemAge_3Pl, size_ratios_3pl_MinMax, 50, 'green', 'filled')
set(gca, "YScale", "log","FontSize",25);
xlabel('Stellar age','FontSize',30);
ylabel('Planet radius ratio','FontSize',30);
legend('Closest planet/Mid planet', 'Mid planet/Furthest planet', 'Closest planet/Furthest planet',fontsize=15);
%title('Planet Radius Ratio of 3 Planet Systems vs Stellar Age',FontSize=30);
hold off;

figure();
scatter(SystemAge_3Pl, orb_per_ratios_3pl_MinMid, 50, 'blue', 'filled')
hold on;
scatter(SystemAge_3Pl, orb_per_ratios_3pl_MidMax, 50, 'red', 'filled')
hold on;
scatter(SystemAge_3Pl, orb_per_ratios_3pl_MinMax, 50, 'green', 'filled')
set(gca, "YScale", "log","FontSize",25);
xlabel('Stellar age','FontSize',30);
ylabel('Orbital period ratio','FontSize',30);
legend('Closest planet/Mid planet', 'Mid planet/Furthest planet', 'Closest planet/Furthest planet',fontsize=15);
%title('Orbital Period Ratio of 3 Planet Systems vs Stellar Age',FontSize=30);
hold off;





%basic tests when the dependence is on stellar MASS:
figure;
scatter(st_mass, pl_radius, 50, 'b', 'filled');
set(gca, "YScale", "log","FontSize",25);
xlabel('Stellar mass','FontSize',30);
ylabel('Planet radius','FontSize',30);
%title('Planet Radius vs Stellar Mass',FontSize=30);

figure;
[n,x] = hist(num_planets(st_mass>median_st_mass));
bar(x, n/sum(n));
hold on;
[n2,x] = hist(num_planets(st_mass<median_st_mass));
bar(x, n2/sum(n2));
xlabel('Number of Planets in System');
ylabel('Amount');
title('Distribution of Number of Planets');
legend('Massive Stars', 'Less Massive Stars');
hold off;

figure();
scatter(SystemMass_2Pl, size_ratios_2pl, 50, 'b', 'filled')
set(gca, "YScale", "log","FontSize",30);
xlabel('Stellar mass','FontSize',30);
ylabel('Planet radius ratio','FontSize',30);
%title('Planet Radius Ratio of 2 Planet Systems vs Stellar Mass',FontSize=30);

figure();
scatter(SystemMass_2Pl, orb_per_ratios_2pl, 50, 'b', 'filled')
set(gca, "YScale", "log","FontSize",25);
xlabel('Stellar mass','FontSize',30);
ylabel('Orbital period ratio','FontSize',30);
%title('Orbital Period Ratio of 2 Planet Systems vs Stellar Mass',FontSize=30);


figure();
scatter(SystemMass_3Pl, size_ratios_3pl_MinMid, 50, 'blue', 'filled')
hold on;
scatter(SystemMass_3Pl, size_ratios_3pl_MidMax,50, 'red', 'filled')
hold on;
scatter(SystemMass_3Pl, size_ratios_3pl_MinMax, 50, 'green', 'filled')
set(gca, "YScale", "log","FontSize",25);
xlabel('Stellar mass','FontSize',30);
ylabel('Planet radius ratio','FontSize',30);
legend('Closest planet/Mid planet', 'Mid planet/Furthest planet', 'Closest planet/Furthest planet',fontsize=15);
%title('Planet Radius Ratio of 3 Planet Systems vs Stellar Mass',FontSize=30);
hold off;

figure();
scatter(SystemMass_3Pl, orb_per_ratios_3pl_MinMid, 50, 'blue', 'filled')
hold on;
scatter(SystemMass_3Pl, orb_per_ratios_3pl_MidMax, 50, 'red', 'filled')
hold on;
scatter(SystemMass_3Pl, orb_per_ratios_3pl_MinMax, 50, 'green', 'filled')
set(gca, "YScale", "log","FontSize",25);
xlabel('Stellar mass',FontSize=30);
ylabel('Orbital period ratio','FontSize',30);
legend('Closest planet/Mid planet', 'Mid planet/Furthest planet', 'Closest planet/Furthest planet',fontsize=15);
%title('Orbital Period Ratio of 3 Planet Systems vs Stellar Mass',FontSize=30);
hold off;






%some 4 variable test:

%controling all of the figures:
set(0,'DefaultFigureVisible','off');


%all threshold conditions:
metal_rich_condition = st_met > 0;
metal_poor_condition = st_met < 0;
median_st_age = median(st_age(~isnan(st_age)));
median_st_mass = median(st_mass(~isnan(st_mass)));



% 4 variubles tests with dependence on metallicity:
scatter(st_met(metal_rich_condition & st_age>median_st_age), pl_radius(metal_rich_condition & st_age>median_st_age), 'blue', 'filled')
hold on
scatter(st_met(metal_rich_condition & st_age<median_st_age), pl_radius(metal_rich_condition & st_age<median_st_age), 'red', 'filled')
hold on
scatter(st_met(metal_poor_condition & st_age>median_st_age), pl_radius(metal_poor_condition & st_age>median_st_age), 'green', 'filled')
hold on
scatter(st_met(metal_poor_condition & st_age<median_st_age), pl_radius(metal_poor_condition & st_age<median_st_age), 'black', 'filled')
xlabel('Stellar Metallicity');
ylabel('Planet Radius');
legend('Metal rich and old', 'Metal rich and young', 'Metal poor and old', 'Metal poor and young');
title('Planet radius vs Stellar Metallicity AND Stellar Age');





%Some histograms:

%Depedence on stellar metallicity:
ParHistCon(0.188, pl_radius, st_met, 'Planet radius', 'Fraction', 'Stellar metallicity > 0.188', 'Stellar metallicity < 0.188')
ParHistDis(0.11, num_planets, st_met, 'Number of planets', 'Fraction', 'Stellar metallicity > 0.11', 'Stellar metallicity < 0.11')
ParHistCon(-0.079, size_ratios_2pl, SystemMet_2Pl, 'Planet radius ratios', 'Fraction', 'Stellar metallicity > -0.079', 'Stellar metallicity < -0.079')
ParHistCon(0.033, orb_per_ratios_2pl, SystemMet_2Pl, 'Planet orbital period ratios', 'Fraction', 'Stellar metallicity > 0.033', 'Stellar metallicity < 0.033')


%Depedence on stellar mass:
ParHistCon(0.865, pl_radius, st_mass, 'Planet radius', 'Fraction', 'Stellar mass > 0.865', 'Stellar mass < 0.865')
ParHistDis(0.825, num_planets, st_mass, 'Number of planets', 'Fraction', 'Stellar mass > 0.825', 'Stellar mass < 0.825')
ParHistCon(0.865, size_ratios_2pl, SystemMass_2Pl, 'Planet radius ratios', 'Fraction', 'Stellar mass > 0.865', 'Stellar mass < 0.865')
ParHistCon(1.054, orb_per_ratios_2pl, SystemMass_2Pl, 'Planet orbital period ratios', 'Fraction', 'Stellar mass > 1.054', 'Stellar mass < 1.054')



%Depedence on stellar age:
ParHistCon(0.784, pl_radius, st_age, 'Planet radius', 'Fraction', 'Stellar age > 0.784', 'Stellar age < 0.784')
ParHistDis(2.314, num_planets, st_age, 'Number of planets', 'Fraction', 'Stellar age > 2.314', 'Stellar age < 2.314')
ParHistCon(3.846, size_ratios_2pl, SystemAge_2Pl, 'Planet radius ratios', 'Fraction', 'Stellar age > 3.846', 'Stellar age < 3.846')
ParHistCon(1.664, orb_per_ratios_2pl, SystemAge_2Pl, 'Planet orbital period ratios', 'Fraction', 'Stellar age > 1.664', 'Stellar age < 1.664')

set(0,'DefaultFigureVisible','off');


figure
scatter(st_age, pl_orb_per)
set(gca, "YScale","log", "XScale", "log");
ylabel("Planet orbital period");
xlabel("Stellar age");





crit_met = 0.11;
figure
hold on
[n,x]= hist(num_planets(st_age>2.314 & st_met>crit_met),1:1:8);
[n2,~] = hist(num_planets(st_age>2.314 & st_met<crit_met),1:1:8);
[n3,~] = hist(num_planets(st_age<2.314 & st_met>crit_met),1:1:8);
[n4,~] = hist(num_planets(st_age<2.314 & st_met<crit_met),1:1:8);
bar(x, n/sum(n))
bar(x, n2/sum(n2), 'BarWidth',0.5)
bar(x, n3/sum(n3), BarWidth=0.3)
bar(x, n4/sum(n4), BarWidth=0.1)
set(gca, "FontSize", 25)
xlabel("Number of planets", "FontSize", 30)
ylabel("Fraction","FontSize", 30)
legend("High age high met", "High age low met", "Low age high met", "Low age low met", FontSize=25)



function ParHistCon(ThreshVal, PlPar, StPar, Xtitle, Ytitle, LegendWriting1, LegendWriting2)
    figure
    hold on
    [n,x]= hist(PlPar(StPar>ThreshVal),30);
    [n2,~] = hist(PlPar(StPar<ThreshVal),30);
    bar(x,n/sum(n))
    bar(x,n2/sum(n2), 'BarWidth',0.4)
    set(gca, "FontSize", 25)
    xlabel(Xtitle,"FontSize", 30)
    ylabel(Ytitle, "FontSize", 30)
    legend(LegendWriting1, LegendWriting2, FontSize=25)
end

function ParHistDis(ThreshVal, PlPar, StPar, Xtitle, Ytitle, LegendWriting1, LegendWriting2)
    figure
    hold on
    [n,x]= hist(PlPar(StPar>ThreshVal),1:1:8);
    [n2,~] = hist(PlPar(StPar<ThreshVal),1:1:8);
    bar(x, n/sum(n))
    bar(x, n2/sum(n2), 'BarWidth',0.4)
    set(gca, "FontSize", 25)
    xlabel(Xtitle, "FontSize", 30)
    ylabel(Ytitle,"FontSize", 30)
    legend(LegendWriting1, LegendWriting2, FontSize=25)
end

function ksPlot1var(PlParDis, stPar, Xtitle, Ytitle)
    figure
    plot(linspace(min(stPar), max(stPar), 1000), (PlParDis),LineWidth=8)
    set(gca, "YScale", "log","FontSize",25);
    xlabel(Xtitle,"FontSize",30);
    ylabel(Ytitle,FontSize=30);
end

function ksPlot3var(PlPar_MinMid, PlPar_MidMax, PlPar_MinMax, stPar, Xtitle, Ytitle)
    figure
    hold on
    plot(linspace(min(stPar), max(stPar), 1000), PlPar_MinMid,'blue',LineWidth=8)
    plot(linspace(min(stPar), max(stPar), 1000), PlPar_MidMax,'red',LineWidth=8)
    plot(linspace(min(stPar), max(stPar), 1000), PlPar_MinMax,'green',LineWidth=8)
    set(gca, "YScale", "log", "FontSize", 25);
    box on
    xlabel(Xtitle,'FontSize',30);
    ylabel(Ytitle,'FontSize',30);
    legend('Closest planet/Mid planet', 'Mid planet/Furthest planet', 'Closest planet/Furthest planet',fontsize=30);
end