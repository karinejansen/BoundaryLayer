% load('L4H12H08-n10.txt');
load('L4H12H08-n20.txt');
load('L4H12H08-n40.txt');
load('L4H12H08-n80.txt');
load('L4H12H08-n160.txt');
load('L4H12H08-n320.txt');

add = [0 0 0.5 0 0 100];

N10 = dlmread('L4H12H08-n10.txt');
N10 = [add;N10];
N20 = dlmread('L4H12H08-n20.txt');
N20 = [add;N20];
N40 = dlmread('L4H12H08-n40.txt');
N40 = [add;N40];
N80 = dlmread('L4H12H08-n80.txt');
N80 = [add;N80];
N160 = dlmread('L4H12H08-n160.txt');
N160 = [add;N160];
N320 = dlmread('L4H12H08-n320.txt');
N320 = [add;N320];


figure (1)
plot(N10(:,2),N10(:,5), '-s', 'MarkerSize',1)
hold on
plot(N20(:,2),N20(:,5), '-s', 'MarkerSize',1)
plot(N40(:,2),N40(:,5), '-s', 'MarkerSize',1)
plot(N80(:,2),N80(:,5), '-s', 'MarkerSize',1)
plot(N160(:,2),N160(:,5), '-s', 'MarkerSize',1)
plot(N320(:,2),N320(:,5),'-s', 'MarkerSize',1)
hold off

title('varying number of discretisation (L=4.0 h1=1.2 h2=0.8)')
xlabel('x-location [m]') 
ylabel('displacement thickness [m]')
%ylim([0 2])
legend('n = 10','n = 20', 'n = 40', 'n = 80', 'n = 160', 'n = 320')

% TIME TAKEN
nUniform = [10 20 40 80 81 160 320];
TimeUniform = [0 0 0.016 0.085 0.069 0.116 0.234];
nNonUniform = [10 20 40 80];
TimeNonUniform = [0 0 0.015 0.085];
labels = cellstr( num2str(nUniform') );

figure (2)
plot(nUniform, TimeUniform, 'o')
title('time taken by CPU')
xlabel('number of discretisations')
ylabel('time taken by CPU')
text(nUniform, TimeUniform, labels,'VerticalAlignment','middle', 'HorizontalAlignment','left')

figure (3)
plot(nUniform, TimeUniform, 'o')
hold on
plot(nNonUniform, TimeNonUniform, 'd')
hold off
title('time taken by CPU')
xlabel('number of discretisations')
ylabel('time taken by CPU')
legend('uniform distribution', 'non-uniform distribution')
text(nUniform, TimeUniform, labels,'VerticalAlignment','middle', 'HorizontalAlignment','left')

% SPEED VARIATION
load('L4H12H08-n80-U0.1.txt');
load('L4H12H08-n80-U1.0.txt');
load('L4H12H08-n80-U10.0.txt');
load('L4H12H08-n80-U100.0.txt');
U01 = dlmread('L4H12H08-n80-U0.1.txt');
U01 = [add;U01];
U1 = dlmread('L4H12H08-n80-U1.0.txt');
U1 = [add;U1];
U10 = dlmread('L4H12H08-n80-U10.0.txt');
U10 = [add;U10];
U100 = dlmread('L4H12H08-n80-U100.0.txt');
U100 = [add;U100];
figure (4) % speed variation displacement thickness
plot(U01(:,2),U01(:,5))
hold on
plot(N80(:,2),N80(:,5))
plot(U1(:,2),U1(:,5))
plot(U10(:,2),U10(:,5))
plot(U100(:,2),U100(:,5))
hold off
title('varying speed (L=4.0 h1=1.2 h2=0.8)')
xlabel('x-location [m]') 
ylabel('displacement thickness [m]')
%ylim([0 2])
legend('U = 0.1','U = 0.5','U = 1.0', 'U = 10.0', 'U = 100.0')

figure (5) % speed variation - friction stress
plot(U01(:,2),U01(:,6))
hold on
plot(N80(:,2),N80(:,6))
plot(U1(:,2),U1(:,6))
plot(U10(:,2),U10(:,6))
plot(U100(:,2),U100(:,6))
hold off
title('varying speed (L=4.0 h1=1.2 h2=0.8)')
xlabel('x-location [m]') 
ylabel('friction stress [N/m^2]')
ylim([0 0.005])
legend('U = 0.1','U = 0.5','U = 1.0', 'U = 10.0', 'U = 100.0')

% H2 VARIATION
load('L4H12H077-n80.txt');
load('L4H12H10-n80.txt');
load('L4H12H12-n80.txt');
load('L4H12H123-n80.txt');
H077 = dlmread('L4H12H077-n80.txt');
H077 = [add;H077];
H10 = dlmread('L4H12H10-n80.txt');
H10 = [add;H10];
H12 = dlmread('L4H12H12-n80.txt');
H12 = [add;H12];
H123 = dlmread('L4H12H123-n80.txt');
H123 = [add;H123];
figure (6) % h2 variation vs. displacement thickness
plot(H077(:,2),H077(:,5))
hold on
plot(N80(:,2),N80(:,5))
plot(H10(:,2),H10(:,5))
plot(H12(:,2),H12(:,5))
plot(H123(:,2),H123(:,5))
hold off
title('varying h2 (L=4.0 h1=1.2)')
xlabel('x-location [m]') 
ylabel('displacement thickness [m]')
%ylim([0 2])
legend('h2 = 0.77','h2 = 0.8','h2 = 1.0', 'h2 = 1.2', 'h2 = 1.23')

figure (7) % speed variation - friction stress
plot(H077(:,2),H077(:,6))
hold on
plot(N80(:,2),N80(:,6))
plot(H10(:,2),H10(:,6))
plot(H12(:,2),H12(:,6))
plot(H123(:,2),H123(:,6))
hold off
title('varying h2 (L=4.0 h1=1.2)')
xlabel('x-location [m]')
ylabel('friction stress [N/m^2]')
ylim([0 0.002])
legend('h2 = 0.77','h2 = 0.8','h2 = 1.0', 'h2 = 1.2', 'h2 = 1.23')

% LENGTH VARIATION
load('L05H12H08-n80.txt');
load('L2H12H08-n80.txt');
load('L8H12H08-n80.txt');
L05 = dlmread('L05H12H08-n80.txt');
L05 = [add;L05];
L2 = dlmread('L2H12H08-n80.txt');
L2 = [add;L2];
L8 = dlmread('L8H12H08-n80.txt');
L8 = [add;L8];

figure (8) % channel length variation vs. displacement thickness
plot((L05(:,2)/0.5),L05(:,5))
hold on
plot((L2(:,2)/2),L2(:,5))
plot((N80(:,2)/4),N80(:,5))
plot((L8(:,2)/8),L8(:,5))
hold off
title('varying channel length L (h1=1.2 h2=0.8)')
xlabel('x/L [-]') 
ylabel('displacement thickness [m]')
%ylim([0 2])
legend('L = 0.5','L = 2.0','L = 4.0', 'L = 8.0')

figure (9) % channel length variation vs. friction stress
plot((L05(:,2)/0.5),L05(:,6))
hold on
plot((L2(:,2)/2),L2(:,6))
plot((N80(:,2)/4),N80(:,6))
plot((L8(:,2)/8),L8(:,6))
hold off
title('varying channel length L (h1=1.2 h2=0.8)')
xlabel('x/L [-]') 
ylabel('friction stress [N/m^2]')
ylim([0 0.002])
legend('L = 0.5','L = 2.0','L = 4.0', 'L = 8.0')

% NON-UNIFORM DISTRIBUTION
in = [0:1:80];
x = 2*cos(in*pi/80 +pi)+2;

figure(10)
plot(in,x)
title('non-uniform distribution')
ylabel('x-location [m]') 
xlabel('n-x [-]')

% PARALLEL CHANNEL - SPEED VARIATION
load('L4H12H12-n80-U5.0.txt');
load('L4H12H12-n80-U10.0.txt');

PU5 = dlmread('L4H12H12-n80-U5.0.txt');
PU5 = [add;PU5];
PU10 = dlmread('L4H12H12-n80-U10.0.txt');
PU10 = [add;PU10];

figure (11) % speed variation displacement thickness
plot(H12(:,2),H12(:,5))
hold on
plot(PU5(:,2),PU5(:,5))
plot(PU10(:,2),PU10(:,5))
hold off
title('PARALLEL - varying speed (L=4.0 h1=1.2 h2=1.2)')
xlabel('x-location [m]') 
ylabel('displacement thickness [m]')
%ylim([0 2])
legend('U = 0.5','U = 5.0','U = 10.0')

figure (12) % speed variation - friction stress
plot(H12(:,2),H12(:,6))
hold on
plot(PU5(:,2),PU5(:,6))
plot(PU10(:,2),PU10(:,6))
hold off
title('PARALLEL - varying speed (L=4.0 h1=1.2 h2=1.2)')
xlabel('x-location [m]') 
ylabel('friction stress [N/m^2]')
ylim([0 0.1])
legend('U = 0.5','U = 5.0','U = 10.0')

% NON-UNIFORM DISTRIBUTION
NU10 = dlmread('nonL4H12H08-n10.txt');
NU10 = [add;NU10];
NU20 = dlmread('nonL4H12H08-n20.txt');
NU20 = [add;NU20];
NU40 = dlmread('nonL4H12H08-n40.txt');
NU40 = [add;NU40];
NU80 = dlmread('nonL4H12H08-n80.txt');
NU80 = [add;NU80];
dummyX = [0 1 2 3 4];
dummyY = [0.0095 0.0095 0.0095 0.0095 0.0095];

figure (13)
plot(dummyX,dummyY, 'w')%dummy
hold on
plot(N10(:,2),N10(:,5), '-s','color',[0 0.4470 0.7410], 'MarkerSize',1)
plot(N20(:,2),N20(:,5), '-s','color',[0.8500 0.3250 0.0980], 'MarkerSize',1)
plot(N40(:,2),N40(:,5), '-s', 'color',[0.9290 0.6940 0.1250], 'MarkerSize',1)
plot(N80(:,2),N80(:,5), '-s', 'color',[0.4940 0.1840 0.5560], 'MarkerSize',1)
plot(N160(:,2),N160(:,5), '-s','color',[0.4660 0.6740 0.1880],  'MarkerSize',1)
plot(dummyX,dummyY, 'w')% dummy
plot(NU10(:,2),NU10(:,5), '-.o','color',[0 0.4470 0.7410],  'MarkerSize',1)
plot(NU20(:,2),NU20(:,5), '-.o','color',[0.8500 0.3250 0.0980], 'MarkerSize',1)
plot(NU40(:,2),NU40(:,5), '-.o', 'color',[0.9290 0.6940 0.1250],'MarkerSize',1)
plot(NU80(:,2),NU80(:,5), '-.o', 'color',[0.4940 0.1840 0.5560],'MarkerSize',1)
hold off

title('uniform vs. non-uniform distribution (L=4.0 h1=1.2 h2=0.8)')
xlabel('x-location [m]') 
ylabel('displacement thickness [m]')
%ylim([0 2])
legend('uniform:', 'n = 10','n = 20', 'n = 40', 'n = 80', 'n = 160', 'non-uniform:','n = 10','n = 20', 'n = 40', 'n = 80')

figure (14)
plot(dummyX,dummyY, 'w')%dummy
hold on
%plot(N10(:,2),N10(:,5), '-s','color',[0 0.4470 0.7410], 'MarkerSize',1)
%plot(N20(:,2),N20(:,5), '-s','color',[0.8500 0.3250 0.0980], 'MarkerSize',1)
%plot(N40(:,2),N40(:,5), '-s', 'color',[0.9290 0.6940 0.1250], 'MarkerSize',1)
plot(N80(:,2),N80(:,5), '-s', 'color',[0.4940 0.1840 0.5560], 'MarkerSize',1)
%plot(N160(:,2),N160(:,5), '-s','color',[0.4660 0.6740 0.1880],  'MarkerSize',1)
plot(dummyX,dummyY, 'w')% dummy
%plot(NU10(:,2),NU10(:,5), '-.o','color',[0 0.4470 0.7410],  'MarkerSize',1)
%plot(NU20(:,2),NU20(:,5), '-.o','color',[0.8500 0.3250 0.0980], 'MarkerSize',1)
plot(NU40(:,2),NU40(:,5), '-.o', 'color',[0.9290 0.6940 0.1250],'MarkerSize',2)
%plot(NU80(:,2),NU80(:,5), '-.o', 'color',[0.4940 0.1840 0.5560],'MarkerSize',1)
hold off

title('uniform vs. non-uniform distribution (L=4.0 h1=1.2 h2=0.8)')
xlabel('x-location [m]') 
ylabel('displacement thickness [m]')
%ylim([0 2])
legend('uniform:', 'n = 80', 'non-uniform:', 'n = 40')

% GEOMETRIC    NON-UNIFORM DISTRIBUTION
G10 = dlmread('geoL4H12H08x0.05-n10.txt');
G10 = [add;G10];
G20 = dlmread('geoL4H12H08x0.05-n20.txt');
G20 = [add;G20];
G40 = dlmread('geoL4H12H08x0.05-n40.txt');
G40 = [add;G40];
G80 = dlmread('geoL4H12H08x0.05-n80.txt');
G80 = [add;G80];

figure (15)
plot(dummyX,dummyY, 'w')%dummy
hold on
plot(N10(:,2),N10(:,5), '-s','color',[0 0.4470 0.7410], 'MarkerSize',1)
plot(N20(:,2),N20(:,5), '-s','color',[0.8500 0.3250 0.0980], 'MarkerSize',1)
plot(N40(:,2),N40(:,5), '-s', 'color',[0.9290 0.6940 0.1250], 'MarkerSize',1)
plot(N80(:,2),N80(:,5), '-s', 'color',[0.4940 0.1840 0.5560], 'MarkerSize',1)
plot(N160(:,2),N160(:,5), '-s','color',[0.4660 0.6740 0.1880],  'MarkerSize',1)
plot(dummyX,dummyY, 'w')% dummy
plot(G10(:,2),G10(:,5), '-.o','color',[0 0.4470 0.7410],  'MarkerSize',1)
plot(G20(:,2),G20(:,5), '-.o','color',[0.8500 0.3250 0.0980], 'MarkerSize',1)
plot(G40(:,2),G40(:,5), '-.o', 'color',[0.9290 0.6940 0.1250],'MarkerSize',1)
plot(G80(:,2),G80(:,5), '-.o', 'color',[0.4940 0.1840 0.5560],'MarkerSize',1)
hold off

title('uniform vs. geometric distribution (L=4.0 h1=1.2 h2=0.8)')
xlabel('x-location [m]') 
ylabel('displacement thickness [m]')
%ylim([0 2])
legend('uniform:', 'n = 10','n = 20', 'n = 40', 'n = 80', 'n = 160', 'geometric:','n = 10','n = 20', 'n = 40', 'n = 80')

figure (16)
plot(dummyX,dummyY, 'w')%dummy
hold on
%plot(N10(:,2),N10(:,5), '-s','color',[0 0.4470 0.7410], 'MarkerSize',1)
%plot(N20(:,2),N20(:,5), '-s','color',[0.8500 0.3250 0.0980], 'MarkerSize',1)
%plot(N40(:,2),N40(:,5), '-s', 'color',[0.9290 0.6940 0.1250], 'MarkerSize',1)
plot(N80(:,2),N80(:,5), '-s', 'color',[0.4940 0.1840 0.5560], 'MarkerSize',2)
%plot(N160(:,2),N160(:,5), '-s','color',[0.4660 0.6740 0.1880],  'MarkerSize',1)
plot(dummyX,dummyY, 'w')% dummy
%plot(G10(:,2),G10(:,5), '-.o','color',[0 0.4470 0.7410],  'MarkerSize',1)
plot(G20(:,2),G20(:,5), '-.o','color',[0.8500 0.3250 0.0980], 'MarkerSize',2)
plot(G40(:,2),G40(:,5), '-.o', 'color',[0.9290 0.6940 0.1250],'MarkerSize',2)
%plot(G80(:,2),G80(:,5), '-.o', 'color',[0.4940 0.1840 0.5560],'MarkerSize',1)
plot(dummyX,dummyY, 'w')% dummy
plot(NU40(:,2),NU40(:,5), '--o', 'color',[0 0.4470 0.7410],'MarkerSize',2)
hold off

title('uniform, cosine and geometric distribution (L=4.0 h1=1.2 h2=0.8)')
xlabel('x-location [m]') 
ylabel('displacement thickness [m]')
%ylim([0 2])
legend('uniform:', 'n = 80', 'geometric:','n = 20', 'n = 40', 'cosine:', 'n = 40')

% GEOMETRIC - TIME TAKEN
nGeometric = [10 20 40 80 160];
TimeGeometric = [0 0 0.068 0.1 0.118];

figure (17)
plot(nUniform, TimeUniform, 'o')
hold on
plot(nNonUniform, TimeNonUniform, 'd')
plot(nGeometric, TimeGeometric, 's')
hold off
title('time taken by CPU')
xlabel('number of discretisations')
ylabel('time taken by CPU')
legend('uniform', 'cosine', 'geometric')
text(nUniform, TimeUniform, labels,'VerticalAlignment','middle', 'HorizontalAlignment','left')

% ITERATIONS
IterUniform = [6.7 5.6 4.2 3.625 3.654321 3.5625 3.703125];
IterNonUniform = [6.8 5.4 4.225 4.8125];
IterGeometric = [5.6 5.05 4.55 4.1 3.725];

figure (18)
plot(nUniform, IterUniform, 'o')
hold on
plot(nNonUniform, IterNonUniform, 'd')
plot(nGeometric, IterGeometric, 's')
hold off
title('average number of iterations')
xlabel('number of discretisations')
ylabel('# of iterations')
legend('uniform', 'cosine', 'geometric')
text(nUniform, TimeUniform, labels,'VerticalAlignment','middle', 'HorizontalAlignment','left')
ylim([3.5 7])
xlim([0 175])

TotalUniform = [67 112 168 290 296 570 1185];
TotalNonUniform = [68 108 169 385 ];
TotalGeometric = [56 101 182 328 596];

figure (19)
plot(nUniform, TotalUniform, 'o')
hold on
plot(nNonUniform, TotalNonUniform, 'd')
plot(nGeometric, TotalGeometric, 's')
hold off
title('total number of iterations')
xlabel('number of discretisations')
ylabel('# of iterations')
legend('uniform', 'cosine', 'geometric')
%text(nUniform, TimeUniform, labels,'VerticalAlignment','middle', 'HorizontalAlignment','left')
%ylim([3.5 7])
xlim([0 175])