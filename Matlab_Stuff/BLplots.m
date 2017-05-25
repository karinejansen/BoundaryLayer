load('L4H12H08-n10.txt');
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
plot(N10(:,2),N10(:,5))
hold on
plot(N20(:,2),N20(:,5))
plot(N40(:,2),N40(:,5))
plot(N80(:,2),N80(:,5))
plot(N160(:,2),N160(:,5))
plot(N320(:,2),N320(:,5))
hold off

title('varying number of discretisation (L=4.0 h1=1.2 h2=0.8)')
xlabel('length [m]') 
ylabel('displacement thickness [m]')
%ylim([0 2])
legend('n = 10','n = 20', 'n = 40', 'n = 80', 'n = 160', 'n = 320')

% TIME TAKEN
nUniform = [10 20 40 80 160 320];
TimeUniform = [0.005 0.013 0.032 0.088 0.128 0.262];
nNonUniform = [10 20 40 80];
TimeNonUniform = [0.005 0.009 0.03 0.085];

figure (2)
plot(nUniform, TimeUniform)

title('time taken by CPU')
xlabel('number of discretisations')
ylabel('time taken by CPU')

figure (3)
plot(nUniform, TimeUniform)
hold on
plot(nNonUniform, TimeNonUniform)
hold off
title('time taken by CPU')
xlabel('number of discretisations')
ylabel('time taken by CPU')
legend('uniform distribution', 'non-uniform distribution')

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
xlabel('length [m]') 
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
xlabel('length [m]') 
ylabel('friction stress [m]')
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
xlabel('length [m]') 
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
xlabel('length [m]')
ylabel('friction stress [m]')
ylim([0 0.002])
legend('h2 = 0.77','h2 = 0.8','h2 = 1.0', 'h2 = 1.2', 'h2 = 1.23')

