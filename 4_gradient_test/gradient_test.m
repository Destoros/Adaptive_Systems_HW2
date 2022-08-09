% Assignment 2, (2.4)
% testing script for stochastic and deterministic gradient algorithms
% Adaptive System UE WS 19/20
% Thomas Wilding, SPSC
clear 
close all
clc

folder_name = 'Figures';

if ~exist(folder_name, 'dir')
    mkdir(folder_name)%create Figures folder
end

%%%% Reference data creation: (1)
% Nx = 1e4; %number of samples
% 
% h  = [0.2 -0.4 1 0.6]';  %impulse response of unknown system
% Nh = length(h);
% Nc  = Nh;                %adaptive filter order is order of unknown system
% 
% rng(1) %set random number generator
% sigma_x = 0.25; %standard deviation of x
% x = sigma_x*randn(1,Nx);
% d = filter(h,1,x);
% 
% save('./data_rng1.mat','x','sigma_x','d','h','Nh','Nc','Nx')


%%%% Reference data creation: (2)
% rng(2) %set random number generator
% sigma_x = 0.25; %standard deviation of x
% x = sigma_x*randn(1,Nx);
% d = filter(h,1,x);
% 
% sigma_w = 0.1;
% d  = d+sigma_w*randn(1,Nx);
% 
% save('./data_rng2.mat','x','sigma_x','d','h','Nh','Nc','Nx')


%% LMS, noise-free, zero initial vector (with reference)
clear
load('./data_rng1.mat')

c0 = zeros(1,Nc);     %initialize with zero vector

mu = 0.01;
alpha = 0;
[~,~,c_lms]  = lms_algorithm(x,d,Nc,mu,alpha,0); %standard LMS

% save('./ref_lms.mat','c_lms')
ref_lms = load('./ref_lms.mat');

figure(1)
subplot(2,1,1), hold on, grid on, box on
plot(c_lms(1,:),'--k','LineWidth',2.5), plot(ref_lms.c_lms(1,:),'k','LineWidth',0.5)
plot(c_lms(2,:),'--b','LineWidth',2.5), plot(ref_lms.c_lms(2,:),'b','LineWidth',0.5)
plot(c_lms(3,:),'--m','LineWidth',2.5), plot(ref_lms.c_lms(3,:),'m','LineWidth',0.5)
plot(c_lms(4,:),'--r','LineWidth',2.5), plot(ref_lms.c_lms(4,:),'r','LineWidth',0.5)
legend('c1[n], LMS','ref','Location','east')
xlim([1 Nx]), ylim([-1.1 1.1])
title('LMS, noise-free, zero initial vector')
subplot(2,1,2), hold on, grid on, box on
plot(d,'--r'), plot(x,'b')
axis tight
legend('d[n]','x[n]')
title(sprintf('signal plots (sigma_x = %2.2f)',sigma_x),'Interpreter','none')

saveas(gcf,'Figures/LMS', 'epsc')


%% NLMS, noise-free, zero initial vector (with reference)
c0 = zeros(1,Nc);     %initialize with zero vector

mu = 0.01;
alpha = 0;
[~,~,c_nlms]  = lms_algorithm(x,d,Nc,mu,alpha,1); %normalized LMS

% save('./ref_nlms.mat','c_nlms')
ref_nlms = load('./ref_nlms.mat');

figure(2)
subplot(2,1,1), hold on, grid on, box on
plot(c_nlms(1,:),'--k','LineWidth',2.5), plot(ref_nlms.c_nlms(1,:),'k','LineWidth',0.5)
plot(c_nlms(2,:),'--b','LineWidth',2.5), plot(ref_nlms.c_nlms(2,:),'b','LineWidth',0.5)
plot(c_nlms(3,:),'--m','LineWidth',2.5), plot(ref_nlms.c_nlms(3,:),'m','LineWidth',0.5)
plot(c_nlms(4,:),'--r','LineWidth',2.5), plot(ref_nlms.c_nlms(4,:),'r','LineWidth',0.5)
legend('c1[n], NLMS','ref','Location','east')
xlim([1 Nx]), ylim([-1.1 1.1])
title('NLMS, noise-free, zero initial vector')
subplot(2,1,2), hold on, grid on, box on
plot(d,'--r'), plot(x,'b')
axis tight
legend('d[n]','x[n]')
title(sprintf('signal plots (sigma_x = %2.2f)',sigma_x),'Interpreter','none')

saveas(gcf,'Figures/NLMS', 'epsc')

%% NLMS, noise-free, zero initial vector (with reference)
mu = 0.01;
alpha = 0;

Rxx = sigma_x^2*eye(Nc); %x is white noise?
p = Rxx*h;

[~,~,c_gs]  = gd_algorithm(x,d,Nc,mu,Rxx,p); %gradient descent

% save('./ref_gd.mat','c_gd')
ref_gd = load('./ref_gd.mat');

figure(3)
subplot(2,1,1), hold on, grid on, box on
plot(c_gs(1,:),'--k','LineWidth',2.5), plot(ref_gd.c_gs(1,:),'k','LineWidth',0.5)
plot(c_gs(2,:),'--b','LineWidth',2.5), plot(ref_gd.c_gs(2,:),'b','LineWidth',0.5)
plot(c_gs(3,:),'--m','LineWidth',2.5), plot(ref_gd.c_gs(3,:),'m','LineWidth',0.5)
plot(c_gs(4,:),'--r','LineWidth',2.5), plot(ref_gd.c_gs(4,:),'r','LineWidth',0.5)
legend('c1[n], GS','ref','Location','east')
xlim([1 Nx]), ylim([-1.1 1.1])
title('gradient descent, noise-free, zero initial vector')
subplot(2,1,2), hold on, grid on, box on
plot(d,'--r'), plot(x,'b')
axis tight
legend('d[n]','x[n]')
title(sprintf('signal plots (sigma_x = %2.2f)',sigma_x),'Interpreter','none')

saveas(gcf,'Figures/GD', 'epsc')

%% algorithm comparison: noisy, random intial vector (no references)
clear
load('./data_rng2.mat')

Rxx = sigma_x^2*eye(Nc);
p = Rxx*h;

c0 = linspace(-2,2,Nc); %pseudo-random initialization

mu = 0.02;
alpha = 0.1;
[~,~,c_lms]  = lms_algorithm(x,d,Nc,mu,alpha,0,c0);
[~,~,c_nlms] = lms_algorithm(x,d,Nc,mu,alpha,1,c0);
[~,~,c_gs]   = gd_algorithm(x,d,Nc,mu,Rxx,p,c0); %gradient search

figure(4), hold on, grid on, box on
hlms = plot(c_lms.','k','LineWidth',0.75);
hnlms = plot(c_nlms.','b','LineWidth',0.75);
hgs = plot(c_gs.','r');
htrue = scatter(Nx*ones(Nh,1),h,'g','o','filled');
axis tight
legend([hlms(1),hnlms(1),hgs(1),htrue],{'LMS','NLMS','GS','h'})
title('gradient algorithms, noisy convergence comparison')

saveas(gcf,'Figures/comparison', 'epsc')

% PLOT DESCRIPTION
% 
% LMS:
% The additional noise causes the LMS to flucatute a lot it looks like it
% will never truly converge to h. In the problem class we derived, that it
% converges to h (if it converges in the first place) for noise free
% enviroment. BUT if there is some noise, e[n] wont be 0 in the optimum
% point, hence it never truly converges
% 
% NMLS:
% Due to the independet step size to signal energy the NMLS overshoots the
% LMS and GS, especially in the beginning and like the LMS it looks like it
% can never truly converge if there is some noise.
% 
% GS:
% Gradient search creates a smooth looking curve with no fluctuation what so
% ever and converges towards h really nicely. It looks like the noise does
% not effect the algorithm. The downside of GS is the fact that it needs the 
% values of cross corellation p and Autocorrelationmatrix Rxx, which are 
% often not known.
% 
% LMS and NMLS algorithms assume ergodic process, which allow to replace
% expectation operator by an time average of one WSS realization for mean
% and variance


 function saveas(~,~,~)
    disp('Figure not saved')
 end

