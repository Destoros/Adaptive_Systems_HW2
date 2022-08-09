close all
clear all
clc

folder_name = 'Figures';

if ~exist(folder_name, 'dir')
    mkdir(folder_name)%create Figures folder
end

load data.mat

%analytical solution
h = [2, -0.5, 4, -2, -1, 2].';

A = 1;
theta = pi/2;
sigma_u = sqrt(4);

rxx_ref = @(k) A^2/2 * cos(theta*k) + sigma_u^2 * kroneckerDelta(sym(k));
Rxx_ref = double(toeplitz(rxx_ref(0:length(h)-1))); %double convertes the sym data to double precision numbers

p_ref = Rxx_ref * h; %c_MSE = h = Rxx^-1 * p


%estimated solution
N = length(h);

for ii = 1:1% amount of different realizations to be plottet
    
 RNG = round(rand()*size(X,1)); %for each iteration choose a random realization

 for P = [10, 50, 100, 500, 1000, 5000, 8000, 10000]   
     
     
     [rxx, mxx] = cross_correlation(X(RNG,:),X(RNG,:),P,N);
     [rdx, mdx] = cross_correlation(D(RNG,:),X(RNG,:),P,N);

     
     figure
        subplot(2,1,1)
            plot(mxx, rxx)
            grid on
            xlabel('m_{xx}')
            title( ['Autocorrelation r_{xx} for P = ' num2str(P)] )
            hold on
            plot(mxx, rxx_ref(mxx), '--r', 'LineWidth', 2);
            hold off
            legend('estimated r_{xx}', 'analytical r_{xx}')
            
        subplot(2,1,2)
            plot(mdx, rdx) %rdx = p
            grid on
            xlabel('m_{xx}')
            title( ['Cross correlation r_{xx} for P = ' num2str(P)] )
            hold on
            plot(mdx, p_ref, '--r', 'LineWidth', 2);
            hold off
            legend('estimated r_{xx}', 'analytical r_{xx}')
            
            
         if ii == 1   %only save for first realization
            saveas(gcf,['Figures/P=' num2str(P)], 'epsc')
         end
            
 end
end  
 
 %With increasing window size P, the estimated values for rxx and rdx get
 %closer and closer to the analytical solution
 
 function saveas(~,~,~)
    disp('Figure not saved')
 end
 