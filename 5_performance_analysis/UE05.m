close all
clear all
clc


folder_name = 'Figures';

if ~exist(folder_name, 'dir')
    mkdir(folder_name)%create Figures folder
end


load data.mat


% Task a) LMS and Gradient Descent for different mu

N = length(h);

mu_list = [0.0001, 0.001, 0.01, 1];
% mu_list = [];

OPTS = 0; %for standard LMS
alpha = 0;

v = zeros(N,size(X,2),size(X,1));
e = zeros(size(X,1), size(X,2));

for mu = mu_list

    for r = 1:size(X,1) %r...realization index
        x = X(r,:);
        d = D(r,:);

        [~,e(r,:),c] = lms_algorithm(x,d,N,mu,alpha,OPTS); %compute e for each realization 

        v(:,:,r) = c - h; %compute v for each realization 

    end


    %compute the expectation by averaging over different realizations
    e_expect = mean(e,1); %mean along the 1 dimension (rows get added)
    v_expect = mean(v,3); %mean along the third dimension (2D matrixes with row and col get added)


    v_plot = log(abs(v_expect)./abs(v_expect(:,1)));
    e_plot = log(e_expect.^2/e_expect(1).^2);


    n = 0:size(X,2)-1;

    figure
        plot(n,v_plot)
        title(['(A) LMS misalignment \mu =' num2str(mu)])  
        legend('v_1[n]', 'v_2[n]', 'v_3[n]', 'v_4[n]', 'v_5[n]', 'v_6[n]')
        xlabel('n')
        grid on
        
        saveas(gcf,['Figures/' own_strrep(gca)], 'epsc')


    figure
        plot(n,e_plot)
        title(['(A) LMS error \mu =' num2str(mu)])  
        legend('e[n]')
        xlabel('n')
        grid on  
        
        saveas(gcf,['Figures/' own_strrep(gca)], 'epsc')

end %end for mu


    
%----------------------
%Gradient Search 


%get Rxx and p
h = [2, -0.5, 4, -2, -1, 2].';

A = 1;
theta = pi/2;
sigma_u = sqrt(4);

rxx_ref = @(k) A^2/2 * cos(theta*k) + sigma_u^2 * kroneckerDelta(sym(k));
Rxx_ref = double(toeplitz(rxx_ref(0:length(h)-1))); %double convertes the sym data to double precision numbers

p_ref = Rxx_ref * h; %c_MSE = h = Rxx^-1 * p


v = zeros(N,size(X,2),size(X,1));
e = zeros(size(X,1), size(X,2));

for mu = mu_list
    
    for r = 1:size(X,1) %r...realization index
        x = X(r,:);
        d = D(r,:);
        
        [~,e(r,:),c] = gd_algorithm(x,d,N,mu,Rxx_ref,p_ref); 

        v(:,:,r) = c - h;
    end

    %compute the expectation by averaging over different realizations
    e_expect = mean(e,1); %mean along the 1 dimension (rows get added)
    v_expect = mean(v,3); %mean along the third dimension (2D matrixes with row and col get added)


    v_plot = log(abs(v_expect)./abs(v_expect(:,1)));
    e_plot = log(e_expect.^2/e_expect(1).^2);


    n = 0:size(X,2)-1;

    figure
        plot(n,v_plot)
        title(['(A) GD misalignment \mu =' num2str(mu)])  
        legend('v_1[n]', 'v_2[n]', 'v_3[n]', 'v_4[n]', 'v_5[n]', 'v_6[n]')
        xlabel('n')
        grid on
        
        saveas(gcf,['Figures/' own_strrep(gca)], 'epsc')


    figure
        plot(n,e_plot)
        title(['(A) GD error \mu =' num2str(mu)])    
        legend('e[n]')
        xlabel('n')
        grid on   
        
        saveas(gcf,['Figures/' own_strrep(gca)], 'epsc')

end %end for mu
    
% -----------------------------------------------------------------------------
% Task b)

mu = 0.0005;
P_list = [10, 100, 1000];
% P_list = [];

r = 1; %used realization
N = length(h);
    
for P = P_list
    
    [rxx, mxx] = cross_correlation(X(r,:),X(r,:),P,N);
    [rdx, mdx] = cross_correlation(D(r,:),X(r,:),P,N);
     
    Rxx = toeplitz(rxx);
    p = rdx;
    
    
    for r = 1:size(X,1) %r...realization index
        x = X(r,:);
        d = D(r,:);
       
        %use estimated statistics
        [~,e(r,:),c] = gd_algorithm(x,d,N,mu,Rxx,p); 
        v(:,:,r) = c - h;
        
        %use analytical statistics
        [~,e_an(r,:),c_an] = gd_algorithm(x,d,N,mu,Rxx_ref,p_ref); 
        v_an(:,:,r) = c_an - h;
    end

    %compute the expectation by averaging over different realizations
    e_expect = mean(e,1); %mean along the 1 dimension (rows get added)
    v_expect = mean(v,3); %mean along the third dimension (2D matrixes with row and col get added)
    
    e_expect_an = mean(e_an,1); %mean along the 1 dimension (rows get added)
    v_expect_an = mean(v_an,3); %mean along the third dimension (2D matrixes with row and col get added)


    v_plot = log(abs(v_expect)./abs(v_expect(:,1)));
    e_plot = log(e_expect.^2/e_expect(1).^2);
    
    v_plot_an = log(abs(v_expect_an)./abs(v_expect_an(:,1)));
    e_plot_an = log(e_expect_an.^2/e_expect_an(1).^2);


    n = 0:size(X,2)-1;

    figure
        plot(n,v_plot)
        hold on
        plot(n,v_plot_an)
        title(['(B) GD misalignment P=' num2str(P)])  
        legend('v_1[n]', 'v_2[n]', 'v_3[n]', 'v_4[n]', 'v_5[n]', 'v_6[n]', ...
               'v_1[n] analytical', 'v_2[n] analytical', 'v_3[n] analytical', 'v_4[n] analytical', 'v_5[n] analytical', 'v_6[n] analytical')
        xlabel('n')
        grid on
        
        saveas(gcf,['Figures/' own_strrep(gca)], 'epsc')


    figure
        plot(n,e_plot)
        hold on
        plot(n,e_plot_an)
        title(['(B) GD error P=' num2str(P)])    
        legend('e[n] estimate', 'e[n] analytical')
        xlabel('n')
        grid on 
        
        saveas(gcf,['Figures/' own_strrep(gca)], 'epsc')

     

end


%Seen from plots: with increasing window size, the gradient descend
%convergers better and better, the misalignment gets lower


% -----------------------------------------------------------------------------
% Task c)


mu = 0.001;
alpha = 0;
OPTS = 0;

N_list = [10, 100];
% N_list = [];

for N = N_list
    
    %reset
    v_LMS = [];
    v_GD = [];
    
    %random initialization for c0
    c0 = rand(N,1); 
    
    A = 1;
    theta = pi/2;
    sigma_u = sqrt(4);

    rxx_ref = @(k) A^2/2 * cos(theta*k) + sigma_u^2 * kroneckerDelta(sym(k));
    Rxx_ref = double(toeplitz(rxx_ref(0:N-1))); %double convertes the sym data to double precision numbers
    
    h_calc = [h; zeros(N-length(h),1)];
    
    p_ref = Rxx_ref * h_calc; %c_MSE = h = Rxx^-1 * p
    
    for r = 1:size(X,1) %r...realization index
        x = X(r,:);
        d = D(r,:);

        [~,e_LMS(r,:),c_LMS] = lms_algorithm(x,d,N,mu,alpha,OPTS,c0); %compute e for each realization 
        [~,e_GD(r,:),c_GD] = gd_algorithm(x,d,N,mu,Rxx_ref,p_ref); 

        v_LMS(:,:,r) = c_LMS - h_calc; %compute v for each realization 
        v_GD(:,:,r) = c_GD - h_calc; %compute v for each realization 

    end
    
    
    %compute the expectation by averaging over different realizations
    e_expect_LMS = mean(e_LMS,1); %mean along the 1 dimension (rows get added)
    v_expect_LMS = mean(v_LMS,3); %mean along the third dimension (2D matrixes with row and col get added)
    
    e_expect_GD = mean(e_GD,1); %mean along the 1 dimension (rows get added)
    v_expect_GD = mean(v_GD,3); %mean along the third dimension (2D matrixes with row and col get added)


    v_plot_LMS = log(abs(v_expect_LMS)./abs(v_expect_LMS(:,1)));
    e_plot_LMS = log(e_expect_LMS.^2/e_expect_LMS(1).^2);
    
    v_plot_GD = log(abs(v_expect_GD)./abs(v_expect_GD(:,1)));
    e_plot_GD = log(e_expect_GD.^2/e_expect_GD(1).^2);


    n = 0:size(X,2)-1;

    figure
        plot(n,v_plot_LMS)
        hold on
        plot(n,v_plot_GD)
        title(['(C) LMS and GD misalignment N=' num2str(N)])  
        %legend not useful, too many coefficients
        xlabel('n')
        grid on
        
        saveas(gcf,['Figures/' own_strrep(gca)], 'epsc')


    figure
        plot(n,e_plot_LMS)
        hold on
        plot(n,e_plot_GD)
        title(['(C) LMS and GD error N=' num2str(N)])    
        legend('e[n] LMS', 'e[n] GD')
        xlabel('n')
        grid on 
        
        saveas(gcf,['Figures/' own_strrep(gca)], 'epsc')
    
    
    
end

%seen from plot: increasing order number has an influence (its overmodelled
%in this case)
%The LMS gets way worse convergence wise, the GD seems less impacted by
%this





%create a placeholder function to overwrite the saveas function
function saveas(~, ~, ~)
    disp('Figure not saved')
end




