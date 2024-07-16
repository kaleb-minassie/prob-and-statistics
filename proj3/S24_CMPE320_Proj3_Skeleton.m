%S21 CMPE320 Project3 Solns
% Modified for S21, EFCL 4/15/2021
% Modified for S24, EFCL 3/23/2024

close all;
clear all;  

Ntrials = 100000; % recommendation, make your own choice
kplot=0;  % this will keep track of MATLAB figure numbers

% This next line is what MATLAB calls an anonymous function, which is a
% function that we use only in one script. the function call is
% (x,mu,sigma^2), and the returned value is a number

fgauss = @(x,mu,sig2) exp(-0.5*((x-mu).^2/sig2))/sqrt(2*pi*sig2);

kplot=1;

%Problem 2.1  Uniformly distributed
Nsum = [2,6,12]; % checked  with S21 assignment
for k=1:length(Nsum)
    
    xd = rand(Nsum(k),Ntrials);  % generate [Nsum(k) by Ntrials] array of random values
    xs = sum(xd); % "sum" adds down the MATLAB columns, thereby giving us the sum of Nsum(k) values
    xmin = 0;
    xmax = Nsum(k); % largest value of the sum
    mu = Nsum(k)*0.5; % 0.5 is the E[X]
    sig2 = Nsum(k)*(1/12); % 1/12 is var[X]
    m = mean(xs); % sample mean
    S = var(xs); % sample variance
    disp(['For ',int2str(Ntrials),' independent trials of the sum of ',int2str(Nsum(k)),' iid rv from U(0,1)']);
    disp(['   the theoretical mean is ',num2str(mu),' and the sample mean is ',num2str(m)]);
    disp(['   the theoretical variance is ',num2str(sig2),' and the sample variance is ',num2str(S)]);
    
    dx=0.01; % fine grain dx for plotting
    x = 0:dx:Nsum(k)+1; % fine grain for fY(y)
    binwidth =  max(xs) / 50; % choose your own
    edges = -binwidth/2:binwidth:Nsum(k)+binwidth/2; % set the bin edges
    
    % create new figure...
    figure(kplot);
    kplot=kplot+1;
    %...and then a new scaled histogram using the values of xs
    H = histogram(xs,'BinEdges',edges,'Normalization','pdf');
    % And unpack the data using unpackHistogram from Project 2
    [Values,Nbins,binCenters]=unpackHistogram(H);
    % And plot the Gaussian pdf on top of the histogram with labels and
    % grids and the other elements of a professional plot
    hold on
       plot(x,fgauss(x,mu,sig2),'r','LineWidth',2);
    hold off
    grid on;
    xlabel('x');
    ylabel('f_X(x)');
    legend('Histogram','CLT Gaussian');
    title(['Project 3, Section 2.1, sum of ',int2str(Nsum(k)),' U(0,1)']);
end
% there are length(Nsum) plots to this point.


%Problem 2.2  Uniformly distributed discrete
Nsum = [2,30,50]; % checked with S24 assignment
Nsides = 12; %   checked with S24 assignment
mu1 = (Nsides + 1) / 2;  % mean of Nsided die
sig2_1 = ((Nsides^2) - 1) / 12; % variance of Nsided die;

% Do the experiment again for a large number of trials (Ntrials) and the
% specified number of terms in the sum (Nsum)

for k=1:length(Nsum)
    
    xd = randi(Nsides,Nsum(k),Ntrials);  % generate [Nsum(k) by Ntrials] array of random values
    xs = sum(xd); % "sum" adds down the MATLAB columns, thereby giving us the sum of Nsum(k) values
    xmin = Nsum(k);
    xmax = Nsum(k) * Nsides; % largest value of the sum
    mu = Nsum(k) * mu1; % set the theoretical mean of the sum
    sig2 = Nsum(k) * sig2_1; % set the theoretical variance of the sum
    m = mean(xs); % sample mean
    S = var(xs); % sample variance
    disp(['For ',int2str(Ntrials),' independent trials of the sum of ',int2str(Nsum(k)),' iid rv from pmf U(1,8)']);
    disp(['   the theoretical mean is ',num2str(mu),' and the sample mean is ',num2str(m)]);
    disp(['   the theoretical variance is ',num2str(sig2),' and the sample variance is ',num2str(S)]);
    
    dx=0.1; % fine grain dx for plotting
    x = 0:dx:Nsum(k)*Nsides+1; % fine grain for fY(y)
    binwidth = 1; % because it's a discrete sum of integers
    edges = -binwidth/2:1:Nsum(k)*Nsides+binwidth/2;
    
    % create new figure...
    figure(kplot);
    kplot=kplot+1;
    %...and then a new scaled histogram using the values of xs
    H = histogram(xs, 'BinEdges', edges, 'Normalization', 'pdf'); % use the correct parameters
    % And unpack the data using unpackHistogram from Project 2
    [Values,Nbins,binCenters]=unpackHistogram(H);
    % And plot the Gaussian pdf on top of the histogram with labels and
    % grids and the other elements of a professional plot
    hold on
       plot(x,fgauss(x,mu,sig2),'r','LineWidth',2);
    hold off
    grid on;
    xlim([min(xs)-5, max(xs)+5]); % set a good x-scale factor
    xlabel('x');
    ylabel('f_X(x)');
    legend('Histogram','CLT Gaussian');
    title(['Project 3, Section 2.2, sum of ',int2str(Nsum(k)),' ',int2str(Nsides),'-sided dice']);
end

%Problem 2.3  Exponentially distributed 
Nsum = [5,40,200]; % checked with S22 assignment 

for k=1:length(Nsum)
    
    lambda = 0.7;  % checked  with S24 assignment
    xd = randx(Nsum(k),Ntrials,lambda); % randx is provide with this assignment.
    %Note: DO NOT use rand2x from Project 1, use randx.

    xs = sum(xd); % compute the sum
    
    % set the appropriate min and max for the x-axis
    xmin = min(xs);
    xmax = max(xs);
    mu = Nsum(k) / lambda; %  the E[sum]
    sig2 = Nsum(k) / (lambda^2); %  var[sum]
    m = mean(xs); % sample mean
    S = var(xs); % sample variance
    disp(['For ',int2str(Ntrials),' independent trials of the sum of ',int2str(Nsum(k)),' iid rv from ',num2str(lambda),'* exp(-',num2str(lambda),'x)']);
    disp(['   the theoretical mean is ',num2str(mu),' and the sample mean is ',num2str(m)]);
    disp(['   the theoretical variance is ',num2str(sig2),' and the sample variance is ',num2str(S)]);
    
    dx = 0.1;
    x=xmin:dx:xmax;
    
    figure(kplot)  %create new plot
    kplot=kplot+1;
   
   % do the histogram
   histogram(xs, 'Normalization', 'pdf');
   % plot the CLT Gaussian approximation on top of it.
    hold on;
    fgauss = @(x, mu, sig2) exp(-0.5*((x - mu).^2/sig2)) / sqrt(2 * pi * sig2);
    plot(x, fgauss(x, mu, sig2), 'r', 'LineWidth', 2);
    hold off;
    grid on;
    xlabel('x');
    ylabel('f_X(x)');
    title(['Project 3, Section 2.3, sum of ', int2str(Nsum(k)), ' iid exponential RVs ']);
    legend('Histogram', 'CLT Gaussian');
end

% Problem 2.4 Sum of iid Bernoulli trials
Nsum = [4, 20, 100];
p = 0.6; % probability of success in each Bernoulli trial

for k = 1:length(Nsum)
    n = Nsum(k); 
    mu = n * p; % Theoretical mean of the sum of n Bernoulli trials
    sig2 = n * p * (1 - p); % Theoretical variance of the sum of n Bernoulli trials
    
    % generate trials
    xd=(rand(n, Ntrials) < p);
    xs = sum(xd);
    m = mean(xs); % sample mean
    S = var(xs); % sample variance

     disp(['For ', int2str(Ntrials), ' trials of the sum of ', int2str(n), ' iid Bernoulli(p=', num2str(p), ') variables:']);
    disp(['   the theoretical mean is ',num2str(mu),' and the sample mean is ',num2str(m)]);
    disp(['   the theoretical variance is ',num2str(sig2),' and the sample variance is ',num2str(S)]);

    % plotting the histogram of the trials
    figure(kplot);
    kplot = kplot + 1;
    subplot(2, 1, 1); % First subplot for the histogram and theoretical PMF
    edges = -0.5:1:(n+0.5); % Histogram bin edges
    histogram(xs, edges, 'Normalization', 'probability');
    title(sprintf('Project 3, Section 2.4, PMF for N = %d ', n));
    xlabel('x');
    ylabel('f_X(x)');
    hold on;
    
    % plot the theoretical PMF using the binomial formula
    y = 0:n;
    pmf = exp(gammaln(n + 1) - gammaln(y + 1) - gammaln(n - y + 1) + y * log(p) + (n - y) * log(1 - p));
    stem(y, pmf, 'ro-');
    
    % label the plot
    legend('Histogram', 'Theoretical PMF');
    hold off;

    subplot(2, 1, 2); % Second subplot for the CLT Gaussian approximation
    histogram(xs, edges, 'Normalization', 'probability');
    hold on;

    % plot the Gaussian approximation
    xi = min(xs):0.1:max(xs);
    normal = exp(-0.5 * ((xi - mu).^2 / sig2)) / sqrt(2 * pi * sig2);
    plot(xi, normal, 'r-');
    
     % label the plot
    title(sprintf('Project 3, Section 2.4, CLT Gaussian Approximation for N = %d', n));
    xlabel('x');
    ylabel('f_X(x)');
    legend('Histogram', 'CLT Gaussian');
    hold off;
end