% S23 CMPE320 Project 1 Skeleton
%    Histograms, PDFs and PMFs
%    EFCL   1/14/2022
%    Revised 2/4/2023
%    Revised 1/31/2024

%These are standard calls that you should make at the beginning of MATLAB
%scripts, but no at the beginning of MATLAB  functions.  They "clean up"
%old work from the workspace and plots.

close all  % remove all existing figures (very useful to avoid confusion)
clear      % remove all existing variables, but not existing break points in scripts or functions

% 3.1 PMF for a single die

%bin_edges contain the upper and lower edges for the histogram that we will
%create. Therefore, it always has one more column than the number of bins,
%with the upper edge of one bin corresponding to the lower edge of the
%next.

Nsides = 8;  % use an 8-sided die (F24)
binwidth = 1;% the width of the bins for the histogram
bin_edges = binwidth/2:binwidth:Nsides+binwidth/2; %  start half bin less than 1 and finish half bin more than Nsides
% Question for thought, does a smaller bin size (say binwidth = 0.25) change the result?
Ntrials = [240, 2400 24000 240000]; % per the assignment
disp(' ');
disp(['Section 3.1 PMF of a single fair, ',int2str(Nsides),'-sided die']);

for ktrials = 1:length(Ntrials)  %loop on the number of trials
    
    %  set a new figure (or figures) for this number of trials
    figure;

    rolls = randi(Nsides, Ntrials(ktrials), 1); % do the random trials;  lookup MATLAB function randi
    
    sample_mean = mean(rolls); % compute the sample mean using MATLAB function mean( )
    sample_var = var(rolls); % compute the sample variance using MATLAB function var( )
    
    %Each call to histogram will produce a plot.  You might consider
    %putting these in separate figures or subplots
    
    % These calls are correct, but I won't do it again.  You're CMPEs; you
    % can figure it out from here
    
    %Raw histogram (just a count of values in each bin)
    histogram(rolls, 'BinEdges', bin_edges);
    title(sprintf('Section 3.1: Raw Histogram for N=%d', Ntrials(ktrials)));
    xlabel('Number of spots');
    ylabel('Count');    
    %Normalized histogram, note calling parameters (for discrete random
    %variables, 'normalization','probability' gives what you want
    %New call creates new plot overwriting the plot from hist_raw.
    %If you want to save the raw histogram plot, add a call to figure
    %before the next statement

    histogram(rolls, 'BinEdges', bin_edges, 'Normalization', 'probability');
    title(sprintf('Section 3.1: Normalized Histogram (PMF) for N=%d', Ntrials(ktrials)));
    xlabel('Number of spots');
    ylabel('Probability');
    
    hold on; % freeze the plot so you can plot something else
    % the theoretical pmf is pk = 1/Nsides, k=1,2,...,Nsides because the die is fair
    stem(1:Nsides, (1/Nsides)*ones(1,Nsides), 'r', 'LineWidth', 2 ); % plot the theoretical
    hold off; %to release the figure hold
    
    % Professional quality plots always have axes labels
    xlabel('Number of spots');
    ylabel('Probability');
    axis([0 Nsides+1 0 ceil(10/Nsides)/10]);  %...and always have an appropriate scale
    
    % ...and always have a title
    
    %...and almost always have a grid
    grid on
    
    %...and a legend.
    legend('Scaled histogram','Prob Mass Fnc (PMF)');
    disp(['For ',int2str(Ntrials(ktrials)),' sample mean: ',num2str(sample_mean),' sample variance: ',num2str(sample_var)]);
end % loop on the trials

mean_th = (Nsides+1)/2; % compute the theoretical mean or average using form given in assignment
var_th = ((Nsides^2 - 1)/12); % compute the theoretical variance using form given in assignment.

disp(['Theoretical mean = ',num2str(mean_th),' theoretical variance: ',num2str(var_th)]);
disp('-----------');
disp(' ');

% account for the fact that you want separate plots for each section

% Section 3.2 PMF for Binary Strings
disp('Section 3.2 PMF for Binary Strings');
n = 100;
ntrials = [100, 1000, 10000];
p1 = [0.5, 0.2, 0.7];
bins_edges = 0.5:1:(n + 0.5);

population_means = zeros(1, length(p1));
population_vars = zeros(1, length(p1));

% Loop over different probability values
for kp1 = 1:length(p1)
    % Calculate population mean and variance for current p1
    population_means(kp1) = 1 / p1(kp1);
    population_vars(kp1) = (1 - p1(kp1)) / (p1(kp1))^2;

    % Create a new figure for each p1
    figure;
    for ktrials = 1:length(ntrials)
        % Create a subplot for each number of trials
        subplot(length(ntrials), 1, ktrials);
        rand_nums = rand(ntrials(ktrials), n);
        work = (rand_nums <= p1(kp1));
        data = zeros(1, ntrials(ktrials));

        % Find the index of the first '1' in each string
        for kn = 1:ntrials(ktrials)
            i1 = find(work(kn,:) == 1);
            if isempty(i1)
                data(kn) = n + 1;
            else
                data(kn) = min(i1);
            end
        end

        % Define the index and PMF values for the stem plot
        firstOneIndex = 1:ceil(5 / p1(kp1));
        pmfValues = p1(kp1) * (1 - p1(kp1)).^(firstOneIndex - 1);

        % Plot histogram and theoretical PMF
        histogram(data, 'BinEdges', bins_edges, 'Normalization', 'probability');
        hold on;
        stem(firstOneIndex, pmfValues, 'r', 'LineWidth', 1);
        hold off;

        % Set labels, title, and legend
        xlabel('Index of the first 1 in the string');
        ylabel('Probability');
        axis([0, ceil(5 / p1(kp1)), 0, p1(kp1) + 0.1]);
        title(['Section 3.2: ', int2str(ntrials(ktrials)), ' Strings of 100 Binary Values, p1 = ', num2str(p1(kp1))]);
        grid on;
        legend(['Histogram, p1 = ', num2str(p1(kp1))], 'Probability Mass Func. (PMF)');

        % Calculate sample mean and variance
        sample_mean1 = mean(data);
        sample_var1 = var(data);

        % Display results
        disp(['For p1 = ', num2str(p1(kp1)), ', N = ', int2str(ntrials(ktrials)), ':']);
        disp(['    Sample Mean = ', num2str(sample_mean1), ', Sample Variance = ', num2str(sample_var1)]);
        disp(['    Population Mean = ', num2str(population_means(kp1)), ', Population Variance = ', num2str(population_vars(kp1))]);
    end
end

disp('-----------');
disp(' ');




% Section 3.3 Laplace distributed using rand2x provided 
disp(' ');
disp('Section 3.3 Laplace distributed');

% Set up new plots as necessary.  Remember, you need ALL of the plots
Jtrials = [100, 1000, 100000];
lambda = 0.7;
pop_mean2 = 0;
pop_var2 = (2/((0.7)^2));
for jtrials = 1:length(Jtrials)
    
    if (Jtrials(jtrials) == 100)
        dx = .8;
        bin_edges1 = -8:dx:8;
    end 
    if (Jtrials(jtrials) == 1000)
        dx = .5;
        bin_edges1 = -12:dx:12;
    end 
    if (Jtrials(jtrials) == 100000)
        dx = .2;
        bin_edges1 = -20:dx:20;
    end 
    
    data1 = rand2x(1, Jtrials(jtrials), lambda);
    x_min = min(data1);
    x_max = max(data1);
    x = x_min:(1/Jtrials(jtrials)):x_max;
     figure;
    hist_raw2 = histogram(data1, 'BinEdges', bin_edges1);
    %figure;
    hist_norm2 = histogram(data1, 'Normalization', 'pdf', 'BinEdges', bin_edges1);
    pdf = (0.5 * lambda * exp(-lambda * abs(x)));
    hold on;
    plot(x, pdf, 'r', 'LineWidth', 1);
    hold off;
    xlabel('Value of the Random Variable')
    ylabel('Probability Density');
    axis([x_min - (1/2), x_max + (1/2), 0, 0.35]);
    title(['Section 3.3: ','N = ',int2str(Jtrials(jtrials)), ' trials of f(x) = 0.5𝜆e^{-𝜆|x|}, 𝜆 = ', num2str(lambda)]);
    grid on;
    legend('Scaled histogram', 'Prob Dens Func. (PDF)');
   
    sample_mean2 = mean(data1);
    sample_var2 = var(data1);
    pop_mean2 = 0;
    pop_var2 = (2/(.7^2));
    disp(['For N = ', num2str(Ntrials(ktrials)), ':']);
      disp(['    Sample Mean = ', num2str(sample_mean2), ', Sample Variance = ', num2str(sample_var2)]);
      disp(['    Population Mean = ', num2str(pop_mean2) ', Population Variance = ', num2str(pop_var2)]);
end
disp('-----------');
disp(' ');


% Section 3.4  N(0,1)
disp(' ');
disp('Section 3.4 Samples from N(0,1)');

Mtrials = [100, 1000, 100000];
pop_mean3 = 0;
pop_var3 = 1;
for mtrials = 1:length(Mtrials)
    
    if (Mtrials(mtrials) == 100)
        dx = .5;
        bin_edges2 = -4:dx:4;
    end 
    if (Mtrials(mtrials) == 1000)
        dx = .25;
        bin_edges2 = -4:dx:4;
    end 
    if (Mtrials(mtrials) == 100000)
        dx = .1;
        bin_edges2 = -4:dx:4;
    end 
    data2 = randn(1, Mtrials(mtrials));
    x_min1 = min(data2);
    x_max1 = max(data2);
    x1 = x_min1:1/(Mtrials(mtrials)):x_max1;
    figure;
    
    hist_raw3 = histogram(data2, 'BinEdges', bin_edges2);
    %figure;
    hist_norm3 = histogram(data2, 'Normalization', 'pdf', 'BinEdges', bin_edges2);
    
    PDF = (1/sqrt(2*pi)) * (exp((-(x1).^2)/2));
    
    hold on;
    plot(x1, PDF, 'r', 'LineWidth', 1);
    hold off;
    xlabel('Value of the Random Variable')
    ylabel('Probability Density');
    axis([x_min1 - (1/2), x_max1 + (1/2), 0, .5]);
    title(['Section 3.4: ', int2str(Mtrials(mtrials)), ' trials of N(0,1)']);
    grid on;
    legend('Scaled histogram', 'Prob Dens Func. (PDF)');
    
    sample_mean3 = mean(data2);
    sample_var3 = var(data2);
    disp(['For N = ', num2str(Ntrials(ktrials)), ':']);
      disp(['    Sample Mean = ', num2str(sample_var3) ', Sample Variance = ', num2str(sample_mean3)]);
      disp(['    Population Mean = ', num2str(pop_mean3), ', Population Variance = ', num2str(pop_var3)]);
end
disp('-----------');
disp(' ');

% Section 3.5  N(2,6.5) per S24 assignment
disp(' ');
disp('Section 3.5 Samples from N(2,6.25)');
pop_mean4 = 2;
pop_var4 = 6.25;
Ptrials = [100, 1000, 100000];
for ptrials = 1:length(Ptrials)
     if (Ptrials(ptrials) == 100)
        dx = 1;
        bin_edges3 = -20:dx:20;
    end 
    if (Ptrials(ptrials) == 1000)
        dx = .5;
        bin_edges3 = -20:dx:20;
    end 
    if (Ptrials(ptrials) == 100000)
        dx = .1;
        bin_edges3 = -20:dx:20;
    end 
    data3 = (sqrt(pop_var4)) * (randn(1, Ptrials(ptrials))) + pop_mean4;
    x_min2 = min(data3);
    x_max2 = max(data3);
    x2 = x_min2:1/(Ptrials(ptrials)):x_max2;
     figure;
    
    hist_raw4 = histogram(data3, 'BinEdges', bin_edges3);
    %figure;
    hist_norm4 = histogram(data3, 'Normalization', 'pdf', 'BinEdges', bin_edges3);
    
    PDF = (1/sqrt(2*pi*pop_var4)) * (exp((-(x2 - pop_mean4).^2)/(2*pop_var4)));
    
    hold on;
    plot(x2, PDF,'r', 'LineWidth', 1);
    hold off;
    xlabel('Value of the Random Variable')
    ylabel('Probability Density');
    axis([x_min2 - (1/2), x_max2 + (1/2), 0, .25]);
    title(['Section 3.5: ', 'N = ',int2str(Ptrials(ptrials)), ' trials of N(2,6.25)']);
    grid on;
    legend('Scaled histogram', 'Prob Dens Func. (PDF)');
   
    sample_mean4 = mean(data3);
    sample_var4 = var(data3);
disp(['For N = ', num2str(Ntrials(ktrials)), ':']);
      disp(['    Sample Mean = ', num2str(sample_mean4), ', Sample Variance = ', num2str(sample_var4)]);
      disp(['    Population Mean = ', num2str(pop_mean4), ', Population Variance = ', num2str(pop_var4)]);
end
disp('-----------');


% 3.6 

% Estimate the requested probability using the raw histogram from 3.5
%   This is a sum followed by a division
count = sum(data >= -0.5 & data < 4.5); % Count the number of trials between -0.5 and +4.5
prob_raw = count / Ntrials(ktrials); % Scale to a probability
disp(['Probability using raw histogram: ', num2str(prob_raw)]);

% Estimate the required probability using the normalized histogram from 3.5
% above.  This is equivalent to a Riemann sum from Calc II
bin_width = bin_edges(2) - bin_edges(1); % Width of the bins
bins = bin_edges(1:end-1) + bin_width / 2; % Bin centers
counts = histcounts(data, bin_edges); % Counts in each bin
prob_norm = sum(counts(bins >= -0.5 & bins < 4.5)) * bin_width; % Riemann sum

disp(['Probability using normalized histogram: ', num2str(prob_norm)]);
% Estimate the required probability by integration.  Hint: you could use
% the Q function!  or the erfc equivalent.
% Using the Q-function for the normal distribution
prob_int = 0.5 * (erf((4.5 - pop_mean4) / sqrt(2 * pop_var4)) - erf((-0.5 - pop_mean4) / sqrt(2 * pop_var4))); % Q function or erf function

disp(['Probability using numerical integration: ', num2str(prob_int)]);
% That's it, you're done!
