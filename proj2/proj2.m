% 2.1 

% Given values
A = 1.5; % Known amplitude of the signal
sigma2 = 0.225; % Given noise variance 
sigma = sqrt(sigma2); % Standard deviation of the noise
N_trials = 10000; % Number of trials

% Generate random signal X with values +A or -A with equal probability
X = (rand(N_trials, 1) > 0.5) * 2 * A - A; % Generates +A or -A

% Generate noise N, which is Gaussian distributed with mean 0 and variance sigma2
N = randn(N_trials, 1) * sigma;

% The received signal R is the sum of X and N
R = X + N;

% Scatterplot of R
figure;
scatter(1:N_trials, R, '.');
title('Section 2.1: Scatterplot of R');
xlabel('Number of trials');
ylabel('R');

% Histogram of R
figure;
histogram(R, 'Normalization', 'pdf');
title('Section 2.1: Histogram of R');
xlabel('R');
ylabel('f_R(r)');

% Calculation of the analytical PDF of R using the law of total probability
r_values = linspace(min(R), max(R), 1000); % Values of r for PDF calculation

% Defining the Gaussian PDFs conditioned on X = A and X = -A
pdf_RA = (1 / (sigma * sqrt(2 * pi))) * exp(-(r_values - A).^2 / (2 * sigma^2));
pdf_RmA = (1 / (sigma * sqrt(2 * pi))) * exp(-(r_values + A).^2 / (2 * sigma^2));

% Total PDF of R combining the conditioned PDFs, since X = A and X = -A are equally likely
f_R = 0.5 * pdf_RA + 0.5 * pdf_RmA; 

% Plot the analytical PDF over the histogram
hold on; % Hold the histogram to overlay the analytical PDF
plot(r_values, f_R, 'LineWidth', 2, 'Color', 'red');
legend('Simulated PDF', 'Analytical PDF');
grid on;
hold off;


% 2.2 Method 1

% 2.2.2 

k = 4; % Given scaling constant for the perfect diode detector

% Apply the scaled perfect diode detector to R
% S is kR if R >= 0, otherwise S is 0
S = max(R, 0) * k;

% Scatterplot of S vs. R
figure;
scatter(R, S, '.');
title('Section 2.2: Scatterplot of S vs. R');
xlabel('R');
ylabel('S');

% Histogram of S to visualize the simulated PDF of S
figure;
histogram(S, 'Normalization', 'pdf');
title('Section 2.2: Histogram of S');
xlabel('S');
ylabel('f_S(s)');

% Ensuring s_values is a row vector for proper plotting
s_values = linspace(min(S), max(S), 1000);

% Calculating the analytical PDF using the expression derived in 2.2.1
% Using the variables A, sigma, and k already defined above
f_S_analytical = (1 / (8 * sqrt(2 * pi * sigma2))) * ...
                  (exp(-((s_values - 4 * A).^2 / (2 * 4^2 * sigma2))) + ...
                   exp(-((s_values + 4 * A).^2 / (2 * 4^2 * sigma2))));

% Overlaying the analytical PDF on the histogram
hold on;
plot(s_values, f_S_analytical, 'LineWidth', 2, 'Color', 'red');
legend('Simulated PDF', 'Analytical PDF');
grid on;
hold off;

% Calculate the expected values (means) of S and R samples
E_S = mean(S); % Expected value of S
E_R = mean(R); % Expected value of R

% Calculate g(E[R]) assuming g(R) = 4R for demonstration purposes
% Adjust this calculation based on your specific method's function
g_E_R = 4 * E_R;

% Display the results
disp('Expected value of S (E[R]):');
disp(E_R);
disp('Expected value of S (E[S]):');
disp(E_S);
disp('g(E[R]): ');
disp(g_E_R);

% 2.3 Method 2

% 2.3.2

% Transform R to S = k| R |
S = k * abs(R);

% Scatterplot of S vs. R
figure;
scatter(R, S, '.');
title('Section 2.3: Scatterplot of S vs. R');
xlabel('R');
ylabel('S');

% Histogram of S to visualize the simulated PDF of S
figure;
histogram(S, 'Normalization', 'pdf');
hold on; % Hold the figure for overlaying the analytical PDF

% Generating values for s for the analytical PDF plot
s_values = linspace(0, max(S), 100); % for a smoother curve

% this is a placeholder for the actual analytical expression
% we can replace it with the correct formula from 2.3.1
analytical_pdf_S = (1/(k*sqrt(2*pi*sigma2))) .* exp(-((s_values/k - A).^2 / (2*sigma2))) + ...
                   exp(-((s_values/k + A).^2 / (2*sigma2)));

% Plot the analytical PDF over the histogram
plot(s_values, analytical_pdf_S, 'LineWidth', 2, 'Color', 'red');
title('Section 2.3: Histogram of S');
xlabel('S');
ylabel('f_S(s)');
legend('Simulated PDF', 'Analytical PDF');
grid on;
hold off; % Release the figure


% Calculate the expected values (means) of S and R samples
E_S = mean(S); % Expected value of S
E_R = mean(R); % Expected value of R

% Function g applied to the expected value of R
% For Method 2, g(R) = k * |R|, so we compute g(E[R]) as k * the absolute value of E[R]
g_E_R = k * abs(E_R);

% Display the expected values
disp('Expected value of S (E[R]):');
disp(E_R);
disp('Expected value of S (E[S]):');
disp(E_S);
disp('g(E[R]): ');
disp(g_E_R);

% 2.4 Method 3

% 2.4.2 

% 2.4.2 Method 3 - Square Law Detector

% Given scaling constant for the square law detector (as previously defined k=4)
k = 4; 

% Using the simulated data from Section 2.1 to simulate the action of the perfect square law detector
% S is kR^2 if R >= 0, otherwise S is 0
% Since R can take negative values and we are squaring it, we do not need to use max(R,0)
S = k * R.^2;

% Scatterplot of S vs R
figure;
scatter(R, S, '.');
title('Section 2.4: Scatterplot of S vs. R');
xlabel('R');
ylabel('S');

% Histogram of S to visualize the simulated PDF of S
figure;
histogram(S, 'Normalization', 'pdf');
title('Section 2.4: Histogram of S');
xlabel('S');
ylabel('f_S(s)');

% Define the range of s for the analytical PDF of S
s_values = linspace(0, max(S), 1000);


f_S_analytical = (1./(4.*sqrt(2.*pi.*sigma2))).*(...
    exp(-((sqrt(s_values./k) - A).^2 ./ (2.*sigma2))) + ...
    exp(-((sqrt(s_values./k) + A).^2 ./ (2.*sigma2))) ...
    );

% Overlaying the analytical PDF on the histogram
hold on;
plot(s_values, f_S_analytical, 'LineWidth', 2, 'Color', 'red');
legend('Simulated PDF', 'Analytical PDF');
grid on;
hold off;

% Calculate the expected values (means) of S and R samples
E_S = mean(S); % Expected value of S
E_R = mean(R); % Expected value of R

% Function g applied to the expected value of R for Method 3
% Assuming the form of g(R) is k*R^2
g_E_R = k * E_R.^2;

% Display the expected values
disp('Expected value of S (E[S]):');
disp(E_S);
disp('g(E[R]): ');
disp(g_E_R);

E_S = mean(S); % Expected value of S
E_R = mean(R); % Expected value of R
g_E_R = k * E_R.^2; 

% Display the results
disp('Expected value of R (E[R]):');
disp(E_R);
disp('Expected value of S (E[S]):');
disp(E_S);
disp('g(E[R]): ');
disp(g_E_R);