%   Skeleton 
%
close all;
clear;

disp('CMPE320 Spring 2024 Project 4: BASK');

% Constants for the project
A = 3; % given in 2024 instructions
p0all = [0.5, 0.25]; % 2024 instructions, the project includes two values of p0

% Remember that p0 = 0.5 is Maximum Likelihood (ML)
% p0 # 0.5 is Maximum A Posteriori (MAP)

% 2.1  Plot the MAP threshold curve as a function of p0
p0 = 0.01:0.01:0.99; % per instructions
gamma_dB = 10; % in decibels, gamma_dB = 10*log10(A^2/sigma^2)
gamma = 10^(gamma_dB / 10); % convert to power ratio using the definition of gamma and gamma_dB

% Compute the equivalent value for sigma2 = sigma^2 using A and gamma
sigma2 = A^2 / gamma;

tauMAP = A/2 * log((1-p0) ./ (p0)) * sqrt(sigma2); % use the result from your derivation of the threshold
figure(1);
plot(p0, tauMAP, 'LineWidth', 2);
xlabel('p_0 = Pr[b_k = 0]');
ylabel('\tau_{MAP}');
grid on;
title(['Section 2.1: \tau_{MAP} with \gamma = ', num2str(gamma_dB), ' dB, A = ', num2str(A), ' and \sigma^2 = ', num2str(sigma2)]);


% 2.2

figure(2); % Create a new figure
p0 = p0all(2); % from 2024 Assignment, p0 # 0.5 is MAP
gamma_dB = 10; % in decibels
gamma = 10^(gamma_dB / 10); % power ratio
sigma2 = A^2 / gamma; % variance value per definition
Ntrials = 500000; % sufficient trials for this exercise
thld_MAP = A/2 * log( (1 - p0)/ (p0)) * sqrt(sigma2); % use your derivation of the MAP threshold

% Generate B, the binary information
B = (rand(1, Ntrials) >= p0); % 0 = binary zero, 1 = binary one

% Convert B to Messages, map 0 to +A, 1 to -A
% B = {0,1}, (0.5-B) = {0.5, -0.5}, (0.5-B)*2*A = {+A,-A}
M = (0.5 - B) * 2 * A; % map 1 to -A, 0 to +A
N = randn(1, Ntrials) * sqrt(sigma2); % Gaussian with zero mean and proper variance

% Use the signal model from the assignment
R = M + N; % The received signal
Rhist = histogram(R, 'Normalization', 'pdf'); % Create and plot the histogram as a pdf

% Fine-grained values of R
r = -2 * A:0.01:2 * A;

% Derivation of fR(r)
fRr = p0 * (1/(sqrt(2 * pi * sigma2)) * exp(-(r - A).^2 / (2 * sigma2))) + ...
       (1 - p0) * (1/(sqrt(2 * pi * sigma2)) * exp(-(r + A).^2 / (2 * sigma2)));

% Plot the analytical PDF on top of the histogram
hold on;
plot(r, fRr, 'LineWidth', 2);
xline(thld_MAP, '--r', 'LineWidth', 2); % Plot the MAP threshold value
hold off;

% Professional labels, note use of subscripts and Greek symbols
grid on;
xlabel('r');
ylabel('f_R(r)');
title(['Section 2.2: f_R(r) for p0 = ', num2str(p0), ' Ntrials = ', int2str(Ntrials), ' \gamma_{dB} = ', num2str(gamma_dB)]);
text(-0.5, 0.26, ['\tau_{MAP} = ', num2str(thld_MAP)]);

%2.3
% Section 2.3.3 - Simulate ML Detector
% Define the Q function
Q = @(x) 0.5 * erfc(x / sqrt(2));

% Set parameters
gamma_dB = [[1:0.5:10] [10:0.25:14]]; % gamma = A^2/sigma^2 in dB
gamma = 10.^(gamma_dB / 10); % as power ratio, same size as gamma_dB
A = 3; % Signal amplitude
sigma = sqrt(A^2 ./ gamma); % Noise standard deviations
Nbits = 5000000; % Number of bits

% Section 2.3.3 - Simulate ML Detector
tau_ML = 0; % ML threshold

% Initialize result arrays
results_ML = zeros(1, length(gamma));
pBT_ML = zeros(1, length(gamma));

% Loop over SNR values
for kSNR = 1:length(gamma)
    % Generate message bits (0 or 1)
    b = randi([0 1], 1, Nbits);
    % Map to BASK symbols (+A or -A)
    m = A * (2 * b - 1);
    % Generate noise
    n = sigma(kSNR) * randn(1, Nbits);
    % Received signal
    r = m + n;
    % ML Detection
    b_hat = r >= tau_ML;
    % Calculate bit errors
    errors = mod(b_hat - b, 2);
    results_ML(kSNR) = sum(errors) / Nbits;
    % Theoretical Probability of Error (ML)
    pBT_ML(kSNR) = Q(A / sigma(kSNR));
end

% Plot results (ML Detector)
figure;
h = semilogy(gamma_dB, pBT_ML, 'k-', gamma_dB, results_ML, 'ro');
set(h, 'LineWidth', 1.5);
legend('Theoretical', 'Simulated');
xlabel('\gamma_{dB} (dB)');
ylabel('Probability of Bit Error');
title('Section 2.3: ML Detector Performance (p_0 = 0.5)');
grid on;

% Section 2.4.2 - Simulate MAP Detector
p0 = 0.25; % Prior probability of binary 0
tau_MAP = (sigma.^2 / (2 * A)) .* log((1-p0) / ( p0)); % MAP threshold

% Initialize result arrays
results_MAP = zeros(1, length(gamma));
pBT_MAP = zeros(1, length(gamma));

% Loop over SNR values
for kSNR = 1:length(gamma)
    % Generate message bits (0 or 1)
    b = randi([0 1], 1, Nbits);
    % Map to BASK symbols (+A or -A)
    m = A * (2 * b - 1);
    % Generate noise
    n = sigma(kSNR) * randn(1, Nbits);
    % Received signal
    r = m + n;
    % MAP Detection
    b_hat = r >= tau_MAP(kSNR);
    % Calculate bit errors
    errors = mod(b_hat - b, 2);
    results_MAP(kSNR) = sum(errors) / Nbits;
    % Theoretical Probability of Error (MAP)
    pBT_given0 = Q((A - tau_MAP(kSNR)) / sigma(kSNR));
    pBT_given1 = Q((tau_MAP(kSNR) + A) / sigma(kSNR));
    pBT_MAP(kSNR) = p0 * pBT_given0 + (1 - p0) * pBT_given1;
end

% Plot results (MAP Detector)
figure;
h = semilogy(gamma_dB, pBT_MAP, 'k-', gamma_dB, results_MAP, 'ro');
set(h, 'LineWidth', 1.5);
legend('Theoretical', 'Simulated');
xlabel('\gamma_{dB} (dB)');
ylabel('Probability of Bit Error');
title(['Section 2.4: MAP Detector Performance (p_0 = ', num2str(p0), ')']);
grid on;

% Section 2.5 - Compare MAP and ML Detector Performance
figure;
h = semilogy(gamma_dB, pBT_ML, 'b-', gamma_dB, pBT_MAP, 'r-');
set(h, 'LineWidth', 1.5);
legend('ML Detector (p_0 = 0.5)', 'MAP Detector (p_0 = 0.25)');
xlabel('\gamma_{dB} (dB)');
ylabel('Probability of Bit Error');
title('Section 2.5: ML vs MAP Detectors');
grid on;

% Compute and plot ratio
ratio = pBT_MAP ./ pBT_ML;
figure;
plot(gamma_dB, ratio, 'k-', 'LineWidth', 1.5);
xlabel('\gamma_{dB} (dB)');
ylabel('\rho (MAP/ML)');
title('Section 2.5: \rho (MAP/ML) vs \gamma_{dB} (dB)');
grid on;
