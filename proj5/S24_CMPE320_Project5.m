%S21 CMPE320 Project 5
%  EFCL 4/28/2021
% Modified for S22
% Modified for S23, students will run the whole script and comment on
% results
% Modified for S24, April 22, 2024


close all
clear
disp(['This version does not require student programming! (You''welcome)']);
disp(['This takes a while!']);
N=10000; % you might have to adjust this for memory. Make it as big as you can.
Nt=1001;
Ntd2 =(Nt+mod(Nt,2))/2; %We'll need this to access the "middle" of output of xcorr


% Initialize Rxx to zeros;  columns are times, each row is a different trial
Rxx = zeros(N,Nt);

% % Initial x to zeros; columns are times, each row is a different trial 
%   x  will store the sample functions
x = zeros(N,Nt);


% Use a loop to generate the N different random sample functions
%     each sample function is 1 x Nt array from randn  N(0,1)

for k = 1:N
    
    % Generate sample function and store in k-th row of x
    x(k,:)=randn(1,Nt);
    
    % Temporarily store the output of xcorr(   )/Nt using the k-th row of x as both inputs 
    junk = xcorr(x(k,:),x(k,:))/Nt; %Matlab cross-correlation creates 2*Nt-1 points
    Rxx(k,:) = junk(Ntd2:Ntd2+Nt-1);% just save the middle Nt;
end;

R_XX = mean(Rxx); % "expected value" over Nrows of Rxx (down the columns)
% fine the max of Rxx

%New Figure
figure(1);

% time array from 0 to Nt
t = [0:Nt-1];

% tau array from [-Ndt/2+1:Ntd/2]

% create and title 4 subplots 
%  1) t vs x  (all the sample functions; across the rows)
%  2) t vs mean(x)  (down the columns)
%  3) tau vs ensemble Rxx (all the functions), xcorr in time as an estimate
%  4) tau vs mean of Rxx computed above, down the columns

subplot(4,1,1);
plot(t,x)
title(['ensemble X, N = ',int2str(N)])
xlabel('time, t');
subplot(4,1,2);
plot(t,mean(x),'b');
title(['ensemble mean'])
xlabel('time, t')
subplot(4,1,3);
tau = [-Ntd2+1:Ntd2-1];
plot(tau,Rxx);
title(['ensemble R_x_x']);
xlabel('\tau = t_1-t_2 = t_2-t_1')
subplot(4,1,4);
plot(tau,R_XX);
title(['R_X_X']);
xlabel('\tau = t_1-t_2 = t_2-t_1')


% Now do the entire thing three more times using a sliding window filter

nfig=1;

% Set the array of FIR filter lengths
L=[16 25 64]; %; Spring 2024 modification.

% Loop on  the filter lengths
for kL=1:length(L)
LL= L(kL);

% Set the current length
thisLength = LL;

% Set the coefficients for this FIR filter for MATLAB function filter (trust me!)

b=ones(1,thisLength)/thisLength; % L point sliding window

% Initialize Ryy and y, as Rxx and x

Ryy = zeros(N,Nt); % 600 points to accommodate transient
y = zeros(N,Nt);

%Loop on  the sample functions
for k = 1:N
    
    %Generate xin as the iid Gaussian, as above, but this time with
    %  Nt+thisLength columns (extra columns)
    xin = randn(1,Nt+LL); % iid Gaussian variance 1 mean zero
    
    % create a temporary output for the filter output
    ytemp = filter(b,1,xin); % filter with the sliding window
    
    % Save into the k-th row of y, but only save the LL+1:end columns of
    % ytemp.  This remove the filter transient from beginning of the filter
    % output
    y(k,:) = ytemp(LL+1:end);
    
    % Create the temporary output of xcorr using the k-th row of y for both
    % inputs; then scale by Nt as before
    junk = xcorr(y(k,:),y(k,:))/Nt;
    
    % Store this output in k-th row of Ryy  save the middle Nt samples as
    % before
    Ryy(k,:) = junk(Ntd2:Ntd2+Nt-1);% just save the middle Nt;

end;  % loop on sample functions

R_YY=mean(Ryy); % create the mean down the  columns, as  before
% find the  max R_YY, as before, deleted for Spring 2023


nfig=nfig+1;
figure(nfig);
t = [0:Nt-1];

% repeat the four previous plots, using y and mean autocorrelation
% make sure to use a new figure each time.
subplot(4,1,1);
plot(t,y);
title(['ensemble Y, L = ',int2str(LL),' N = ',int2str(N)])
subplot(4,1,2);
xlabel('time, t')
plot(t,mean(y));
title('mean Y');
xlabel('time, t')
subplot(4,1,3);
tau = [-Ntd2+1:Ntd2-1];
plot(tau,Ryy);
title(['ensemble Ryy']);
xlabel('\tau = t_1-t_2 = t_2-t_1')
subplot(4,1,4);
plot(tau,R_YY,'.');
title(['R_Y_Y']);
xlabel('\tau = t_1-t_2 = t_2-t_1')

% display the ratio of the peak autocorrelations RXX/RYY
% deleted for Spring 2023

end; %Loop on another length




    
    






