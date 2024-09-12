function [R2,lags2,conf2] = whiteness_test(e,lagmax,t)
% [R,lags,conf] = whiteness_test(e,lagmax,t)
% This function gives several attributes of a signal to determine if it is approximately a white noise:
% i) autocorrelation values, ii) spectrum (it has been approximately a flat spectrum), iii) mean value
% (Ee ~= 0), % iv) Ljung-Box test with significant level equal to 0.05 (Pass = insufficient evidence exists
% to reject the null hypothesis that the signal is a white noise through LAGMAX lags)
%
% e         Signal data e(t)
% lagmax    Maximum lag in e(t-tau). Optional (by default, lagmax = 50)
% t         Independent variable. Optional (by default the sample time is 1)
% lags      Lag vector
% R         Autocorrelation values R(lag) for every lag
% conf      Limit of confidence interval (99%)
%
% Examples:
% N = 200000; t = 0:N-1; e = randn(N,1); whiteness_test(e)
% N = 200000; t = 0:N-1; e1 = rand(N,1); e = e1-mean(e1); whiteness_test(e)
% N = 200000; t = 0:N-1; e = wgn(N,1,1); whiteness_test(e,100)
% e = wgn(2000,1,1); Gd = tf([1 0],[1 0.8],'Ts',1); u = lsim(Gd,e); whiteness_test(u)
% N = 20000; u = idinput([N, 1, 1], 'PRBS', [0, 1], [-1,1]);
%
% Carlos Mario Vélez S.
% Universidad EAFIT
% Medellín, Antioquia, Colombia
% % https://sis-control.blogspot.com/

N = length(e);
if size(e,2)==1
    e = e';
end
switch nargin
    case 1
        lagmax = 50;
        t = 0:N-1;
    case 2
        t = 0:N-1;
end
[R,lags] = xcorr(e,lagmax,'coeff'); % Autocorrelation values R(lag) for every lag
conf = sqrt(2)*erfcinv(0.01)/sqrt(N); % Limit of confidence interval (99%)

% Ljung-Box Test. See: https://bit.ly/3oqnkIR
R1 = R(lagmax+2:2*lagmax+1);
lags2 = lags(lagmax+2:2*lagmax+1)';
idx = N-lags2;
Q = N*(N+2)*sum(R1.^2./idx);
pvalue = 1-chi2cdf(Q,lagmax);
if pvalue <= 0.05
    lbtest = 'Not pass with α = 0.05';
else
    lbtest = 'Pass with α = 0.05'; %
end

switch nargout
    case 0
        subplot(2,2,1) % Plot of signal, mean value (Ee ~= 0) and Ljung-Box Test (lbtest = 0)
        plot(t,e)
        xlabel('Time')
        ylabel('e(t)')
        title({['Signal, μ_e = ' num2str(mean(e))]; ['Ljung-Box Test = ' lbtest ' and ' num2str(lagmax) ' lags']})
        ylim([-max(abs(e))-0.1  max(abs(e))+0.1])

        subplot(2,2,2) % Signal spectrum (approximately a flat frequency spectrum)
        pwelch(e)

        subplot(2,2,3) % Plot of autocorrelation values for several lags: Ee(t)e(t-lag)
        plot(lags,R,'b.')
        line([lags(1) lags(end)],conf*[1 1],'color','red','linestyle','--','linewidth',1)
        line([lags(1) lags(end)],-conf*[1 1],'color','red','linestyle','--','linewidth',1)
        line([0 0],[0 max(R)],'color','blue','linestyle','--','linewidth',0.1)
        xlabel('Lag')
        ylabel('Amplitude')
        title({'Autocorrelation';'(confidence interval = 99%)'})
        ylim([ min([R,-conf])-0.002    max([R(1:ceil(lagmax)),conf])+0.002 ])

        subplot(2,2,4)
        histogram(e)
        xlabel('Values')
        ylabel('Number of cases')
        title('Signal histogram')
    case 1
        R2 = R;
    case 2
        R2 = R; lags2 = lags;
    case 3
        R2 = R; lags2 = lags; conf2 = conf;
end

