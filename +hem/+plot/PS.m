function [f, P, LF, HF, LFn, HFn] = PS(t, x, make_plot)
%Estimates the power spectrum using an autoregressive model generated with
%the Yule-Walker algorithm. But, it is easy to change to another spectral
%density estimation method. LF and HF are also calculated, along with
%normalized values LFn and HFn. Input data should be a short period of
%time... 5 to 10 minutes usually.

Nf = 1000;
f = linspace(0,0.4,Nf);
Fs = 1/(t(2)-t(1));

% [P,Prob] = lomb(t,x-mean(x),f);
% [P,F] = peig(x,5,f,Fs);
% [Pxx,F] = pwelch(x, round(2^11), [], f, Fs);
% [Pxx,F] = pwelch(x, length(x), 0, f, Fs);
% P = log(Pxx);
[P, F] = pyulear(x-mean(x), 12, f, Fs);
LF = sum(P((0.04 <= f) & (f <= 0.15)));
HF = sum(P((0.15 <= f) & (f <= 0.5)));
LFn = LF/(LF+HF);
HFn = HF/(LF+HF);
% keyboard
P = 10*log10(P); % Convert to dB

if make_plot
    plot(f, P, 'k', 'LineWidth', 3);
%     title(sprintf('LF: %g; HF: %g', LF, HF));
    xlabel('Frequency (Hz)');
    ylabel('Power/frequency (dB/Hz)');
end
