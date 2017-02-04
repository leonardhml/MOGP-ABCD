Fs = 1000;            % Sampling frequency
T = 1/Fs;             % Sampling period
L = 1000;             % Length of signal
t = (0:L-1)*T;        % Time vector

x1 = cos(2*pi*50*t);
Y = fft(x1)
