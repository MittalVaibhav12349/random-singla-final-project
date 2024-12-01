function [output, error] = subband_af(x, di, sysorder, L, mu, sigma)

% Initialization
N = length(x); % the number of Iterations
b = fir1(sysorder-1,0.5); % FIR system to be identified
n = 0.01*randn(N,1); % Uncorrelated noise signal
d = filter(b,1,di)+n; % Desired signal = output of H + Uncorrelated noise signal
w = zeros(sysorder, L); % Initially filter weights are zero
x_sb = zeros(N, L);
y_sb = zeros(N, L);
ERLE = zeros(N, 1);
SNR = zeros(N, 1);

% Parameters for ERLE calculation
window_size = 100; % Moving average window size for ERLE
start_point = 500; % Start calculating after 500 samples to avoid initial transients
min_error_power = 1e-10; % Minimum error power to avoid extreme ERLE values

for n = sysorder : N
    % Analysis filter bank
    for k = 1:L
        x_sb(n, k) = x(n) * cos((2 * pi * k * n) / L);
    end
    
    % Subband processing
    for k = 1:L
        % Filter the subband signal
        y_sb(n, k) = filter(w(:, k), 1, x_sb(n, k));
        
        % Calculate the error signal
        e_sb(n, k) = d(n) - y_sb(n, k);
        
        % Update the filter coefficients
        w(:, k) = w(:, k) + mu * e_sb(n, k) * x_sb(n, k);
    end
    
    % Synthesis filter bank
    y(n) = 0;
    for k = 1:L
        y(n) = y(n) + y_sb(n, k) * cos((2 * pi * k * n) / L);
    end
    
    % Calculate the overall error signal
    e(n) = d(n) - y(n);
    
    % Calculate ERLE and SNR, avoiding division by very small numbers
    if n > start_point
        error_power = sum(e(max(1,n-window_size+1):n).^2);
        error_power = max(error_power, min_error_power); % Avoid division by zero
        ERLE(n) = 10 * log10(sum(d(max(1,n-window_size+1):n).^2) / error_power);
    end
    
    % Calculate SNR
    SNR(n) = 20 * log10(rms(d(1:n)) / rms(e(1:n)));
end

mean_ERLE = mean(ERLE(start_point:end));
mean_SNR = mean(SNR(start_point:end));

% Plots
figure();
hold on;
plot(ERLE, 'g');
title('System performance(ERLE)');
xlabel('Samples');
ylabel('Power(dB)');
legend('ERLE');

figure();
xc = xcorr(y, x);
plot(xc);
title('Cross Correlation');
xlabel('Samples');
ylabel('Magnitude');

figure();
hold on;
plot(x, 'g');
plot(y, 'r');
semilogy(abs(e), 'm');
title('System output');
xlabel('Samples');
ylabel('Magnitude');
legend('Desired', 'Output', 'Error');

% Return output
output = y;
error = e;

end

