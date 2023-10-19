
% ////////////////////////////////////////////
% Question 1
N=1024;
NQ1 = 2049; %samples
omega = linspace(0, pi, NQ1);
std1 = sqrt(1/26);
std2 = sqrt(0.51);
spectrum1 = std1^2 * (26 + 18 * cos(omega) - 8 * cos(2 * omega));
spectrum2 = std2^2 ./ (1.49 - 1.4 * cos(omega));

figure(1);
plot(omega, spectrum1);
title('Analytic Spectrum for X1');
xlabel('\omega');
ylabel('X1(\omega)');

figure(2);
plot(omega, spectrum2);
title('Analytic Spectrum for X2');
xlabel('\omega');
ylabel('X2(\omega)');
% ////////////////////////////////////////////

% Question 2A1
NQ2A = 1024;
NQ2B = 2048;

rng('default')
s=rng;
w1 = normrnd(0,std1,[1,NQ2A]); % random w1[n] - 1024 samples
rng(s);
w2 = normrnd(0,std2,[1,NQ2B]); % random w2[n] - 2048 samples

w1_0 = normrnd(0,std1); % for the use of x1(1) and x1(2)
x1 = zeros(size(w1));
x1(1) = w1(1)-3*w1_0-4*normrnd(0,std1);
x1(2) = w1(2)-3*w1(1) - 4*w1_0;

for n = 3:length(w1)
    x1(n) = w1(n)-3*w1(n-1)-4*w1(n-2);
end

x2_Q21=zeros(size(w2)); % Q21 for the Question 2.1
x2_0_Q21 = normrnd(0,std2);
x2_Q21(1) = 0.7*x2_0_Q21+w2(1);

for n = 2:length(w2)
    x2_Q21(n) = 0.7*x2_Q21(n-1)+w2(n);
end

% ////////////////////////////////////////////
% Question 2A2
w2_2A2 = normrnd(0,std2,[1,NQ2A]);
x2_Q22=zeros(size(w2_2A2)); % Q22 for the Question 2.2
x2_0_Q22 = normrnd(0,1);
x2_Q22(1) = 0.7*x2_0_Q22+w2_2A2(1);

for n = 2:length(w2_2A2)
    x2_Q22(n) = 0.7*x2_Q22(n-1)+w2_2A2(n);
end
% ////////////////////////////////////////////

% Question 2B1
M = 4096;
x1_2B1 = [x1, zeros(1, M-numel(x1))]; % zero padding
X1_2B1 = fft(x1_2B1);
periSpectFunc_x1 = 1/N*X1_2B1(1:M/2+1).*conj(X1_2B1(1:M/2+1));

% ////////////////////////////////////////////

% Question 2C1
Rx1_2C1 = xcorr(x1);
Rx1_2C1 = Rx1_2C1 / N;
Rx1_2C1 = [Rx1_2C1(N:2*N-1) zeros(1, M/2+1) Rx1_2C1(1:N-1)];
% ////////////////////////////////////////////

corrSpectFunc_x1 = fft(Rx1_2C1);

% ////////////////////////////////////////////

% Question 2D Pre - making x2 signals
% Periodogram
x2_2DP = [x2_Q22, zeros(1, M-numel(x2_Q22))]; % zero padding
X2_2DP = fft(x2_2DP);
periSpectFunc_x2 = 1/N*abs(X2_2DP(1:M/2+1)).^2;

% Correlogram
Rx2_2DP =xcorr(x2_Q22); % Rx2(0) in the middle
Rx2_2DP = 1/N*[Rx2_2DP(length(x2_Q22):length(Rx2_2DP)) zeros(1,length(Rx2_2DP)+2) Rx2_2DP(1:length(x2_Q22)-1)];
corrSpectFunc_x2 = fft(Rx2_2DP);
% ////////////////////////////////////////////

% Question 2D1
figure(3)
subplot(1,2,1)
hold on
plot(omega, abs(periSpectFunc_x1(1:(M/2)+1)))
plot(omega,abs(corrSpectFunc_x1(1:M/2+1)))
plot(omega, spectrum1, color = [0,0,0])
hold off
title('X1 Periodogram Vs Correlogram')
xlabel('\omega')
ylabel('S(\omega)')
legend('Periodogram', 'Correlogram', 'Analytic')
hold off;

% Question 2D2
figure(3)
subplot(1,2,2)
hold on
plot(omega, abs(periSpectFunc_x2(1:(M/2)+1)))
plot(omega,abs(corrSpectFunc_x2(1:M/2+1)))
plot(omega, spectrum2, color = [0,0,0])
hold off
title('X2 Periodogram Vs Correlogram')
xlabel('\omega')
ylabel('S(\omega)')
legend('Periodogram', 'Correlogram', 'Analytic')
hold off;
% ////////////////////////////////////////////


% % Section 2

% Question 3

%****PARAMETERS****
%Additional Global Parameters
std1 = sqrt(1/26);
Mc=100;
NW=2049;
M=4096;
omega = linspace(0, pi, NW);

%Periodogram Parameters
spectrum_periodogram = zeros(1,M,Mc);

%bartlett Parameters
K1_bartlett=16;
L1_bartlett=N/K1_bartlett;
K2_bartlett=64;
L2_bartlett=N/K2_bartlett;
spectrum_bartlett_A = zeros(1,M,Mc);
spectrum_bartlett_B = zeros(1,M,Mc);

%Welch Parameters
K1_welch = 61;
L1_welch = 64;
D1_welch = L1_welch-48;
K2_welch = 253;
L2_welch = 16;
D2_welch = L2_welch-12;
spectrum_welch_A = zeros(1,M,Mc);
spectrum_welch_B = zeros(1,M,Mc);

%Blackman-Tukey Parameters
l1_bt = 4;
l2_bt = 2;
spectrum_blackman_A = zeros(1,M,Mc);
spectrum_blackman_B = zeros(1,M,Mc);

% ////////////////////////////////////////////

for monte = 1:Mc
    % New x1 construction
    wS2Q1 = normrnd(0,std1,[1,N]);
    wS2Q1_0 = normrnd(0,std1); % for the use of x1(1) and x1(2)
    x1_S2Q1 = zeros(size(wS2Q1));
    x1_S2Q1(1) = wS2Q1(1)-3*wS2Q1_0-4*normrnd(0,std1);
    x1_S2Q1(2) = wS2Q1(2)-3*wS2Q1(1) - 4*wS2Q1_0;
    for n = 3:length(wS2Q1)
        x1_S2Q1(n) = wS2Q1(n)-3*wS2Q1(n-1)-4*wS2Q1(n-2);
    end
    x1_S2Q1 = [x1_S2Q1, zeros(1, M-numel(x1_S2Q1))];
    % ////////////////////////////////////////////

    % Periodogram
    X1S2Q1_periodogram = fft(x1_S2Q1);
    spectrum_periodogram(1,:,monte) = 1/N*abs(X1S2Q1_periodogram).^2;
    % ////////////////////////////////////////////
    
    % bartlett
    
    x1_bartlett_A = zeros(K1_bartlett, L1_bartlett);
    X1_bartlett_A = zeros(1, M, K1_bartlett);
    x1_bartlett_B = zeros(K2_bartlett, L2_bartlett);
    X1_bartlett_B = zeros(1, M, K2_bartlett);

    % Fourier transform of x1 for bartlett part A
    for k=1:K1_bartlett
        for n=1:L1_bartlett
            x1_bartlett_A(k,n) = x1_S2Q1((k-1)*L1_bartlett+n);
        end
    end

    for k=1:K1_bartlett
        X1_bartlett_A(1,:,k) = fft(x1_bartlett_A(k,:),M);
    end

    
    % bartlett part A
    for k = 1:K1_bartlett
        spectrum_bartlett_A(1, :, monte) = spectrum_bartlett_A(1, :, monte) + 1/L1_bartlett * abs(X1_bartlett_A(1, :, k)).^2;
    end
    spectrum_bartlett_A(1, :, monte) = 1/K1_bartlett * spectrum_bartlett_A(1, :, monte);


    % Fourier transform of x1 for bartlett part B
    for k=1:K2_bartlett
        for n=1:L2_bartlett
            x1_bartlett_B(k,n) = x1_S2Q1((k-1)*L2_bartlett+n);
        end
    end

    for k=1:K2_bartlett
        X1_bartlett_B(1,:,k) = fft(x1_bartlett_B(k,:),M);
    end
    
    % bartlett part B
    for k = 1:K2_bartlett
        spectrum_bartlett_B(1, :, monte) = spectrum_bartlett_B(1, :, monte) + 1/L2_bartlett * abs(X1_bartlett_B(1, :, k)).^2;
    end
    spectrum_bartlett_B(1, :, monte) = 1/K2_bartlett * spectrum_bartlett_B(1, :, monte);
    % ////////////////////////////////////////////

    % Welch
    
    % Welch part A
    x1_welch_A = zeros(K1_welch,L1_welch);
    for k = 1:K1_welch
        for n = 1:L1_welch
            x1_welch_A(k,n) = x1_S2Q1((k-1)*D1_welch+n);
        end
    end

    X1_welch_A = zeros(1,M,K1_welch);
    for k = 1:K1_welch
        X1_welch_A(1,:,k) = fft(x1_welch_A(k,:),M);
    end

    for k=1:K1_welch
        spectrum_welch_A(1,:,monte) = spectrum_welch_A(1,:,monte) + 1/L1_welch*abs(X1_welch_A(1,:,k)).^2;
    end

    spectrum_welch_A(1,:,monte) = spectrum_welch_A(1,:,monte) / K1_welch;

    % Welch part B
    x1_welch_B = zeros(K2_welch,L2_welch);
    for k = 1:K2_welch
        for n = 1:L2_welch
            x1_welch_B(k,n) = x1_S2Q1((k-1)*D2_welch+n);
        end
    end

    X1_welch_B = zeros(1,M,K2_welch);
    for k = 1:K2_welch
        X1_welch_B(1,:,k) = fft(x1_welch_B(k,:),M);
    end

    for k=1:K2_welch
        spectrum_welch_B(1,:,monte) = spectrum_welch_B(1,:,monte) + 1/L2_welch*abs(X1_welch_B(1,:,k)).^2;
    end

    spectrum_welch_B(1,:,monte) = spectrum_welch_B(1,:,monte) / K2_welch;
%     ////////////////////////////////////////////

    % Blackman-Tukey
    RX1_bt = xcorr(x1_S2Q1(1:N));
    RX1_bt = RX1_bt/N;

    % BT Part A
    RX1_bt_A = [RX1_bt(N+1:N+1+l1_bt) zeros(1, M-2*l1_bt-1) RX1_bt(N+1-l1_bt:N)];
    spectrum_blackman_A(1,:,monte) = abs(fft(RX1_bt_A));
    
    % BT Part B
    RX1_bt_B = [RX1_bt(N+1:N+1+l2_bt) zeros(1, M-2*l2_bt-1) RX1_bt(N+1-l2_bt:N)];
    spectrum_blackman_B(1,:,monte) = abs(fft(RX1_bt_B));
%   //////////////////////////////////////////////

end

% Average
spectrum_periodogram_mean = zeros (1,M);
spectrum_bartlett_A_mean = zeros (1,M);
spectrum_bartlett_B_mean = zeros (1,M);
spectrum_welch_A_mean = zeros (1,M);
spectrum_welch_B_mean = zeros (1,M);
spectrum_blackman_A_mean = zeros (1,M);
spectrum_blackman_B_mean = zeros (1,M);

for mean=1:Mc
    spectrum_periodogram_mean = spectrum_periodogram_mean + spectrum_periodogram(1,:,mean);
    spectrum_bartlett_A_mean = spectrum_bartlett_A_mean + spectrum_bartlett_A(1,:,mean);
    spectrum_bartlett_B_mean = spectrum_bartlett_B_mean + spectrum_bartlett_B(1,:,mean);
    spectrum_welch_A_mean = spectrum_welch_A_mean + spectrum_welch_A(1,:,mean);
    spectrum_welch_B_mean = spectrum_welch_B_mean + spectrum_welch_B(1,:,mean);
    spectrum_blackman_A_mean = spectrum_blackman_A_mean + spectrum_blackman_A(1,:,mean);
    spectrum_blackman_B_mean = spectrum_blackman_B_mean + spectrum_blackman_B(1,:,mean);
end
spectrum_periodogram_mean = spectrum_periodogram_mean / Mc;
spectrum_bartlett_A_mean = spectrum_bartlett_A_mean / Mc;
spectrum_bartlett_B_mean = spectrum_bartlett_B_mean / Mc;
spectrum_welch_A_mean = spectrum_welch_A_mean / Mc;
spectrum_welch_B_mean = spectrum_welch_B_mean / Mc;
spectrum_blackman_A_mean = spectrum_blackman_A_mean / Mc;
spectrum_blackman_B_mean = spectrum_blackman_B_mean / Mc;

% Bias
bias_periodogram = spectrum_periodogram_mean(1:NW) - spectrum1;
bias_bartlett_A = spectrum_bartlett_A_mean(1:NW) - spectrum1;
bias_bartlett_B = spectrum_bartlett_B_mean(1:NW) - spectrum1;
bias_welch_A = spectrum_welch_A_mean(1:NW) - spectrum1;
bias_welch_B = spectrum_welch_B_mean(1:NW) - spectrum1;
bias_blackman_A = spectrum_blackman_A_mean(1:NW) - spectrum1;
bias_blackman_B = spectrum_blackman_B_mean(1:NW) - spectrum1;

% Estimator Variance
var_periodogram = zeros(1,NW);
var_bartlett_A = zeros(1,NW);
var_bartlett_B = zeros(1,NW);
var_welch_A = zeros(1,NW);
var_welch_B = zeros(1,NW);
var_blackman_A = zeros(1,NW);
var_blackman_B = zeros(1,NW);
for m=1:Mc
    var_periodogram = var_periodogram + abs(spectrum_periodogram(1,1:NW,m) - spectrum1).^2;
    var_bartlett_A = var_bartlett_A + abs(spectrum_bartlett_A(1,1:NW,m) - spectrum1).^2;
    var_bartlett_B = var_bartlett_B + abs(spectrum_bartlett_B(1,1:NW,m) - spectrum1).^2;
    var_welch_A = var_welch_A + abs(spectrum_welch_A(1,1:NW,m) - spectrum1).^2;
    var_welch_B = var_welch_B + abs(spectrum_welch_B(1,1:NW,m) - spectrum1).^2;
    var_blackman_A = var_blackman_A + abs(spectrum_blackman_A(1,1:NW,m) - spectrum1).^2;
    var_blackman_B = var_blackman_B + abs(spectrum_blackman_B(1,1:NW,m) - spectrum1).^2;
end

var_periodogram = var_periodogram / Mc;
var_bartlett_A = var_bartlett_A / Mc;
var_bartlett_B = var_bartlett_B / Mc;
var_welch_A = var_welch_A / Mc;
var_welch_B = var_welch_B / Mc;
var_blackman_A = var_blackman_A / Mc;
var_blackman_B = var_blackman_B / Mc;

% MSE - Mean Square Error
mse_periodogram = var_periodogram + bias_periodogram.^2;
mse_bartlett_A = var_bartlett_A + bias_bartlett_A.^2;
mse_bartlett_B = var_bartlett_B + bias_bartlett_B.^2;
mse_welch_A = var_welch_A + bias_welch_A.^2;
mse_welch_B = var_welch_B + bias_welch_B.^2;
mse_blackman_A = var_blackman_A + bias_blackman_A.^2;
mse_blackman_B = var_blackman_B + bias_blackman_B.^2;


figure(5)
subplot(2,2,1)
hold on
plot (omega, (spectrum_periodogram_mean(1:M/2+1)))
plot (omega, spectrum_bartlett_A_mean(1:M/2+1))
plot (omega, spectrum_bartlett_B_mean(1:M/2+1))
plot (omega, spectrum_welch_A_mean(1:M/2+1))
plot (omega, spectrum_welch_B_mean(1:M/2+1))
plot (omega, spectrum_blackman_A_mean(1:M/2+1))
plot (omega, spectrum_blackman_B_mean(1:M/2+1))
plot(omega, spectrum1, color = [0 0 0], LineWidth=2)
hold off
xlabel('\omega')
ylabel('S(\omega)')
legend ('Periodogram', 'bartlett A', 'bartlett B', 'Welch A', 'Welch B', 'BT A', 'BT B', 'Analytic')
title('x1 Average Graph')

figure(5)
subplot(2,2,2)
hold on
plot (omega, bias_periodogram)
plot (omega, bias_bartlett_A)
plot (omega, bias_bartlett_B)
plot (omega, bias_welch_A)
plot (omega, bias_welch_B)
plot (omega, bias_blackman_A)
plot (omega, bias_blackman_B)
hold off
xlabel('\omega')
ylabel('S(\omega)')
legend ('Periodogram', 'bartlett A', 'bartlett B', 'Welch A', 'Welch B', 'Blackman-Tukey A', 'Blackman-Tukey B')
title('x1 Bias Graph')

figure(5)
subplot(2,2,3)
hold on
plot (omega, var_periodogram)
plot (omega, var_bartlett_A)
plot (omega, var_bartlett_B)
plot (omega, var_welch_A)
plot (omega, var_welch_B)
plot (omega, var_blackman_A)
plot (omega, var_blackman_B)
hold off
xlabel('\omega')
ylabel('S(\omega)')
legend ('Periodogram', 'bartlett A', 'bartlett B', 'Welch A', 'Welch B', 'Blackman-Tukey A', 'Blackman-Tukey B')
title('x1 Variance Graph')

figure(5)
subplot(2,2,4)
hold on
plot (omega, mse_periodogram)
plot (omega, mse_bartlett_A)
plot (omega, mse_bartlett_B)
plot (omega, mse_welch_A)
plot (omega, mse_welch_B)
plot (omega, mse_blackman_A)
plot (omega, mse_blackman_B)
hold off
xlabel('\omega')
ylabel('S(\omega)')
legend ('Periodogram', 'bartlett A', 'bartlett B', 'Welch A', 'Welch B', 'Blackman-Tukey A', 'Blackman-Tukey B')
title('x1 MSE Graph')
% ////////////////////////////////////////////////////////////////////////

%****PARAMETERS****
%Global Parameters
std2 = sqrt(0.51);

%Periodogram Parameters
spectrum_periodogram_x2 = zeros(1,M,Mc);

%bartlett Parameters
spectrum_bartlett_A_x2 = zeros(1,M,Mc);
spectrum_bartlett_B_x2 = zeros(1,M,Mc);

%Welch ParametersK1_welch = 61;
spectrum_welch_A_x2 = zeros(1,M,Mc);
spectrum_welch_B_x2 = zeros(1,M,Mc);

%Blackman-Tukey Parameters
spectrum_blackman_A_x2 = zeros(1,M,Mc);
spectrum_blackman_B_x2 = zeros(1,M,Mc);

% ////////////////////////////////////////////

for monte = 1:Mc
    % New x2 construction
    w2_S2Q1 = normrnd(0,std2,[1,NQ2A]);
    x2_S2Q1=zeros(size(w2_S2Q1)); % Q22 for the Question 2.2
    x2_0_S2Q1 = normrnd(0,1);
    x2_S2Q1(1) = 0.7*x2_0_S2Q1+w2_S2Q1(1);
    
    for n = 2:length(w2_S2Q1)
        x2_S2Q1(n) = 0.7*x2_S2Q1(n-1)+w2_S2Q1(n);
    end
    x2_S2Q1 = [x2_S2Q1 zeros(1, M-length(x2_S2Q1))];
    % ////////////////////////////////////////////
    
    % Periodogram
    X2S2Q1_periodogram = fft(x2_S2Q1);
    spectrum_periodogram_x2(1,:,monte) = 1/N*abs(X2S2Q1_periodogram).^2;
    % ////////////////////////////////////////////
    
    % bartlett
    
    x2_bartlett_A = zeros(K1_bartlett, L1_bartlett);
    X2_bartlett_A = zeros(1, M, K1_bartlett);
    x2_bartlett_B = zeros(K2_bartlett, L2_bartlett);
    X2_bartlett_B = zeros(1, M, K2_bartlett);

    % Fourier transform of x2 for bartlett part A
    for k=1:K1_bartlett
        for n=1:L1_bartlett
            x2_bartlett_A(k,n) = x2_S2Q1((k-1)*L1_bartlett+n);
        end
    end

    for k=1:K1_bartlett
        X2_bartlett_A(1,:,k) = fft(x2_bartlett_A(k,:),M);
    end

    
    % bartlett part A
    for k = 1:K1_bartlett
        spectrum_bartlett_A_x2(1, :, monte) = spectrum_bartlett_A_x2(1, :, monte) + 1/L1_bartlett * abs(X2_bartlett_A(1, :, k)).^2;
    end
    spectrum_bartlett_A_x2(1, :, monte) = 1/K1_bartlett * spectrum_bartlett_A_x2(1, :, monte);


    % Fourier transform of x2 for bartlett part B
    for k=1:K2_bartlett
        for n=1:L2_bartlett
            x2_bartlett_B(k,n) = x2_S2Q1((k-1)*L2_bartlett+n);
        end
    end

    for k=1:K2_bartlett
        X2_bartlett_B(1,:,k) = fft(x2_bartlett_B(k,:),M);
    end
    
    % bartlett part B
    for k = 1:K2_bartlett
        spectrum_bartlett_B_x2(1, :, monte) = spectrum_bartlett_B_x2(1, :, monte) + 1/L2_bartlett * abs(X2_bartlett_B(1, :, k)).^2;
    end
    spectrum_bartlett_B_x2(1, :, monte) = 1/K2_bartlett * spectrum_bartlett_B_x2(1, :, monte);
    % ////////////////////////////////////////////

    % Welch
    
    % Welch part A
    x2_welch_A = zeros(K1_welch,L1_welch);
    for k = 1:K1_welch
        for n = 1:L1_welch
            x2_welch_A(k,n) = x2_S2Q1((k-1)*D1_welch+n);
        end
    end

    X2_welch_A = zeros(1,M,K1_welch);
    for k = 1:K1_welch
        X2_welch_A(1,:,k) = fft(x2_welch_A(k,:),M);
    end

    for k=1:K1_welch
        spectrum_welch_A_x2(1,:,monte) = spectrum_welch_A_x2(1,:,monte) + 1/L1_welch*abs(X2_welch_A(1,:,k)).^2;
    end

    spectrum_welch_A_x2(1,:,monte) = spectrum_welch_A_x2(1,:,monte) / K1_welch;

    % Welch part B
    x2_welch_B = zeros(K2_welch,L2_welch);
    for k = 1:K2_welch
        for n = 1:L2_welch
            x2_welch_B(k,n) = x2_S2Q1((k-1)*D2_welch+n);
        end
    end

    X2_welch_B = zeros(1,M,K2_welch);
    for k = 1:K2_welch
        X2_welch_B(1,:,k) = fft(x2_welch_B(k,:),M);
    end

    for k=1:K2_welch
        spectrum_welch_B_x2(1,:,monte) = spectrum_welch_B_x2(1,:,monte) + 1/L2_welch*abs(X2_welch_B(1,:,k)).^2;
    end

    spectrum_welch_B_x2(1,:,monte) = spectrum_welch_B_x2(1,:,monte) / K2_welch;
%     ////////////////////////////////////////////

    % Blackman-Tukey
    RX2_bt = xcorr(x2_S2Q1(1:N));
    RX2_bt = RX2_bt/N;

    % BT Part A
    RX2_bt_A = [RX2_bt(N+1:N+1+l1_bt) zeros(1, M-2*l1_bt-1) RX2_bt(N+1-l1_bt:N)];
    spectrum_blackman_A_x2(1,:,monte) = abs(fft(RX2_bt_A));
    
    % BT Part B
    RX2_bt_B = [RX2_bt(N+1:N+1+l2_bt) zeros(1, M-2*l2_bt-1) RX2_bt(N+1-l2_bt:N)];
    spectrum_blackman_B_x2(1,:,monte) = abs(fft(RX2_bt_B));
%   //////////////////////////////////////////////

end

% Average
spectrum_periodogram_mean_x2 = zeros (1,M);
spectrum_bartlett_A_mean_x2 = zeros (1,M);
spectrum_bartlett_B_mean_x2 = zeros (1,M);
spectrum_welch_A_mean_x2 = zeros (1,M);
spectrum_welch_B_mean_x2 = zeros (1,M);
spectrum_blackman_A_mean_x2 = zeros (1,M);
spectrum_blackman_B_mean_x2 = zeros (1,M);

for mean=1:Mc
    spectrum_periodogram_mean_x2 = spectrum_periodogram_mean_x2 + spectrum_periodogram_x2(1,:,mean);
    spectrum_bartlett_A_mean_x2 = spectrum_bartlett_A_mean_x2 + spectrum_bartlett_A_x2(1,:,mean);
    spectrum_bartlett_B_mean_x2 = spectrum_bartlett_B_mean_x2 + spectrum_bartlett_B_x2(1,:,mean);
    spectrum_welch_A_mean_x2 = spectrum_welch_A_mean_x2 + spectrum_welch_A_x2(1,:,mean);
    spectrum_welch_B_mean_x2 = spectrum_welch_B_mean_x2 + spectrum_welch_B_x2(1,:,mean);
    spectrum_blackman_A_mean_x2 = spectrum_blackman_A_mean_x2 + spectrum_blackman_A_x2(1,:,mean);
    spectrum_blackman_B_mean_x2 = spectrum_blackman_B_mean_x2 + spectrum_blackman_B_x2(1,:,mean);
end
spectrum_periodogram_mean_x2 = spectrum_periodogram_mean_x2 / Mc;
spectrum_bartlett_A_mean_x2 = spectrum_bartlett_A_mean_x2 / Mc;
spectrum_bartlett_B_mean_x2 = spectrum_bartlett_B_mean_x2 / Mc;
spectrum_welch_A_mean_x2 = spectrum_welch_A_mean_x2 / Mc;
spectrum_welch_B_mean_x2 = spectrum_welch_B_mean_x2 / Mc;
spectrum_blackman_A_mean_x2 = spectrum_blackman_A_mean_x2 / Mc;
spectrum_blackman_B_mean_x2 = spectrum_blackman_B_mean_x2 / Mc;

% Bias
bias_periodogram_x2 = spectrum_periodogram_mean_x2(1:NW) - spectrum2;
bias_bartlett_A_x2 = spectrum_bartlett_A_mean_x2(1:NW) - spectrum2;
bias_bartlett_B_x2 = spectrum_bartlett_B_mean_x2(1:NW) - spectrum2;
bias_welch_A_x2 = spectrum_welch_A_mean_x2(1:NW) - spectrum2;
bias_welch_B_x2 = spectrum_welch_B_mean_x2(1:NW) - spectrum2;
bias_blackman_A_x2 = spectrum_blackman_A_mean_x2(1:NW) - spectrum2;
bias_blackman_B_x2 = spectrum_blackman_B_mean_x2(1:NW) - spectrum2;

% Estimator Variance
var_periodogram_x2 = zeros(1,NW);
var_bartlett_A_x2 = zeros(1,NW);
var_bartlett_B_x2 = zeros(1,NW);
var_welch_A_x2 = zeros(1,NW);
var_welch_B_x2 = zeros(1,NW);
var_blackman_A_x2 = zeros(1,NW);
var_blackman_B_x2 = zeros(1,NW);
for m=1:Mc
    var_periodogram_x2 = var_periodogram_x2 + abs(spectrum_periodogram_x2(1,1:NW,m) - spectrum2).^2;
    var_bartlett_A_x2 = var_bartlett_A_x2 + abs(spectrum_bartlett_A_x2(1,1:NW,m) - spectrum2).^2;
    var_bartlett_B_x2 = var_bartlett_B_x2 + abs(spectrum_bartlett_B_x2(1,1:NW,m) - spectrum2).^2;
    var_welch_A_x2 = var_welch_A_x2 + abs(spectrum_welch_A_x2(1,1:NW,m) - spectrum2).^2;
    var_welch_B_x2 = var_welch_B_x2 + abs(spectrum_welch_B_x2(1,1:NW,m) - spectrum2).^2;
    var_blackman_A_x2 = var_blackman_A_x2 + abs(spectrum_blackman_A_x2(1,1:NW,m) - spectrum2).^2;
    var_blackman_B_x2 = var_blackman_B_x2 + abs(spectrum_blackman_B_x2(1,1:NW,m) - spectrum2).^2;
end

var_periodogram_x2 = var_periodogram_x2 / Mc;
var_bartlett_A_x2 = var_bartlett_A_x2 / Mc;
var_bartlett_B_x2 = var_bartlett_B_x2 / Mc;
var_welch_A_x2 = var_welch_A_x2 / Mc;
var_welch_B_x2 = var_welch_B_x2 / Mc;
var_blackman_A_x2 = var_blackman_A_x2 / Mc;
var_blackman_B_x2 = var_blackman_B_x2 / Mc;

% MSE - Mean Square Error
mse_periodogram_x2 = var_periodogram_x2 + bias_periodogram_x2.^2;
mse_bartlett_A_x2 = var_bartlett_A_x2 + bias_bartlett_A_x2.^2;
mse_bartlett_B_x2 = var_bartlett_B_x2 + bias_bartlett_B_x2.^2;
mse_welch_A_x2 = var_welch_A_x2 + bias_welch_A_x2.^2;
mse_welch_B_x2 = var_welch_B_x2 + bias_welch_B_x2.^2;
mse_blackman_A_x2 = var_blackman_A_x2 + bias_blackman_A_x2.^2;
mse_blackman_B_x2 = var_blackman_B_x2 + bias_blackman_B_x2.^2;


figure(6)
subplot(2,2,1)
hold on
plot (omega, (spectrum_periodogram_mean_x2(1:M/2+1)))
plot (omega, spectrum_bartlett_A_mean_x2(1:M/2+1))
plot (omega, spectrum_bartlett_B_mean_x2(1:M/2+1))
plot (omega, spectrum_welch_A_mean_x2(1:M/2+1))
plot (omega, spectrum_welch_B_mean_x2(1:M/2+1))
plot (omega, spectrum_blackman_A_mean_x2(1:M/2+1))
plot (omega, spectrum_blackman_B_mean_x2(1:M/2+1))
plot(omega, spectrum2, color = [0 0 0], LineWidth=2)
hold off
xlabel('\omega')
ylabel('S(\omega)')
legend ('Periodogram', 'bartlett A', 'bartlett B', 'Welch A', 'Welch B', 'BT A', 'BT B', 'Analytic')
title('x2 Average Graph')

figure(6)
subplot(2,2,2)
hold on
plot (omega, bias_periodogram_x2)
plot (omega, bias_bartlett_A_x2)
plot (omega, bias_bartlett_B_x2)
plot (omega, bias_welch_A_x2)
plot (omega, bias_welch_B_x2)
plot (omega, bias_blackman_A_x2)
plot (omega, bias_blackman_B_x2)
hold off
xlabel('\omega')
ylabel('S(\omega)')
legend ('Periodogram', 'bartlett A', 'bartlett B', 'Welch A', 'Welch B', 'Blackman-Tukey A', 'Blackman-Tukey B')
title('x2 Bias Graph')

figure(6)
subplot(2,2,3)
hold on
plot (omega, var_periodogram_x2)
plot (omega, var_bartlett_A_x2)
plot (omega, var_bartlett_B_x2)
plot (omega, var_welch_A_x2)
plot (omega, var_welch_B_x2)
plot (omega, var_blackman_A_x2)
plot (omega, var_blackman_B_x2)
hold off
xlabel('\omega')
ylabel('S(\omega)')
legend ('Periodogram', 'bartlett A', 'bartlett B', 'Welch A', 'Welch B', 'Blackman-Tukey A', 'Blackman-Tukey B')
title('x2 Variance Graph')

figure(6)
subplot(2,2,4)
hold on
plot (omega, mse_periodogram_x2)
plot (omega, mse_bartlett_A_x2)
plot (omega, mse_bartlett_B_x2)
plot (omega, mse_welch_A_x2)
plot (omega, mse_welch_B_x2)
plot (omega, mse_blackman_A_x2)
plot (omega, mse_blackman_B_x2)
hold off
xlabel('\omega')
ylabel('S(\omega)')
legend ('Periodogram', 'bartlett A', 'bartlett B', 'Welch A', 'Welch B', 'Blackman-Tukey A', 'Blackman-Tukey B')
title('x2 MSE Graph')
% ////////////////////////////////////////////////////////////////////////

% % Question 4

% signal x1
bias_periodogram_square = 1/NW * sum(bias_periodogram.^2,2);
bias_bartlett_A_square = 1/NW * sum(bias_bartlett_A.^2,2);
bias_bartlett_B_square = 1/NW * sum(bias_bartlett_B.^2,2);
bias_welch_A_square = 1/NW * sum(bias_welch_A.^2,2);
bias_welch_B_square = 1/NW * sum(bias_welch_B.^2,2);
bias_blackman_A_square = 1/NW * sum(bias_blackman_A.^2,2);
bias_blackman_B_square = 1/NW * sum(bias_blackman_B.^2,2);


var_periodogram_av = 1/NW * sum(var_periodogram, 2);
var_bartlett_A_av = 1/NW * sum(var_bartlett_A, 2);
var_bartlett_B_av = 1/NW * sum(var_bartlett_B, 2);
var_welch_A_av = 1/NW * sum(var_welch_A, 2);
var_welch_B_av = 1/NW * sum(var_welch_B, 2);
var_blackman_A_av = 1/NW * sum(var_blackman_A, 2);
var_blackman_B_av = 1/NW * sum(var_blackman_B, 2);

mse_periodogram_av = 1/NW * sum(mse_periodogram, 2);
mse_bartlett_A_av = 1/NW * sum(mse_bartlett_A, 2);
mse_bartlett_B_av = 1/NW * sum(mse_bartlett_B, 2);
mse_welch_A_av = 1/NW * sum(mse_welch_A, 2);
mse_welch_B_av = 1/NW * sum(mse_welch_B, 2);
mse_blackman_A_av = 1/NW * sum(mse_blackman_A, 2);
mse_blackman_B_av = 1/NW * sum(mse_blackman_B, 2);

% signal x2
bias_periodogram_square_x2 = 1/NW * sum(bias_periodogram_x2.^2,2);
bias_bartlett_A_square_x2 = 1/NW * sum(bias_bartlett_A_x2.^2,2);
bias_bartlett_B_square_x2 = 1/NW * sum(bias_bartlett_B_x2.^2,2);
bias_welch_A_square_x2 = 1/NW * sum(bias_welch_A_x2.^2,2);
bias_welch_B_square_x2 = 1/NW * sum(bias_welch_B_x2.^2,2);
bias_blackman_A_square_x2 = 1/NW * sum(bias_blackman_A_x2.^2,2);
bias_blackman_B_square_x2 = 1/NW * sum(bias_blackman_B_x2.^2,2);


var_periodogram_av_x2 = 1/NW * sum(var_periodogram_x2, 2);
var_bartlett_A_av_x2 = 1/NW * sum(var_bartlett_A_x2, 2);
var_bartlett_B_av_x2 = 1/NW * sum(var_bartlett_B_x2, 2);
var_welch_A_av_x2 = 1/NW * sum(var_welch_A_x2, 2);
var_welch_B_av_x2 = 1/NW * sum(var_welch_B_x2, 2);
var_blackman_A_av_x2 = 1/NW * sum(var_blackman_A_x2, 2);
var_blackman_B_av_x2 = 1/NW * sum(var_blackman_B_x2, 2);

mse_periodogram_av_x2 = 1/NW * sum(mse_periodogram_x2, 2);
mse_bartlett_A_av_x2 = 1/NW * sum(mse_bartlett_A_x2, 2);
mse_bartlett_B_av_x2 = 1/NW * sum(mse_bartlett_B_x2, 2);
mse_welch_A_av_x2 = 1/NW * sum(mse_welch_A_x2, 2);
mse_welch_B_av_x2 = 1/NW * sum(mse_welch_B_x2, 2);
mse_blackman_A_av_x2 = 1/NW * sum(mse_blackman_A_x2, 2);
mse_blackman_B_av_x2 = 1/NW * sum(mse_blackman_B_x2, 2);





