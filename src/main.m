src/main.m

%% section 1
clear all; close all; clc;

% Parameters
w0 = 2*pi*250;     
f = -300:0.1:300;  
jw2pif = 1j*2*pi*f;

% Fourier Transform 
Xf = (1/(2j))*( 1./(-2 + 1j*(w0 + 2*pi*f)) - 1./(-2 + 1j*(-w0 + 2*pi*f)) );

% Plot
figure(1);
plot(f, abs(Xf), 'LineWidth', 1.5);
xlabel('f [Hz]');
ylabel('|X^{F}_{\omega_0}(f)|');
title('|X_{ω₀}ᶠ(f)| for ω₀ = 2π · 250 [rad/sec]');
grid on;

%% Section 2

Fs = 16000;           
T_total = 4;          % 4 secs 
N = Fs * T_total;     

t = (0:N-1) / Fs;     % time vector

% angular frequencies
w1 = 2*pi*250;        % [rad/sec]
w2 = 2*pi*315;
w3 = 2*pi*375;
w4 = 2*pi*500;

% x_w(t) = sin(wt) * exp(-2t) * u(t)
xw = @(w, tt) sin(w.*tt) .* exp(-2*tt) .* (tt >= 0);

% allocate vector for x(t)
x = zeros(1, N);

% indices for each 1 second interval
idx1 = (t >= 0   & t < 1);
idx2 = (t >= 1   & t < 2);
idx3 = (t >= 2   & t < 3);
idx4 = (t >= 3   & t < 4);

% building x(t)
x(idx1) = xw(w1, t(idx1));          % 0 <= t < 1
x(idx2) = xw(w2, t(idx2) - 1);      % 1 <= t < 2
x(idx3) = xw(w3, t(idx3) - 2);      % 2 <= t < 3
x(idx4) = xw(w4, t(idx4) - 3);      % 3 <= t < 4


soundsc(x, Fs);

% plot
figure(2);
plot(t, x);
xlabel('t [sec]');
ylabel('x(t)');
title('Sampled signal x(t) at 16kHz');
grid on;

%% Section 3

% the time interval
t_min = 2.99;
t_max = 3.01;


idx = find(t >= t_min & t <= t_max);

% Extract time and signal in the zoom range
t_zoom = t(idx);
x_zoom = x(idx);


figure(4);
plot(t_zoom, x_zoom, ':', 'LineWidth', 1.5);
xlabel('t [sec]');
ylabel('x(t)');
title('Zoomed signal around t = 3 sec');
grid on;
hold on;   


figure(5);
plot(t_zoom, x_zoom, ':', 'LineWidth', 1.5);
xlabel('t [sec]');
ylabel('x(t)');
title('Zoomed signal around t = 3 sec');
grid on;
hold on;  



%% Section 4 

Fs_high = 16000;      
Fs_low  = 2000;       % new sampling rate
M = Fs_high / Fs_low; 

x_n = downsample(x, M);   % length 8000

n   = 0:length(x_n)-1;    % 0 ... 7999
t_n = n / Fs_low;         % time in seconds for each sample

% zoom interval in time
t_min = 2.99;
t_max = 3.01;

% find indices of samples in the desired time range
idx_zoom_n = find(t_n >= t_min & t_n <= t_max);

t_zoom_n = t_n(idx_zoom_n);
x_zoom_n = x_n(idx_zoom_n);

figure(4);
stem(t_zoom_n, x_zoom_n, 'filled');
xlabel('t [sec]');
ylabel('x[n]');
title('Discrete-time samples x[n] around t = 3 sec (Fs = 2 kHz)');
grid on;
hold on;   

%% Section 5 
Fs = 2000;               % sampling frequency [Hz]
N  = length(x_n);        % N = 8000

% DFT
Xk = fft(x_n, N);

% shift zero-frequency to the center
Xk_shift = fftshift(Xk);

% frequency axis: -Fs/2 ... Fs/2 - Fs/N
f = (-N/2 : N/2-1) * (Fs/N);

figure(2);
plot(f, abs(Xk_shift), 'LineWidth', 1.5);
xlabel('f [Hz]');
ylabel('|X[k]|');
title('Magnitude of DFT of x[n]');
grid on;

%% Section 7 

Fs_low  = 2000;      
Fs_high = 16000;     
L = Fs_high / Fs_low;    


x_u = upsample(x_n, L);      % length 8*8000 = 64000

% sinc interpolation kernel (512 taps: -256 ... 255)
n_h   = -256:255;
h_sinc = sinc(n_h / L);      % MATLAB: sinc(x) = sin(pi x)/(pi x)


x_ideal = conv(x_u, h_sinc, 'same');    % length = 64000

% time axis
t_ideal = (0:length(x_ideal)-1) / Fs_high;

% zoom interval
t_min = 2.99;
t_max = 3.01;

idx_zoom_ideal = find(t_ideal >= t_min & t_ideal <= t_max);


plot(t_ideal(idx_zoom_ideal), x_ideal(idx_zoom_ideal), 'LineWidth', 1.5);
xlabel('t [sec]');
ylabel('x(t)');
title('x(t), samples x[n] and reconstructed x_{ideal}(t) around t = 3 sec');
grid on;
hold on;


%% Section 8 

Fs_low  = 2000;      
Fs_high = 16000;     
L = Fs_high / Fs_low;    


x_u = upsample(x_n, L);      % length 64000

% ZOH interpolation kernel: rectangle pulse
h_zoh = ones(1, L);          

% 3. convolvulotion
x_zoh = conv(x_u, h_zoh, 'same');

% time axis
t_zoh = (0:length(x_zoh)-1) / Fs_high;

% 4. zoom to 2.99 <= t <= 3.01 
t_min = 2.99;
t_max = 3.01;
idx_zoom_zoh = find(t_zoh >= t_min & t_zoh <= t_max);

figure(4);
plot(t_zoh(idx_zoom_zoh), x_zoh(idx_zoom_zoh), 'LineWidth', 1.5);
xlabel('t [sec]');
ylabel('x(t)');
title('x(t), ideal reconstruction and ZOH reconstruction around t = 3 sec');
grid on;
hold on;  
legend('x(t)','samples x[n]','x_{ideal}(t)','x_{zoh}(t)','Location','best');


%% Section 9 
Fs_low  = 2000;          
Fs_high = 16000;         
L = Fs_high / Fs_low;   

x_u = upsample(x_n, L);      

% FOH interpolation kernel discrete triangle
n_h    = -L:L;
h_foh  = max(1 - abs(n_h)/L, 0);   % triangular kernel

% 3. convolvulotion
x_foh = conv(x_u, h_foh, 'same');

% 4. time axis
t_foh = (0:length(x_foh)-1) / Fs_high;

% 5. zoom to 2.99 <= t <= 3.01
t_min = 2.99;
t_max = 3.01;
idx_zoom_foh = find(t_foh >= t_min & t_foh <= t_max);

figure(4);
plot(t_foh(idx_zoom_foh), x_foh(idx_zoom_foh), 'LineWidth', 1.5);
xlabel('t [sec]');
ylabel('x(t)');
title('x(t) and reconstructed signals (ideal, ZOH, FOH) around t = 3 sec');
grid on;
hold on;


%% Section 11 

Fs_high = 16000;         
Fs_low2 = 800;            % new sampling rate
M2 = Fs_high / Fs_low2;   

x_n2 = downsample(x, M2);     % length 4*800 = 3200
n2   = 0:length(x_n2)-1;
t_n2 = n2 / Fs_low2;          % time axis of samples



% zoom interval
t_min = 2.99;
t_max = 3.01;

idx_zoom_n2 = find(t_n2 >= t_min & t_n2 <= t_max);

figure(5);
stem(t_n2(idx_zoom_n2), x_n2(idx_zoom_n2), 'filled');
xlabel('t [sec]');
ylabel('x[n]  (Fs = 800 Hz)');
title('Samples of x(t) around t = 3 sec  (Fs = 800 Hz)');
grid on;
hold on;



N2  = length(x_n2);          % 3200 samples
Xk2 = fft(x_n2, N2);
Xk2_shift = fftshift(Xk2);

f2 = (-N2/2 : N2/2-1) * (Fs_low2/N2);   % freq axis: -Fs/2...Fs/2

figure(3);
plot(f2, abs(Xk2_shift), 'LineWidth', 1.5);
xlabel('f [Hz]');
ylabel('|X_2[k]|');
title('Magnitude of DFT of x[n]  (Fs = 800 Hz)');
grid on;



L2 = Fs_high / Fs_low2;      

x_u2 = upsample(x_n2, L2);   

% sinc kernel length 512
n_h2    = -256:255;
h_sinc2 = sinc(n_h2 / L2);

x_ideal2 = conv(x_u2, h_sinc2, 'same');
t_ideal2 = (0:length(x_ideal2)-1) / Fs_high;


idx_zoom_ideal2 = find(t_ideal2 >= t_min & t_ideal2 <= t_max);

figure(5);
plot(t_ideal2(idx_zoom_ideal2), x_ideal2(idx_zoom_ideal2), 'LineWidth', 1.5);



h_zoh2 = ones(1, L2);          % rectangular kernel
x_zoh2 = conv(x_u2, h_zoh2, 'same');
t_zoh2 = (0:length(x_zoh2)-1) / Fs_high;

idx_zoom_zoh2 = find(t_zoh2 >= t_min & t_zoh2 <= t_max);

figure(5);
plot(t_zoh2(idx_zoom_zoh2), x_zoh2(idx_zoom_zoh2), 'LineWidth', 1.5);



n_hfoh2 = -L2:L2;
h_foh2  = max(1 - abs(n_hfoh2)/L2, 0);  % triangular kernel

x_foh2 = conv(x_u2, h_foh2, 'same');
t_foh2 = (0:length(x_foh2)-1) / Fs_high;

idx_zoom_foh2 = find(t_foh2 >= t_min & t_foh2 <= t_max);

figure(5);
plot(t_foh2(idx_zoom_foh2), x_foh2(idx_zoom_foh2), 'LineWidth', 1.5);

xlabel('t [sec]');
ylabel('amplitude');
title('Reconstruction of x(t) from Fs = 800 Hz around t = 3 sec');
grid on;
hold on;

legend('samples x[n]','x_{ideal,800}(t)','x_{zoh,800}(t)','x_{foh,800}(t)','Location','best');



%% All the graphs of FT amplitudes 

w0  = 2*pi*250;
f1  = -300:0.1:300;
Xf1 = (1/(2j))*( 1./(-2 + 1j*(w0 + 2*pi*f1)) ...
               - 1./(-2 + 1j*(-w0 + 2*pi*f1)) );

% 2) DFT for Fs = 2 kHz
Fs1 = 2000;
N1  = length(x_n);
f_dft1 = (-N1/2 : N1/2-1) * (Fs1/N1);
Xk1    = fftshift( fft(x_n, N1) );

% 3) DFT for Fs = 800 Hz
Fs2 = 800;
N2  = length(x_n2);
f_dft2 = (-N2/2 : N2/2-1) * (Fs2/N2);
Xk2s   = fftshift( fft(x_n2, N2) );

% Plot all three in one figure
figure(100); clf;

subplot(3,1,1);
plot(f1, abs(Xf1), 'LineWidth', 1.5);
xlabel('f [Hz]');
ylabel('|X^{F}_{\omega_0}(f)|');
title('Analytic |X_{ω_0}^F(f)|');
grid on;

subplot(3,1,2);
plot(f_dft1, abs(Xk1), 'LineWidth', 1.5);
xlabel('f [Hz]');
ylabel('|X[k]|');
title('DFT of x[n], F_s = 2 kHz');
grid on;

subplot(3,1,3);
plot(f_dft2, abs(Xk2s), 'LineWidth', 1.5);
xlabel('f [Hz]');
ylabel('|X_2[k]|');
title('DFT of x_2[n], F_s = 800 Hz');
grid on;






%% Summary figure – time zoom + reconstructions (Fs = 2 kHz and 800 Hz)

figure(101); clf;

% ---------- Top subplot: Fs = 2 kHz ----------
subplot(2,1,1);
plot(t_zoom, x_zoom, ':', 'LineWidth', 1.5);       % original x(t) zoomed
hold on;
stem(t_zoom_n, x_zoom_n, 'filled');                % samples x[n]
plot(t_ideal(idx_zoom_ideal), x_ideal(idx_zoom_ideal), 'LineWidth', 1.5);
plot(t_zoh(idx_zoom_zoh),     x_zoh(idx_zoom_zoh), 'LineWidth', 1.5);
plot(t_foh(idx_zoom_foh),     x_foh(idx_zoom_foh), 'LineWidth', 1.5);

xlabel('t [sec]');
ylabel('x(t)');
title('F_s = 2 kHz – x(t), samples and reconstructions');
grid on;
legend('initial','samples x[n]','ideal','ZOH','FOH','Location','best');

% ---------- Bottom subplot: Fs = 800 Hz ----------
subplot(2,1,2);
plot(t_zoom, x_zoom, ':', 'LineWidth', 1.5);       % same original x(t) zoom
hold on;
stem(t_n2(idx_zoom_n2), x_n2(idx_zoom_n2), 'filled');  % samples at 800 Hz
plot(t_ideal2(idx_zoom_ideal2), x_ideal2(idx_zoom_ideal2), 'LineWidth', 1.5);
plot(t_zoh2(idx_zoom_zoh2),     x_zoh2(idx_zoom_zoh2), 'LineWidth', 1.5);
plot(t_foh2(idx_zoom_foh2),     x_foh2(idx_zoom_foh2), 'LineWidth', 1.5);

xlabel('t [sec]');
ylabel('x(t)');
title('F_s = 800 Hz – x(t), samples and reconstructions');
grid on;
legend('initial','samples x_2[n]','ideal 800','ZOH 800','FOH 800','Location','best');
