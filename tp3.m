close all; clear all; clc;

%%

M = 20;

h = [1 0.6 0.2 -0.2 0.1 0.05];
var_V = 0.005;

h_pad = zeros(1, M);
h_pad(1:min(numel(h),M)) = h(1:min(numel(h),M));

h_flip = fliplr(h_pad);

R_Y = conv(h_pad, h_flip)(M:2*M-1);
R_Y(1) += var_V;
RR_Y = toeplitz(R_Y);

R_YX = [h(1) zeros(1, M-1)];
RR_YX = R_YX';

w0 = RR_Y \ RR_YX;

%%

rand('seed', 42);

N = 100;

B = rand(1, N) < 0.5;
X = 2*B - 1;
V = normrnd(0, var_V, 1, N);
Y = filter(h, 1, X) + V;
Z = filter(w0, 1, Y);

function plot_seqRS(X)
  nfft = 2048;
  freqrange = linspace(-1, 1, nfft);
  freqticks = {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'};

  N = numel(X);
  nn = 0:N-1;

  figure();
  subplot(2, 2, [1 2]);
  stairs(nn, X);
  grid;
  xlim([0 N-1]);
  ylim([-2.5 2.5]);
  xlabel('Tiempo');
  title(['Secuencia ' inputname(1)]);
  set(findall(gcf,'type','line'),'LineWidth',2);
  set(gca, 'FontSize', 16);

  [R_Xest kk] = xcorr(X);
  S_Xest = 1/N * abs(fftshift(fft(X, nfft))).^2;

  subplot(2, 2, 3);
  stem(kk, R_Xest);
  grid;
  xlabel('Retraso');
  ylim([-50 150]);
  title(['Autocorrelación de ' inputname(1)]);
  set(findall(gcf,'type','line'),'LineWidth',2);
  set(gca, 'FontSize', 16);

  subplot(2, 2, 4);
  plot(freqrange, mag2db(S_Xest));
  grid;
  xticklabels(freqticks);
  ylim([-60 40]);
  xlabel('Frecuencia normalizada');
  ylabel('Densidad espectral de potencia [dB]');
  title(['PSD de ' inputname(1)]);
  set(findall(gcf,'type','line'),'LineWidth',2);
  set(gca, 'FontSize', 16);
endfunction

plot_seqRS(X)
plot_seqRS(Y)
plot_seqRS(Z)

mm = 0:M-1;

figure();
stem(mm',[h_pad' w0]);
grid;
xlim([0 M-1]);
ylim([-0.7 1.1]);
xlabel('Tiempo');
title('Coeficientes de h(n) y w(n)');
legend('h(n)', 'w(n)');
set(findall(gcf,'type','line'),'LineWidth',2);
set(gca, 'FontSize', 16);

%%

for M = [2 5 10 30]
  h_pad = zeros(1, M);
  h_pad(1:min(numel(h),M)) = h(1:min(numel(h),M));

  h_flip = fliplr(h_pad);

  R_Y = conv(h_pad, h_flip)(M:2*M-1);
  R_Y(1) += var_V;
  RR_Y = toeplitz(R_Y);

  R_YX = [h(1) zeros(1, M-1)];
  RR_YX = R_YX';

  w0 = RR_Y \ RR_YX;

  Z = filter(w0, 1, Y);
  plot_seqRS(Z);

  mm = 0:M-1;

  figure();
  subplot(2, 2, 1);
  stem(mm,w0);
  grid;
  xlim([0 M-1]);
  ylim([-0.7 1.1]);
  xlabel('Tiempo');
  title('Coeficientes de w(n)');
  set(findall(gcf,'type','line'),'LineWidth',2);
  set(gca, 'FontSize', 16);

  hconvw = conv(h, w0);

  subplot(2, 2, 2);
  stem(0:numel(hconvw)-1, hconvw);
  grid;
  xlim([0 M-1]);
  ylim([-0.7 1.1]);
  xlabel('Tiempo');
  title('Convolución entre h(n) y w(n)');
  set(findall(gcf,'type','line'),'LineWidth',2);
  set(gca, 'FontSize', 16);

  freqrange = linspace(-1, 1, 2048);
  freqticks = {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'};

  H = freqz(h, 1, freqrange);
  W = freqz(w0, 1, freqrange);

  subplot(2, 2, [3 4]);
  plot(freqrange, (abs(H)));
  hold on;
  plot(freqrange, (abs(W)));
  plot(freqrange, (abs(H).*abs(W)));
  grid;
  xticklabels(freqticks);
  xlabel('Frecuencia normalizada');
  ylabel('Magnitud [dB]');
  title('Respuesta en frecuencia de |H(\omega)|, |W(\omega)| y |H(\omega)|\cdot|W(\omega)|');
  legend('|H(\omega)|', '|W(\omega)|', '|H(\omega)|\cdot|W(\omega)|');
  set(findall(gcf,'type','line'),'LineWidth',2);
  set(gca, 'FontSize', 16);
endfor
