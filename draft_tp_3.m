

clear all;
close all;
clc

h = [1, 0.6, 0.2, -0.2, 0.1, 0.05];
sigma_v_2 = 0.005;
M = 20;

% Hay que calcular los coeficientes del equalizador.
% En primer lugar tengo que hallar R_Y que se forma a partir de la
% autocorrelacion de Y(n). Entonces hay que obtener esta autocorrelación
% primero. Se sabía que era la convolución de h con su flipped + sigma_v_2
% en k = 0.

L = length(h);

% Invierto h
h_flip = fliplr(h);
figure;
stem(h_flip);
% Convoluciono h con h moño
R_y = conv(h, h_flip);
% me quedo con la parte de R_y que me interesa, osea de 0 a 2L-1 que seria 
Ry_k = R_y(L:2*L-1);


if length(Ry_k) < M
    Ry_k = [Ry_k, zeros(1, M - length(Ry_k))];
end


% le sumo la varianza de v
Ry_k(1) = Ry_k(1) + sigma_v_2;



R_Y = toeplitz(Ry_k);

R_YX = [zeros(1, M - length(h)), h_flip];
figure;
stem(R_YX);
R_YX = fliplr(R_YX);
w_opt = inv(R_Y)*R_YX';

N = 100;

B = randi([0, 1], 1, N);
X = 2*B - 1;
[Rx, k] = xcorr(X, 'biased');
Nfft = 2048;
Sx = fftshift(fft(Rx, Nfft));
w = linspace(-pi, pi, Nfft);

figure;
subplot(3, 1, 1);
stairs(X);
subplot(3,1,2);
plot(k, Rx);
subplot(3,1,3);
plot(w, abs(Sx));

Y_1 = filter(h,1,X);
V = normrnd(0, 0.005, size(Y_1));
Y = Y_1 + V;
[Ry, k_y] = xcorr(Y, 'biased');
Sy = fftshift(fft(Ry, Nfft));

figure;
subplot(3, 1, 1);
stairs(Y);
hold on;
stairs(X);
grid on;
subplot(3,1,2);
plot(k_y, Ry);
subplot(3,1,3);
plot(w, abs(Sy));

Z = filter(w_opt, 1, Y);

[Rz, k_z] = xcorr(Z, 'biased');
Sz = fftshift(fft(Rz, Nfft));

figure;
subplot(3, 1, 1);
stairs(X);
hold on;
stairs(Y);
stairs(Z);
legend('x', 'y', 'z');
subplot(3,1,2);
plot(k_z, Rz);
subplot(3,1,3);
plot(w, abs(Sz));

figure;
stem(h);
hold on;
stem(w_opt);
grid on;


% Ejercicio 4

EMES = [2, 5,  10, 30];
H = fftshift(fft(h, Nfft));




for M=EMES
    if length(Ry_k) < M
        Ry_k = [Ry_k, zeros(1, M - length(Ry_k))];
    elseif length(Ry_k) > M
        Ry_k = Ry_k(1:M);
    end


    le sumo la varianza de v
    Ry_k(1) = Ry_k(1) + sigma_v_2;

    R_Y = toeplitz(Ry_k);

    if M > length(h_flip)
        R_YX = [zeros(1, M - length(h)), h_flip];
    elseif M < length(h_flip)
        R_YX = h_flip(length(h_flip)-M+1:length(h_flip))
    end
    R_YX = fliplr(R_YX);
    w_opt = inv(R_Y)*R_YX';

    W = fftshift(fft(w, Nfft));

    figure;
    subplot(3, 1, 1);
    stem(w_opt);
    grid on;
    subplot(3,1,2);
    stem(conv(h, w_opt));
    grid on;
    subplot(3,1,3);
    plot(w, abs(H));
    hold on;
    plot(w, abs(W));
    plot(w, abs(H.*W));
    legend('H', 'W', 'H*W');

    grid on;

end











