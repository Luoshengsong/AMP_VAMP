

clc; clear all;

addpath('./AMP_Yu');
addpath('./VAMP_Ph');
addpath('./AMP_Ph');
addpath('./main')
addpath('./stateEvo')

load("A_CS.mat");

M = 1;

epsilon = 100/256;

SNR_dB = 20;
sigma2 = 1 ./ 10^(SNR_dB/10);

% A = A_CS;
A = sqrt(0.5/60) * (randn(60,256) + 1i*randn(60,256));

[U, S, V] = svd(A);

S = max(S, 0.01);
A = U * S * V';
D = sum(abs(A).^2, 1);



[L, N] = size(A);

K = N * epsilon;

path_loss = 1e-14 + (1e-6 - 1e-14) * rand(N,1);

supp = randperm(N);
h=zeros(N,M);
for n=1:N
    h(n,:) = sqrt(1/2)*(randn(1,M) + sqrt(-1)*randn(1,M))* path_loss(n);
end

x = zeros(N,M);
for m=1:M
    x(supp(1:K),:) = h(supp(1:K),:);
end

w = sqrt(0.5 * sigma2) * (randn(L,M) + 1j*randn(L,M));

y = A * x + w;

y_ = U' * y;

[xnoise,xhat,mse,tau_real,tau_est] = noisyCAMPmmseforKLS(A, N, M, L, y, x, 50, K/N, path_loss, sigma2);

figure; semilogy(1:length(mse), tau_est, 'b-s'); grid on;