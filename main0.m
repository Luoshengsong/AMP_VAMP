clc;clear all;

addpath('./AMP_Yu');
addpath('./VAMP_Ph');
addpath('./AMP_Ph');
addpath('./main')
addpath('./stateEvo')

L = 100;
N = 2000;
M = 20;
K = 100;
coverage = 'Rayleigh';

SNR_dB = 80;
sigma2 = 1 ./ 10^(SNR_dB/10);

% power_Tx_dBm = 0.0003;
% noise_psd_dBm = -169; % / Hz
% bandWidth = 1e6;
% power = 10^(power_Tx_dBm/10)*10^(-3);
% noise_power = 10^(noise_psd_dBm/10)*10^(-3);
% noise = noise_power * bandWidth;
% sigma2 = noise/power/L;

path_loss = channel_beta(N, coverage);

% unmatch_dB = -unifrnd(0,10,[N,1]); 
% unmatch = 10.^(unmatch_dB / 10);
% path_loss = path_loss .* unmatch;

supp = randperm(N);
h=zeros(N,M);
for n=1:N
    h(n,:) = sqrt(1/2)*(randn(1,M) + sqrt(-1)*randn(1,M))*sqrt(path_loss(n));
end
x = zeros(N,M);
for m=1:M
    x(supp(1:K),:) = h(supp(1:K),:);
end

A = sqrt(0.5 / L) * (randn(L,N) + 1j*randn(L,N));
w = sqrt(0.5 * sigma2) * (randn(L,M) + 1j*randn(L,M));

y = A * x + w;

[A1,Ah1,vampEstFin]  = VAMP_pro(y, A, x, K/N, path_loss, sigma2);
vampNMSEdB_ = vampEstFin.err; 
vampNit = vampEstFin.nit;
vampNMSE = 10.^(vampNMSEdB_ / 10);

X_2 = diag( sum(abs(x).^2,1) );
vampNMSE = sum( X_2 * vampNMSE , 1) / trace(X_2);

Afro2 = N;
Aamp = FxnhandleLinTrans(L,N,A1,Ah1,Afro2/(M*N));
clear optAMP;
optAMP = AmpOpt();
tstart = tic;
[~,optFin,ampEstHist] = ampEst_(y, A, optAMP, M, path_loss, K/N);

%[~,optAMPfin,ampEstHist] = ampEst(EstimIn,y,Aamp,optAMP);
time_amp = toc(tstart);
ampNit = length(ampEstHist.it);
ampNMSEdB_ = nan(M,ampNit);
for m=1:M
ampNMSEdB_(m,:) = 10*log10(sum(abs(ampEstHist.xhat((m-1)*N+[1:N],:)-x(:,m)*ones(1,ampNit)).^2,1)/norm(x(:,m))^2);
end
ampNMSE = 10.^(ampNMSEdB_ / 10);
ampNMSE = sum( X_2 * ampNMSE , 1) / trace(X_2);


[xnoise,xhat,mse,tau_real,tau_est] = noisyCAMPmmseforKLS(A,N,M,L,y,x,50,K/N,path_loss,sigma2);
mse = mse * (N*M) / norm(x,'fro')^2;

figure(2)
handy = semilogy(1:vampNit,vampNMSE,'b.-');
set(handy(1),'Displayname',['-VAMP']);
hold on;
handy = [handy, semilogy(1:ampNit,ampNMSE,'r-o')]; 
set(handy(1,end),'Displayname',['-AMP_Ph']);
hold on;
handy = [handy, semilogy(1:50,mse,'k--')]; 
set(handy(1,end),'Displayname',['-AMP_Yu']);
ylabel('NMSE')
xlabel('iterations')
grid on


%handy = [handy, semilogx(1;50,10.*log10(tau_est),'k--')]; 
%figure(2);semilogy(1:50,10.*log10(tau_est),'b-o');xlabel('No.iter');ylabel('\tau');

