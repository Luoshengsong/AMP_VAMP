clc;clear all;
warning('off')

poolobj = gcp('nocreate'); delete(poolobj);  % > Parallel operation related
%rng('default');

% handle random seed
if verLessThan('matlab','7.14')
  defaultStream = RandStream.getDefaultStream;
else
  defaultStream = RandStream.getGlobalStream;
end
% if 0 % new RANDOM trial
  savedState = defaultStream.State;
  save random_state.mat savedState;
% else % repeat last trial
%   load random_state.mat
% end
defaultStream.State = savedState;

addpath('./AMP_Yu');
addpath('./VAMP_Ph');
addpath('./AMP_Ph');
addpath('./main');
addpath('./stateEvo');
addpath('./SensingMat');

L = 31;
Nset = 1024; %[128, 1024, 2^13, 2^14];
M = 1;
K = 8;
coverage = 'Rayleigh';  % ['area','linear', 'Rayleigh']
MonteCarlo = 250;
SNR_dB = [5,10,15,20,25];
isCmplx = true;
damp = 0.8;
Potential_Q = 2 * K;  % 4 -> find 8 potential users to check; 8 -> 16

NMSE_MMSE = nan(length(Nset), length(SNR_dB));
NMSE_VAMP = nan(length(Nset), length(SNR_dB));
NMSE_AMP = nan(length(Nset), length(SNR_dB));

P_MD_VAMP = nan(length(Nset), length(SNR_dB));
P_MD_AMP = nan(1, length(SNR_dB));

p = parpool(4);     % > Parallel operation related

for nIdx = 1: length(Nset)
    N = Nset(nIdx);
    % fix A
    %% generate linear transform
    if N == 2^10
        A=struct2array(load('N=1024.mat'));
    elseif N == 2^7
        A=struct2array(load('N=128.mat'));
    elseif N == 2^15
        A=struct2array(load('N=2^15.mat'));
    elseif N == 2^11
        A=struct2array(load('N=2^11.mat'));
    elseif N == 2^12
        A=struct2array(load('N=2^12.mat'));
    elseif N == 2^13
        A=struct2array(load('N=2^13.mat'));
    elseif N == 2^14
        A=struct2array(load('N=2^14.mat'));
    end

    A=normalize(A);
    %A = sqrt(0.5 / L) * (randn(L,N) + 1j*randn(L,N));
    % VAMP-preprocess by Ph
    linearStage = 'exact';%'lsqr'; % VAMP linear stage in {'exact','cg','lsqr','gd'} [lsqr]
    [A_f,vampOpt] = VAMP_Prepro(A,damp, K/N,  linearStage, isCmplx);
    % whether consider pathloss: is isPathloss = True
    if strcmp(coverage,'Rayleigh') 
        isPathloss = false; 
        xmean1 = 0;  % prior mean of nonzero coefficient
        xvar1 = 1;   % prior variance of nonzero coefficient
        if isCmplx
            EstimIn = SparseScaEstim(CAwgnEstimIn(xmean1,xvar1),K/N);  %  K/N: probability of  non-zero coef
          else
            EstimIn = SparseScaEstim(AwgnEstimIn(xmean1,xvar1),K/N);
        end
    else
        isPathloss = true;  
        EstimIn = [];
    end
   

    for snrIdx = 1: length (SNR_dB)
        SNRdB = SNR_dB(snrIdx);
        sigma2 = 1 ./ 10^(SNRdB/10);
        MSE_MMSE = nan(MonteCarlo,1); 
        MSE_VAMP = nan(MonteCarlo,1);  MSE_AMP = nan(MonteCarlo,1);
        MD_VAMP = nan(MonteCarlo,1);   MD_AMP = nan(MonteCarlo,1);
        %MD_singleUser = nan(MonteCarlo,1);

        tic,
        display(['No. user = ', num2str(N), ', where current SNR = ', num2str(SNRdB), 'dB']);
        parfor mt = 1: MonteCarlo  % > Parallel operation related: parfor
            pathloss = channel_beta(N, coverage);
            supp = randperm(N);
            Active_Set = supp(1:K);
            A_active = A(:,Active_Set);
            x = zeros(N,M);
            x(Active_Set,:) = diag(sqrt(pathloss(Active_Set))) *  sqrt(1/2)*(randn(K,M) + sqrt(-1)*randn(K,M));
            x_active = x(Active_Set,:);
            w = sqrt(0.5 * sigma2) * (randn(L,M) + 1j*randn(L,M));
            y = A * x + w;

            MSE_MMSE(mt) = LMMSE_pro(A_active, x_active, eye(K), y, sigma2);

            % VAMP processing 
            [vampEstFin] = VAMP_pro1(linearStage, isPathloss, N, pathloss, y,x, sigma2, A_f, vampOpt, EstimIn);
            x_hat_vamp = vampEstFin.x2;  % final estimate 
            vampNit = vampEstFin.nit;
            MSE_VAMP(mt) = norm(x_hat_vamp-x, 'fro')^2;
            MD_VAMP(mt) =  checkActive(x_hat_vamp, Active_Set, Potential_Q);

            % AMP processing by Yu
            [xnoise,xhat_amp,mse,tau_real,tau_est] = noisyCAMPmmseforKLS(A,N,M,L,y,x,50,K/N,pathloss,sigma2);
            MSE_AMP(mt) = mse(end) * (N*M) ;
            MD_AMP(mt) =  checkActive(xhat_amp, Active_Set, Potential_Q);
        end
        toc,
        NMSE_MMSE(nIdx, snrIdx) = sum(MSE_MMSE) / MonteCarlo;
        NMSE_VAMP(nIdx, snrIdx) = sum(MSE_VAMP) / K / MonteCarlo;
        NMSE_AMP(nIdx, snrIdx) = sum(MSE_AMP) / K / MonteCarlo;

        P_MD_VAMP(nIdx, snrIdx) = sum(MD_VAMP) / K / MonteCarlo;
        P_MD_AMP(nIdx, snrIdx) = sum(MD_AMP) / K / MonteCarlo;
    end
end
delete(p);   % > Parallel operation related

colors = ['b','m', 'g', 'r', 'c', 'y'];
figure
handy = [semilogy(SNR_dB, NMSE_MMSE(1,:), join(['k','-s']) )];
set(handy(1,end),'Displayname',join(['MMSE:K=',num2str(K)]) );
hold on;
for nn = 1: length(Nset)
    handy = [handy, semilogy(SNR_dB, NMSE_VAMP(nn,:), join([colors(nn),'-o']) )];
    set(handy(1,end),'Displayname',join(['VAMP:N=2^{',num2str(log2(Nset(nn))),'}']) );
    hold on;

    handy = [handy, semilogy(SNR_dB, NMSE_AMP(nn,:), join([colors(nn),'-.^']))]; 
    set(handy(1,end),'Displayname',join([' AMP:N=2^{',num2str(log2(Nset(nn))),'}']) );
    if nn < length(Nset)
        hold on;
    end
end

ylabel('NMSE')
xlabel('SNR(dB)')
title(join(["NMSE(K=",num2str(K), ' active users)']));
set(gca,'yscale','log')
grid on
legend(handy(1,:)); 


% single usr bound
P_MD_singleUser = SingleUserBound_MD(SNR_dB, 1e4);

figure
 handy = [semilogy(SNR_dB, P_MD_singleUser, join(['k','-*']) )];
 set(handy(1,end),'Displayname','SingleUserBound' );
hold on;
for nn = 1: length(Nset)
    handy = [ handy, semilogy(SNR_dB, P_MD_VAMP(nn,:), join([colors(nn),'-o']) )];
    set(handy(1,end),'Displayname',join(['VAMP:N=2^{',num2str(log2(Nset(nn))),'}']) );
    hold on;

    handy = [handy, semilogy(SNR_dB, P_MD_AMP(nn,:), join([colors(nn),'-.v']))]; 
    set(handy(1,end),'Displayname',join(['AMP:N=2^{',num2str(log2(Nset(nn))),'}']) );
    if nn < length(Nset)
        hold on;
    end
end

ylabel('Miss Detection Probability')
xlabel('SNR(dB)')
title(join(["Miss Detection Probability(K=",num2str(K), ' active users)']));
set(gca,'yscale','log')
grid on
legend(handy(1,:)); 



