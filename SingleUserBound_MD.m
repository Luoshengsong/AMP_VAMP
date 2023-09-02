function MD_singleUser = SingleUserBound_MD(SNR_dB, Monte)
% always active (u = 1) -> check MD
MD_singleUser = zeros(1,length(SNR_dB));
for snrIdx = 1: length (SNR_dB)
    SNRdB = SNR_dB(snrIdx);
    sigma2 = 1 ./ 10^(SNRdB/10);
    
    cnt_error = 0;
    for mt = 1: Monte
        uh = sqrt(0.5 * 1) * (randn(1,1) + 1j*randn(1,1));
        r = uh + sqrt(0.5 * sigma2) * (randn(1,1) + 1j*randn(1,1));
        if abs(r)^2 < sigma2 * (1+ sigma2) * log( 1 + 1/sigma2 )
            error = 1;
            cnt_error = cnt_error +1;
        else
            error = 0;
        end
    end
    MD_singleUser(snrIdx) = cnt_error / Monte;
end
end