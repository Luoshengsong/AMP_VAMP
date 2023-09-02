function beta = channel_beta(N, coverage)

switch coverage
    case 'area'
        x = 2000 * rand(1e6,1) - 1000;
        y = 2000 * rand(1e6,1) - 1000;
        rx = x( ((x.^2 + y.^2) > 50^2) & ((x.^2 + y.^2) < 1000^2));
        ry = y( ((x.^2 + y.^2) > 50^2) & ((x.^2 + y.^2) < 1000^2));
        distance = sqrt(rx.^2 + ry.^2);
        if length(distance) < N
            error('N is too large for this setting!');
        else
            distance = distance(1:N);
            beta_dB = -128.1 - 36.7 * log10(distance./1000);
            beta = 10 .^(beta_dB / 10);
        end
    case 'linear'
        distance = 0.05 + (1-0.05) * rand(N,1);
        beta_dB = -128.1 - 36.7 * log10(distance);
        beta = db2mag(beta_dB).^2;
    case 'Rayleigh'
        beta = ones(N,1);
    otherwise
        error('unknown type of coverage! It shoud be chosen from [area, linear, Rayleigh].');
end

end