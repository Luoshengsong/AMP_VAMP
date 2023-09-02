function NMSE = LMMSE_pro(A_active, h_active, cov_h, y, sigma2)
    h_LMMSE = cov_h * A_active' * inv(A_active*cov_h*A_active'+sigma2*eye(length(y))) * y;
    
    NMSE = norm(h_active - h_LMMSE)^2 / length(h_active);
end