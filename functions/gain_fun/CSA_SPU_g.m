function gain = CSA_SPU_g(xi_hat,gamma_k,qk)
% Calculate the gain factor for the 
% short-time Cosine Spectral Amplitude MMSE + Speech Presence
% Uncertainty estimator(CSA + SPU)
% Gamma speech priori/Gaussian noise
% double-precision 
    xi_hat = sqrt(3)./xi_hat./2; %faster
%     xi_hat = 1./xi_hat; 
    delta_hat = sqrt(gamma_k); % note this  is different with Laplacian priori model

    M_pls = xi_hat+delta_hat;
    M_mns = abs(xi_hat-delta_hat);
    Msgn = sign(xi_hat-delta_hat);

%     M_p2 = M_pls.^2/4;
%     M_m2 = M_mns.^2/4;
    M_p2 = (M_pls/2).^2;
    M_m2 = (M_mns/2).^2;

    S_pls = sqrt(M_pls);
    S_mns = sqrt(M_mns);

    C_pls = sqrt(M_pls).^3;
    C_mns = sqrt(M_mns).^3;

    I_m_n_quat = besseli(-1/4,M_m2);
    I_m_p_quat = besseli(1/4,M_m2);

    Ik_p_quat = besselk(1/4,M_p2);

    c = pi/sqrt(2);
%==========================================================================    
    qkr=qk/(1-qk);
    C_hat = exp(-M_p2);
    D_mns = exp(-xi_hat.*delta_hat);

    N = C_pls.*(besselk(3/4,M_p2)+...
        -Ik_p_quat)+...
        c.*D_mns.*C_mns.*(besseli(-3/4,M_m2)+...
        +I_m_p_quat...
        -Msgn.*(besseli(3/4,M_m2)+I_m_n_quat));

    D = 2*delta_hat.*(S_pls.*Ik_p_quat+...
        c.*D_mns.*S_mns.*(I_m_n_quat...
        -Msgn.*I_m_p_quat)+2*sqrt(2*pi)*C_hat./sqrt(xi_hat)*qkr); %double check on this
    
%     D(D == 0) = eps;  
    gain = N./D;
%==========================================================================
%     N = exp(M_p2).*C_pls.*(besselk(3/4,M_p2)-Ik_p_quat)+...
%         c.*exp(M_m2).*C_mns.*(besseli(-3/4,M_m2)...
%         +I_m_p_quat-Msgn.*(besseli(3/4,M_m2)+I_m_n_quat));
%     
%     D =  2*delta_hat.*(exp(M_p2).*S_pls.*Ik_p_quat+...
%         c.*exp(M_m2).*S_mns.*(I_m_n_quat-Msgn.*I_m_p_quat));
%     
%     gain  = real(N./D);
%     
%     qkr=(1-qk)/qk;
%     lambda = qkr.*sqrt(xi_hat)./(2*sqrt(2*pi)).*(exp(M_p2).*S_pls.*Ik_p_quat...
%         +c.*exp(M_m2).*S_mns.*(I_m_n_quat-Msgn.*I_m_p_quat));
%     
%     gain = gain.*lambda./(lambda+1);
end