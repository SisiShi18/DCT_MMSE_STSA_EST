function gain = CSA_SPU_l(xi_hat,gamma_k,qk)
% Calculate the gain factor for the 
% short-time Cosine Spectral Amplitude MMSE + Speech Presence
% Uncertainty estimator(CSA + SPU)
% Laplacian speech priori/Gaussian noise
% double-precision 
%     Xi = Xi./(1-qk);
%     qkr=qk/(1-qk);
%     xi_hat = 1./sqrt(Xi);
%     delta_hat = sqrt(gamma_k/2);
% 
%     E_plus = xi_hat+delta_hat;
%     E_minus= xi_hat-delta_hat;
%     
% %     D_plus = exp(4.*(xi_hat.*delta_hat));
% %     I_plus = D_plus.*erfc(E_plus);
% %     I_minus = erfc(E_minus);
% %     
% %     C_hat = exp(-E_minus.^2);
% %     c = C_hat.*2/sqrt(pi);
% %     
% %     gain = (c...
% %         -(E_plus.*I_plus+E_minus.*I_minus))./(delta_hat.*...
% %         (I_plus+I_minus+c.*sqrt(Xi).*qkr));
%     
%     C_hat = exp(-E_plus.^2);
%     c = C_hat*2./sqrt(pi);
%     
%     D_minus = exp(-4.*(xi_hat.*delta_hat));
%     I_plus = erfc(real(E_plus));
%     I_minus = D_minus.*erfc(real(E_minus));    
%     
%     N = (c...
%         -(E_plus.*I_plus+E_minus.*I_minus));
%     
%     D = (delta_hat.*...
%          (I_plus+I_minus+c.*sqrt(Xi).*qkr));
%      
%     D(D == 0) = eps;  
%     gain = N./D;         
%     xi_hat = 1./sqrt(Xi);
    xi_hat = 1./xi_hat;
    delta_hat = sqrt(gamma_k/2);
    
    E_plus = xi_hat+delta_hat;
    E_minus= xi_hat-delta_hat;
%===============================================    
%     C_hat = exp(-E_minus.^2);
% 
%     c =C_hat.*2./sqrt(pi);
%     D_plus = exp(4.*(xi_hat.*delta_hat));
%     I_plus = D_plus.*erfc(E_plus);
%     I_minus = erfc(E_minus);
% 
%     gain = (c...
%         -(E_plus.*I_plus+E_minus.*I_minus))./(delta_hat.*(I_plus+I_minus));
%     gain(isnan(gain)) = eps;  
%     
%     qkr=qk/(1-qk);
%     Lambda=1+qkr*c.*xi_hat./(I_plus+I_minus);
%     gain= gain./Lambda;  
%================================================
    I_plus = exp(E_plus.^2).*erfc(E_plus);
    I_minus = exp(E_minus.^2).*erfc(E_minus);
    
    gain = (2/sqrt(pi)-(E_plus.*I_plus+E_minus.*I_minus))./(delta_hat.*(I_plus+I_minus));

    qkr=(1-qk)/qk;
    
    lambda = qkr.*xi_hat.*(sqrt(pi)/2).*(I_plus+I_minus);
     
    gain = gain.*lambda./(lambda+1);
end



