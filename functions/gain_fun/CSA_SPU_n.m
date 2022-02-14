function gain = CSA_SPU_n(xi,gamma_k,qk) 
% Calculate the gain factor for the 
% short-time Discrte Cosine Transform (DCT) 
% Spectral Amplitude MMSE + Speech Presence Uncertainty estimator(CSA +
% SPU).
% Gaussian speech priori/Gaussian noise
% double-precision 
%     Xi = Xi./(1-qk);
%     vk = gamma_k.*xi./(1+xi); 
%     qkr=qk/(1-qk);
% %     qkr=(1-qk)/qk;
%     Lambda=1+qkr*sqrt(1+xi).*exp(-vk/2);
%     
% %     gain = (erf(sqrt(vk/2)) ...
% %         + sqrt(2/pi./vk./exp(vk))).*Xi./(1+Xi);    
%     gain = (erf(sqrt(vk/2)) ...
%         + sqrt(2/pi)./sqrt(vk)./exp(vk/2)).*xi./(1+xi);
%     gain= gain./Lambda; 
    
    vk=xi.*gamma_k./(1+xi);
    qkr=(1-qk)/qk;
    c=sqrt(2/pi);
    C=exp(-0.5*vk);
    A=c./sqrt(vk).*C;
    B = erf(sqrt(0.5*vk));
    
    Lambda=qkr*exp(0.5*vk)./sqrt(1+xi);
   
    gain = (A+B).*xi./(1+xi);
    pSAP=Lambda./(1+Lambda);
    gain=gain.*pSAP;   
end
