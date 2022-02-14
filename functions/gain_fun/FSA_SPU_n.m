function gain = FSA_SPU_n(Xi,gamma_k,qk)
% Calculate the gain factor for the 
% short-time Discrete Fourier Transform (DFT)
% Spectral Amplitude MMSE + Speech Presence
% Uncertainty estimator( FSA + SPU)
% Gaussian speech priori/Gaussian noise
% double-precision 
    qkr=(1-qk)/qk;
%     qkr=qk/(1-qk);
%     vk = gamma_k.*Xi./(1+Xi);
% 
%     pSAP= 1+qkr*(1+Xi).*exp(-vk);
% 
%     A=(sqrt(pi)/2)*sqrt(vk).*exp(-vk/2)./gamma_k;
%     B=(1+vk).*(besseli(0,vk/2))+vk.*(besseli(1,vk/2));
% 
%     gain= A.*B./pSAP;
    
    vk=Xi.*gamma_k./(1+Xi);
    j0=besseli(0,vk/2);
    j1=besseli(1,vk/2);
    c=sqrt(pi)/2;

    C=exp(-0.5*vk);
    A=((c*(vk.^0.5)).*C)./gamma_k;
    B=(1+vk).*j0+vk.*j1;
    hw=A.*B;
    % --- estimate speech presence probability
    evk=exp(vk);
    Lambda=qkr*evk./(1+Xi);
    pSAP=Lambda./(1+Lambda);
    gain=hw.*pSAP;
end