function [xfinal,sig_mod,n_mod]= NBLG_t(varargin)
% Description:
%   Linear Bilateral Laplactian Gain Estimator (LBLG)
%--------------------------------------------------------------------------
%-------------------------- Input validation ------------------------------
%--------------------------------------------------------------------------
in = inputParser;
addParameter(in,'noisy',@(x) isnumeric(x)); % noisy speech signal
addOptional(in,'noise',@(x) isnumeric(x));  % noise signal for ideal variance estimate
addOptional(in,'clean',@(x) isnumeric(x));  % clean signal for ideal variance estimate

% define signal parameters and short-time parameters
addParameter(in,'Fs',@(x) isnumeric(x) && x>0 ); % Sampling frequency (samples/sec)
addParameter(in,'frame_len',20,@(x) isnumeric(x) && x>0); % frame length in samples
addParameter(in,'shift_len',10,@(x) isnumeric(x) && x>0); % frame shift in samples
addParameter(in,'xi_min',10^(-35/10),@(x) isnumeric(x) && x>0); %in the text book is -15dB. p219
addParameter(in,'gain_min',eps,@(x) isnumeric(x) && x>0); %in the text book is -15dB. p219

% synthesis window type
default_opt = 'modified';
valid_opt = {'modified','hamming','rect'};
check_opt = @(x) any(validatestring(x,valid_opt));
addParameter(in,'synWinType',default_opt,check_opt);

% a-priori SNR estimation method
default_opt = 'DD'; % Decision-directed
valid_opt = {'ideal','DD','noise_ideal'};
check_opt = @(x) any(validatestring(x,valid_opt));
addParameter(in,'xi_type',default_opt,check_opt);

default_opt = 'normal';
valid_opt = {'norm','lap','gamma','laplace','normal','gauss','gaussian','laplacian'};
check_opt = @(x) any(validatestring(x,valid_opt));
addParameter(in,'speech_priori',default_opt,check_opt);

in.parse(varargin{:});

noisy = in.Results.noisy;
noise = in.Results.noise;
clean = in.Results.clean;
xi_min = in.Results.xi_min;
gain_min = in.Results.gain_min;

Fs = in.Results.Fs;
frame_len = in.Results.frame_len;
shift_len = in.Results.shift_len;
inputs.synWinType = in.Results.synWinType;

speech_priori = in.Results.speech_priori;
xi_type = in.Results.xi_type;
%%
f_inputs.frameLen = frame_len; % frame length in ms
f_inputs.shiftLen = shift_len; % frame shift in ms
f_noisy = feature_Class(noisy,Fs,f_inputs);
f_noise = feature_Class(noise,Fs,f_inputs);
f_clean = feature_Class(clean,Fs,f_inputs);
%%  constant
inputs.Ts = f_noisy.Ts; %window-shift length in discrete time
noisy_mag = f_noisy.abs; % get noisy STDCT magnitude spectrum
noisy_pow = noisy_mag.^2;
spec = f_noisy.spec_dct;
[Nframes,~] = size(noisy_mag);
%--------------------------------------------------------------------------
mag_noise = f_noise.abs;
mag_clean = f_clean.abs;
spec_noise = f_noise.spec_dct;
spec_clean = f_clean.spec_dct;

%%
switch speech_priori
    
    case {'norm','normal','gaussian'}
        
        G_p = @(xi,gamma) xi./(xi+1) + sqrt(xi./(xi+1)).*sqrt(2/pi./gamma).*...
            (exp(-gamma./xi./(1+xi)/2)-exp(-gamma.*xi./(1+xi)/2))...
            ./(erf(sqrt(gamma./xi./(1+xi)/2))+erf(sqrt(gamma.*xi./(1+xi)/2)));
        
        G_m = @(xi,gamma) xi./(xi+1) - sqrt(xi./(xi+1)).*sqrt(2/pi./gamma).*...
            (exp(-gamma./xi./(1+xi)/2) - exp(-gamma.*xi./(1+xi)/2))...
            ./(erfc(sqrt(gamma./xi./(1+xi)/2)) + erfc(sqrt(gamma.*xi./(1+xi)/2)));
        
    case {'lap','laplace','laplacian'}
        
        G_p = @(xi,gamma) sqrt(xi./gamma/2)./(1 - sqrt(xi)) - ...
            exp(sqrt(2*gamma))./(exp(sqrt(2*gamma./xi))-exp(sqrt(2*gamma)));
        
        G_m = @(xi,gamma) sqrt(xi./gamma/2)./(1 + sqrt(xi)) .* ...
            (exp(-sqrt(2*gamma./xi)) - exp(-sqrt(2*gamma))) ./ ...
            (exp(-sqrt(2*gamma./xi)) + exp(-sqrt(2*gamma))) +  ...
            exp(-sqrt(2*gamma./xi))./ ...
             (exp(-sqrt(2*gamma./xi)) + exp(-sqrt(2*gamma)));     
end

order = '2';
M = 9; % number of neibouring frames

%%
switch xi_type
    case 'ideal' % get ideal xi
        noise_pow = mag_noise.^2;
        clean_pow = mag_clean.^2;
        
        xi_frames = clean_pow./noise_pow;
        xi_frames = max(xi_min,xi_frames);
        gamma_k= max(min(noisy_pow./(noise_pow),1000),0.001);

        for n=1:Nframes
          
            if n >= M
                noisy_Bins = noisy_mag(n-(M-1):n,:);
                est = LSF(order, (1:M)*0.01, noisy_Bins);
            else
                noisy_Bins =  noisy_mag(1:M,:);
                est = LSF(order, (1:M)*0.01, noisy_Bins);
            end
            
            estXk = abs(est(:,end))';
        
            sign_est = ( noisy_mag(n,:) >= estXk );
            
            spec(n,:)  = sign_est.*spec(n,:).*G_p(xi_frames(n,:),gamma_k(n,:)) + ~sign_est.*spec(n,:).*G_m(xi_frames(n,:),gamma_k(n,:));      
            
            noisy_mag(n,:) = abs(spec(n,:));
        end
        polar_noisy = sign(spec);
    case {'DD','noise_ideal'}
        
        alpha = 0.98;
        clean_est_frame = [];
        
        if strcmp(xi_type , 'noise_ideal')
            noise_pow = mag_noise.^2;
        else
            %Initialize noise sample mean vector using first 6 frames- which assumed to be noise/silence.
            noise_pow = init_noise(noisy_mag,frame_len,shift_len);
            qp=struct;
            [noise_pow]= estnoiseg_dct(noisy_pow,shift_len.*1e-3,qp,noise_pow);
        end
        
        for n=1:Nframes
            %% ====================== Compute the MMSE Gain ==========================  
            [xi,gamma_k] = est_xi_norm(noisy_pow(n,:),noise_pow(n,:),xi_min,alpha,n,clean_est_frame);
                        
            if n >= M
                noisy_Bins = noisy_mag(n-(M-1):n,:);
                est = LSF(order, (1:M)*0.01, noisy_Bins);
            else
                noisy_Bins =  noisy_mag(1:M,:);
                est = LSF(order, (1:M)*0.01, noisy_Bins);
            end
            
            estXk = abs(est(:,end))';
            
            sign_est = ( noisy_mag(n,:) >= estXk );
          
            spec(n,:)  = sign_est.*spec(n,:).*max(G_p(xi,gamma_k),gain_min)+ ~sign_est.*spec(n,:).*max(G_m(xi,gamma_k),gain_min);
            
            clean_est_frame = abs(spec(n,:)).^2;
            
            noisy_mag(n,:) = abs(spec(n,:));   
 %--------------------------------------------------------------------------
            spec_clean(n,:) = sign_est.*spec_clean(n,:).*max(G_p(xi,gamma_k),gain_min) + ~sign_est.*spec_clean(n,:).*max(G_m(xi,gamma_k),gain_min);
            mag_clean(n,:) = abs(spec_clean(n,:)); 
            spec_noise(n,:) = sign_est.*spec_noise(n,:).*max(G_p(xi,gamma_k),gain_min) + ~sign_est.*spec_noise(n,:).*max(G_m(xi,gamma_k),gain_min);
            mag_noise(n,:) = abs(spec_noise(n,:));
            %--------------------------------------------------------------------------
        end
        polar_noisy = sign(spec);
        polar_sig = sign(spec_clean);
        polar_noise = sign(spec_noise);
end
%--------------------------------------------------------------------------
xfinal = f_noisy.ISTDCT(noisy_mag,polar_noisy,f_noisy.idx_mat,inputs);
sig_mod = f_noisy.ISTDCT(mag_clean,polar_sig,f_noisy.idx_mat,inputs);
n_mod = f_noisy.ISTDCT(mag_noise,polar_noise,f_noisy.idx_mat,inputs);
end
%================================E.O.F.====================================