function [xfinal,sig_mod,n_mod]= CSC_MMSE_t(varargin)
% Description:
%   short-time (Cosine) Spectral Amplitude estimator with
%   real Gaussian speech/noise priori.
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
%%  constant
inputs.Ts = f_noisy.Ts; %window-shift length in discrete time
noisy_mag = f_noisy.abs; % get noisy STDCT magnitude spectrum
noisy_pow = noisy_mag.^2;
spec = f_noisy.spec_dct;
polar_noisy = f_noisy.polar;
[Nframes,~] = size(noisy_mag);
%--------------------------------------------------------------------------
f_noise = feature_Class(noise,Fs,f_inputs);
f_clean = feature_Class(clean,Fs,f_inputs);
mag_noise = f_noise.abs;
mag_clean = f_clean.abs;
spec_noise = f_noise.spec_dct;
spec_clean = f_clean.spec_dct;
      
polar_noise = f_noise.polar;
polar_clean = f_clean.polar;
% ---------------------------
%%
switch xi_type
    case 'ideal' % get ideal xi
        noise_pow = mag_noise.^2;
        clean_pow = mag_clean.^2;
        
        xi_frames = clean_pow./noise_pow;
        xi_frames = max(xi_min,xi_frames);
        
        gamma_k=min(noisy_pow./noise_pow,40);     
        
        switch speech_priori
            case {'norm','normal','gaussian'}
                
               gain = Wiener(xi_frames);
                
            case {'lap','laplace','laplacian'}
                
               gain = CSC_lap(sqrt(xi_frames),sqrt(gamma_k));
        end
        
        spec  = spec.*gain;
        noisy_mag = abs(spec);
        
    case {'DD','noise_ideal'}
        
        switch speech_priori
            case {'norm','normal','gaussian'}
                alpha1 = 0.98;
            case {'lap','laplace','laplacian'}
                alpha2 = 0.9;
                xi_min = 10^(-25/10);
        end
        
        if strcmp(xi_type , 'noise_ideal')
            noise_pow = mag_noise.^2;
        else
            %Initialize noise sample mean vector using first 6 frames- which assumed to be noise/silence.
            noise_pow = init_noise(noisy_mag,frame_len,shift_len);
            
            qp=struct;
            
            [noise_pow]=estnoiseg_dct(noisy_pow,shift_len.*1e-3,qp,noise_pow);
        end

        clean_est_frame = [];
        
        for n=1:Nframes
            %% ====================== Process speech frames ==========================
            switch speech_priori
                case {'norm','normal','gaussian'}
                    [xi,~] = est_xi_norm(noisy_pow(n,:),noise_pow(n,:),xi_min,alpha1,n,clean_est_frame);
                    
                    gain = Wiener(xi);

                    gain = max(gain,gain_min); 

                    spec(n,:)  = spec(n,:).*gain;
                    
                    noisy_mag(n,:) = abs(spec(n,:));
                    
                    clean_est_frame = noisy_mag(n,:).^2;
                    
                case {'lap','laplace','laplacian'}
                    
                    [xi_hat,gamma_k] = est_xi_lap_ref(noisy_pow(n,:),noise_pow(n,:), xi_min, alpha2 ,n,clean_est_frame);
                    
                    gain = CSC_lap(xi_hat,sqrt(gamma_k));
                   
                    gain = max(gain,gain_min);

                    spec(n,:)  = spec(n,:).*gain;
                    
                    noisy_mag(n,:) = abs(spec(n,:));
                    
                    clean_est_frame = noisy_mag(n,:);     
            end
            spec_clean(n,:)  = spec_clean(n,:).*gain;
            mag_clean(n,:) = abs(spec_clean(n,:));
            spec_noise(n,:)  = spec_noise(n,:).*gain;
            mag_noise(n,:) = abs(spec_noise(n,:));
        end
end

xfinal = f_noisy.ISTDCT(noisy_mag,polar_noisy,f_noisy.idx_mat,inputs);
sig_mod = f_noisy.ISTDCT(mag_clean,polar_clean,f_noisy.idx_mat,inputs);
n_mod = f_noisy.ISTDCT(mag_noise,polar_noise,f_noisy.idx_mat,inputs);

end
%%-------------------------------------------------------------------------
%-------------------------------Function-----------------------------------
%--------------------------------------------------------------------------
function gain = CSC_lap(xi_hat,gamma_hat)
%
% Cosine Spectral coefficient estimator with laplacian speech prior.
%
% reference :
% 'Speech enhancement using an MMSE short time DCT coefficients estimator
% with supergaussian speech modeling', 2007, Zou and Zhang
f1 = (1./xi_hat + gamma_hat)/sqrt(2);
f2 = (1./xi_hat - gamma_hat)/sqrt(2);

N1 = f1*sqrt(2).*exp(f1.^2).*erfc(f1);
N2 = f2*sqrt(2).*exp(f2.^2).*erfc(f2);
D = exp(f1.^2).*erfc(f1) + exp(f2.^2).*erfc(f2); 

gain = (N1-N2)./D./gamma_hat;

end
%================================E.O.F.====================================