function [xfinal,sig_mod,n_mod]= CSA_MMSE_SPU_t(varargin)
% Description:
%   short-time (Cosine) Spectral Amplitude estimator with
%   real Gaussian speech/noise priori and Speech Presence Uncertainty.
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
addParameter(in,'qk',0.2,@(x) isnumeric(x)); % noisy speech signal
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
%--------------------------------------------------------------------------
f_noisy = feature_Class(noisy,Fs,f_inputs);
f_noise = feature_Class(noise,Fs,f_inputs);
f_clean = feature_Class(clean,Fs,f_inputs);
%--------------------------------------------------------------------------
inputs.Ts = f_noisy.Ts; %window-shift length in discrete time
%% =============== Transform data / Set up Constants ======================
noisy_mag = f_noisy.abs; % get noisy STDCT magnitude spectrum
noisy_pow = noisy_mag.^2;
polar = f_noisy.polar;
%--------------------------------------------------------------------------
noise_mag = f_noise.abs;
clean_mag = f_clean.abs;
polar_noise = f_noise.polar;
polar_sig = f_clean.polar;
%--------------------------------------------------------------------------
switch xi_type
    case 'ideal' % get ideal xi
        noise_pow = noise_mag.^2;
        clean_pow = clean_mag.^2;
        
        gamma_k= max(min(noisy_pow./(noise_pow),1000),0.001);
        
        switch speech_priori
            case {'norm','normal','gaussian'}
                xi_frames = max(xi_min,clean_pow./noise_pow);
                gain = CSA_SPU_n(xi_frames,gamma_k,0.2);
            case {'lap','laplace','laplacian'}
%                 xi_min = 10^(-25/10);
                xi_frames = max(xi_min,clean_pow./noise_pow);
                gain = CSA_SPU_l(sqrt(xi_frames),gamma_k,0.1);
            case {'gamma'}
%                 xi_min = 10^(-30/10);
                xi_frames = max(xi_min,clean_pow./noise_pow);
                gain = CSA_SPU_g(sqrt(xi_frames),gamma_k,0.1);
        end
        noisy_mag= noisy_mag.*gain;
        polar = f_noisy.polar;
    case {'DD','noise_ideal'}
        switch speech_priori
            case {'norm','normal','gaussian'}
                alpha1 = 0.98;
                qk1 = 0.3;
            case {'lap','laplace','laplacian'}
%                 xi_min = 10^(-32/10);
                alpha2 = 0.92;
                qk2 = 0.2;
            case {'gamma'}
%                 xi_min = 10^(-30/10);
                alpha3 = 0.92;
                qk3 = 0.2;
        end
        %%
        [Nframes] = size(noisy_mag);
        clean_est_frame = [];
        %--------------------------------------------------------------------------
        if strcmp(xi_type , 'noise_ideal')
            noise_pow = noise_mag.^2;
            
            for n=1:Nframes
                %% ====================== Compute the MMSE Gain ==========================
                switch speech_priori
                    case {'norm','normal','gaussian'}
                        
                        [xi,gamma_k] = est_xi_norm(noisy_pow(n,:),noise_pow(n,:),xi_min,alpha1,n,clean_est_frame);
                        
                        gain = CSA_SPU_n(xi,gamma_k,qk1);
                        
                        gain = max(gain,gain_min); 
                        
                        noisy_mag(n,:) = noisy_mag(n,:).*gain;
                        
                        clean_est_frame = noisy_mag(n,:).^2;
                        
                    case {'lap','laplace','laplacian'}
                        
                        [xi_hat,gamma_k] = est_xi_lap(noisy_pow(n,:),noise_pow(n,:), xi_min, alpha2 ,n,clean_est_frame);
                        
                        gain = CSA_SPU_l(xi_hat,gamma_k,qk2);
                        
                        gain = max(gain,gain_min);
                        
                        noisy_mag(n,:) = noisy_mag(n,:).*gain;
                        clean_est_frame = noisy_mag(n,:);
                        
                    case {'gamma'}
                        
                        [xi_hat, gamma_k] = est_xi_gamma(noisy_pow(n,:),noise_pow(n,:), xi_min, alpha3 ,n,clean_est_frame);
                        
                        gain = CSA_SPU_g(xi_hat,gamma_k,qk3);
                        
                        gain = max(gain,gain_min);
                      
                        noisy_mag(n,:) = noisy_mag(n,:).*gain;
                        clean_est_frame = noisy_mag(n,:);
                        
                end
                %--------------------------------------------------------------------------
                noise_mag(n,:) = noise_mag(n,:).*gain;
                clean_mag(n,:) = clean_mag(n,:).*gain;
                %--------------------------------------------------------------------------
            end
        else
            noise_pow = init_noise(noisy_mag,frame_len,shift_len);
            %             noise_pow = repmat(noise_pow,size(noisy_mag,1),1);
            %--------------------------------------------------------------------------
            %             qp=struct;
            %             [noise_pow] = estnoiseg_dct(noisy_pow,shift_len.*1e-3,qp,noise_pow);
            %--------------------------------------------------------------------------
            q = 0.5; % a priori probability of speech presence:
            PH1mean = 0.5;
            alphaPH1mean = 0.95;
            alphaPSD = 0.95;
            priorFact  = q./(1-q);
            xiOptDb    = 15.6; % optimal fixed a priori SNR for SPP estimation
            xiOpt      = 10.^(xiOptDb./10);
            logGLRFact = log(1./sqrt(1+xiOpt));
            GLRexp     = xiOpt./(1+xiOpt)/2;
            for n=1:Nframes
                %% ====================== Compute the MMSE Gain ==========================
                switch speech_priori
                    case {'norm','normal','gaussian'}
                        GLR  = priorFact .* exp(min(logGLRFact + GLRexp.*(noisy_pow(n,:)./noise_pow),200));
                        PH1   = GLR./(1+GLR); % a posteriori speech presence probability
                        
                        PH1mean = alphaPH1mean * PH1mean + (1-alphaPH1mean) * PH1;
                        stuckInd = PH1mean > 0.99;
                        PH1(stuckInd) = min(PH1(stuckInd),0.99);
                        estimate =  PH1 .* noise_pow + (1-PH1) .* noisy_pow(n,:) ;
                        noise_pow = alphaPSD *noise_pow+(1-alphaPSD)*estimate;
                        
                        %                     [xi,gamma_k] = est_xi_norm(noisy_pow(n,:),noise_pow(n,:),xi_min,alpha1,n,clean_est_frame);
                        [xi,gamma_k] = est_xi_norm(noisy_pow(n,:),noise_pow,xi_min,alpha1,n,clean_est_frame);
                        
                        gain = CSA_SPU_n(xi,gamma_k,qk1);
                        
                        %                     M = [xi.',(gamma_k).'];
                        %                     gain = Cn_SPU(M)';
                        %                     gain = max(gain,MIN_GAIN);
                        gain = max(gain,gain_min);
                        
                        noisy_mag(n,:) = noisy_mag(n,:).*gain;
                        
                        clean_est_frame = noisy_mag(n,:).^2;
                        
                    case {'lap','laplace','laplacian'}
                        
                        GLR  = priorFact .* exp(min(logGLRFact + GLRexp.*(noisy_pow(n,:)./noise_pow),200));
                        PH1   = GLR./(1+GLR); % a posteriori speech presence probability
                        
                        PH1mean = alphaPH1mean * PH1mean + (1-alphaPH1mean) * PH1;
                        stuckInd = PH1mean > 0.99;
                        PH1(stuckInd) = min(PH1(stuckInd),0.99);
                        estimate =  PH1 .* noise_pow + (1-PH1) .* noisy_pow(n,:) ;
                        noise_pow = alphaPSD *noise_pow+(1-alphaPSD)*estimate;
                        
                        [xi_hat,gamma_k] = est_xi_lap(noisy_pow(n,:),noise_pow, xi_min, alpha2 ,n,clean_est_frame);

                        gain = CSA_SPU_l(xi_hat,gamma_k,qk2);
                        
                        gain = max(gain,gain_min);
                        
                        noisy_mag(n,:) = noisy_mag(n,:).*gain;
                        clean_est_frame = noisy_mag(n,:);

                    case {'gamma'}
                        
                        GLR  = priorFact .* exp(min(logGLRFact + GLRexp.*(noisy_pow(n,:)./noise_pow),200));
                        PH1   = GLR./(1+GLR); % a posteriori speech presence probability
                        
                        PH1mean = alphaPH1mean * PH1mean + (1-alphaPH1mean) * PH1;
                        stuckInd = PH1mean > 0.99;
                        PH1(stuckInd) = min(PH1(stuckInd),0.99);
                        estimate =  PH1 .* noise_pow + (1-PH1) .* noisy_pow(n,:) ;
                        noise_pow = alphaPSD *noise_pow+(1-alphaPSD)*estimate;
                        
                        [xi_hat,gamma_k] = est_xi_gamma(noisy_pow(n,:),noise_pow, xi_min, alpha3 ,n,clean_est_frame);
                        
                        gain = CSA_SPU_g(xi_hat,gamma_k,qk3);
                        
                        gain = max(gain,gain_min);
                        
                        noisy_mag(n,:) = noisy_mag(n,:).*gain;
                        clean_est_frame = noisy_mag(n,:);
                        
                end
                %--------------------------------------------------------------------------
                noise_mag(n,:) = noise_mag(n,:).*gain;
                clean_mag(n,:) = clean_mag(n,:).*gain;
                %--------------------------------------------------------------------------
            end
        end
end

xfinal = f_noisy.ISTDCT(noisy_mag,polar,f_noisy.idx_mat,inputs);
sig_mod = f_noisy.ISTDCT(clean_mag,polar_sig,f_noisy.idx_mat,inputs);
n_mod = f_noisy.ISTDCT(noise_mag,polar_noise,f_noisy.idx_mat,inputs);
end