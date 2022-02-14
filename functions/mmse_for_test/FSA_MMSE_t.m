function [xfinal,sig_mod,n_mod]= FSA_MMSE_t(varargin)
% Description:
%   short-time (Fourier) Spectral Amplitude estimator with
%   complex Gaussian speech/noise priori.
%--------------------------------------------------------------------------
%-------------------------- Input validation ------------------------------
%--------------------------------------------------------------------------
in = inputParser;
addParameter(in,'noisy',@(x) isnumeric(x)); % noisy speech signal
addOptional(in,'noise',@(x) isnumeric(x));  % noise signal for ideal variance estimate
addOptional(in,'clean',@(x) isnumeric(x));  % clean signal for ideal variance estimate

% define signal parameters and short-time parameters
addParameter(in,'Fs',@(x) isnumeric(x) && x>0 ); % Sampling frequency (samples/sec)
addParameter(in,'frame_len',32,@(x) isnumeric(x) && x>0); % frame length in samples
addParameter(in,'shift_len',8,@(x) isnumeric(x) && x>0); % frame shift in samples
addParameter(in,'xi_min',10^(-30/10),@(x) isnumeric(x) && x>0); %in the text book is -15dB. p219
addParameter(in,'gain_min',eps,@(x) isnumeric(x) && x>0); %in the text book is -15dB. p219

% synthesis window type
default_opt = 'modified';
valid_opt = {'modified','hamming','rect','hanning'};
check_opt = @(x) any(validatestring(x,valid_opt));
addParameter(in,'synWinType',default_opt,check_opt);

% a-priori SNR estimation method
default_opt = 'DD'; % Decision-directed
valid_opt = {'ideal','DD','noise_ideal'};
check_opt = @(x) any(validatestring(x,valid_opt));
addParameter(in,'xi_type',default_opt,check_opt);

default_opt = 'normal';
valid_opt = {'norm','lap','laplace','gamma','normal','gauss','gaussian','laplacian'};
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
%--------------------------------------------------------------------------
clean_mag = f_clean.mag;
noise_mag = f_noise.mag;
% get noisy STDFT spectra
noisy_mag = f_noisy.mag; % short-time magnitude spectrum
noisy_pow = noisy_mag.^2;
pha_noise = f_noise.pha;
pha_sig = f_clean.pha;
pha = f_noisy.pha;
%%
inputs.Ts = f_noisy.Ts; %window-shift length in discrete time
%%
global Gmag

switch xi_type
    case 'ideal'
        noise_pow = noise_mag.^2;
        clean_pow = clean_mag.^2;
        xi_frames = clean_pow./noise_pow;
        xi_frames = max(xi_min,xi_frames);

        gamma_k=min(noisy_pow./noise_pow,40);     
        
        M = num2cell(xi_frames,2);
        N = num2cell(gamma_k,2);
        L = cat(2,M,N);
       
        switch speech_priori
            case {'norm','normal','gaussian'}
                gain = stsa_dft_n(xi_frames,gamma_k);
            case {'lap','laplace','laplacian'}
                for i = 1:length(L)
                    gain(i,:) = lookup_gain_in_table(Gmag,L{i,2}.',L{i,1}.',-40:1:50,-40:1:50,1);
                end
        end
        noisy_mag = noisy_mag.*gain;
    case {'DD','noise_ideal'}
        %Set parameters for decistion-directed xi-estimate in STDFT domain
        alpha=0.98; 
        
        [Nframes,~] = size(noisy_mag);
        %% Initialize noise sample mean vector using first 6 frames- which assumed to be noise/silence.
        % initialize noise variance (approx) for decision-directed method
        if strcmp(xi_type , 'noise_ideal')
            noise_pow = noise_mag.^2;
        else       
            noise_pow = init_noise(noisy_mag,frame_len,shift_len);
            qp=struct;
            [noise_pow]=estnoiseg_dft(noisy_pow,shift_len.*1e-3,qp,noise_pow);
        end
        
        for n=1:Nframes
            %% ========================Calculate posteriori SNR========================
            % posteriori SNR : gamma_k = R_k^2/lambda_D;
            % where R_k is the noisy magnitude, lambda_D is the noise variance
            gamma_k= max(min(noisy_pow(n,:)./noise_pow(n,:),1000),0.001);
            %% =============decision-direct estimate of a priori SNR===================
            % alpha control the trade-off between speech distortion and residual noise.
            if n==1  % initial condition for recursive computation
                xi=alpha+(1-alpha)*max(gamma_k-1,0);
            else           
                xi=alpha*Xk_prev./noise_pow(n,:) + (1-alpha)*max(gamma_k-1,0);
                xi=max(xi_min,xi);  % limit ksi to -25 dB, but in the text book is -15dB. p219
            end   
             %% update/estimate the noise spectrum using VAD algorithm or noise-estimation alghorithm
             switch speech_priori
                 case {'norm','normal','gaussian'}
                    gain = stsa_dft_n(xi,gamma_k);
                    
                 case {'lap','laplace','laplacian'}
                     
                    gain = lookup_gain_in_table(Gmag,gamma_k,xi,-40:1:40,-40:1:40,1);
                    
             end
             
             gain = max(gain,gain_min);
             
             noisy_mag(n,:) = noisy_mag(n,:).*gain;
             
             Xk_prev = noisy_mag(n,:).^2;
             %--------------------------------------------------------------------------
             noise_mag(n,:) = noise_mag(n,:).*gain;
             clean_mag(n,:) = clean_mag(n,:).*gain;
             %--------------------------------------------------------------------------
        end   
end

xfinal = f_noisy.ISTDFT(noisy_mag,pha,f_noisy.idx_mat,f_noisy.N_analysis,inputs);
sig_mod = f_noisy.ISTDFT(clean_mag,pha_sig,f_noisy.idx_mat,f_noisy.N_analysis,inputs);
n_mod = f_noisy.ISTDFT(noise_mag,pha_noise,f_noisy.idx_mat,f_noisy.N_analysis,inputs);
end
%--------------------------------------------------------------------------
function xhat = stsa_dft_l_double(ksi,gammak,Yk,Yk2) % also slow 190s

f= @(x,v)(x.^v).*exp((x.^2)*(-gammak./Yk2)-x.*sqrt(2*gammak./ksi)./Yk).*besseli(0,x.*2.*gammak./Yk);

n = integral(@(x)f(x,2),0,80,'ArrayValued',true);
d = integral(@(x)f(x,1),0,80,'ArrayValued',true);

n(isnan(n)) = 0;
d(isnan(d)) = 1;
xhat = n./d;
end
%==========================================================================
function gain = stsa_dft_g_double(xi,gamma_k,v)
% Minimum Mean-Square Error Estimation of Discrete Fourier Coefficients
% With Generalized Gamma Priors (2007);
% gamma = 1; v = 1; for Laplacian speech priori
    mu = sqrt(v*(v+1));
    f = @(y,x,g,v,mu)(y.^v).*exp(-y.^2./(4*g)-mu*y./(2*sqrt(x.*g))).*besseli(0,y); %eq.(14)
    n = integral(@(y)f(y,xi,gamma_k,v,mu),0,80,'ArrayValued',true);
    d = integral(@(y)f(y,xi,gamma_k,v-1,mu),0,80,'ArrayValued',true);
    gain = n./d./gamma_k/2;
end
%==========================================================================
function gain = stsa_dft_gAprx_double(xi,gamma_k,v,K)
    gain1 = A1_S(xi,gamma_k,v,K);
    gain2 = A1_L(xi,gamma_k,v);
    gain = max(gain1,gain2);
end
%==========================================================================
% function [gain] = A1(xi,gamma_k,v)
%     mu = sqrt(v*(v+1));
%     f = @(y,x,g,v,mu)(y.^v).*exp(-y.^2./(4*g)-mu*y./(2*sqrt(x.*g))).*besseli(0,y);
% %     gain1 = zeros(size(xi));
% %     parfor i = 1: length(xi) %slower
% %         x =xi(i);
% %         g = gamma_k(i);
% %         gain1(i) = integral(@(y)f1(y,x,g,v,mu),0,50);
% %     end
%     n = integral(@(y)f(y,xi,gamma_k,v,mu),0,80,'ArrayValued',true);
%     d = integral(@(y)f(y,xi,gamma_k,v-1,mu),0,80,'ArrayValued',true);
%     gain = n./d./gamma_k/2;
% end
%==========================================================================
function gain = A1_S(Xi,gamma_k,v,K)
% eq.(16): for low SNRs
% Minimum Mean-Square Error Estimation of Discrete Fourier Coefficients
% With Generalized Gamma Priors (2007);
% gamma = 1; v = 1; for Laplacian speech priori
D1= @(p,z) pu(-(p+1/2),z);

K = 0:K;
mu = sqrt(v*(v+1));
q = sqrt(2*gamma_k);
p = mu./sqrt(2.*Xi);
N = size(gamma_k,2);
nu = zeros(1,N);
de = zeros(1,N);

G = repmat(gamma_k,size(K,2),1)';
P = repmat(p,1,N);

parfor i = 1:N
    nu(i) = trapz(K,(G(i,:)/2).^K.*gamma(v+2*K+1).*D1(-(v+2*K+1),P(i))'./(factorial(K).^2));
    de(i) = trapz(K,(G(i,:)/2).^K.*gamma(v+2*K).*D1(-(v+2*K),P(i))'./(factorial(K).^2));
end

gain = nu./de./q;
end
%==========================================================================
function gain = A1_L(Xi,gamma_k,v)
% eq.(18): for large SNRs
% Minimum Mean-Square Error Estimation of Discrete Fourier Coefficients
% With Generalized Gamma Priors
% gamma = 1; v = 1; for Laplacian speech priori
%--------------------------------------------------------------------------
global parabolic_x_neg3_2 parabolic_x_neg1_2 parabolic_y_neg3_2 parabolic_y_neg1_2
    h = v-1/2;
    mu = sqrt(v*(v+1));
    
    q = sqrt(2*gamma_k);
    p = mu./sqrt(2.*Xi)-q;
    
    nu = interp1(parabolic_x_neg3_2,parabolic_y_neg3_2,p,'spline');
    de = interp1(parabolic_x_neg1_2,parabolic_y_neg1_2,p,'spline');    
    gain = h*nu./de./q;
%--------------------------------------------------------------------------    
% D1= @(p,z) pu(-(p+1/2),z);
% 
% l = v+1/2;
% h = v-1/2;
% mu = sqrt(v*(v+1));
% 
% q = sqrt(2*gamma_k);
% p = mu./sqrt(2.*Xi)-q;
% N = size(gamma_k,2);
% 
% nu = zeros(1,N);
% de = zeros(1,N);
% 
% P = repmat(p,1,N);
% 
% parfor i = 1:N
%     nu(i) = real(D1(-l,P(i)));
%     de(i) = real(D1(-h,P(i)));
% end
% 
% gain = h*nu./de./q;
end
%==========================================================================
function gain = stsa_dft_b_double(Xi,gamma_k)
% beta order, 2015: eq.(28)
% Speech enhancement based on beta-order MMSE estimation of STSA and
% Laplacian speech modeling
%
%     W = @(l,m,z) (gamma(-2*m)./gamma((1/2)-m-l)).*whittakerM(l,m,z)+...
%         (gamma(2*m)./gamma((1/2)+m-l)).*whittakerM(l,-m,z);
%     D0 = @(p,z) 2.^(1/4+p/2).*W((1/4+p/2),-1/4,(z.^2)/2)./sqrt(z);
%     D0= @(p,z) 2.^(p/2).*exp(-z.^2/4).*...
%         (sqrt(pi)./gamma((1-p)/2).*hypergeom(-p/2,1/2,(z.^2)/2)-...
%         sqrt(2*pi).*z./gamma(-p/2).*hypergeom((1-p)/2,3/2,(z.^2)/2));
% D0&D2 is same as D1, but slower
% array as input
%     D1= @(p,z) 2.^(1/4+p/2).*whittakerW((1/4+p/2),-1/4,(z.^2)/2)./sqrt(z);
% scalar as input
%     D1= @(p,z) 2^(1/4+p/2)*whittakerW((1/4+p/2),-1/4,(z^2)/2)/sqrt(z);
%==========================================================================
D1= @(p,z) pu(-(p+1/2),z);

q = sqrt(2*gamma_k);
p = 2*sqrt(2./Xi)./pi;
z = p-q;

nu = zeros(size(gamma_k));
de = zeros(size(gamma_k));

%     NU = zeros(size(gamma_k));
%     DE = zeros(size(gamma_k));

%     G = repmat(gamma_k,size(m,2),1)';

parfor i = 1:length(z)
    nu(i) = real(D1(-5/2,z(i)));
    de(i) = real(D1(-3/2,z(i)));
    %         NU(i) = trapz(m,(G(i,:)/2).^m.*gamma(2*m+3).*D1(-(2*m+3),p(i))'./(factorial(m).^2));
    %         DE(i) = trapz(m,(G(i,:)/2).^m.*gamma(2*m+2).*D1(-(2*m+2),p(i))'./(factorial(m).^2));
end

%     gain_1 = NU./DE./q;
%     gain_2 = (3/2).*nu./de./q;
%     gain = max(gain_1,gain_2);
gain = (3/2).*nu./de./q;
%     gain = NU./DE./q;
end
% --------------------------------------------------------------------------
% function xhat = stsa_dft_l_double(ksi,gammak,Yk,Yk2) % also slow 119s
% %usage:  ksi:aprioriSNR; gammak:aposteriorSNR
% %Yk: vector of DFT coefficients of noisy spectrum
% %infi: a value used to approximate infinity (here infi=40)
%     infi = 40;
%     x= 0:0.01:infi;
%     theta=0:0.01:pi/4;
%
%     aa = (-sqrt(2)*gammak./ksi./Yk2);
%
%     a = (x.^2)'*(-gammak./Yk2);
%     b = x'*(2*gammak./Yk);
%
%     fn2= (x.^2)'.*exp(a).*besseli(0,b);
%     fn3= x'.*exp(a).*besseli(0,b);
%
%     fn2(isnan(fn2)==1)=0; %workaround for the case that inf*infsmall
%     fn3(isnan(fn3)==1)=0; %it’s set to a very small value because
%
%     fn1 = zeros(size(fn2));
% tic
%     parfor i = 1:length(x)
%         fn1(i,:) = trapz(theta,exp(cos(theta').*x(i)*aa),1);
%     end
% toc
% tic
%     A = trapz(bsxfun(@times, fn1, fn2),1);
%     C = trapz(bsxfun(@times, fn1, fn3),1);
% toc
%     xhat = A./C;
% end
% --------------------------------------------------------------------------
% function xhat = stsa_dft_l_double(ksi,gammak,Yk,Yk2) % slow 222s
%  %usage:  ksi:aprioriSNR; gammak:aposteriorSNR
%  %Yk: vector of DFT coefficients of noisy spectrum
%  %infi: a value used to approximate infinity (here infi=40)
%     infi = 40;
%     x= 0:0.01:infi;
%     theta=0:0.01:pi/4;
%
%     aa = (-sqrt(2)*gammak./ksi./Yk2);
%
%     [X,T,AA] = meshgrid(x,theta,aa);
%
%     a = (x.^2)'*(-gammak./Yk2);
%     b = x'*(2*gammak./Yk);
%
%     fn2= (x.^2)'.*exp(a).*besseli(0,b);
%     fn3= x'.*exp(a).*besseli(0,b);
%
%     fn2(isnan(fn2)==1)=0; %workaround for the case that inf*infsmall
%     fn3(isnan(fn3)==1)=0; %it’s set to a very small value because
%
%     fn1 = trapz(theta,exp(AA.*cos(T).*X),1);
%     fn1 = squeeze(fn1);
%
%     A = trapz(bsxfun(@times, fn1, fn2),1);
%     C = trapz(bsxfun(@times, fn1, fn3),1);
%
%     xhat = A./C;
% end
% % --------------------------------------------------------------------------
function xhat = stsa_dft_l_double1(ksi,gammak,Yk,Yk2) % also slow 190s
%usage:  ksi:aprioriSNR; gammak:aposteriorSNR
%Yk: vector of DFT coefficients of noisy spectrum
%infi: a value used to approximate infinity (here infi=40)
infi = 40;
x= 0:0.1:infi;
theta=0:0.01:pi/4;

aa = -sqrt(2)*gammak./ksi./Yk2;

a = (x.^2)'*(-gammak./Yk2);
b = x'*(2*gammak./Yk);
%     a = (x.^2)'*(-gammak./Yk2)-(x'*4.*sqrt(gammak./ksi)./Yk./pi);
%     b = x'*(2*gammak./Yk);
%
fn2= (x.^2)'.*exp(a).*besseli(0,b);
fn3= x'.*exp(a).*besseli(0,b);
%

fn2(isnan(fn2)==1)=eps*1e-1; %workaround for the case that inf*infsmall
fn3(isnan(fn3)==1)=eps*1e-1; %it’s set to a very small value because

fn1 = zeros(size(fn2));

parfor i = 1: length(aa)
    fn1(:,i) = trapz(theta,exp(aa(i).*cos(theta').*x),1);
end
%     A = trapz(fn2,1);
%     C = trapz(fn3,1);
A = trapz(bsxfun(@times, fn1, fn2),1);
C = trapz(bsxfun(@times, fn1, fn3),1);
xhat = A./C;
end
% %--------------------------------------------------------------------------
% function xhat = stsa_dft_l_double(ksi,gammak,Yk,Yk2) % slow 133s/180
%  %usage:  ksi:aprioriSNR; gammak:aposteriorSNR
%  %Yk: vector of DFT coefficients of noisy spectrum
%  %infi: a value used to approximate infinity (here infi=40)
%     infi = 40;
%     x= 0:0.1:infi;
%
%     aa = -sqrt(2*gammak./ksi)./Yk2;
%     a = (x.^2)'*(-gammak./Yk2);
%     b = x'*(2*gammak./Yk);
%
%     fn2= (x.^2)'.*exp(a).*besseli(0,b);
%     fn3= x'.*exp(a).*besseli(0,b);
%     fn2(isnan(fn2)==1)=0; %workaround for the case that inf*infsmall
%     fn3(isnan(fn3)==1)=0; %it’s set to a very small value because
%
%     fn1 = zeros(size(x,2),size(gammak,2));
%
%      parfor i = 1:size(x,2)
%         fn1(i,:) = integral(@(theta)exp(aa.*x(i).*cos(theta))...
%                 ,0,pi/4,'ArrayValued',true,'RelTol', 1e-4, 'AbsTol', 1e-6);
%      end
%
%     A = trapz(bsxfun(@times, fn1, fn2),1);
%     C = trapz(bsxfun(@times, fn1, fn3),1);
%     xhat = A./C;
% end
%--------------------------------------------------------------------------
function gain = stsa_dft_l_double2(ksi,gamma_k,Yk)

m= 0:20;
D1= @(p,z) pu(-(p+1/2),z);

q = sqrt(2*gamma_k);
p = sqrt(1./ksi);
z1 = sqrt(gamma_k);
z2 = p-q;

nu = zeros(size(gamma_k));
de = zeros(size(gamma_k));
NU = zeros(size(gamma_k));
DE = zeros(size(gamma_k));

G = repmat(gamma_k./Yk,size(m,2),1)';
K = repmat(2./Yk.*(sqrt(gamma_k./ksi)),size(m,2),1)';

parfor i = 1:length(gamma_k)
    %         nu(i) = real(D1(-5/2,z2(i)));
    %         de(i) = real(D1(-3/2,z2(i)));
    NU(i) = trapz(m,G(i,:).^(2.*m).*K(i,:).^(-m)...
        .*gamma(2*m+3).*exp(1./(4.*ksi(i))).*D1(-(2*m+3),z1(i))'./(factorial(m).^2));
    DE(i) = trapz(m,G(i,:).^(2.*m).*K(i,:).^(-m)...
        .*gamma(2*m+2).*exp(1./(4.*ksi(i))).*D1(-(2*m+2),z1(i))'./(factorial(m).^2));
end

%       gain_1 = NU./DE./q;
%       gain_2 = (3/2).*nu./de./q;
%       gain = max(gain_1,gain_2);
%       gain = (3/2).*nu./de./q;
gain = NU./DE./q;
end