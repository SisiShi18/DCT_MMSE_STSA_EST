%% Script Description:
%   Implement DCT Short-Time Spectral Amplitude MMSE estimators (CSA);
%   Output enhanced speech files in the audios folder.
%   Output plots of objective quality measures in the fig folder.
% %--------------------------------------------------------------------------
clear all; close all;
clc;
%--------------------------------------------------------------------------
%%                          Progress indication
%--------------------------------------------------------------------------
progress = {'-','\','\','|','/'}; P = length(progress); pp = 0;
%--------------------------------------------------------------------------
%%                 Define working directory and file path
%--------------------------------------------------------------------------
fprintf('Setting working directory and file paths...\n')
currentFolder = pwd;
addpath(genpath(currentFolder));
exp_dir = [currentFolder filesep];
aud_dir = [exp_dir 'audios' filesep];
fig_dir = [exp_dir 'figs' filesep];
LUT_dir = [exp_dir 'gain_LUT' filesep];
%--------------------------------------------------------------------------
%%                    Create clean and noise file paths
%--------------------------------------------------------------------------
clean_files = dir([aud_dir 'clean' filesep '*.wav']); % path to clean speech test files.
noise_files = dir([aud_dir 'noise' filesep '*.wav']); % path to noise test files.
%--------------------------------------------------------------------------
%%                           Define Globals
%--------------------------------------------------------------------------
global Cl Cg min_db max_db Gmag
min_db =-10; max_db = 35;
%--------------------------------------------------------------------------
%%                       Load Look Up Tables (LUT)
%--------------------------------------------------------------------------
% don't have to load every time; comment out 'clear all;' and this section
%--------------------------------------------------------------------------
fprintf('Loading Gain Look Up Table...\n');
savefile = strcat(LUT_dir,'CSA_gain_LUT.mat'); load(savefile);
%--------------------------------------------------------------------------
%%                      Gain function Interpolation
%--------------------------------------------------------------------------
fprintf('Interpolating gain functions...\n');
[~,Gmag,~]=Tabulate_gain_functions(1,sqrt(2));
%--------------------------------------------------------------------------
Cl = griddedInterpolant({xi_CSA,(gammaK_CSA')},(CSA_l_gain),'spline','nearest');
Cg = griddedInterpolant({xi_CSA,(gammaK_CSA')},(CSA_g_gain),'spline','nearest');
% --------------------------------------------------------------------------
%%          Define MMSE estimator types/parameters & functions
%--------------------------------------------------------------------------
% CSC MMSE : Cosine Spectral Coefficients/Components MMSE estimator
%
% CSA/FSA MMSE:
%   short-time (Cosine/Fourier) Spectral Amplitude MMSE estimator
%
% CSA/FSA MMSE + SPU:
%   short-time (Cosine/Fourier) Spectral Amplitude MMSE estimator
%   with Speech Presence Uncertainty (SPU).
% LBLG :  Linear bilateral gain estimator (DCT)
% LBLG normal : or dual gain Wiener filter (DGW)
% NBLB :  Non-linear bilateral gain estimator (DCT)
%--------------------------------------------------------------------------
MMSE_types = {'CSA MMSE','FSA MMSE','CSA MMSE+SPU','FSA MMSE+SPU',...
    'LBG','NBLG','CSC MMSE','CSC MMSE + SPU'};
%--------------------------------------------------------------------------
%           (1)        (2)          (3)            (4)
myfun = {@CSA_MMSE_t, @FSA_MMSE_t, @CSA_MMSE_SPU_t, @FSA_MMSE_SPU_t, ...
    @LBG_t, @NBLG_t, @CSC_MMSE_t, @CSC_MMSE_SPU_t} ;
%          (5)    (6)     (7)          (8)
%--------------------------------------------------------------------------
% CSA normal/laplace/gamma:
%   short-time (Cosine/Fourier) Spectral Amplitude MMSE estimator
%   with Gaussian/Laplacian/Gamma speech priori and Gaussian noise priori.
%--------------------------------------------------------------------------
% methods = {'Unprocessed','Wiener','CSC laplace',...
%             'DGW','LBLG laplace',...
%             'NBLG normal','NBLG laplace',...
%             'CSA normal','CSA laplace', 'SEE RANK CSA gamma',...
%             'FSA normal','FSA laplace'};
%--------------------------------------------------------------------------
%%                  Define experiment parameter index
plot_spect = true;
spu =  'spu_'; % uncomment this for SPU weighting 
% spu = ''; % uncomment this 
if ~strcmp(spu,'spu_')
    %     methods = {'Unprocessed','FSA normal','FSA laplace',...
    %         'Wiener','CSC laplace','DGW','NBLG',...
    %         'CSA normal','CSA laplace','CSA gamma'};
    %     para_idx = [[0,0,0];[2,1,2];[2,2,2];...
    %         [7,1,7];[7,2,7];[5,1,5];[6,2,6];...
    %         [1,1,1];[1,2,1];[1,3,1]];
    %
    %     methods = {'Unprocessed','FSA normal','FSA laplace',...
    %         'CSC laplace','DGW','NBLG',...
    %         'CSA normal','CSA laplace','CSA gamma'};
    methods = {'Unprocessed','$G_{EM}$','$G_{L-FSA}$',...
        '$G_{L-CSC}$','DGW','NBLG',...
        '$G_N$','$G_L$','$G_G$'};
    para_idx = [[0,0,0];[2,1,2];[2,2,2];...
        [7,2,7];[5,1,5];[6,2,6];...
        [1,1,1];[1,2,1];[1,3,1]];
else
    methods = {'Unprocessed', 'FSA normal+SPU',...
        'Wiener+SPU', 'CSC laplace+SPU',...
        'CSA normal+SPU','CSA laplace+SPU','CSA gamma+SPU'};
    para_idx = [[0,0,0];[4,1,4];...
        [8,1,8];[8,2,8];...
        [3,1,3];[3,2,3];[3,3,3]];
end
%--------------------------------------------------------------------------
mmse_inputs.xi_min = 10^(-25/10);
mmse_inputs.gain_min = 0.1;
%--------------------------------------------------------------------------
%%             Define speech/noise variance estimation method
%--------------------------------------------------------------------------
% xi_type = {'noise_ideal'}; %ideal noise psd estimate

xi_type = {'DD'}; % Decision-direct a priori SNR estimation
%--------------------------------------------------------------------------
%%                 Define speech quality measure types
%--------------------------------------------------------------------------
measure = {'pesq','stoi','seg','seg_S','seg_N','sti'};
m_idx = [1,2,3]; % measure index
%--------------------------------------------------------------------------
%%                  Define SNR range and noise types
%--------------------------------------------------------------------------
% plt.SNR_arr = 0:5:15; % SNR range
plt.SNR_arr = 5; % SNR range
x = plt.SNR_arr;

noise_type = [19];
% noise types--------------------------------------------------------------
% 2  : 'Pink noise';
% 3  : 'white';
% 6  : 'Speech noise';
% 19 : 'Voice Babble';
% 20 : 'F-16 two-seat';
% 21 : 'Car Factory electrical welding';
% 22 : 'Car Factory production hall';
% 23 : 'Car Volvo-340 asphalt road';
%--------------------------------------------------------------------------
%%              Set signal analysis/synthesis parameters
%--------------------------------------------------------------------------
Fs = 16e3;
synWinType = {'hamming'};
% synWinType = {'hamming','modified'};
speech_priori = {'norm','lap','gamma'};
% frame_len = [20,20,32]; % duration in time (ms)
% shift_len = [10,5,8];  % 75% over-lapping
frame_len = [20]; % duration in time (ms)
shift_len = [5];  % 50% over-lapping
Tw_arr= round(Fs.*frame_len./1000);
%--------------------------------------------------------------------------
%%             Set the parameters for save the results
%--------------------------------------------------------------------------
plt.methods = methods;
plt.pesq_mode ='wb';
plt.save_plots = false; % save the plots
plt.save_data = false; % save the outputs
plt.plotType = 'bar';
plt.save_dir = fig_dir;

% plt.xTicks = {'0dB','5dB','10dB','15dB'};

plt.xTicks = {'5 dB'};

plt.unit = 'inch';
plt.width= 7;
plt.height= plt.width/((1+sqrt(5))/2);
plt.fontType = 'Times';
plt.fontSize = 10;
plt.legendLoc = 'northeastoutside';
plt.STI_method = 'ansi'; % sti method (i.e. 'ansi' or 'payton' or 'drullman')
plt.improv = true;
%--------------------------------------------------------------------------
MM = length(methods);RR = length(plt.SNR_arr);WW = length(synWinType);
NN = length(noise_type);SS = length(shift_len);XX = length(xi_type);
%--------------------------------------------------------------------------
latex.horiz =  plt.xTicks;
latex.vert = methods(2:end);
latex.hline =  [0,1,NaN];
latex.vline =  [1];
latex.format = '%.3f';
%--------------------------------------------------------------------------
%%                      Start to process speech files
%--------------------------------------------------------------------------
for ww = 1:WW %length(synWinType)
    mmse_inputs.synWinType = synWinType{ww};
    plt.synWinType = synWinType{ww};
    %---------------------------------------------------------------
    for ss = 1: SS % length(shift_len)
        
        mmse_inputs.frame_len = frame_len(ss);
        mmse_inputs.shift_len = shift_len(ss);
        plt.frame_len = frame_len(ss);
        plt.shift_len = shift_len(ss);
        Tw = Tw_arr(ss);
        %---------------------------------------------------------------
        for n = 1:XX % length(xi_type)
            
            mmse_inputs.xi_type = xi_type{n};
            
            plt.xi_type = xi_type{n};
            
            plt.save_note = [mmse_inputs.xi_type,'_',mmse_inputs.synWinType,...
                '_Syn_',num2str(mmse_inputs.shift_len),'ms_Shif']; % differentiate the tests
            
            output(ss).xi_type =  plt.xi_type;
            %---------------------------------------------------------------
            for nn = 1:NN % length(noise_type)
                %% Read noise file
                noise_no = noise_type(nn);
                
                [noise_init, ~] = audioread([noise_files.folder,...
                    '/', noise_files.name]); % noise waveform.
                
                plt.noise_name = get_noise_name(noise_no);
                output(ss).noise(nn).noise_name = plt.noise_name;
                %---------------------------------------------------------
                
                %% Read clean speech signal
                [clean, Fs] = audioread([clean_files(1).folder, filesep, clean_files(1).name]); % clean waveform.
                %                     clean(clean == 0) = eps;
                
                output(ss).noise(nn).audio.clean = clean;
                
                %% Create noisy signal according to global SNR
                for rr = 1:RR
                    [noisy,noise]= get_noisy(clean,noise_init,plt.SNR_arr(rr),Tw,'global');
                    
                    mmse_inputs.noisy= noisy;
                    mmse_inputs.noise = noise;
                    mmse_inputs.clean = clean;
                    mmse_inputs.Fs = Fs;
                    
                    S = length(noisy);
                    
                    output(ss).noise(nn).audio.noisy(rr,:) = noisy;
                    
                    clean = clean./norm(clean);
                    noisy = noisy./norm(noise);
                    
                    est_file = ['audios' filesep clean_files(1).name(1:end-4),'_',...
                        'clean.wav'];
                    
                    est_file = regexprep(est_file, ' ', '_');
                    audiowrite([exp_dir est_file],clean./max(abs(clean)),Fs,'BitsPerSample',16);
                    
                    est_file = ['audios' filesep clean_files(1).name(1:end-4),'_',...
                        'noisy','_',plt.noise_name,'_',...
                        num2str(plt.SNR_arr(rr)), 'dB.wav'];
                    
                    est_file = regexprep(est_file, ' ', '_');
                    audiowrite([exp_dir est_file],noisy./max(abs(noisy)),Fs,'BitsPerSample',16);
                    
                end
                
                %% ========================== Apply MMSE estimator ==========================
                for m = 2:MM
                    %---------------------------------------------------------------
                    output(ss).noise(nn).mmse(m-1).name = MMSE_types{para_idx(m,1)};
                    
                    mmse_inputs.speech_priori = speech_priori{para_idx(m,2)};
                    
                    output(ss).noise(nn).mmse(m-1).priori = mmse_inputs.speech_priori;
                    
                    plt.priori{m} = mmse_inputs.speech_priori;
                    
                    [mod_noisy,mod_clean,mod_noise] = myfun{para_idx(m,3)}(mmse_inputs);
                    
                    output(ss).noise(nn).audio.mod_sig(rr,m-1,:) = mod_noisy;
                    
                    S = length(mod_noisy);
                    
                    est_file = ['audios' filesep MMSE_types{para_idx(m,1)},'_',...
                        mmse_inputs.speech_priori,'_', plt.xi_type,'_',plt.noise_name,'_',...
                        num2str(plt.SNR_arr(rr)), 'dB.wav'];
                 
                    est_file = regexprep(est_file, ' ', '_');
                    audiowrite([exp_dir est_file],mod_noisy./max(abs(mod_noisy)),Fs,'BitsPerSample',16);                    
                    
                end
            end
        end
        
    end
end
%--------------------------------------------------------------------------
%                     Plot spectrogram of outputs
%--------------------------------------------------------------------------
if plot_spect == 1
    a = {'(a) ','(b) ','(c) ','(d) ','(e) ','(f) ','(g) ','(h) ','(i) ','(j) ','(k) ','(l) ','(m) ','(n) ','(o) ','(p) '};
    
    sl = length(output.noise.mmse)+2;
    si = ceil(sl/2);
    
    figure;
    p = panel();
    p.pack(si, 2);
    
    p(1,1).select();
    p.fontname = 'Times Roman';
    p.fontsize = 12;
    p.fontweight = 'bold'; 
    
    myspectrogram(output.noise.audio(1).clean,Fs);
    
    title('(a) Clean speech','Interpreter','latex');
    set(gca,'XTick',[0 0.5 1 1.5 2 2.5], 'YTick',[0 4000 8000],'ylim',[0 8000]);
    ylabel('Freq.[kHz]','fontsize',6,'FontWeight','bold');
    
    set(gca,'xticklabel',{[]},'yticklabel',{'0' '4' '8'}) ;
    
    p(1,2).select();
    
    myspectrogram(output.noise.audio(1).noisy(1,:),Fs);
    
    title(['(b) Noisy speech, , ' plt.noise_name ' SNR = ' num2str(plt.SNR_arr(1)) ' dB'],'Interpreter','latex');
    set(gca,'XTick',[0 0.5 1 1.5 2 2.5], 'YTick',[0 4000 8000],'ylim',[0 8000]);
    ylabel('Freq.[kHz]');
    set(gca,'xticklabel',{[]},'yticklabel',{'0' '4' '8'});
    
    k = 3;
    for r = 2:si
        for c = 1:2
            p(r, c).select();
            
            myspectrogram(squeeze(output.noise.audio(1).mod_sig(1,k-2,:)),Fs);
            
            title([a{k} methods{k-1}],'Interpreter','latex');
            set(gca,'XTick',[0 0.5 1 1.5 2 2.5], 'YTick',[0 4000 8000],'ylim',[0 8000]);
            
            ylabel('Freq.[kHz]','fontsize',6,'FontWeight','bold');
            set(gca,'xticklabel',{[]},'yticklabel',{'0' '4' '8'}) ;
            
            if k == sl-1|| k== sl
                set(gca,'xticklabel',{[0 0.5 1 1.5 2 2.5]},'yticklabel',{[]}) ;
                xlabel('Time[s]');
            end
            
            k=k+1;
            
            if k == sl+1
                break
            end
        end
    end
    %     end
    save_name = strcat('spect_MMSE_compare_', plt.noise_name,'_',num2str(plt.SNR_arr(1)),'_dB','.eps');
    
    print('-depsc', '-r600',sprintf(strcat(plt.save_dir,save_name)));
    
end
%% ============================== Functions ===============================
function noise_t = get_noise_name(noise_no)
switch noise_no
    case 1
        noise_t = 'Sinusoid';
    case 2
        noise_t = 'Pink noise';
    case 3
        noise_t = 'White noise';
    case 4
        noise_t = 'White noise_n6dB';
    case 5
        noise_t = 'White noise_n12dB';
    case 6
        noise_t = 'Speech noise';
    case 7
        noise_t = 'M 109';
    case 8
        noise_t = 'Buccaneer 190 Knots 1000 Feet';
    case 9
        noise_t = 'Leopard 2';
    case 10
        noise_t = 'Wheel carrier';
    case 11
        noise_t = 'Buccaneer 450 Knots 300 Feet';
    case 12
        noise_t = 'Lynx';
    case 13
        noise_t = 'Leopard 1';
    case 14
        noise_t = 'Operation room of destroyer';
    case 15
        noise_t = 'Engine room of destroyer';
    case 16
        noise_t = 'Machine gun repeated';
    case 17
        noise_t = 'HF radio channel';
    case 18
        noise_t =  'STITEL test signal';
    case 19
        noise_t = 'Voice babble';
    case 20
        noise_t = 'F-16 two-seat';
    case 21
        noise_t = 'Car Factory electrical welding';
    case 22
        noise_t = 'Car Factory production hall';
    case 23
        noise_t = 'Car Volvo-340 asphalt road';
    case 24
        noise_t = 'Car Volvo-340 brick road';
    otherwise
        noise_t = 'unknown';
end
end
%=========================================================================
%--------------------------------E.O.F.------------------------------------
%==========================================================================