# DCT_MMSE_STSA_EST

Here contains the scripts for implementing the DCT-based MMSE spectral amplitude estimators.
Output of the program:
1. Spectrograms of enhanced speech using different methods, saved in the 'figs' folder
2. Enhanced speech wave files saved in the 'audios' folder.

The following methods are included:

FSA: DFT short-time Spectral Amplidue Estimator

CSA: DCT short-time Spectral Amplidue Estimator (Proposed)

CSC : DCT short-time Spectral Coefficient Estimator

CSC normal : also called Wiener filter, using Guassian PDF

LBLG: Linear bilateral gain estimator (DCT)

LBLG normal : or dual gain Wiener filter, using Gaussian PDF (DGW)

NBLB : Non-linear bilateral gain estimator (DCT)

'FSA normal', [1],

'FSA laplacian', [2],
            
'CSA normal', proposed ,

'CSA laplacian', proposed ,

'CSA gamma', proposed ,

'FSA normal+SPU', [1] ,

'CSA normal+SPU', proposed , 

'CSA laplacian+SPU', proposed ,

'CSA gamma+SPU', proposed ,

'CSC laplacian', [7] ,

'CSC laplacian + SPU', [7] ,

'LBLG normal', [3] ,

'LBLG laplacian', [5] ,

'NBLG normal', [4] ,

'NBLG laplacian', [5] ,

'Wiener filter', (DCT), [8] ,

'Wiener filter + SPU' (DCT) [8] .

Two tests:

1.  The noise psd is estimated as the power of short-time amplitude : |N|^2 
    
2.  The noise psd is estimated by using the noise tracker given in [6]

In both tests, the 'a priori SNR' is estimate using the 'Decision Direct' method in [1]


[1] Y. Ephraim and D. Malah, “Speech enhancement using a minimum-mean square error short-time spectral amplitude estimator,” IEEE Transactions on acoustics, speech, and signal processing , vol. 32, no. 6, pp. 1109–1121, 1984

[2] J. S. Erkelens, R. C. Hendriks, R. Heusdens, and J. Jensen, “Minimum mean-square error estimation of discrete fourier coefficients with general-ized gamma priors,”IEEE Transactions on Audio, Speech, and Language Processing , vol. 15, no. 6, pp. 1741–1752, 2007

[3] I. Soon and S. Koh, “Low distortion speech enhancement,” IEEE Proceedings-Vision, Image and Signal Processing, vol. 147, no. 3, pp.247–253, 2000.

[4] T. Hasan and M. K. Hasan, “Mmse estimator for speech enhancement considering the constructive and destructive interference of noise,” IET signal processing, vol. 4, no. 1, pp. 1–11, 2010

[5] . M. Mahmmod, A. R. Ramli, S. H. Abdulhussian, S. A. R. Al-Haddad, and W. A. Jassim, “Low-distortion mmse speech enhancement estimator based on laplacian prior,” IEEE Access, vol. 5, pp. 9866–9881, 2017.

[6] T. Gerkmann and R. C. Hendriks, “Unbiased mmse-based noise power estimation with low complexity and low tracking delay,” IEEE Transactions on Audio, Speech, and Language Processing, vol. 20, no. 4, pp. 1383–1393, 2011

[7] Zou Xia, Zhang Xiongwei, "Speech enhancement using an MMSE short time DCT coefficients estimator with supergaussian speech modeling," Journal of Electronics, 2007

[8] Y. Soon, S. N. Koh, and C. K. Yeo, “Noisy speech enhancement using discrete cosine transform,” Speech communication, vol. 24, no. 3, pp. 249–257, 1998

