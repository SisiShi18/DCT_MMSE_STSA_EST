# DCT_MMSE_STSA_EST

Here contains the scripts for implementing the DCT-based MMSE spectral amplitude estimators.
Output of the program:
1. Spectrograms of enhanced speech using different methods, saved in the 'figs' folder
2. Enhanced speech wave files saved in the 'audios' folder.

The following methods are included:

FSA: DFT short-time Spectral Amplidue Estimator
CSA: DCT short-time Spectral Amplidue Estimator (Proposed)
LBLG: Linear bilateral gain estimator (DCT)
LBLG normal : or dual gain Wiener filter (DGW)
NBLB : Non-linear bilateral gain estimator (DCT)

'FSA normal', [1],

'FSA laplace', [2],
            
'CSA normal', proposed ,

'CSA laplace', proposed ,

'CSA gamma', proposed ,

'FSA normal+SPU', [1] ,

'CSA normal+SPU', proposed , 

'CSA laplace+SPU', proposed ,

'CSA gamma+SPU', proposed ,

'LBLG normal', [3] ,

'LBLG laplace', [5] ,

'NBLG normal', [4] ,

'NBLG laplace', [5] ,


Two tests:

1.  The noise psd is estimated as the power of short-time amplitude : |N|^2 
    
2.  The noise psd is estimated by using the noise tracker given in [6]

In both tests, the 'a priori SNR' is estimate using the 'Decision Direct' method in [1]


[1] 'Speech enhancement using a minimum-mean square error short-time spectral amplitude estimator' , 1984

[2] 'Minimum mean-square error estimation of discrete Fourier coefficients with generalized Gamma priors' , 2007

[3] 'Low distortion speech enhancement' , 2000

[4] 'MMSE estimator for speech enhancement considering the constructive and destructive interference of noise' , 2010

[5] 'Low-distortion MMSE speech enhancement estimator based on laplacian prior' , 2017

[6] 'Unbiased MMSE-based noise power estimation with low complexity and low tracking delay' , 2011

