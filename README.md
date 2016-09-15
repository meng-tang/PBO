###Introduction ###

This software library implements the Parametric Pseudo-Bound Cuts algorithm in the 
Pseudo-Bound Optimization framework as is described in

	"Pseudo-Bound Optimization for Binary Energies."
	Meng Tang, Ismail Ben Ayed, Yuri Boykov.
	In European Conference on Computer Vision (ECCV), Zurich, Switzerland, September, 2014.

If you use this software for research purposes, please cite
the aforementioned paper in any resulting publication.
For any technical issue concerning the code, you can contact Meng Tang. (mtang73@csd.uwo.ca)
   
###How to compile the code ###

You can open solution file PBO.sln with Visual Studio 2010. Add where you store the openSEG 
library as an include directory. Part of the code also depends on Matlab library. So set include directory
and library directory as follows if necessary and compile under Release or Debug mode.

D:\Program Files (x86)\MATLAB\R2009a\extern\include (replace with your Matlab path)
D:\Program Files (x86)\MATLAB\R2009a\extern\lib\win32\microsoft (replace with your Matlab path)

###How to use the code for specific applications ###

The code can be used for the following applications discussed in the paper:

A. Entropy-based image foreground / background segmentation.
   (supervised with bounding box or Unsupervised)
B. Matching target distribution (KL divergence or Bhattacharyya coefficient minimization)
C. Nonsubmodular pairwise energy minimization 
   (curvature regularization, binary image deconvolution)
   
In main.cpp, uncomment function EntropyMinimize(argc, argv) for application A, matchdist(argc,argv) for B 
and qpprogrammain(argc, argv) for C. You also need to set variable string rootdir for each function.

For application A (entropy minimization), input command line is of the format:
pbo.exe + image name + color bin number per channel + weight of smoothness term
Example: pbo.exe 227092 16 10.0 
Three methods are implemented, BCD as is in GrabCut (histogram instead of GMM), PPBC with and without bounding box.
Resulting segmentation will be saved in 'images' directory. 
The console will show energies and running time of BCD and supervised / unsupervised pPBC.

For application B (matching distribution), input command line is of the format:
pbo.exe + image name + color bin number per channel + 'KL' or 'BHA' + weight of smoothness term + weight of distribution matching term
Example: pbo.exe 153093 32 KL 1.0 1000
Three algorithms are implemented, Fast Trust Region, Auxiliary Cut and pPBC (Parametric Pseudo-Bound Cuts 
as in PBO framework). Each algorithm gives segmentation saved in 'images' directory.
The console will show energies of segmentation obtained by each method.

For application C (nonsubmudolar pairwise), input command line is of the format:
pbo.exe + path to quadratic program mat file + image width + image height
The mat file contains three variables 'C' 'U' and 'M'. For variable s_p = 1 or 0, the energy is defined as,
E(S) = C + \sum_i{U_i*s_i} + \sum{pq\in \mathcal{N}}{M_pq*s_p*s_q}
'C' is a double constant. 'U' is a vector where each element U_i is the weight of unary term.
'M' is a matrix where each row gives a pairwise term (either positive or negative weight).
Each row contains three doubles. The first and second are indexes of variables s_p and s_q. The third is the weight M_pq for pairwise term between the two variables.
Note that indexes of pixels are column-wise, i.e., index = pixel_column_number*image_height + pixel_row_number.

Example of command: pbo.exe D:/repositories/PBO_public_code/pairwise/qpdata.mat 120 104
If your energy is not based on image, then set image width to 1 and image height to number of variables. 
Two specific applications of pairwise nonsubmodular energies are included:
C1: stenosis segmentation.
C2: binary image deconvolution.

To test C1, open runStenosis.m and modify exepath and qpdatapath. You should see resulting energies showing up like 54648.3 55191.2 56121.2 56766.4 ... (for different curvature weight)
To test C2, open runDeblur.m and modify exepath and qpdatapath. In our paper, we played with variable noiseSigma from 0.05, 0,1, 0.15 to 0.2. The matlab file randomly adds noise for 10 times
and run different algorithms. Average energy is shown.

qpproram.h implements three algorithms, fast trust region, auxiliary cut and Pseudo-bound optimization.
There are three versions of PBO, PPBC-B, L or T, as is described in the paper.
To switch mode, call function PPBCqp.setmode(mode) where mode can be PPBCT (recommended) or PPBCT or PPBCL.
Also set parameters of PPBC for different applications by calling PPBC.setpara(...)

CopyRight 2013-2014 Meng Tang


