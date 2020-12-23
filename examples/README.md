## Contents

Similar to the parent repository affinewarp by Alex H. Williams, these jupyter notebooks demonstrate the basic functionality of this package and reproduce the essential figures in our manuscript in the repository.

1) [**`Fig1_Landscape.ipynb`**](/examples/Fig1_Landscape.ipynb) - demonstrates the objective value landscape. Its complexity confirmed that we should only use brute force greedy search to find the solution of warping to the problem.
2) [**`Fig2_SimulateData.ipynb`**](/examples/Fig2_SimulateData.ipynb) - simulates data, and demonstrates the full range of drift correction models (shift-only, linear, piecewise linear).
3) [**`Fig3_RealData.ipynb`**](/examples/Fig3_RealData.ipynb/) - demonstrates drift correction on real data, with time averaging window 1 second (data provided by Nick Steinmetz (https://www.biorxiv.org/content/10.1101/2020.10.27.358291v1)). Note that part of this notebook needs to work with Kilosort MATLAB code (Please download Kilosort MATLAB repo (https://github.com/MouseLand/Kilosort) and then add the code from the folder [**`/examples/kilosort_code_forFig3/`**](/examples/kilosort_code_forFig3/).
4) [**`Fig3_RealData_tensofasecond`**](/examples/Fig3_RealData_tensofasecond/) - demonstrates drift correction on real data, with time averaging window 0.1 second (data provided by Nick Steinmetz (https://www.biorxiv.org/content/10.1101/2020.10.27.358291v1)). Note that part of this notebook needs to work with MATLAB.
5) [**`Fig4_shift_warp_PCA.ipynb`**](/examples/Fig4_shift_warp_PCA.ipynb) - demonstrates drift correction on real data (data provided by Nick Steinmetz (https://www.biorxiv.org/content/10.1101/2020.10.27.358291v1)).
