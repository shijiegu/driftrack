## Contents

Similar to the parent repository affinewarp by Alex H. Williams, these jupyter notebooks demonstrate the basic functionality of this package and reproduce the essential figures in our manuscript in the repository.

0) [**`Fig0_shift_warp_PCA.ipynb`**](/examples/Fig4_shift_warp_PCA.ipynb) - Implements the 2D matrix version of "Time-Shifted Tensor Decomposition" as in this [preprint](https://www.biorxiv.org/content/10.1101/2020.03.02.974014v2.full.pdf). After applying to real data, we found that the second PC is significantly smaller than the first, and so we directly enforce rank=1 in the method development.
1) [**`Fig1_Landscape.ipynb`**](/examples/Fig1_Landscape.ipynb) - demonstrates the objective value landscape. Its complexity confirmed that we should only use brute force greedy search to find the solution of warping to the problem.
2) [**`Fig2_SimulateData.ipynb`**](/examples/Fig2_SimulateData.ipynb) - simulates data, and demonstrates the full range of drift correction models (shift-only, linear, piecewise linear).
3) [**`Fig3_RealData.ipynb`**](/examples/Fig3_RealData.ipynb/) - demonstrates drift correction on real data, with time averaging window 1 second (data provided by Nick Steinmetz in their Neuropixel [publication](https://www.biorxiv.org/content/10.1101/2020.10.27.358291v1). Note that part of this notebook needs to work with Kilosort MATLAB code (Please download Kilosort MATLAB repo (https://github.com/MouseLand/Kilosort) and then add the code from the folder [**`/examples/kilosort_code_forFig3/`**](/examples/kilosort_code_forFig3/).
4) [**`Fig3_RealData_tensofasecond`**](/examples/Fig3_RealData_tensofasecond/) - Similar to [**`Fig3_RealData.ipynb`**], but the time average window is 0.1 second.
