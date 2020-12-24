# Drift-tracking for NeuralPixel-a-like large scale ephys data

This repo contains code to do electrophysiological data drift-correction for large scale recordings such as that obtained by Neuralpixel Recordings. The method is similar to that in the recent [Neuralpixel 2.0 paper](https://www.biorxiv.org/content/10.1101/2020.03.02.974014v2.full), but we pursue a version based on 1-d time series with extensions in modeling abrupt global amplitude change in spiking data. 

For details on the problem as well as the methodology, please refer to the [manuscript](https://github.com/shijiegu/driftrack/blob/master/Drift_Tracking%20(2).pdf) included in this repository. The manuscipt demonstrates that our method achieves the same level of performance as the image registration method. This project was developed for UC Berkeley [EE290-001](https://people.eecs.berkeley.edu/~yima/courses/EE290-Fall2019/EE290-2019-syllabus.pdf) course project. The machinary of solving the problem largely comes from [Discovering precise temporal patterns in large-scale neural recordings through robust and interpretable time warping](https://doi.org/10.1016/j.neuron.2019.10.020), with minor changes.


<img width="708" alt="Screen Shot 2020-12-23 at 10 11 41 AM" src="https://user-images.githubusercontent.com/29357775/103011078-52bae600-4507-11eb-9bbe-58d44242a942.png">

## Getting started

After installing (see below), check out the demos in the [`examples/`](https://github.com/shijiegu/driftrack/tree/master/examples) folder.

Either download or clone the repo:

```
git clone https://github.com/shijiegu/driftrack/
```

Then navigate to the downloaded folder:

```
cd /path/to/driftrack
```

Install the package and requirements:

```
pip install .
pip install -r requirements.txt
```

You will need to repeat these steps if we update the code.


## Contact

shijiegu@berkeley.edu (or open an issue here).
