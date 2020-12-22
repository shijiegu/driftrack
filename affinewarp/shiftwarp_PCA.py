import numpy as np
import numba
from tqdm import trange
import scipy as sci
from sklearn.utils.validation import check_is_fitted
import scipy.optimize

from .spikedata import SpikeData
from .utils import check_dimensions
from ._optimizers import _diff_gramian, PoissonObjective

from .bmat import nnls_solveh_banded

class ShiftWarping_PCA(object):
    """
    Time-Shifted PCA in Williams 2020 (a PCA version of Tensor version as in https://www.biorxiv.org/content/10.1101/2020.03.02.974014v2.full)
    Models translations + PCA in time across a collection of multi-dimensional time
    series. Does not model stretching or compression of time across trials.

    Attributes
    ----------
    shifts : ndarray
        Number of time bins shifted on each trial.
    fractional_shifts : ndarray
        Time shifts expressed as a fraction of trial length.
    loss_hist : list
        History of objective function over optimization.
    """

    def __init__(self, maxlag=.5, warp_reg_scale=0, smoothness_reg_scale=0,
                 l2_reg_scale=1e-7, loss='quadratic', center_shifts=False,
                 nonneg=False,num_components=5):
        """Initializes ShiftWarping object with hyperparameters.

        Parameters
        ----------
        maxlag : float
            Maximal allowable shift.
        warp_reg_scale : float
            Penalty strength on the magnitude of the shifts.
        smoothness_reg_scale : int or float
            Penalty strength on L2 norm of second temporal derivatives of the
            warping templates.
        l2_reg_scale : int or float
            Penalty strength on L2 norm of the warping template.
        loss : str
            Specifies loss function, either 'quadratic' or 'poisson'
        center_shifts : bool
            If True, shifts are mean-centered on each iteration. Defaults
            to False.
        nonneg : bool
            If True, forces template to be nonnegative. Defaults to False.
        """

        if (maxlag < 0) or (maxlag > .5):
            raise ValueError('Value of maxlag must be between 0 and 0.5')

        self.maxlag = maxlag
        self.loss_hist = []
        self.warp_reg_scale = warp_reg_scale
        self.smoothness_reg_scale = smoothness_reg_scale
        self.l2_reg_scale = l2_reg_scale
        self.loss = loss
        self.center_shifts = center_shifts
        self.nonneg = nonneg
        self.num_components= num_components

        if loss == 'quadratic':
            self._shifted_loss = _compute_shifted_quad_loss
            self._eval_loss = _eval_quad_loss
        elif loss == 'poisson':
            self._shifted_loss = _compute_shifted_poiss_loss
            self._eval_loss = _eval_poiss_loss
        else:
            raise ValueError(
                "'loss' parameter should be one of ('quadratic', 'poisson').")

    def fit(self, data, iterations=20, verbose=True, warp_iterations=None,
            trial_idx=slice(None), neuron_idx=slice(None)):
        """
        Fit shift warping to data.

        Parameters
        ----------
        data : ndarray
            Collection of time series with shape (num_trials, num_timepoints,
            num_units).
        iterations : int
            Number of iterations to run.
        verbose : bool
            If True, prints progress bar. Otherwise, no output is displayed.
        warp_iterations : int or None
            This is ignored, used to provide a consistent API with
            PiecewiseWarping.
        trial_idx : indices
            Indices along first dimension of `data` array. If provided, the
            warping template is only fit to these trials. By default, all
            trials are used for fitting.
        neuron_idx : indices
            Indices along last dimension of `data` array. If provided, the
            warping functions are only fit to these units. By default, all
            units are used for fitting.
        """

        # Check input.
        if isinstance(data, SpikeData):
            raise ValueError("Expected input 'data' to be a 3d numpy "
                             "holding binned spike counts.")

        # data dimensions:
        #   K = number of trials
        #   T = number of timepoints
        #   N = number of features/units
        #   L = number of time lags/shifts to search over in each direction.
        K, T, N = data.shape
        L = int(self.maxlag * T)

        # Warps are penalized based on distance from identity. These penalties
        # can be pre-computed up front.
        self._warp_penalty = self.warp_reg_scale * T * N * \
            np.abs(np.linspace(-self.maxlag, self.maxlag, 2 * L + 1)[None, :])

        # Initialize shifts and model template.
        self.shifts = np.zeros(K, dtype=int)
        self.amps = np.zeros((K,self.num_components))
        self.template = np.zeros((T,N,self.num_components))
        self._fit_template(
            data[trial_idx, :, :], self.shifts[trial_idx],self.amps[trial_idx,:])

        # Initialize loss history
        if self.loss == "quadratic":
            resid = self.template[None, :, :] - data
            self.loss_hist = [np.mean(resid ** 2)]
        else:
            z1 = np.exp(self.template)[None, :, :]
            z2 = data * self.template[None, :, :]
            self.loss_hist = [np.mean(z1 - z2)]

        # progress bar
        pbar = trange(iterations) if verbose else range(iterations)

        # main loop
        for i in pbar:
            # Update parameters and compute loss.
            self._fit_warps(data[:, :, neuron_idx])
            self._fit_template(
                data[trial_idx, :, :], self.shifts[trial_idx], self.amps[trial_idx])
            self._record_loss(data)

            # update loss display
            if verbose:
                pbar.set_description(
                    'Loss: {0:.2f}'.format(self.loss_hist[-1]))

        # compute shifts as a fraction of trial length
        self.fractional_shifts = self.shifts / T
        self._losses = None

    def _fit_warps(self, data):
        """Updates shift parameters."""

        # Data dimensions.
        K, T, N = data.shape
        L = int(self.maxlag * T)

        # Compute reconstruction errors for each shift.
        losses = np.zeros((K, 2 * L + 1))
        self._shifted_loss(data, self.template, losses)

        # Compute total objective and find optimal shifts.
        obj = losses + self._warp_penalty
        self.shifts = np.argmin(obj, axis=1) - L

        if self.center_shifts:
            self.shifts = \
                self.shifts - round(float(np.mean(self.shifts)))

    def _fit_template_r(self, r, data, shifts, amps):
        """Updates template firing rates."""

        # Data dimensions.
        K, T, N = data.shape
        assert len(shifts) == K

        if self.loss == 'quadratic':
            DtD = _diff_gramian(T, self.smoothness_reg_scale * K, self.l2_reg_scale * K)
            WtW = np.zeros((3, T))
            WtX = np.zeros((T, N))
            _fill_WtW(shifts, amps**2, WtW[-1])
            res=np.squeeze(data)-self.predict_wo_r(r)
            _fill_WtX(res, shifts, amps, WtX)
            if self.nonneg:
                self.template[:,:,r] = nnls_solveh_banded((WtW + DtD), WtX, self.template)
            else:
                self.template[:,:,r] = sci.linalg.solveh_banded((WtW + DtD), WtX)

        elif self.loss == 'poisson':
            # Get initial parameters.
            try:
                x0 = self.template.ravel()
            except AttributeError:
                x0 = np.zeros(T * N)

            # Set up optimization problem.
            obj = PoissonObjective(data, self.smoothness_reg_scale,
                                   self.l2_reg_scale, shifts=shifts)

            # Fit template.
            opt = scipy.optimize.minimize(obj, x0,
                                          jac=True,  # hessp=obj.hessp,
                                          method='L-BFGS-B',  # 'newton-cg',
                                          options=dict(maxiter=4))
            self.template = (opt.x).reshape(T, N)

    def _fit_template(self, data, shifts, amps):
        for r in range(self.num_components):
            self._fit_template_r(r, data, shifts, amps[:,r])

    def _fit_amps(self, r, data, shifts, template):
        # Data dimensions.
        K, T, N = data.shape
        assert len(shifts) == K

        if self.loss == 'quadratic':
            for r in range(self.num_components):
                DtD = _diff_gramian(T, self.smoothness_reg_scale * K, self.l2_reg_scale * K)
                WtX = np.zeros((T, N))
                WtW_t = np.zeros((T,K))
                _fill_WtW_t(shifts, WtW_t)
                template_shifted = np.zeros((K, T,1))
                _predict(template[:,r], self.shifts, template_shifted)

                res=np.squeeze(data)-self.predict_wo_r(r)
                for k in range(K):
                    amp_nominator=np.dot(res[k,:].reshape((-1,1)),template_shifted[k,:,1].reshape((-1,1)))
                    amp_denom=template[:,r].T*np.diag(WtW_t[:,k])*template[:,r]                
                    self.amps[k,r]=amp_nominator/amp_denom

        else:
            raise ValueError(
                "Have not implemented thgis yet.")


    def _record_loss(self, data):
        """Computes loss on all data."""
        loss = self._eval_loss(data, self.template, self.shifts)
        self.loss_hist.append(loss)

    def argsort_warps(self):
        """
        Returns an ordering of the trials based on the learned shifts.
        """
        self.assert_fitted()
        return np.argsort(self.shifts)

    def predict(self):
        """
        Returns model prediction (warped version of template on each trial).
        """
        self.assert_fitted()

        # Allocate space for prediction.
        K = len(self.shifts)
        T, N, R = self.template.shape
        pred = np.zeros((K, T, R))

        # Compute prediction in JIT-compiled function.
        _predict(self.template, self.shifts, pred)
        return np.sum(pred,2)

    def predict_wo_r(self,r):
        """
        Returns model prediction without r component (warped version of template on each trial).
        """

        # Allocate space for prediction.
        K = len(self.shifts)
        T, N, R = self.template.shape
        pred = np.zeros((K, T, R-1))

        all_components=np.arange(self.num_components)
        ind=np.setdiff1d(all_components,r)

        # Compute prediction in JIT-compiled function.
        _predict(self.template[:,:,ind], self.shifts, pred)
        print(np.shape(pred))

        return np.sum(pred*np.expand_dims(self.amps[:,ind],axis=1),2)

    def transform(self, data):
        """
        Applies inverse warping functions to align raw data across trials.
        """
        self.assert_fitted()
        data, is_spikes = check_dimensions(self, data)

        # For SpikeData objects.
        if is_spikes:
            return data.shift_each_trial_by_fraction(self.fractional_shifts)

        # For dense data (trials x timebins x units).
        else:
            # warp dense data
            out = np.empty_like(data)
            _warp_data(data, self.shifts, out)
            return out

    def event_transform(self, trials, frac_times):
        """
        Time warp events by applying inverse warping functions.

        Parameters
        ----------
        trials : array-like
            Vector of ints holding trial indices for each event.
        frac_times : array-like
            Vector of floats holding fractional times for each event
            (``frac_times[i] == 0`` means that event ``i`` occured at
            trial start; ``frac_times[i] == 1`` means that event ``i``
            occured at trial end).

        Returns
        -------
        aligned_times : array-like
            Transformed fractional event times.
        """
        self.assert_fitted()
        return np.asarray(frac_times)[trials] - self.fractional_shifts[trials]

    def copy_fit(self, model):
        """
        Applies inverse warping functions to align raw data across trials.
        """
        if not isinstance(model, ShiftWarping):
            raise ValueError("Can only copy another ShiftWarping instance.")
        model.assert_fitted()
        self.template = model.template.copy()
        self.shifts = model.shifts
        self.fractional_shifts = model.fractional_shifts
        return self

    def assert_fitted(self):
        check_is_fitted(self, 'shifts')


#@numba.jit(nopython=True)
def _predict(template, shifts, out):
    """
    Produces model prediction. Applies shifts to template on each trial.

    Parameters
    ----------
    template : array_like
        Warping template, shape: (time x features)
    shifts : array_like
        Learned shifts on each trial, shape: (trials)
    out : array_like
        Storage for model prediction, shape: (trials x time x features)

    Notes
    -----
    This private function is just-in-time compiled by numba. See
    ShiftWarping.predict(...) wraps this function.
    """
    K = len(shifts)
    T, N, R = template.shape
    for r in range(R):
        template_=template[:,:,r]
        for k in range(K):
            i = -shifts[k]
            t = 0
            while i < 0:
                out[k, t, r] = template_[0]
                t += 1
                i += 1
            while (i < T) and (t < T):
                out[k, t, r] = template_[i]
                t += 1
                i += 1
            while t < T:
                out[k, t, r] = template_[-1]
                t += 1


@numba.jit(nopython=True)
def _fill_WtW(shifts, amps, out):
    T = len(out)
    for s,a in zip(shifts,amps):
        if s < 0:
            out[-1] += a*(1 - s)
            for i in range(2, T + s + 1):
                out[-i] += a
        else:
            out[0] += a*(1 + s)
            for i in range(1, T - s):
                out[i] += a

def _fill_WtW_t(shifts, out):
    T = len(out)
    for t in range(np.shape(shifts)[0]):
        s=shifts[t]
        if s < 0:
            out[-1,t] += (1 - s)
            for i in range(2, T + s + 1):
                out[-i,t] += 1
        else:
            out[0,t] += (1 + s)
            for i in range(1, T - s):
                out[i,t] += 1


@numba.jit(nopython=True)
def _fill_WtX(data, shifts, amps, out):
    K, T, N = data.shape
    for k in range(K):
        i = -shifts[k]
        t = 0

        for t in range(T):
            if i < 0:
                out[0] += data[k, t]*amps[k]
            elif i >= T:
                out[-1] += data[k, t]*amps[k]
            else:
                out[i] += data[k, t]*amps[k]
            i += 1


@numba.jit(nopython=True)
def _warp_data(data, shifts, out):
    K, T, N = data.shape
    for k in range(K):
        i = shifts[k]
        t = 0
        while i < 0:
            out[k, t] = data[k, 0]
            t += 1
            i += 1
        while (i < T) and (t < T):
            out[k, t] = data[k, i]
            t += 1
            i += 1
        while t < T:
            out[k, t] = data[k, -1]
            t += 1


@numba.jit(nopython=True)
def _eval_quad_loss(data, template, shifts):
    K, T, N = data.shape
    total_loss = 0.0

    for k, s in zip(range(K), shifts):
        for t in range(T):

            # shifted index
            i = t - s
            if i < 0:
                i = 0
            elif i >= T:
                i = T-1

            # add loss for each neuron
            for n in range(N):
                total_loss += ((data[k, t, n] - template[i, n]) ** 2)

    return total_loss / data.size


@numba.jit(nopython=True)
def _eval_poiss_loss(data, template, shifts):
    K, T, N = data.shape
    total_loss = 0.0
    exp_template = np.exp(template)

    for k, s in zip(range(K), shifts):
        for t in range(T):

            # shifted index
            i = t - s
            if i < 0:
                i = 0
            elif i >= T:
                i = T-1

            # add loss for each neuron
            for n in range(N):
                total_loss += (exp_template[i, n] - template[i, n] * data[k, t, n])

    return total_loss / data.size


@numba.jit(nopython=True, parallel=True)
def _compute_shifted_quad_loss(data, template, losses):

    K, T, N = data.shape
    L = losses.shape[1] // 2

    for k in numba.prange(K):
        for t in range(T):
            for l in range(-L, L+1):

                # shifted index
                i = t - l
                if i < 0:
                    i = 0
                elif i >= T:
                    i = T-1

                # quadratic loss
                for n in range(N):
                    losses[k, l+L] += ((data[k, t, n] - template[i, n]) ** 2)


@numba.jit(nopython=True, parallel=True)
def _compute_shifted_poiss_loss(data, template, losses):

    K, T, N = data.shape
    L = losses.shape[1] // 2

    exp_template = np.exp(template)

    for k in numba.prange(K):
        for t in range(T):
            for l in range(-L, L+1):

                # shifted index
                i = t - l
                if i < 0:
                    i = 0
                elif i >= T:
                    i = T-1

                # poisson loss
                for n in range(N):
                    losses[k, l+L] += (exp_template[i, n] - template[i, n] * data[k, t, n])
