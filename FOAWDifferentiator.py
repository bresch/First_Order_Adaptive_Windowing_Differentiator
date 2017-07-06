##################################################################################
# File: FOAWDifferentiator.py
# Brief :A class that implements a first order adaptive windowing differentiator 
# Author: Mathieu Bresciani <brescianimathieu@gmail.com>
# From:  Discrete-Time Adaptive Windowing for Velocity Estimation
# Farrokh Janabi-Sharifi, Vincent Hayward, and Chung-Shin J. Chen
##################################################################################

class FOAWDifferentiator:

    def __init__(self, dt, noise_level):
        self.set_noise_level(noise_level)
        self.set_sample_time(dt)
        self._max_window_size = 14
        self.reset()

    def set_noise_level(self, noise_level):
        self._delta = noise_level

    def set_sample_time(self, dt):
        self._dt = dt

    def get_last_window_size(self):
        return self._last_window_size

    def reset(self):
        self._buffer = [0]*(self._max_window_size+1)
        self._nb_samples = 0
        self._last_window_size = 0

    def add_sample(self, sample):
        if self._nb_samples <= self._max_window_size:
            self._nb_samples = self._nb_samples + 1
        else:
            self.shift_buffer()
            self._nb_samples = self._max_window_size + 1

        self._buffer[self._nb_samples-1] = sample

    def shift_buffer(self):
        for i in range(0,self._nb_samples-1):
            self._buffer[i] = self._buffer[i+1]

    def best_fit_FAOW(self, window_size):
        sum1 = 0.0
        sum2 = 0.0
        last_sample_pos = self._nb_samples - 1

        for i in range(0, window_size+1):
            sum1 = sum1 + self._buffer[last_sample_pos - i]
            sum2 = sum2 + self._buffer[last_sample_pos - i]*i

        y_mean = sum1/(window_size+1)
        sum1 = sum1*window_size
        sum2 = sum2*2
        
        den = self._dt*window_size*(window_size+1.0)*(window_size+2.0)/6.0
        self._a = (sum1 - sum2)/den
        self._b = y_mean + self._a * window_size * self._dt/2.0
        self._y0 = y_mean - self._b * window_size * self._dt/2.0

    def fit(self):

        last_sample_pos = self._nb_samples - 1
        window_size = 1
        passed = 0
        result = 0

        self.best_fit_FAOW(window_size)
        slope = self._a
        result = slope
        self._last_window_size = window_size

        if last_sample_pos == 0:
            self._last_window_size = 0
            return 0.0
        
        for window_size in range(2, self._nb_samples-1):
            self.best_fit_FAOW(window_size)
            slope = self._a

            for j in range(1,window_size):
                pos = self._b - slope*j*self._dt

                max_bound = pos + self._delta
                min_bound = pos - self._delta

                sample_to_check = self._buffer[last_sample_pos-j]

                if sample_to_check <= max_bound and sample_to_check >= min_bound:
                    passed = passed + 1
                else:
                    break

            if passed == (window_size - 1):
                result = slope
                passed = 0
                self._last_window_size = window_size
            else:
                break

        return result

    def apply(self, sample):
        self.add_sample(sample)
        derivative = self.fit()

        return derivative

    
