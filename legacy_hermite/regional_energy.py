"""Regional energy accounting; no assumption that errors decay from the clamp."""
import numpy as np


class RegionalEnergy:
    def __init__(self, beam, exclude_clamp_lengths=()):
        self.z0 = float(beam.grid_points[:, 2].min())
        self.z1 = float(beam.grid_points[:, 2].max())
        self.length = self.z1-self.z0
        lengths = np.asarray(list(exclude_clamp_lengths), dtype=float)
        if not np.all(np.isfinite(lengths)) or np.any(lengths < 0) or np.any(lengths >= self.length):
            raise ValueError('Excluded clamp distances must be finite, >=0, and less than beam length')
        self.lengths = sorted(set(lengths.tolist()))
        self.cutoffs = [self.z0+d for d in self.lengths]
        self.total = np.zeros(2)
        self.retained = np.zeros((len(self.cutoffs), 2))
        self.excluded = np.zeros((len(self.cutoffs), 2))

    def split(self, low, high):
        return [low] + [z for z in self.cutoffs if low < z < high] + [high]

    def add(self, difference_sq, reference_sq, midpoint_z):
        values = np.array([difference_sq, reference_sq])
        self.total += values
        for i, cutoff in enumerate(self.cutoffs):
            target = self.retained if midpoint_z >= cutoff else self.excluded
            target[i] += values

    def _metrics(self, values, volume_percent):
        difference_sq, reference_sq = np.maximum(values, 0.)
        return {'difference_energy_norm': float(np.sqrt(difference_sq)),
                'difference_energy_norm_squared_lbf_in': float(difference_sq),
                'reference_energy_norm_squared_lbf_in': float(reference_sq),
                'reference_energy_norm': float(np.sqrt(reference_sq)),
                'relative_energy_difference_percent':
                    float(100*np.sqrt(difference_sq/reference_sq)) if reference_sq > 0 else None,
                'volume_percent': float(volume_percent),
                'share_of_global_squared_difference_percent':
                    float(100*difference_sq/self.total[0]) if self.total[0] > 0 else None}

    def report(self):
        return [{'excluded_clamp_length_in': d,
                 'retained_z_range_in': [self.z0+d, self.z1],
                 'retained': self._metrics(self.retained[i], 100*(1-d/self.length)),
                 'excluded': self._metrics(self.excluded[i], 100*d/self.length)}
                for i, d in enumerate(self.lengths)]
