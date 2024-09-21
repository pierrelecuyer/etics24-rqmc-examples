package etics24qmc;

import umontreal.ssj.mcqmctools.MonteCarloModelDouble;
import umontreal.ssj.probdist.JohnsonSUDist;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.probdist.NormalDist;

// Experiments for WSC 2023 paper

public class RidgeJohnsonSU implements MonteCarloModelDouble {

	int s;
	double xi = 0.0, lambda = 1.0, gamma = 1.0, delta = 1.0;
	double sum, val;
	double invsqrts, mean;

	// Constructor.
	public RidgeJohnsonSU(int s) {
		this.s = s;
		invsqrts = 1.0 / Math.sqrt ((double)s);
		mean = xi - lambda * Math.exp (0.5 / (delta * delta)) * Math.sinh (gamma/delta);
	}

	public void simulate (RandomStream stream) {
		sum = 0.0;
		for (int j = 0; j < s; j++) {
			sum += NormalDist.inverseF01 (stream.nextDouble());
		}
		sum *= invsqrts;
		val = JohnsonSUDist.inverseF (gamma, delta, xi, lambda, NormalDist.cdf01(sum)) - mean;	
	}

	public double getPerformance () {
		return val;
	}

	// Descriptor of this model.
	@Override
	public String toString () {
		return "RidgeJohnsonSU: Ridge for the sum of s SumJohnsonSU";
	}

	// Short descriptor (tag) for this model.
	public String getTag () {
		return "RidgeJohnsonSU";
	}
}
