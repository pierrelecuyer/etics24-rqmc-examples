package etics24qmc;

import umontreal.ssj.mcqmctools.MonteCarloModelDouble;
import umontreal.ssj.rng.RandomStream;

// Experiments for WSC 2023 paper.

public class Gaussian implements MonteCarloModelDouble {

	int s;
	double sum;

	// Constructor.
	public Gaussian(int s) {
		this.s = s;
	}

	// Generates and returns X, without IS.
	public void simulate (RandomStream stream) {
		sum = 0.0;
		double u;
		for (int j = 0; j < s; j++) {
		   u = stream.nextDouble();
			sum += u * u;
		}
	}

	// Generates and returns X, without IS.
	public double getPerformance () {
		return Math.exp(sum);
	}

	// Descriptor of this model.
	@Override
	public String toString () {
		return "Gaussian: exponential of the sum of squares";
	}

	// Short descriptor (tag) for this model.
	public String getTag () {
		return "Gaussian";
	}
}