package etics24qmc;

import java.io.*;

import umontreal.ssj.mcqmctools.MonteCarloModelDouble;
import umontreal.ssj.mcqmctools.RQMCExperiment;
import umontreal.ssj.randvar.NormalGen;
import umontreal.ssj.rng.LFSR113;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stochprocess.BrownianMotion;
import umontreal.ssj.stochprocess.GeometricBrownianMotion;
import umontreal.ssj.util.Chrono;


// Generate and store RQMC replicates for Asian option model, for WSC 2023 paper

public class AsianRQMC extends RQMCExperiment {
	
	public static void main(String[] args) throws IOException {

		int m = 10000;                // Number of RQMC randomizations.
		int s = 12;
		double T1 = 1.0 / s;
		double T = 1.0;
		double strike = 100.0;
		double s0 = 100.0;
		double r = 0.05;
		double sigma = 0.5;
		RandomStream noise = new LFSR113();
		NormalGen gen = new NormalGen(noise);
		//MonteCarloModelDoubleTag asian = new AsianOption(r, s, T1, T, strike);
		AsianOption asian = new AsianOption(r, s, T1, T, strike);
		asian.setProcess(new GeometricBrownianMotion(s0, r, sigma,
				new BrownianMotion(0, 0, 1, gen)));

		Chrono timerTotal = new Chrono();
		for (int k = 6; k <= 14; k = k+2) {
		   	// RQMCSamples.simulRepsAllTypes (asian, s, k, m);
		   	// RQMCSamplesEtics24.simulRepsLMS (asian, s, k, m);
			}

        System.out.println ("Total time for everything: " 
		  + timerTotal.format() + "\n================================== \n");
	}
}
