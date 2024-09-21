package etics24qmc;

import umontreal.ssj.hups.*;
import umontreal.ssj.mcqmctools.*;
import umontreal.ssj.randvar.NormalGen;
import umontreal.ssj.rng.*;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.stochprocess.BrownianMotion;
import umontreal.ssj.stochprocess.GeometricBrownianMotion;
import umontreal.ssj.util.Chrono;

// Program for QMC and RQMC examples in the introduction of my book, in 2 and 12 dimensions.
public class TestAsianOptionRQMC {

   double secondsMC;
   double varianceMC;

   public TestAsianOptionRQMC(AsianOption asian) {
      // this.asian = asian;
   }


   // Make an RQMC experiment and compare with MC, for each point set.
   // Sobol point set will have 2^k points.
   // Korobov (LCG) point set will have n1 points.
   // Multiplier for LCG point set will be a1.
   // This can be used for an arbitrary `model`.
   public void experOneDimRQMC(MonteCarloModelDouble model, int dim, int k, int n1, int a1, int m, RandomStream noise) {
      PointSetRandomization dShift = new LMScrambleShift(noise);
      PointSetRandomization rShift = new RandomShift(noise);
      Tally statRQMC = new Tally("Stats on payoff with RQMC");
      Chrono timer = new Chrono();
      
      // Latin hypercube sampling with 2^k points.
      PointSet pLH = new LatinHypercube(1 << k, dim);  // 2^k intervals
      System.out.println(RQMCExperiment.simulReplicatesRQMCDefaultReportCompare(model, pLH, rShift, m, statRQMC,
            varianceMC, secondsMC));
      
      // We do stratification only in up to 4 dimensions 
      PointSet pStrat = null;
      if (dim <= 4) {
         // numInt is the number of intervals for each coordinate, to get approx. 2^k points.
         int numInt = (int)Math.round (Math.pow ((double)(1 << k), 1.0 / (double)dim));
         System.out.println("Stratification with " + numInt + " intervals");
         pStrat = new StratifiedUnitCube(numInt, dim);
         System.out.println(RQMCExperiment.simulReplicatesRQMCDefaultReportCompare(model, pStrat, rShift, m, statRQMC,
               varianceMC, secondsMC));
      }

      // Korobov lattice, alone and then with baker's transformation.
      PointSet pLCG = new LCGPointSet(n1, a1);
      System.out.println(RQMCExperiment.simulReplicatesRQMCDefaultReportCompare(model, pLCG, rShift, m, statRQMC,
            varianceMC, secondsMC));   
      BakerTransformedPointSet pLCGBaker = new BakerTransformedPointSet(pLCG);
      System.out.println(RQMCExperiment.simulReplicatesRQMCDefaultReportCompare
            (model, pLCGBaker, rShift, m, statRQMC, varianceMC, secondsMC));
      
      // We create a Sobol' sequence in d-1 dimensions, the we add
      // a first coordinate equal to i/n for the point i.
      // In two dimensions, this gives the Hammersley point set.
      DigitalSequenceBase2 p0 = new SobolSequence(k, 31, dim - 1);
      PointSet pSobol = p0.toNetShiftCj();      
      System.out.println(RQMCExperiment.simulReplicatesRQMCDefaultReportCompare(model, pSobol, dShift, m, statRQMC,
            varianceMC, secondsMC));
      System.out.println("Total CPU time:      " + timer.format() + "\n");
   }

   // Main program: QMC and RQMC experiment with Asian option.
   public static void main(String[] args) {
      int numObsTimes = 4;
      double T1 = 1.0 / numObsTimes;
      double T = 1.0;
      double strike = 100.0;
      double s0 = 100.0;
      double r = 0.05;
      double sigma = 0.5;
      RandomStream noise = new LFSR113();
      NormalGen gen = new NormalGen(noise);
      AsianOption asian = new AsianOption(r, numObsTimes, T1, T, strike);
      asian.setProcess(new GeometricBrownianMotion(s0, r, sigma, new BrownianMotion(0, 0, 1, gen)));
      TestAsianOptionRQMC test = new TestAsianOptionRQMC(asian);
      System.out.println("\n******************************************");
      System.out.println("Asian option with sequential sampling \n");

      // We do the experiment with MC only once, and save the variance and average time per run.
      Tally statValueMC = new Tally("Stats on payoff with crude MC");
      int n = 10000000; // 10 million runs for Monte Carlo.
      Chrono timer = new Chrono();
      System.out.println("We first perform crude MC for comparison \n");
      System.out.println(MonteCarloExperiment.simulateRunsDefaultReportStudent
            (asian, n, noise, statValueMC, 0.95, 4, timer));
      // We memorize CPU time and variance to compare with RQMC.
      test.secondsMC = timer.getSeconds() / n;
      test.varianceMC = statValueMC.variance();

      // RQMC experiments.
      int m = 100;  // Number of RQMC replicates.
      // test.experOneDimRQMC(asian, numObsTimes, 7, 101, 12, m, noise);
      test.experOneDimRQMC(asian, numObsTimes, 10, 1024, 115, m, noise);
      test.experOneDimRQMC(asian, numObsTimes, 16, 65521, 944, m, noise);
      test.experOneDimRQMC(asian, numObsTimes, 18, 262139, 21876, m, noise);
   }

}
