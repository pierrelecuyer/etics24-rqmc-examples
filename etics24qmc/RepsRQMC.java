package etics24qmc;

import java.io.*;

//import umontreal.ssj.hups.PointSet;
//import umontreal.ssj.hups.PointSetRandomization;
import umontreal.ssj.mcqmctools.MonteCarloModelDouble;
import umontreal.ssj.mcqmctools.RQMCExperiment;
//import umontreal.ssj.stat.Tally;
import umontreal.ssj.util.Chrono;

// Generate and store RQMC replicates for various models, s, and k.

public class RepsRQMC extends RQMCExperiment {

   public static void main(String[] args) throws IOException {

      int m = 100; // Number of RQMC randomizations.
      MonteCarloModelDouble model;
      Chrono timerTotal = new Chrono();

      for (int s = 4; s <= 32; s *= 2) {
         for (int k = 10; k <= 14; k = k + 2) {
            // Uncomment the model you want below. ***
            // model = new SumUeU(s);
            // model = new IndBox(s);
            // model = new PieceLinGauss(s);
            // model = new IndSumNormal(s);
            // model = new SmoothGauss(s);
            // model = new SumUniforms(s);
            // model = new MC2(s);
            model = new RidgeJohnsonSU(s);
            model = new ProductExpCosRQMC(s, 1.0, 1.0);
            
            RQMCSamplesEtics24.simulRepsManyTypes(model, s, k, m);
            // RQMCSamples.simulRepsNUS (model, s, k, m);
         }
      }

      // model = new MC2(8);
      // RQMCSamples.simulRepsAllTypes (model, 8, 6, m);

      System.out
            .println("Total time for everything: " + timerTotal.format() + "\n================================== \n");
   }
}
