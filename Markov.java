
import java.lang.Math;
import java.util.*;
import java.io.PrintWriter;


class Markov{

  // ========================================================================================
  // ========================================================================================
  // ========================================================================================

  public static  double getTransProb(int s1, int s2, int numStates){
    Turtle t = new Turtle(s1, numStates, true, false);
    t.setEqualSSProb();
    t.setDiscreteNbrs();
    return t.getTransProbTo(s2);
  }

  // ========================================================================================
  // ========================================================================================
  // ========================================================================================

  public static double getSejProb(int s1, int s2, int numStates, double TS){
    int n = 1000000;
    // int n = 100;
    int[] count = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    Turtle t = new Turtle(s1, numStates, true, false);
    t.setEqualSSProb();
    t.setDiscreteNbrs();

    for (int i = 0; i < n; i++){
        // set mover back to start state
        int currentState = s1;
        t.setState(currentState);

        for (int j = 0; j < TS; j++){
          // tower sample for TS steps
          t.towerSample();
          currentState = t.getState();
        }

      // increment counter for whatever state we are at after TS time steps
      count[currentState] ++;
    }
    return (double)count[s2] / n;
  }

  // ========================================================================================
  // ========================================================================================
  // ========================================================================================


  public static double getBiasTransProb(int s1, int s2, double[] ssprob){
    Turtle t = new Turtle(s1, 9, false, false);
    t.inputSSProb(ssprob);
    t.setDiscreteNbrs();
    return t.getTransProbTo(s2);
  }

  // ========================================================================================
  // ========================================================================================
  // ========================================================================================


  public static double  getContTransProb(int s1, int s2, double[] rates){
    Turtle t = new Turtle(s1, 3, false, true);
    t.setContTransitions(rates);
    return t.getTransProbTo(s2);
  }

  // ========================================================================================
  // ========================================================================================
  // ========================================================================================


  public static double getContSejProb(int s1, int s2, double[] rates, double TSC){
      int n = 1000000;
      // int n = 100;
      int[] count = { 0, 0, 0, 0 };

      Turtle t = new Turtle(s1, 3, false, true);
      t.setContTransitions(rates);

      for (int i = 0; i < n; i++){
          int currentState = s1;
          t.setState(currentState);
          double time = 0.0;
          double deltaT = 0.0;

          while (time < TSC){

              // continuously increment time by new deltaT
              deltaT = t.calculateDeltaT();
              time += deltaT;

              // tower sample while we are within the desired time range
              // if we have already surpassed TSC in the first deltaT for example, then
              // the state the turtle is in after TSC is the same state it started in
              if (time < TSC) { t.towerSample(); }
              currentState = t.getState();
          }
        // increment counter for the state we are in after TSC time
        count[currentState] ++;
      }
      double estimatedProb = (double)count[s2]/n;
      return estimatedProb;
  }

}
