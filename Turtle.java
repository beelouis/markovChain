import java.lang.Math;
import java.util.*;

class Turtle{

  private int numStates;
  private int currentState;
  private double[] ssprobs;
  private double ssprobCurrentState;

  private HashMap<Integer, HashMap<Integer, Double>> totalTable;
  private HashMap<Integer, Double> nbrsTransProb;

  private HashMap<Integer, HashMap<Integer, Double>> totalContRates;
  private HashMap<Integer, Double> nbrsContRates;
  private boolean continuous;
  private double[] rates;

  private Random ran;

  // ========================================================================================
  // ======================== Constructor ===================================================
  // ========================================================================================

  public Turtle(int currentState, int numStates, boolean equalSSProb, boolean continuous){
    ran = new Random();
    this.currentState = currentState;
    this.numStates = numStates;
    this.continuous = continuous;
  }

  // ========================================================================================
  // ========= Create Complete neighbour table for DTMC  ====================================
  // ========================================================================================


  public void setDiscreteNbrs(){
    // this table never changes in all of n iterations of moving through the chain
    // after changing state, the turtle simply has to get its new neighbours
    totalTable = new HashMap<>();
    for (int i = 1; i <= numStates; i++){

      // the localTable is a mapping of neighbour states and the probability of moving to that state
      // it is local to the ith state in the loop, where i is not actually in the map
      HashMap<Integer, Double> localTable = new HashMap<>();

      // Using a HashSet because it will remove duplicate entries of ghost states which simplifies things
      HashSet<Integer> nbrs = calculateNbrs(i);

      // add the neighbours to the local table, with proposal probability of 1 / directions
      for (Integer state : nbrs) localTable.put(state, (1.0 / 4));

      for (Integer state : localTable.keySet()){
        // get transition probabilities from proposal probabilities and steady state probabilities
        double ssprobTo = calcSSProbTo(state);
        double pacc = calcPAcc(ssprobTo, ssprobs[i-1]);
        double finalTrans = pacc * localTable.get(state);

        // update the local table with new transition probabilities
        // now, the transition probability of moving to a ghost state will be 0, instead of 1/4
        // this wont contribute to the sum, so if there are any ghost states then (sum of probs != 1)
        localTable.put(state, finalTrans);
      }

      // sum all transition probabilities; if it is less than 1, add a self transition
      double sum = localTable.values().stream().reduce(0.0, (p, c) -> p + c);
      if (sum < 1.0) localTable.put(i, 1.0 - sum);

      // remove ghost state from the list of potential neighbours:
      localTable.remove(0);

      // add an entry into the hashmap for the ith state, to its neighbours with their probabilities:
      totalTable.put(i, localTable);
    }
    // get the current state's neighbour table for immediate use
    nbrsTransProb = totalTable.get(currentState);
  }

  // =========== Subroutine for getting all the neighbours for the markov chain =============

  private HashSet<Integer> calculateNbrs(int fromState){
    Integer east, north, west, south;
    int sqrt = (int)Math.sqrt(numStates);

    // check if the necessary move would take you outside the bounds of the grid
    // being assigned a state of 0 represents a "ghost state"
    east  = (fromState % sqrt != 0)         ? fromState + 1     : 0;
    north = (fromState + sqrt <= numStates) ? fromState + sqrt  : 0;
    west  = ((fromState - 1) % sqrt != 0)   ? fromState - 1     : 0;
    south = (fromState - sqrt > 0)          ? fromState - sqrt  : 0;

    return new HashSet<>(Arrays.asList(east, north, west, south));
  }

  // ========================================================================================
  // =========== Create a local neighbour table for current state in CTMC ==================
  // ========================================================================================

  public void setContTransitions(double[] rates){
    totalTable = new HashMap<>();
    totalContRates = new HashMap<>();
    this.rates = rates;

        for (int fromState = 1; fromState <= 3; fromState++){
            HashMap<Integer, Double> localContTable = new HashMap<>();
            HashMap<Integer, Double> localTransProb = new HashMap<>();

              int index = (numStates - 1) * (fromState - 1);
              // index is the element in the rates array that I begin pulling rates from,
              // and it will break after 2 loops since we only have 2 transitions per state

              // "numstates - 1" accounts for skipping self transitions.
              // 9 possible combinations of transitions (n^2), but only 6 (n^2-n) without self transition
              // "fromState - 1" means that it starts the index of the rates array:
              //  at index 0 for state 1  : (2 * (1-1 = 0) = 0),
              //  at index 2 for state 2  : (2 * (2-1 = 1) = 2),
              //  at index 4 for state 3  : (2 * (3-1 = 2) = 4)

              for (int toState = 1; toState <= 3; toState++){
                if (fromState == toState)  continue; // skip self transitions
                else {
                  localContTable.put(toState, rates[index]);
                  index++;  // only increment through the rates array when we actually use one
                }
              }

              double sumRates = localContTable.values().stream().reduce(0.0, (p, c) -> p + c);

              for (Integer state : localContTable.keySet()){
                // transition probability of i -> j = (rate of ij / sum of all rates from i)
                double contTransProb = localContTable.get(state) / sumRates;
                localTransProb.put(state, contTransProb);
              }

          // enter rates for all neighbours to ith state into local map
          totalContRates.put(fromState, localContTable);
          // enter this map as the value to the ith state's key in the total table
          totalTable.put(fromState, localTransProb);
        }

    // get currentState's neighbour tables for immediate use
    nbrsTransProb = totalTable.get(currentState);
    nbrsContRates = totalContRates.get(currentState);
    // printContTransitionDetails();
  }

  // ========================================================================================
  // ================ Calculate DeltaT in CTMC ==============================================
  // ========================================================================================


  public double calculateDeltaT(){
    double r = ran.nextDouble();
    double propensity = nbrsContRates.values().stream().reduce(0.0, (p, c) -> p + c);
    double deltaT = (-1 / propensity) * Math.log(r);
    return deltaT;
  }

  // ========================================================================================
  // ======== Perform tower sampling to change state ========================================
  // ========================================================================================


  public void towerSample(){
    // put all transition probabilitiies in arraylist, then order it highest to lowest
    ArrayList<Double> preT = new ArrayList<>(nbrsTransProb.values());
    Collections.sort(preT, Collections.reverseOrder());

    double[] t = new double[preT.size() + 1];
    t[0] = 0.0;
    for (int i = 1; i < t.length; i++){
      // enter the ordered probabilities, after the 0.0 in index 0, into array t
      t[i] = preT.get(i-1);
    }

    int nextState = -1;
    double r = ran.nextDouble();
    double sum = 0.0;

    for (int i = 0; i < t.length; i++){
      sum += t[i];
      if (sum > r){
        nextState = decideNextState(t[i]);
        break;
      }
    }
    // printTowerSampleDetails(r, nextState);
    changeState(nextState);
  }

  // ========= Subroutine for Tower sampling ================================================

  private int decideNextState(double prob){
      // Randomly picks between states with the same probability
      ArrayList<Integer> possibleStates = new ArrayList<>();
      for (Integer state : nbrsTransProb.keySet()){
        // if the probability of moving to this neighbour is the same as the one chosen in tower sampling, it is an option:
        if (nbrsTransProb.get(state) == prob){
          possibleStates.add(state);
        }
      }
      int r = ran.nextInt(possibleStates.size());
      return possibleStates.get(r);
  }

  // ========= Subroutine for Tower sampling ================================================

  private void changeState(int nextState){
      currentState = nextState;
      nbrsTransProb = totalTable.get(currentState);
      if (!continuous)  ssprobCurrentState = ssprobs[currentState - 1];
      else              nbrsContRates = totalContRates.get(currentState);
  }


  // ========================================================================================
  // ========= getters, setters, and debuggers  =============================================
  // ========================================================================================

  public double getTransProbTo(int nextState){
    // if the proposed state is not a neighbour, return probability 0
    // in CTMC, states are not neighbours to themselves, so this handles that as well
    return nbrsTransProb.containsKey(nextState) ? nbrsTransProb.get(nextState) : 0.0;
  }

  private double calcSSProbTo(int state){
    // State 0 denotes a ghost state; Steady state probability of leaving the chain is 0
    return state == 0 ? 0.0 : ssprobs[state-1];
  }

  private double calcPAcc(double ssprobTo, double ssprobFrom){
    return Math.min(1, (ssprobTo / ssprobFrom) );
  }

  public int getState(){
    return currentState;
  }

  public void setState(int state){
    changeState(state);
  }

  public void setEqualSSProb(){
    double[] ssprobs = new double[numStates];
    for (int i = 0; i < ssprobs.length; i++) ssprobs[i] = (1.0 / numStates);
    setSSProb(ssprobs);
  }

  public void inputSSProb(double[] ssprobs){
    setSSProb(ssprobs);
  }

  private void setSSProb(double[] ssprobs){
    this.ssprobs = ssprobs;
    this.ssprobCurrentState = ssprobs[currentState - 1];
  }


  private void printTowerSampleDetails(double r, int nextState){
    System.out.println("=====================================");
    System.out.println("Current State: " + currentState);
    nbrsTransProb.entrySet().forEach(e -> {
      System.out.println("To State: " + e.getKey());
      System.out.println("Prob: " + e.getValue());
    });
    System.out.println("R: " + r);
    System.out.println("Changing to state: " + nextState);
  }

  private void printContTransitionDetails(){
    nbrsContRates.entrySet().forEach(e -> {
      System.out.println("state: " + e.getKey());
      System.out.println("rate: " + e.getValue());
    });
    nbrsTransProb.entrySet().forEach(e -> {
      System.out.println("state: " + e.getKey());
      System.out.println("prob: " + e.getValue());
    });
  }
}
