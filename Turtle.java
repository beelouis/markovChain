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

  /**
    * Constructor:
    * holds a complete transition table from each state to all others
    * has a current state, as well as the local neighbour table to that state
    * must have either inputted steady state probabilities, or all equal
    * must be set to either continuous or discrete time
  */

  public Turtle(int currentState, int numStates, boolean equalSSProb, boolean continuous){
    ran = new Random();
    this.currentState = currentState;
    this.numStates = numStates;
    this.continuous = continuous;
  }

  // set equal steady state probabilities

  public void setEqualSSProb(){
    double[] ssprobs = new double[numStates];
    for (int i = 0; i < ssprobs.length; i++) ssprobs[i] = (1.0 / numStates);
    setSSProb(ssprobs);
  }

  // set the steady state probabilities to what the user specifies

  public void inputSSProb(double[] ssprobs){
    setSSProb(ssprobs);
  }

  // method for setting the ssprob

  private void setSSProb(double[] ssprobs){
    this.ssprobs = ssprobs;
    this.ssprobCurrentState = ssprobs[currentState - 1];
  }

  /**
    * Create a complete neighbour table for Discrete Time MC, and assign to turtle the current state's table:
    * --------------------------------------------------------------------------------------------------
    * this table never changes in all of n iterations of moving through the chain
    * after changing state, the turtle simply has to get its new neighbours
    * HashMap<>() localTable = mapping of neihbour states and probability of moving to that state
    *                          for each ith state - i itself is not the state in the map
    * then this local table is put into the complete hash table, with the ith state as the key
  */

  public void setDiscreteNbrs(){
    totalTable = new HashMap<>();
    for (int i = 1; i <= numStates; i++){

      HashMap<Integer, Double> localTable = new HashMap<>();
      HashSet<Integer> nbrs = calculateNbrs(i);

      // add the neighbours to the local table, with proposal probability of 1 / directions
      for (Integer state : nbrs) localTable.put(state, (1.0 / 4));

      for (Integer state : localTable.keySet()){
        // get transition probabilities from proposal probabilities and steady state probabilities
        double ssprobTo = calcSSProbTo(state);
        double pacc = calcPAcc(ssprobTo, ssprobs[i-1]);
        double finalTrans = pacc * localTable.get(state);

        localTable.put(state, finalTrans);
      }

      // sum all transition probabilities; if it is less than 1, add a self transition
      double sum = localTable.values().stream().reduce(0.0, (p, c) -> p + c);
      if (sum < 1.0) localTable.put(i, 1.0 - sum);

      // remove ghost state from the list of potential neighbours:
      localTable.remove(0);
      totalTable.put(i, localTable);
    }
    // get the current state's neighbour table for immediate use
    nbrsTransProb = totalTable.get(currentState);
  }

  /**
    * Subroutine called from setDiscreteNbrs() which is assigned to HashMap<>() localTable:
    * ---------------------------------------------------------------------------------------
    * check if the necessary move would take you outside the bounds of the grid
    * being assigned a state of 0 represents a "ghost state"
  */

  private HashSet<Integer> calculateNbrs(int fromState){
    Integer east, north, west, south;
    int sqrt = (int)Math.sqrt(numStates);

    east  = (fromState % sqrt != 0)         ? fromState + 1     : 0;
    north = (fromState + sqrt <= numStates) ? fromState + sqrt  : 0;
    west  = ((fromState - 1) % sqrt != 0)   ? fromState - 1     : 0;
    south = (fromState - sqrt > 0)          ? fromState - sqrt  : 0;

    return new HashSet<>(Arrays.asList(east, north, west, south));
  }

  // Subroutine from setDiscreteNbrs()

  private double calcSSProbTo(int state){
    // State 0 denotes a ghost state; Steady state probability of leaving the chain is 0
    return state == 0 ? 0.0 : ssprobs[state-1];
  }

  // Subroutine from setDiscreteNbrs()

  private double calcPAcc(double ssprobTo, double ssprobFrom){
    return Math.min(1, (ssprobTo / ssprobFrom) );
  }


  /**
    * Create a total neighbour table for Continuous Time MC, and assign to turtle the current states' table:
    * ----------------------------------------------------------------------------------------------------
    * index is the element in the rates array that I begin the loop pulling rates from,
    * and it will break after 2 loops since we only have 2 transitions per state
    * "numstates - 1" accounts for skipping self transitions.
    * 9 possible combinations of transitions (n^2), but only 6 (n^2-n) without self transition
    * "fromState - 1" means that it starts the index of the rates array:
    * at index 0 for state 1  : (2 * (1-1 = 0) = 0),
    * at index 2 for state 2  : (2 * (2-1 = 1) = 2),
    * at index 4 for state 3  : (2 * (3-1 = 2) = 4)
  */

  public void setContTransitions(double[] rates){
    totalTable = new HashMap<>();
    totalContRates = new HashMap<>();
    this.rates = rates;

        for (int fromState = 1; fromState <= 3; fromState++){
            HashMap<Integer, Double> localContTable = new HashMap<>();
            HashMap<Integer, Double> localTransProb = new HashMap<>();

              int index = (numStates - 1) * (fromState - 1);
              for (int toState = 1; toState <= 3; toState++){
                if (fromState == toState) { continue; } // skip self transitions
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
  }

  // Calculate DeltaT in CTMC: the amount of time that passes between steps

  public double calculateDeltaT(){
    double r = ran.nextDouble();
    double propensity = nbrsContRates.values().stream().reduce(0.0, (p, c) -> p + c);
    double deltaT = (-1 / propensity) * Math.log(r);
    return deltaT;
  }


  // Tower sample over possible states to stochastically move in the next time step

  public void towerSample(){
    // put all transition probabilitiies in arraylist, then order it highest to lowest
    ArrayList<Double> preT = new ArrayList<>(nbrsTransProb.values());
    Collections.sort(preT, Collections.reverseOrder());

    double[] t = new double[preT.size() + 1];
    t[0] = 0.0;
    for (int i = 1; i < t.length; i++){
      // enter the ordered probabilities, AFTER the 0.0 in index 0, into array t
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
    changeState(nextState);
  }

  // Subroutine from towerSample() which chooses between states with equal probabilities

  private int decideNextState(double prob){
      ArrayList<Integer> possibleStates = new ArrayList<>();
      for (Integer state : nbrsTransProb.keySet()){
        if (nbrsTransProb.get(state) == prob){
          possibleStates.add(state);
        }
      }
      int r = ran.nextInt(possibleStates.size());
      return possibleStates.get(r);
  }

  // Subroutine from towerSample() to change the turtle's state

  private void changeState(int nextState){
      currentState = nextState;
      nbrsTransProb = totalTable.get(currentState);
      if (!continuous)  ssprobCurrentState = ssprobs[currentState - 1];
      else              nbrsContRates = totalContRates.get(currentState);
  }

  /**
    * return the probability of moving to a specified state:
    * if the proposed state is not a neighbour, return probability 0
    * in CTMC, states are not neighbours to themselves, so this handles that as well

  */

  public double getTransProbTo(int nextState){
    return nbrsTransProb.containsKey(nextState) ? nbrsTransProb.get(nextState) : 0.0;
  }

  public int getState(){
    return currentState;
  }

  public void setState(int state){
    changeState(state);
  }

}
