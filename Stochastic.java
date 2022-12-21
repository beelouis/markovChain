//This is the brief for the Stochastic Assessment
///
//_______________________________________________________
//
//         Edit this file at your own peril.
//
//         I recommend only changing numbers.
//
//         Do not add anything. During marking I will use a different version of Stochastic.java
//
//         Any code you add here will ne lost and ignored.
//
//______________________________________________________

import java.lang.Math;
import java.util.*;


class Stochastic {
		public static void main(String[] args){
				//Do not delete/alter the next line
		    long startT=System.currentTimeMillis();


//Q1
//Assume a square 3x3 grid on which a turtle performs a random walk in discrete time.
// The cells in the grid should be numbered as follows:
//	7 8 9
//	4 5 6
//	1 2 3
//The turtle can walk east, west, south, north.
//Translate the grid into a Markov chain.
//Assume that the desired steady state probabilities of all states are equal.
// Return  as A1 the  transition probability that the turtle walks from state s1 to state s2
//You must use the method we discussed in the lectures. I do not accept other answers (even if correct).



//numStates .... number of states of the Markov state (constrained to be k*k, where k=1,2,3,4,...).
//s1, s2 ...  the states of the Markov chain.

int  numStates =9;
int s1 =1;
int s2 = 2;

System.out.println("Question 1 ===========================");
double A1 = Markov.getTransProb(s1,s2,numStates);


//Q2
//Same as in the first problem but now for a grid of size 2x2
//This grid should be numbered as follows:
//
//	3 4
//	1 2
//
numStates = 4;
 s1 =1;
 s2 = 1;
 System.out.println("Question 2 ===========================");

 double A2 = Markov.getTransProb(s1,s2,numStates);


//Q3
//Create a Markov chain as in Q1 and Q2 and return as A3 the estimated probability that the turtle  is in state s2 at time step TS when it started in state s1.
//The Markov chain should represent a 3x3 grid (numbered as above) on which the turtle makes transitions up/down/left/right
//For this question, I will not change the number of states during marking
//All other values may be changed.



numStates =9; // fixed for this question

s1=1;
s2=2;
int TS=3;

double A3 = Markov.getSejProb(s1,s2,numStates,TS);



//Q4
//Create a Markov chain such that the steady state probabilities for state 1, 2, etc.. are as in the array ssprob.
//The Markov chain should represent a 3x3 grid as in Q1 - Q3.
// A4 should return the transition probability to go from state s1 to state s2
// such that in steady state the probabilities to be in a state are as desired
//I will not change the number of states of the chain for this question (but I will change the probabilities and the value of s1,s2)
//Note that it could be that s1=s2
//You must use the method we discussed in the lectures. I do not accept other answers (even if correct).

//ssprob[0] ....steady state probability to be in state 1
//ssprob[1] ... steady state probability to be in state 2
//etc....

double[] ssprob  = {
	0.1,
	0.1,
	0.1,
	0.2,
	0.05,
	0.2,
	0.1,
	0.05,
	0.1};
s1=6;
s2=5;

System.out.println("Question 4 ===========================");

double A4 = Markov.getBiasTransProb(s1,s2,ssprob);


//Q5
//We now  switch to a continuous time Markov chain, consisting of 3 states.
//This no longer models a grid with a turtle, but an arbitrary CTMC.
//The rates are encoded in the array:

double[] rates = {10.0,5.0,1.0,2.0,5.0,2.0};
//This means the following transition rates:
//
//
//   State 1 ---> State 2:   10.0
//   State 1   ----> State 3: 5.0
//   State 2  -->  State 1 : 1.0
//   State 2  ---> State 3 : 2.0
//   State 3 --->  State 1 : 5.0
//   State 3 ---> Sytate 2:  2.0
// However, note that I will change the rates for marking (but not the number of states)
s1 = 1;
s2 = 2;
//Return as A5 the transition probability from state s1 to state s2
//Make sure that the function deals with the case when s1 and s2 are the same.
System.out.println("Question 5 ===========================");

double A5 = Markov.getContTransProb(s1,s2,rates);

//Q6
// Using the rates as defined above (albeit with possibly different numerical values), return as A6 the (estimated) probability that at time TSC the Markov chain is in state s2 when it was started in state s1.
//For marking, I will change the numerical values, but keep the number of states constant.

double TSC=.02;

rates[0]=10.0;
rates[1]=5.0;
rates[2]=1.0;
rates[3]=1.0;
rates[4]=2.0;
rates[5]=3.0;

s1=1;
s2=1;

System.out.println("Question 6 ===========================");

double A6 = Markov.getContSejProb(s1,s2,rates,TSC);




//Below, we just print out the results of your attempt



System.out.println("Your answers to the questions were: " + " \n A1  "+A1 + "\n A2  "+A2+ "\n A3  "+A3+ "\n A4  "+A4+ "\n A5  "+A5+ "\n A6  "+A6);




		        long endT= System.currentTimeMillis();
			System.out.println("Total execution time was: " +  ((endT - startT)/1000.0) + " seconds");



	  }











}
