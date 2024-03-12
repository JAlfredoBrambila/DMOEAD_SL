

package algorithm.reinforcementlearning;

import java.util.HashMap;
import java.util.Map;
import java.util.Random;

/**
 *
 * @author J. Alfredo Brambila H. <alfredo.brambila@outlook.com>
 */
public class SARSALambdaBase {
    private int numStates;
    private int numActions;
    private double alpha;
    private double gamma;
    private double lambda;
    private int numEpisodes;
    private int numSteps;
    
    private Map<StateActionPair, Double> qTable;
    private Map<StateActionPair, Double> eligibilityTraces;
    //double [] s;
    
   private double epsilon;
    
    public SARSALambdaBase(int numStates, int numActions, double alpha, double gamma, int numEpisodes, int numSteps) {
        this.numStates = numStates;
        this.numActions = numActions;
        this.alpha = alpha;
        this.gamma = gamma;
        this.numEpisodes = numEpisodes;
        this.numSteps = numSteps;
        
        qTable = new HashMap<>();
        eligibilityTraces = new HashMap<>();
    }
    
    public SARSALambdaBase(int numStates, int numActions, double alpha, double gamma, double lambda) {
        this.numStates = numStates;
        this.numActions = numActions;
        this.alpha = alpha;
        this.gamma = gamma;
        this.lambda = lambda;
        
        qTable = new HashMap<>();
        eligibilityTraces = new HashMap<>();
    }
    
    public void initTable() {
        qTable = new HashMap<>();
        eligibilityTraces = new HashMap<>();
    }
    
    public int chooseAction(int state) {
        Random random = new Random();
        if (random.nextDouble() < 0.1) {
            return random.nextInt(this.numActions);
        } else {
            double bestQvalue = Double.NEGATIVE_INFINITY;
            int bestAction = 0;
            
            for(int action = 0; action < this.numActions; action++) {
                double qValue = getQValue(state, action);
                if(qValue > bestQvalue) {
                    bestQvalue = qValue;
                    bestAction = action;
                }
            }
            
            return bestAction;
        }
    }
    
    public int chooseAction(int state, double ep) {
        this.epsilon = ep;
        Random random = new Random();
        if (random.nextDouble() < this.epsilon) {
            //System.out.print("Rand__ ");
            return random.nextInt(this.numActions);
        } else {
            double bestQvalue = Double.NEGATIVE_INFINITY;
            int bestAction = 0;
            
            
            for(int action = 0; action < this.numActions; action++) {
                double qValue = getQValue(state, action);
                System.out.println("QVal act: " + action + " st " + state + " = " + qValue);
                if(qValue > bestQvalue) {
                    bestQvalue = qValue;
                    bestAction = action;
                }
            }
            
            //System.out.println("Table__ ");
            //System.out.println("Select: " + bestAction);
            
            return bestAction;
        }
    }
    
    public double getQValue(int state, int action) {
        StateActionPair pair = new StateActionPair(state, action);
        return qTable.getOrDefault(pair, 0.0);
    }
    
    public void updateEligibilityTraces(int state, int action) {
        StateActionPair pair = new StateActionPair(state,action);
        eligibilityTraces.put(pair, 0.0);
    }
    
    public void updateQTable(double reward, double currentQValue, double nextQValue) {
        double delta = reward + this.gamma * nextQValue - currentQValue;
        for(StateActionPair pair : eligibilityTraces.keySet() ) {
            double oldValue = getQValue(pair.getState(), pair.getAction());
            double eligibilityTrace = eligibilityTraces.get(pair);
            double newValue = oldValue + this.alpha * delta * eligibilityTrace;
            qTable.put(pair, newValue);
        }
    }
    
    public void decayEligibilityTraces() {
        for(StateActionPair pair : eligibilityTraces.keySet() ) {
            double decayedValue = this.gamma * this.lambda * eligibilityTraces.get(pair);
            eligibilityTraces.put(pair, decayedValue);
        }
    }
    
    public double getReward(int state, int action) {
        double reward = 0;
        
        if(state > 0) {
            
            reward = -1;
        }
        
        return reward;
    }
    
    public int getNextState(int state) {
        if(state == 0) {
            return 1;
        } else {
            return 0;
        }
    }
    
}


class StateActionPair {
    private int state;
    private int action;
    
    public StateActionPair(int state, int action) {
        this.state = state;
        this.action = action;
    }
    
    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + action;
        result = prime * result + state;
        
        return result;
    }
    
    public boolean equals(Object obj) {
        
        if(this == obj)
            return true;
        if(obj == null)
            return false;
        if(getClass() != obj.getClass())
            return false;
        
        StateActionPair other = (StateActionPair) obj;
        
        if(action != other.action)
            return false;
        if(state != other.state)
            return false;
        
        return true;
    }
    
    public int getState() {
        return this.state;
    }
    
    public int getAction() {
        return action;
    }
    
    
}