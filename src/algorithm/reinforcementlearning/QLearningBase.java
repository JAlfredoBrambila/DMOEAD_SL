package algorithm.reinforcementlearning;


import java.util.Random;


/**
 *
 * @author J. Alfredo Brambila H. <alfredo.brambila@outlook.com>
 */
public class QLearningBase {

    private int numStates;
    private int numActions;
    private double alpha;
    private double gamma;
    private int numEpisodes;
    private int numSteps;
    
    double [][] qTable;
    double [] s;
    
    public QLearningBase(int numStates, int numActions, double alpha, double gamma, int numEpisodes, int numSteps) {
        this.numStates = numStates;
        this.numActions = numActions;
        this.alpha = alpha;
        this.gamma = gamma;
        this.numEpisodes = numEpisodes;
        this.numSteps = numSteps;
        
        qTable = new double[numStates][numActions];
        initTable();
        s = new double[numStates];
    }
    
    public QLearningBase(int numStates, int numActions, double alpha, double gamma) {
        this.numStates = numStates;
        this.numActions = numActions;
        this.alpha = alpha;
        this.gamma = gamma;
        this.numEpisodes = numEpisodes;
        this.numSteps = numSteps;
        
        qTable = new double[numStates][numActions];
        initTable();
        s = new double[numStates];
    }
    
    /*
    public void run() {
        for(int i=0; i<numEpisodes; i++) {
            
            for(int t=0; t<numSteps; t++) {
                
            }
        }
    }
    */
    
    public int chooseAction(int state) {
        Random random = new Random();
        if(random.nextDouble() < 0.1) {
            return random.nextInt(this.numActions);
        } else {
            double[] stateActions = qTable[state];
            int bestAction = 0;
            for(int i=1; i<this.numActions; i++) {
                if(stateActions[i] > stateActions[bestAction]) {
                    bestAction = i;
                }
            }
            return bestAction;
        }
    }
    
    public double argMax(int state) {
        double[] stateActions = qTable[state];
        double maxQValue = stateActions[0];
        
        for(int i=1; i<this.numActions; i++) {
            if(stateActions[i] > maxQValue) {
                maxQValue = stateActions[i];
            }
        }
        
        return maxQValue;
    }
    
    public void updateQTable(int state, int action, int nextState, double reward) {
        double currentQValue = qTable[state][action];
        // Bellman Equation
        double newQValue = currentQValue + alpha * (reward + gamma * argMax(nextState) - currentQValue);
        qTable[state][action] = newQValue;
    }
    
    public void initTable() {
        for(int i=0; i<qTable.length; i++) {
            for(int j=0; j<qTable[0].length; j++) {
                qTable[i][j] = 0.0;
            }
        }
    }
    
    
    public double getReward(int state, int action) {
        double reward = 0;
        
        if(state > 0) {
            
            reward = -1;
        }
        
        return reward;
    }
    
}
