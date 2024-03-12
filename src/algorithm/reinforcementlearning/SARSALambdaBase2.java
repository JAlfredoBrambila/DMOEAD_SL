package algorithm.reinforcementlearning;


import java.util.Random;
import utils.DoubleUtils;


/**
 *
 * @author J. Alfredo Brambila H. <alfredo.brambila@outlook.com>
 */
public class SARSALambdaBase2 {

    private int numStates;
    private int numActions;
    private double alpha;
    private double gamma;
    private double lambda;
    private int numEpisodes;
    private int numSteps;
    
    double [][] qTable;
    double [] s;
    double [][] E;
    
    public SARSALambdaBase2(int numStates, int numActions, double alpha, double gamma, double lambda, int numEpisodes, int numSteps) {
        this.numStates = numStates;
        this.numActions = numActions;
        this.alpha = alpha;
        this.gamma = gamma;
        this.lambda = lambda;
        this.numEpisodes = numEpisodes;
        this.numSteps = numSteps;
        
        qTable = new double[numStates][numActions];
        initTable();
        E = new double[numStates][numActions];
        initE();
        s = new double[numStates];
    }
    
    public SARSALambdaBase2(int numStates, int numActions, double alpha, double gamma, double lambda) {
        this.numStates = numStates;
        this.numActions = numActions;
        this.alpha = alpha;
        this.gamma = gamma;
        this.lambda = lambda;
        this.numEpisodes = numEpisodes;
        this.numSteps = numSteps;
        
        qTable = new double[numStates][numActions];
        initTable();
        E = new double[numStates][numActions];
        initE();
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
            
            //System.out.println("Select: " + bestAction);
            return bestAction;
        }
    }
    
    public int chooseAction(int state, double epsilon) {
        Random random = new Random();
        if(random.nextDouble() < epsilon) {
            int indice = random.nextInt(this.numActions);
            //System.out.println("Rand: " + indice);
            return indice;
            
        } else {
            double[] stateActions = qTable[state];
            int bestAction = 0;
            for(int i=1; i<this.numActions; i++) {
                //System.out.println(stateActions[i] + " > " + stateActions[bestAction]);
                if(stateActions[i] > stateActions[bestAction]) {
                    
                    bestAction = i;
                    
                    //System.out.println("Best: " + bestAction);
                }
            }
            
            //System.out.println("Select: " + bestAction);
            return bestAction;
        }
    }
    
    /*public double argMax(int state) {
        double[] stateActions = qTable[state];
        double maxQValue = stateActions[0];
        
        for(int i=1; i<this.numActions; i++) {
            if(stateActions[i] > maxQValue) {
                maxQValue = stateActions[i];
            }
        }
        
        return maxQValue;
    }*/
    
    public void updateQTable(int state, int action, int nextState, int nextAction, double reward) {
        double currentQValue = qTable[state][action];
        double nextQValue = qTable[nextState][nextAction];
        // Bellman Equation
        double newQValue = currentQValue + alpha * (reward + gamma * nextQValue - currentQValue);
        qTable[state][action] = newQValue;
    }
    
    public void initTable(double[][] qTable) {
        for(int i=0; i<qTable.length; i++) {
            for(int j=0; j<qTable[0].length; j++) {
                //qTable[i][j] = 0.0;
                qTable[i][j] = DoubleUtils.getRandomNumber0_1();
            }
        }
    }
    
    public void initTable() {
        for(int i=0; i<qTable.length; i++) {
            for(int j=0; j<qTable[0].length; j++) {
                //qTable[i][j] = 0.0;
                qTable[i][j] = DoubleUtils.getRandomNumber0_1();
            }
        }
    }
    
    public void initE() {
        for(int i=0; i<qTable.length; i++) {
            for(int j=0; j<qTable[0].length; j++) {
                E[i][j] = 0.0;
                //qTable[i][j] = DoubleUtils.getRandomNumber0_1();
            }
        }
    }
    
    
    public double getReward(int state, int action) {
        double reward = 0;
        
        if(state > 0) {
            
            reward = 0.1; //0.0001
            //reward = 0.0001; //0.0001
            
        } else {
            //reward = 0.5; //0.5
            reward = 0.7;
        }
        
        return reward;
    }
    
    // Calcula el error temporal δ como:
    /*public double updateDelta(double reward,int state,int action,int nextState,int nextAction) {
        // δ=R+γQ(S',A')−Q(S,A)
        double delta = 0.0;
        delta = reward + this.gamma * this.qTable[nextState][nextAction] - this.qTable[state][action];
        
        return delta;
    }*/
    
    //
    public void updateEligibilityTraces(int state, int action) {
        // E(S,A)=λγE(S,A)+1
        //this.E[state][action] =  this.lambda * this.lambda * this.E[state][action] + 1;
		this.E[state][action] =  this.lambda * this.gamma * this.E[state][action] + 1;
    }
    
    public void updateEligibilityTraces(int state, int action, double val) {
        // E(S,A)=λγE(S,A)+1
        this.E[state][action] =  val;
    }
    
    
    public void updateQTable(double reward,int state,int action,int nextState,int nextAction) {
        // Q(S,A)=Q(S,A)+αδE(S,A)
        double delta = 0.0;
        delta = reward + this.gamma * this.qTable[nextState][nextAction] - this.qTable[state][action];
        
        this.qTable[state][action] = this.qTable[state][action] + this.alpha * delta * this.E[state][action];
    }
    
    public int getNextState(int state) {
        if(state == 0) {
            return 1;
        } else {
            return 0;
        }
    }
}
