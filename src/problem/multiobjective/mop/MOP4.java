
package problem.multiobjective.mop;

import solution.real.Individual;
import problem.multiobjective.Problem;

/**
 *
 * @author J. Alfredo Brambila H. <alfredo.brambila@outlook.com>
 */
public class MOP4 implements Problem {
    //double[] rangeValues;
    double[] lower;
    double[] upper;
    String problemName;
    static int contEval=0;
    int nVars;
    int nObjectives;
    
    public MOP4() {

        problemName = "MOP4";
        nVars = 3;
        nObjectives = 2;
        upper = new double[nVars];
        lower = new double[nVars];

        for(int i=0; i<nVars; i++) {
            lower[i] = -5.0;
            upper[i] = 5.0;
        }
    }
    @Override
    public void costFunction(Individual ind) {
        double a = 0.8;
        int b = 3;
        
        int n = ind.getPosition().length;
        
        double[] sum = new double[2];
        sum[0] = 0.0;
        sum[1] = 0.0;
        for(int i=0; i<n-1; i++) {
            sum[0] += -10 * Math.exp(-0.2 * Math.sqrt(Math.pow(ind.getPosition()[i], 2) + Math.pow(ind.getPosition()[i+1], 2)));
        }
        
        for(int i=0; i<n; i++) {
            sum[1] += Math.pow(Math.abs(ind.getPosition()[i]), a) + 5 * Math.pow(Math.sin(ind.getPosition()[i]), b);
        }
        
        contEval++; 
        
        ind.setCost(sum);
    }
    
    @Override
    public int getContEvals() {
        return contEval;
    }

    @Override
    public double getVarMin() {
        return this.lower[0];
    }

    @Override
    public double getVarMax() {
        return this.upper[0];
    }

    @Override
    public String getProblemName() {
        return this.problemName;
    }

    @Override
    public int getNumVars() {
        return nVars;
    }

    @Override
    public double[] getLowerBound() {
        return this.lower;
    }

    @Override
    public double[] getUpperBound() {
        return this.upper;
    }
    
    @Override
    public int getNumObjectives() {
        return this.nObjectives;
    }

}
