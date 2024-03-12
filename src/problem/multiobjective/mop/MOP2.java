

package problem.multiobjective.mop;

import solution.real.Individual;
import problem.multiobjective.Problem;

/**
 *
 * @author J. Alfredo Brambila H. <alfredo.brambila@outlook.com>
 */
public class MOP2 implements Problem{
    
    //double[] rangeValues;
    double[] lower;
    double[] upper;
    String problemName;
    static int contEval=0;
    int nVars;
    int nObjectives;
    
    public MOP2() {

        problemName = "MOP2";
        nVars = 3; // 10
        nObjectives = 2;
        lower = new double[nVars];
        upper = new double[nVars];

        for(int i=0; i<nVars; i++) {
            lower[i] = -5.0;
            upper[i] = 5.0;
        }
    }

    @Override
    public void costFunction(Individual ind) {
        int n = ind.getPosition().length;
        
        double[] sum = new double[2];
        sum[0] = 0.0;
        sum[1] = 0.0;
        for(int i=0; i<n; i++) {
            sum[0] += (Math.pow(ind.getPosition()[i] - (1/Math.sqrt(n)), 2));
            sum[1] += (Math.pow(ind.getPosition()[i] + (1/Math.sqrt(n)), 2));
        }
        
        sum[0] *= -1;
        sum[1] *= -1;
        
        sum[0] = 1 - Math.exp(sum[0]);
        sum[1] = 1 - Math.exp(sum[1]);
        
        contEval++;
        ind.setCost(sum);
        
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
    public int getContEvals() {
        return contEval;
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
