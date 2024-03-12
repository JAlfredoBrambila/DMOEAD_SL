
package operator.mutation;

import java.util.LinkedList;
import solution.real.Individual;
import utils.DoubleUtils;

/**
 *
 * @author J. Alfredo Brambila H. <alfredo.brambila@outlook.com>
 */
public class PolinomialMutation implements MutationOperator<Individual>{

    double[] lo; 
    double[] up;
    double muP; 
    double distIndex;
    
    public PolinomialMutation() {
        
    }
    
    //double muP, double[] lo, double[] up, double distIndex
    public PolinomialMutation(double[] lower, double[] upper, double muP, double distIndex) {
        this.lo = lower;
        this.up = upper;
        this.muP = muP;
        this.distIndex = distIndex;
    }
    
    @Override
    public Individual execute(Individual ind) {
        
        Individual ind2 = new Individual();
        
        ind2.copy(ind);
        //ind2.setPosition(j, (ind.getPosition()[j] + sigma * Math.random()));
        ind2.setIndex(ind.getIndex());
        double rnd, delta1, delta2, mutPow, deltaq;
        double y, yl, yu, val, xy;
        
        for(int i=0; i<ind.getPosition().length; i++) {
            if(DoubleUtils.getRandomNumber0_1() <= muP) {
                y = ind2.getPosition()[i];
                yl = lo[i];
                yu = up[i];
                if(yl == yu) {
                    y = yl;
                } else {
                    delta1 = (y - yl) / (yu - yl);
                    delta2 = (yu - y) / (yu - yl);
                    rnd = DoubleUtils.getRandomNumber0_1();
                    mutPow = 1.0 / (distIndex + 1.0);
                    
                    if (rnd <= 0.5) {
                        xy = 1.0 - delta1;
                        val = 2.0 * rnd + (1.0 - 2.0 * rnd) * (Math.pow(xy, distIndex + 1.0));
                        deltaq = Math.pow(val, mutPow) - 1.0;
                    } else {
                        xy = 1.0 - delta2;
                        val = 2.0 * (1.0 - rnd) + 2.0 * (rnd - 0.5) * (Math.pow(xy, distIndex + 1.0));
                        deltaq = 1.0 - Math.pow(val, mutPow);
                    }
                    y = y + deltaq * (yu - yl);
                    
                    y = DoubleUtils.repairSolutionVariableValue(y, yl, yu);
                    //y = solutionRepair.repairSolutionVariableValue(y, yl, yu);
                    //System.out.println("X");
                }
                
                ind2.getPosition()[i] = y;
            }
        }
        return ind2;
        
    }

}
