
package operator.mutation;

import java.util.LinkedList;
import solution.real.Individual;
import utils.DoubleUtils;

/**
 *
 * @author J. Alfredo Brambila H. <alfredo.brambila@outlook.com>
 */
public class RealMutation implements MutationOperator<Individual>{

    double[] lo; 
    double[] up;
    double muP; 
    double[] sigmaV;
    
    public RealMutation() {
        
    }
    
    //Individual ind, double mu, double[] sigmaV, double[] lo, double[] up
    //double muP, double[] lo, double[] up, double distIndex
    public RealMutation(double[] lower, double[] upper, double muP, double[] sigma) {
        this.lo = lower;
        this.up = upper;
        this.muP = muP;
        this.sigmaV = sigma;
    }
    
    @Override
    public Individual execute(Individual ind) {
        
        Individual ind2 = new Individual();
        ind2.clone(ind);
        
        int nVars = ind.getPosition().length;
        int nMu = (int) (Math.ceil(muP * nVars));
        
        int j=0;

        j = DoubleUtils.getRandomIntNumber(0,nVars);
        //j = Utils.getRandomIntNumber(nMu,nVars);
        
        
        
        double sigma = 0.0;
        
        if(sigmaV.length>1){
            sigma = sigmaV[j];
        } else {
            sigma = sigmaV[0];
        }
        
        double newVal = (ind.getPosition()[j] + sigma * Math.random());
        if(newVal < lo[j]) {
            newVal = lo[j];
        }
        if(newVal > up[j]) {
            newVal = up[j];
        }
        //sigma = 0.8;
        //System.out.println("Sigma " + sigma);
        ind2.copy(ind);
        ind2.setPosition(j, newVal);
        ind2.setIndex(ind.getIndex());
        
        return ind2;        
    }

}
