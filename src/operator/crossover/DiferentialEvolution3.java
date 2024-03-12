
package operator.crossover;

import java.util.LinkedList;
import solution.real.Individual;
import utils.DoubleUtils;

/**
 *
 * @author J. Alfredo Brambila H. <alfredo.brambila@outlook.com>
 */
public class DiferentialEvolution3 implements CrossoverOperator<Individual> {

    double[] lo; 
    double[] up;
    public DiferentialEvolution3() {
        
    }
    
    public DiferentialEvolution3(double[] lower, double[] upper) {
        this.lo = lower;
        this.up = upper;
    }
    
    @Override
    public LinkedList<Individual> execute(LinkedList<Individual> parents) {
        
        throw new UnsupportedOperationException("This method returns only one child, please use the executeOnlyChild function");
        
    }
    
    public Individual doCrossover(LinkedList<Individual> parents) {
        
        Individual offspring = new Individual();
        int nVars = parents.get(0).getPosition().length;
        //double[] alpha = new double[nVars];
        double CR = 1.0;
        double F = 0.5;
        double K = 0.5;
        
        //offspring.clone(parents.get(0));
        offspring.clone(this.current);
        
        int k = DoubleUtils.getRandomIntNumber(0, nVars);
        
        for(int j=0; j<nVars; j++) {
            if(DoubleUtils.getRandomNumber0_1() < CR || j == k) {
                offspring.getPosition()[j] = this.current.getPosition()[j] + K * (parents.get(0).getPosition()[j]-this.current.getPosition()[j]) + F * (parents.get(1).getPosition()[j]-parents.get(2).getPosition()[j])  + F * (parents.get(3).getPosition()[j]-parents.get(4).getPosition()[j]);
                //offspring.getPosition()[j] = parents.get(1).getPosition()[j] + K * (parents.get(0).getPosition()[j]-parents.get(1).getPosition()[j]) + F * (parents.get(1).getPosition()[j]-parents.get(2).getPosition()[j])  + F * (parents.get(3).getPosition()[j]-parents.get(4).getPosition()[j]);
            } /*else {
                offspring.getPosition()[j] = parents.get(0).getPosition()[j];
            }*/
            
            if(offspring.getPosition()[j] < lo[j]) {
                double r2 = DoubleUtils.getRandomNumber0_1();
                offspring.getPosition()[j] = lo[j] + r2 * (parents.get(0).getPosition()[j] - lo[j]);
                //offspring.getPosition()[j] = lo[j];
            } 
            
            if(offspring.getPosition()[j] > up[j]) {
                double r3 = DoubleUtils.getRandomNumber0_1();
                offspring.getPosition()[j] = up[j] - r3 * (up[j] - parents.get(0).getPosition()[j]);
                //offspring.getPosition()[j] = up[j];
            } 
            
        }
        
        return offspring;
        
    }

    @Override
    public Individual executeOnlyChild(LinkedList<Individual> parents) {
        
        Individual child = this.doCrossover(parents);

        return child;
    }

    Individual current;
    
    @Override
    public void setCurrentIndividual(Individual c) {
        this.current = new Individual();
        this.current.clone(c);
    }

}
