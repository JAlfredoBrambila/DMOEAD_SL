
package operator.crossover;

import java.util.LinkedList;
import solution.real.Individual;
import utils.DoubleUtils;

/**
 *
 * @author J. Alfredo Brambila H. <alfredo.brambila@outlook.com>
 */
public class SBXCrossover implements CrossoverOperator<Individual> {

    double[] lo; 
    double[] up;
    public SBXCrossover() {
        
    }
    
    public SBXCrossover(double[] lower, double[] upper) {
        this.lo = lower;
        this.up = upper;
    }
    
    @Override
    public LinkedList<Individual> execute(LinkedList<Individual> parents) {
        
        //LinkedList<Individual> offspring = new LinkedList<Individual>();
        
        LinkedList<Individual> offspring = this.doCrossover(parents);
        
        /*Individual offspring1 = new Individual();
        offspring1.copyConfigAtributes(parents.get(0));
        Individual offspring2 = new Individual();
        offspring2.copyConfigAtributes(parents.get(1));

        int size = parents.get(0).getPosition().length;
        
        double nc = 1;
	double u = DoubleUtils.getRandomNumber0_1();
        
	double beta = 0;
	if(u <= 0.5) {
	    beta = Math.pow((2*u), (1/(nc+1)));
	} else {
	    beta = Math.pow((1/(2*(1-u))), (1/(nc+1)));
	}
        
        for(int i=0; i<size; i++) {
            offspring1.setPosition(i,DoubleUtils.repairSolutionVariableValue((0.5*((parents.get(0).getPosition()[i]+parents.get(1).getPosition()[i])-beta*Math.abs(parents.get(1).getPosition()[i]-parents.get(0).getPosition()[i]))), lo[i], up[i]));
            offspring2.setPosition(i,DoubleUtils.repairSolutionVariableValue((0.5*((parents.get(0).getPosition()[i]+parents.get(1).getPosition()[i])+beta*Math.abs(parents.get(1).getPosition()[i]-parents.get(0).getPosition()[i]))), lo[i], up[i]));
        }
        
        offspring.add(offspring1);
        offspring.add(offspring2);
        */
        
        
        
        return offspring;
        
    }
    
    public LinkedList<Individual> doCrossover(LinkedList<Individual> parents) {
        
        LinkedList<Individual> offspring = new LinkedList<Individual>();
        
        Individual offspring1 = new Individual();
        offspring1.copyConfigAtributes(parents.get(0));
        Individual offspring2 = new Individual();
        offspring2.copyConfigAtributes(parents.get(1));

        int size = parents.get(0).getPosition().length;
        
        double nc = 1;
	double u = DoubleUtils.getRandomNumber0_1();
        
	double beta = 0;
	if(u <= 0.5) {
	    beta = Math.pow((2*u), (1/(nc+1)));
	} else {
	    beta = Math.pow((1/(2*(1-u))), (1/(nc+1)));
	}
        
        for(int i=0; i<size; i++) {
            offspring1.setPosition(i,DoubleUtils.repairSolutionVariableValue((0.5*((parents.get(0).getPosition()[i]+parents.get(1).getPosition()[i])-beta*Math.abs(parents.get(1).getPosition()[i]-parents.get(0).getPosition()[i]))), lo[i], up[i]));
            offspring2.setPosition(i,DoubleUtils.repairSolutionVariableValue((0.5*((parents.get(0).getPosition()[i]+parents.get(1).getPosition()[i])+beta*Math.abs(parents.get(1).getPosition()[i]-parents.get(0).getPosition()[i]))), lo[i], up[i]));
        }
        
        offspring.add(offspring1);
        offspring.add(offspring2);
        
        return offspring;
        
    }

    @Override
    public Individual executeOnlyChild(LinkedList<Individual> parents) {
        
        Individual child = new Individual();
        
        LinkedList<Individual> offspring = this.doCrossover(parents);
        
        child.clone(offspring.get(0));
        
        for(int i=0; i<child.getPosition().length; i++) {
            if(DoubleUtils.getRandomNumber0_1() <= 0.5) {
                child.getPosition()[i] = offspring.get(1).getPosition()[i];
            } 
        }
        
        return child;
    }

    Individual current;
    @Override
    public void setCurrentIndividual(Individual c) {
        this.current = new Individual();
        this.current.clone(c);
    }

}
