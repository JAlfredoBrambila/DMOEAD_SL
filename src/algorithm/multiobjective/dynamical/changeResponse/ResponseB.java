/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */

package algorithm.multiobjective.dynamical.changeResponse;

import java.util.LinkedList;
import operator.mutation.MutationOperator;
import operator.mutation.PolinomialMutation;
import solution.real.Individual;
import utils.DoubleUtils;
import problem.multiobjective.Problem;
import problem.multiobjective.dynamical.DynamicProblem;

/**
 *
 * @author J. Alfredo Brambila H. <alfredo.brambila@outlook.com>
 */
public class ResponseB implements ChangeResponseMechanism<Individual> {

    private double zeta;
    int nObjectives;
    DynamicProblem problem;
    MutationOperator mutation;
    
    public ResponseB(DynamicProblem problem, double zeta) {
        this.zeta = zeta;
        this.nObjectives = problem.getNumObjectives();
        this.problem = problem;
        this.mutation = new PolinomialMutation(problem.getLowerBound(),problem.getUpperBound(),0.2,20.0);
    }
    @Override
    public void execute(LinkedList<Individual> pop) {
        int nPop = pop.size();
        int nRemp = (int)( nPop * this.zeta);
        
        for(int i=0; i<nRemp; i++) {
            int j = DoubleUtils.getRandomIntNumber(0, nPop);
            
            Individual ind = (Individual)this.mutation.execute(pop.get(j));
            pop.get(j).clone(ind);
            // aplica la funcion de costo y obtiene los nObj costos
            problem.costFunction(pop.get(j));
        }
        
        // set z
        /*for(Individual indUp : pop) {
            for(int jUp=0; jUp<nObjectives; jUp++) {
                z[jUp] = Math.min(z[jUp], indUp.getCost()[jUp]);
            }
        }*/
    }

}
