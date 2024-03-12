/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */

package algorithm.multiobjective.dynamical.changeResponse;

import java.util.LinkedList;
import solution.real.Individual;
import utils.DoubleUtils;
import problem.multiobjective.Problem;
import problem.multiobjective.dynamical.DynamicProblem;

/**
 *
 * @author J. Alfredo Brambila H. <alfredo.brambila@outlook.com>
 */
public class ResponseA implements ChangeResponseMechanism<Individual> {

    private double zeta;
    int nObjectives;
    DynamicProblem problem;
    public ResponseA(DynamicProblem problem, double zeta) {
        this.zeta = zeta;
        this.nObjectives = problem.getNumObjectives();
        this.problem = problem;
    }
    @Override
    public void execute(LinkedList<Individual> pop) {
        int nPop = pop.size();
        int nRemp = (int)( nPop * this.zeta);
        
        for(int i=0; i<nRemp; i++) {
            int j = DoubleUtils.getRandomIntNumber(0, nPop);
            // crea nuevo individuo
            Individual ind = new Individual();
            
            //genera valores reales aletatorios(entre varMin y varMax) para el cromosoma de tamaño nVar
            ind.setPosition(DoubleUtils.getRandomPos(problem.getLowerBound(), problem.getUpperBound(), problem.getNumVars()));
            
            // crea vector de costos de tamaño nObj para el individuo 
            ind.setCost(new double[this.nObjectives]);
            
            // aplica la funcion de costo y obtiene los nObj costos
            problem.costFunction(ind);
            ind.setIndex(pop.get(j).getIndex());
            
            //Agrega el nuevo individuo a la poblacion
            pop.set(j, ind);
            //pop.add(j, ind);
        }
        
        // set z
        /*for(Individual indUp : pop) {
            for(int jUp=0; jUp<nObjectives; jUp++) {
                z[jUp] = Math.min(z[jUp], indUp.getCost()[jUp]);
            }
        }*/
    }

}
