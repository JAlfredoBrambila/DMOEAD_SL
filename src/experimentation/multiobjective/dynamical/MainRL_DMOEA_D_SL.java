/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */

package experimentation.multiobjective.dynamical;

import algorithm.multiobjective.dynamical.DynamicMultiobjectiveAlgorithm;
import algorithm.multiobjective.dynamical.dmoead.DMOEAD_SARSALambda;
import algorithm.multiobjective.dynamical.dmoead.DMOEAD_TAOS;
import operator.crossover.DiferentialEvolution1;
import operator.mutation.PolinomialMutation;
import problem.multiobjective.dynamical.DynamicProblem;
import problem.multiobjective.dynamical.df.DF10;
import problem.multiobjective.dynamical.df.DF12;
import problem.multiobjective.dynamical.df.DF4;
import problem.multiobjective.dynamical.df.DF6;
import problem.multiobjective.dynamical.dmop.DMOP1;
import problem.multiobjective.dynamical.dmop.DMOP2;
import problem.multiobjective.dynamical.fda.FDA1;
import problem.multiobjective.dynamical.fda.FDA3;

/**
 *
 * @author J. Alfredo Brambila H. <alfredo.brambila@outlook.com>
 */
public class MainRL_DMOEA_D_SL {

    public static void main(String[] args) {
        System.setProperty("java.util.Arrays.useLegacyMergeSort", "true");
        
        int NUM_EXP = 30; // number of experiments
        
        // experiment output folder
        String expFolder = "expMarch24";
        
        /*
        problem = new FDA1(tt,nt);
        problem = new FDA3(tt,nt);
        problem = new DMOP1(tt,nt);
        problem = new DMOP2(tt,nt);
        problem = new DF4(tt,nt);
        problem = new DF6(tt,nt);
        problem = new DF10(tt,nt);
        problem = new DF12(tt,nt);
        */
        
        // Dinamic parameters
        int maxItVar = 30; // maximum number of generations
        int tt=maxItVar; // change frequency
        int nt=10; // severity of change
        
        DynamicProblem problem;
        
        // uncomment and comment to use dynamic test problems
        problem = new FDA1(tt,nt);
        //problem = new FDA3(tt,nt);
        //problem = new DMOP1(tt,nt);
        //problem = new DMOP2(tt,nt);
        //problem = new DF4(tt,nt);
        //problem = new DF6(tt,nt);
        //problem = new DF10(tt,nt);
        //problem = new DF12(tt,nt);
        
        
        //DMOEA/D-TAOS
        for (int i = 0; i < NUM_EXP; i++) {
            DynamicMultiobjectiveAlgorithm algoritm = new DMOEAD_TAOS(
                    problem,// Problem
                    maxItVar, // Maximum Number of Iterations
                    100, // Population Size
                    100, // Archive Size
                    0.5, // gamma
                    0.1, // mu
                    0.2, // Zeta
                    //new SBXCrossover(problem.getLowerBound(),problem.getUpperBound()), // Crossover operator
                    new DiferentialEvolution1(problem.getLowerBound(), problem.getUpperBound()), // Default Crossover Operator
                    new PolinomialMutation(problem.getLowerBound(), problem.getUpperBound(), 0.2, 20.0), // Mutation operator
                    true, // Experimentation
                    (i+1) // Num Exp
            );

            algoritm.setExperimentationFolder(expFolder);
            algoritm.execute();
        }
        
        
        
        //DMOEA/D-SARSALambda
        for (int i = 0; i < NUM_EXP; i++) {
            DynamicMultiobjectiveAlgorithm algoritm = new DMOEAD_SARSALambda(
                    problem,// Problem
                    maxItVar, // Maximum Number of Iterations
                    100, // Population Size
                    100, // Archive Size
                    0.5, // gamma
                    0.1, // mu
                    0.2, // Zeta
                    //new SBXCrossover(problem.getLowerBound(),problem.getUpperBound()), // Crossover operator
                    new DiferentialEvolution1(problem.getLowerBound(), problem.getUpperBound()), // Default Crossover Operator
                    new PolinomialMutation(problem.getLowerBound(), problem.getUpperBound(), 0.2, 20.0), // Mutation operator
                    true, // Experimentation
                    (i+1) // Num Exp
            );

            algoritm.setExperimentationFolder(expFolder);
            algoritm.execute();
        }
        
        
        
    }
}
