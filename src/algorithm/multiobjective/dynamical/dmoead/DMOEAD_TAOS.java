/**
 * Java multi-Objective Evolutionary Algorithm Mini Framework (JMOEAMF)
 * Clase: MOEAD.java
 * Paquete: algorithm.moead
 * Info: MOEA/D Multiobjective Evolutionary Algorithm Based on Decomposition
 * @version: 0.7
 * @autor: José Alfredo Brambila Hernámdez <alfredo.brambila@outlook.com>
 * @autor: Miguel Angel Garcia Morales <talivan12@hotmail.com>
 * @autor: Hector Fraire Huacuja <hector.fraire2014@gmail.com>
 * Proyecto Biblioteca de clases para JMOEAMF 
 * Desarrollo: Enero de 2022
 * Actualización: abril de 2022
 * Se permite el uso total o parcial de este código fuente siempre y cuando se le dé el crédito correspondiente a los autores
 *
 */

package algorithm.multiobjective.dynamical.dmoead;


import algorithm.multiobjective.moead.utils.SubProblem;
import algorithm.multiobjective.dynamical.DynamicMultiobjectiveAlgorithm;
//import static algorithm.multiobjective.dynamical.dmoead.DMOEAD.expNum;
import algorithm.multiobjective.utils.Front;
import java.nio.file.Paths;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;
import operator.crossover.CrossoverOperator;
import operator.crossover.DiferentialEvolution1;
import operator.crossover.DiferentialEvolution2;
import operator.crossover.DiferentialEvolution3;
import operator.crossover.DiferentialEvolution4;
import operator.crossover.SBXCrossover;
import operator.mutation.MutationOperator;
import operator.selection.SelectionOperator;
import org.apache.commons.math3.distribution.BetaDistribution;
import problem.multiobjective.Problem;
import problem.multiobjective.dynamical.DynamicProblem;
import problem.multiobjective.dynamical.fda.FDA2;
import solution.real.Individual;
import utils.Dominance;
import utils.DoubleUtils;
import utils.GeneralUtils;
import utils.PlotS;
import utils.adaptative.AdaptativeModule1;
import utils.adaptative.SlidingWindow;
import utils.adaptative.TS;
import utils.adaptative.Window;


/**
 *
 * @author J. Alfredo Brambila H. <alfredo.brambila@outlook.com>
 */
public class DMOEAD_TAOS implements DynamicMultiobjectiveAlgorithm<Individual>{

    int nVar; // tamaño cromosoma
    double[] varMin;
    double[] varMax;
    int nObjectives; // numero de objetivos
    int maxIt; //numero de iteraciones maximo
    int nPop;// tamaño de la poblacion
    int nArchive; 
    int T;
    double gamma;
    LinkedList<SubProblem> sp;
    double[] z;
    double[] sigma;
    double mu;

    
    LinkedList<Individual> pop; // Listado de poblacion POP
    LinkedList<Individual> ndPop;
    DynamicProblem problem; // Problema
    LinkedList<Individual> EP; // 
    
    //Dynamic
    double zeta;
    LinkedList<Individual> lookouts;
    
    static int expNum = 1;
    static int algConf = 4;
    boolean EXPERIMENTATION;
    
    CrossoverOperator crossover;
    MutationOperator mutation;
    
    //Adaptative
    int numStrategies; // Number of operators
    int[] SUC;
    double FEtest;
    double FEapply;
    double k;
    
    int[] usage;
    
    
    String carpetaExperimentos;
    
    int numAlg;


    
    public DMOEAD_TAOS(DynamicProblem problem, int maxIt, int nPOP, int nArchive, double gamma, double mu, double zeta, CrossoverOperator crossover, MutationOperator mutation) {
        this.problem = problem;
        this.nVar = problem.getNumVars();
        this.varMin = problem.getLowerBound();
        this.varMax = problem.getUpperBound();
        this.nObjectives = problem.getNumObjectives();
        this.maxIt = maxIt;
        this.nPop = nPOP;
        this.nArchive = nArchive;
        this.T = (int) Math.max(Math.ceil(0.15 * nPop), 2); // numero de vecinos
        this.T = (int) Math.min(Math.max(T, 2), 15);
        this.gamma = gamma;
        this.mu = mu;
        this.zeta = zeta;

        //this.selection = selection;
        this.crossover = crossover;
        this.mutation = mutation;

        //
        sp = createSubProblems(); // Create Sub-problems
        z = zeros(this.nObjectives); // Initialize Goal Point
        sigma = DoubleUtils.getSigmaValue(0.1, varMin, varMax, nVar);
        pop = new LinkedList<Individual>();
        EP = new LinkedList<Individual>();
        ndPop = new LinkedList<Individual>();

        
        //Adaptative
        numStrategies = 4;
        SUC = new int[numStrategies];
        k = 1.0;
        
        
        EXPERIMENTATION = false;

    }
    
    public DMOEAD_TAOS(DynamicProblem problem, int maxIt, int nPOP, int nArchive, double gamma, double mu, double zeta, CrossoverOperator crossover, MutationOperator mutation, boolean experimentation, int exp) {
        this.problem = problem;
        this.nVar = problem.getNumVars();
        this.varMin = problem.getLowerBound();
        this.varMax = problem.getUpperBound();
        this.nObjectives = problem.getNumObjectives();
        this.maxIt = maxIt;
        this.nPop = nPOP;
        this.nArchive = nArchive;
        this.T = (int) Math.max(Math.ceil(0.15 * nPop), 2); // numero de vecinos
        this.T = (int) Math.min(Math.max(T, 2), 15);
        this.gamma = gamma;
        this.mu = mu;
        this.zeta = zeta;

        //this.selection = selection;
        this.crossover = crossover;
        this.mutation = mutation;

        //
        sp = createSubProblems(); // Create Sub-problems
        z = zeros(this.nObjectives); // Initialize Goal Point
        sigma = DoubleUtils.getSigmaValue(0.1, varMin, varMax, nVar);
        pop = new LinkedList<Individual>();
        EP = new LinkedList<Individual>();
        ndPop = new LinkedList<Individual>();

        
        //Adaptative
        numStrategies = 4;
        SUC = new int[numStrategies];
        k = 1.0;
        usage = new int[numStrategies];
        
        EXPERIMENTATION = experimentation;
        expNum = exp;
        numAlg = 4;
        
        
    }
    
    public void setExperimentationFolder(String folderName) {
        this.carpetaExperimentos = folderName;
    }
    
    public DMOEAD_TAOS() {
        
        // seleccionar problema dinamico [ FDA1 | FDA2 | FDA3 | DMOP1 |DMOP2 ]
        problem = new FDA2();
        
        System.out.println("Selected problem " + problem.getProblemName());
        
        nVar = problem.getNumVars();
        varMin = problem.getLowerBound();
        varMax = problem.getUpperBound();
        nObjectives = problem.getNumObjectives();
        maxIt = 200; //100
        nPop = 100;//50
        nArchive = 100;
        T = (int)Math.max(Math.ceil(0.15*nPop), 2); // numero de vecinos
        T = (int)Math.min(Math.max(T, 2), 15);
        gamma = 0.5;
        sp = createSubProblems(); // Create Sub-problems
        z = zeros(this.nObjectives); // Initialize Goal Point
        
        mu = 0.1;
        sigma = DoubleUtils.getSigmaValue(0.1, varMin, varMax, nVar);
        
        pop = new LinkedList<Individual>();
        EP = new LinkedList<Individual>();
        ndPop = new LinkedList<Individual>();
        
        //Dynamic
        zeta = 0.2; 
        
        EXPERIMENTATION = false;
    }

    /**
     * Funcion principal de la clase MOEA/D
     */
    public void execute() {
        //String carpetaExperimentos="ExpAGO23";
        String nombreAlgoritmo = "DMOEAD_A_TAOS_";
        
        if(this.carpetaExperimentos == null || this.carpetaExperimentos.isEmpty()) {
            this.carpetaExperimentos = "ExpGenericRoot";
        }
        
        problem.update(0);
        // Adaptative
        int strategy_flag;
        for(int i=0; i<numStrategies; i++)
            SUC[i] = 0;
        FEtest = numStrategies * this.nPop;
        FEapply = k * FEtest; // revisar K y k
        int e = 0;
        
        
        
        generatePop();
        
        // set g
        for(int i=0; i<this.pop.size(); i++) {
            pop.get(i).setG(decomposedCost(pop.get(i),z,sp.get(i).getLambda()));
        }
        
         
        determineDomination(pop); // Determine Population Domination Status
        EP = getNonDominated(pop); // Initialize Estimated Pareto Front
        
        // *************************
        // *** Problema DINAMICO ***
        // *************************
        // Se toma el 10% de la población para usar como detectores
        lookouts = getLookouts(pop);
        PlotS plot = new PlotS();
        
        int k1=0;
        int k2=0;
        int j1=0;
        int j2=0;
         
        Individual y;
        
        double C = 100.0;
        LinkedList<TS> pds = createPDSList(this.numStrategies,1.0,1.0);
        
        int contW=1;
        for(int it = 0; it<this.maxIt*5; it++) {
            
            
            // Adaptative
            
            for(int i=0; i<this.nPop; i++) {
                
                
                
                // Adaptative
                if(e % (this.FEtest + this.FEapply) < this.FEtest ) {
                    if(e % (this.FEtest + this.FEapply) == 0) {
                        for(int ii=0; ii<numStrategies; ii++)
                            SUC[ii] = 0;
                    }
                    
                    strategy_flag = (e % this.numStrategies) + 1;
                } else {
                    strategy_flag = argMax(SUC);
                }
                
                //strategy_flag = this.betaSample(pds);
                
                usage[strategy_flag-1] = usage[strategy_flag-1]+1;
                
                y = this.matingEvolution(strategy_flag, i);
                problem.costFunction(y);
                
                //z = Math.min(z, y.get)
                for(int j=0; j<nObjectives; j++) {
                    z[j] = Math.min(z[j], y.getCost()[j]);
                }
                
                int SUI = 0;
                for(double j : sp.get(i).getNeighbors()) {
                    y.setG(this.decomposedCost(y, z, sp.get((int)j).getLambda()));
                    
                     
                    // Adaptative
                    
                    if(y.getG() <= pop.get((int)j).getG()) {
                        pop.get((int)j).clone(y);
                        SUI = 1;
                         
                    }
                    
                    SUC[strategy_flag-1] = SUC[strategy_flag-1] + SUI;
                    e++;
                }
                
                

            }
            
            determineDomination(pop);
            ndPop = getNonDominated(pop);
            
            DoubleUtils.mergeList(EP, ndPop);
            
            determineDomination(EP);
            EP = getNonDominated(EP);
            
            //EP = getNonDominated(EP);

            int valIndice=0;
            if(EP.size() > this.nArchive) {                
                while(EP.size() > this.nArchive) {
                    valIndice = DoubleUtils.getRandomIntNumber(0, EP.size());
                    EP.remove(valIndice);
                }
            }
            
            
            // *************************
            // *** Problema DINAMICO ***
            // *************************
            
            // Actualiza tiempo en cada iteracion
            problem.update(it);
            //System.out.println("UP: " + it);
            
            //DIN
            // Detectores de cambio en funcion objetivo "lookouts"
            // hasChangedObjectiveFunction revisa si hubo cambios en la funcion objetivo.
            // Se comparan los costos de los vigias con una nueva evaluacion de los mismos
            // si los costos son diferentes, el problema cambió
            if (hasChangedObjectiveFunction(lookouts)) {
                /*System.out.println("*******************************************");
                System.out.print("change in objective function occurred ");
                System.out.println("Iteration number: " + it);*/
                
                // Adaptative
                //this.refreshment();
                
                plot.addSerieM(EP, contW);
                
                
                
                /////////
                if (EXPERIMENTATION) {
                    GeneralUtils.SaveExperimentationFile(carpetaExperimentos, problem, contW, nombreAlgoritmo, numAlg, EP, expNum);
                } else {
                    // imprime informacion de poblacion
                    System.out.println("\n*** EP List ***\n");
                    DoubleUtils.printList(EP);
                    System.out.println("Objetivos: ");
                    DoubleUtils.printObj(EP);
                    System.out.println("");
                    System.out.println("");
                }
                /////////
                
                
                
                //sp = createSubProblems(); // Create Sub-problems
                z = zeros(this.nObjectives); // Initialize Goal Point
                
                //System.out.println("POP: " + pop.size());
                // reevaluar POP
                evaluatePop();//
                
                //System.out.println("POP: " + pop.size());
                // Regenerar población inicial con el metodo A
                // reemplaza 10% de la poblacion de forma aleatoria
                this.generarePopA();
                
                //System.out.println("POP: " + pop.size() + " sp: " + sp.size());
                for(int i=0; i<this.pop.size(); i++) {
                    //System.out.println("i: " + i);
                    pop.get(i).setG(decomposedCost(pop.get(i),z,sp.get(i).getLambda()));
                }
                
                
                determineDomination(pop);
                ndPop = getNonDominated(pop);

                EP.clear();
                DoubleUtils.mergeList(EP, ndPop);

                EP = getNonDominated(EP);

                EP = getNonDominated(EP);

                valIndice = 0;
                if (EP.size() > this.nArchive) {
                    while (EP.size() > this.nArchive) {
                        valIndice = DoubleUtils.getRandomIntNumber(0, EP.size());
                        EP.remove(valIndice);
                    }
                }
                
                
                
                // contador de número de cambios
                contW++;
                
                
            }
            
            
        }
        
        
        
        ////////////
        if (EXPERIMENTATION) {
            GeneralUtils.SaveExperimentationFile(carpetaExperimentos, problem, contW, nombreAlgoritmo, numAlg, EP, expNum);
            expNum++;
        } else {
            
            System.out.println("EP ");
            for (Individual ind : EP) {
                System.out.println(ind);
            }

            plot.addSerieM(EP, contW);
            //Plot.plotFront(EP);
            plot.plotSeries();
        }
        ///////////
        
    }
    
    /**
     * Genera la poblacion inicial aleatoria
     */
    public void generatePop() {
        for(int i=0; i<this.nPop; i++) {
            // crea nuevo individuo
            Individual ind = new Individual();
            
            //genera valores reales aletatorios(entre varMin y varMax) para el cromosoma de tamaño nVar
            ind.setPosition(DoubleUtils.getRandomPos(varMin, varMax, nVar));
            
            // crea vector de costos de tamaño nObj para el individuo 
            ind.setCost(new double[this.nObjectives]);
            
            // aplica la funcion de costo y obtiene los nObj costos
            problem.costFunction(ind);
            ind.setIndex(i);
            
            // set z
            for(int j=0; j<nObjectives; j++) {
                z[j] = Math.min(z[j], ind.getCost()[j]);
            }
            
            
            //Agrega el nuevo individuo a la poblacion
            pop.add(ind);
        }
    }
    
    public int argMax(int[] ops) {
        double max = ops[0];
        int op=0;
        for(int i=1; i<ops.length; i++) {
            if(ops[i] > max) {
                max = ops[i];
                op = i;
            }
        }
        
        
        return op+1; // max
    }
    
    // DYTS 
    public int betaSample(LinkedList<TS> pds) {
        int n = pds.size();
        double[] ops = new double[n];
        BetaDistribution beta;
        for(int i=0; i<n; i++) {
            beta = new BetaDistribution(pds.get(i).getA(), pds.get(i).getB());
            ops[i] = beta.inverseCumulativeProbability(Math.random());
        }
        
        int op = 0;
        
        double max = ops[0];
        for(int i=1; i<n; i++) {
            if(ops[i] > max) {
                max = ops[i];
                op = i;
            }
        }
        
        
        return op+1; // max
    }
    
    public LinkedList<TS> createPDSList(int nStrategies, double a, double b) {
        LinkedList<TS> p = new LinkedList<TS>();
        for(int i=0; i<nStrategies; i++) {
            p.add(new TS(a,b));
        }
        
        return p;
    }
    
    public Individual matingEvolution(int strategySelected, int i) {
        Individual child;
        
        
        switch(strategySelected) {
            case 1:
                child = this.operator1(i);
                break;
            case 2:
                child = this.operator2(i);
                break;
            case 3:
                child = this.operator3(i);
                break;
            case 4:
                child = this.operator4(i);
                break;
            default:
                child = this.operator1(i);
                break;
        }
        
        return child;
    }

    public Individual operator1(int i) {

        // DE/rand/1/bin: u = x1 + F * (x2 - x3)
        this.crossover = new DiferentialEvolution1(problem.getLowerBound(),problem.getUpperBound());
        this.crossover.setCurrentIndividual(this.pop.get(i));
        //this.crossover = new SBXCrossover(problem.getLowerBound(),problem.getUpperBound());
        
        int k1 = 0;
        int k2 = 0;
        int k3 = 0;
        int j1 = 0;
        int j2 = 0;
        int j3 = 0;
        Individual p1;
        Individual p2;
        Individual p3;
        Individual y;
        
        k1 = DoubleUtils.getRandomIntNumber(0, T);
        k2 = DoubleUtils.getRandomIntNumber(0, T);
        k3 = DoubleUtils.getRandomIntNumber(0, T);
        
        while (k1 == k2) {
            k2 = DoubleUtils.getRandomIntNumber(0, T);
        }
        
        
        while (k3 == k2 || k3 == k1) {
            k3 = DoubleUtils.getRandomIntNumber(0, T);
        }

        j1 = (int) sp.get(i).getNeighbors()[k1];
        p1 = new Individual();
        p1.clone(pop.get(j1));

        j2 = (int) sp.get(i).getNeighbors()[k2];
        p2 = new Individual();
        p2.clone(pop.get(j2));

        j3 = (int) sp.get(i).getNeighbors()[k3];
        p3 = new Individual();
        p3.clone(pop.get(j3));

        y = new Individual();

        //y.clone(Operator.crossover(p1, p2, gamma, varMin, varMax));
        LinkedList<Individual> parentL = new LinkedList<Individual>();
        parentL.add(p1);
        parentL.add(p2);
        parentL.add(p3);

        y.clone((Individual) this.crossover.executeOnlyChild(parentL));

        y.clone((Individual) this.mutation.execute(y));
        
        return y;
    }
    
    public Individual operator2(int i) {

        //DE/rand/2/bin: u = xi + F * (x1 - x2) + F * (x3 - x4) 
        
        this.crossover = new DiferentialEvolution2(problem.getLowerBound(),problem.getUpperBound());
        this.crossover.setCurrentIndividual(this.pop.get(i));
        
        //this.crossover = new SBXCrossover(problem.getLowerBound(),problem.getUpperBound());
        
        int k1 = 0;
        int k2 = 0;
        int k3 = 0;
        int k4 = 0;
        int j1 = 0;
        int j2 = 0;
        int j3 = 0;
        int j4 = 0;
        Individual p1;
        Individual p2;
        Individual p3;
        Individual p4;
        
        Individual y;
        
        k1 = DoubleUtils.getRandomIntNumber(0, T);
        k2 = DoubleUtils.getRandomIntNumber(0, T);
        k3 = DoubleUtils.getRandomIntNumber(0, T);
        k4 = DoubleUtils.getRandomIntNumber(0, T);
        
        while (k1 == k2) {
            k2 = DoubleUtils.getRandomIntNumber(0, T);
        }
        
        
        while (k3 == k2 || k3 == k1) {
            k3 = DoubleUtils.getRandomIntNumber(0, T);
        }
        
        while (k4 == k2 || k4 == k1 || k4 == k3) {
            k4 = DoubleUtils.getRandomIntNumber(0, T);
        }

        j1 = (int) sp.get(i).getNeighbors()[k1];
        p1 = new Individual();
        p1.clone(pop.get(j1));

        j2 = (int) sp.get(i).getNeighbors()[k2];
        p2 = new Individual();
        p2.clone(pop.get(j2));

        j3 = (int) sp.get(i).getNeighbors()[k3];
        p3 = new Individual();
        p3.clone(pop.get(j3));
        
        j4 = (int) sp.get(i).getNeighbors()[k4];
        p4 = new Individual();
        p4.clone(pop.get(j4));

        y = new Individual();

        //y.clone(Operator.crossover(p1, p2, gamma, varMin, varMax));
        LinkedList<Individual> parentL = new LinkedList<Individual>();
        parentL.add(p1);
        parentL.add(p2);
        parentL.add(p3);
        parentL.add(p4);

        //this.crossover.setCurrentIndividual(this.pop.get(i));
        
        y.clone((Individual) this.crossover.executeOnlyChild(parentL));

        y.clone((Individual) this.mutation.execute(y));

        return y;
    }
    
    public Individual operator3(int i) {

        this.crossover = new DiferentialEvolution3(problem.getLowerBound(),problem.getUpperBound());
        this.crossover.setCurrentIndividual(this.pop.get(i));
        
        int k1 = 0;
        int k2 = 0;
        int k3 = 0;
        int k4 = 0;
        int k5 = 0;
        
        int j1 = 0;
        int j2 = 0;
        int j3 = 0;
        int j4 = 0;
        int j5 = 0;
        
        Individual p1;
        Individual p2;
        Individual p3;
        Individual p4;
        Individual p5;
        
        Individual y;
        
        k1 = DoubleUtils.getRandomIntNumber(0, T);
        k2 = DoubleUtils.getRandomIntNumber(0, T);
        k3 = DoubleUtils.getRandomIntNumber(0, T);
        k4 = DoubleUtils.getRandomIntNumber(0, T);
        k5 = DoubleUtils.getRandomIntNumber(0, T);
        
        while (k1 == k2) {
            k2 = DoubleUtils.getRandomIntNumber(0, T);
        }
        
        
        while (k3 == k2 || k3 == k1) {
            k3 = DoubleUtils.getRandomIntNumber(0, T);
        }
        
        while (k4 == k2 || k4 == k1 || k4 == k3) {
            k4 = DoubleUtils.getRandomIntNumber(0, T);
        }
        
        while (k5 == k2 || k5 == k1 || k5 == k3 || k5 == k4) {
            k5 = DoubleUtils.getRandomIntNumber(0, T);
        }

        j1 = (int) sp.get(i).getNeighbors()[k1];
        p1 = new Individual();
        p1.clone(pop.get(j1));

        j2 = (int) sp.get(i).getNeighbors()[k2];
        p2 = new Individual();
        p2.clone(pop.get(j2));

        j3 = (int) sp.get(i).getNeighbors()[k3];
        p3 = new Individual();
        p3.clone(pop.get(j3));
        
        j4 = (int) sp.get(i).getNeighbors()[k4];
        p4 = new Individual();
        p4.clone(pop.get(j4));
        
        j5 = (int) sp.get(i).getNeighbors()[k5];
        p5 = new Individual();
        p5.clone(pop.get(j5));

        y = new Individual();

        //y.clone(Operator.crossover(p1, p2, gamma, varMin, varMax));
        LinkedList<Individual> parentL = new LinkedList<Individual>();
        parentL.add(p1);
        parentL.add(p2);
        parentL.add(p3);
        parentL.add(p4);
        parentL.add(p5);

        //this.crossover.setCurrentIndividual(this.pop.get(i));
        
        y.clone((Individual) this.crossover.executeOnlyChild(parentL));

        y.clone((Individual) this.mutation.execute(y));

        return y;
    }
    
    public Individual operator4(int i) {

        
        
        this.crossover = new DiferentialEvolution4(problem.getLowerBound(),problem.getUpperBound());
        this.crossover.setCurrentIndividual(this.pop.get(i));
        //this.crossover = new SBXCrossover(problem.getLowerBound(),problem.getUpperBound());
        
        int k1 = 0;
        int k2 = 0;
        int k3 = 0;
        int j1 = 0;
        int j2 = 0;
        int j3 = 0;
        Individual p1;
        Individual p2;
        Individual p3;
        Individual y;
        
        k1 = DoubleUtils.getRandomIntNumber(0, T);
        k2 = DoubleUtils.getRandomIntNumber(0, T);
        k3 = DoubleUtils.getRandomIntNumber(0, T);
        
        while (k1 == k2) {
            k2 = DoubleUtils.getRandomIntNumber(0, T);
        }
        
        
        while (k3 == k2 || k3 == k1) {
            k3 = DoubleUtils.getRandomIntNumber(0, T);
        }

        j1 = (int) sp.get(i).getNeighbors()[k1];
        p1 = new Individual();
        p1.clone(pop.get(j1));

        j2 = (int) sp.get(i).getNeighbors()[k2];
        p2 = new Individual();
        p2.clone(pop.get(j2));

        j3 = (int) sp.get(i).getNeighbors()[k3];
        p3 = new Individual();
        p3.clone(pop.get(j3));

        y = new Individual();

        //y.clone(Operator.crossover(p1, p2, gamma, varMin, varMax));
        LinkedList<Individual> parentL = new LinkedList<Individual>();
        parentL.add(p1);
        parentL.add(p2);
        parentL.add(p3);
        
        //this.crossover.setCurrentIndividual(this.pop.get(i));

        y.clone((Individual) this.crossover.executeOnlyChild(parentL));

        y.clone((Individual) this.mutation.execute(y));
        
        return y;
    }
    
    
    public LinkedList<SubProblem> createSubProblems() {
        LinkedList<SubProblem> subProblems = new LinkedList<SubProblem>();
        SubProblem subProb;
        
        double[] lambda = new double[this.nObjectives];
        double norma = 0.0;
        
        for (int i = 0; i < this.nPop; i++) {
            subProb = new SubProblem();
            lambda = new double[this.nObjectives];
            for (int j = 0; j < this.nObjectives; j++) {
                lambda[j] = DoubleUtils.getRandomNumber0_1();
            }

            // devuelve la norma euclidiana del vector v
            norma = euclideanNorm(lambda);
            //norma = unoNorm(lambda);
            for (int j = 0; j < this.nObjectives; j++) {
                lambda[j] = lambda[j] / norma;
            }
            
            
            
            subProb.setLambda(lambda);
            
            subProblems.add(subProb);
        }
        
       
        
        double[][] D = euclideanDistance(subProblems);
        double[] columna;
        double[] SO;
        for (int i = 0; i < this.nPop; i++) {
            columna = null;
            columna = DoubleUtils.getColumFromMTX(D, i);
            //SO = getSortArrayPos(D[i]);
            //Utils.printArray(columna);
            SO = getSortArrayPos(columna);
            //Utils.printArray(SO);
            subProblems.get(i).setNeighbors(SO);
            //System.out.println("i: " + i+ " > " + subProblems.get(i));
        }
        
        
        return subProblems;
    }
    
    public boolean inList(LinkedList<Individual> l, Individual b) {
        boolean inList = false;
        for(Individual ind : l) {
            if(ind.cmp(b)) {
                return true;
            }
        }
        return inList;
    }
    
    public double[] getSortArrayPos(double[] a) {
        
        double[] posArray = new double[T];
        double[][] arr = new double[a.length][2];
        
        for(int i=0; i<a.length; i++) {
            arr[i][0] = a[i];
            arr[i][1] = i;
            
        }
        
        DoubleUtils.Sort2DArrayBasedOnColumnNumber(arr, 1);
        
        //System.out.println("***");
        for(int i=0; i<T; i++) {
            posArray[i] = arr[i][1];
            
        }
        
        return posArray;
    }
    
    public double euclideanNorm(double[] v) {
        
        double sum = 0.0;
        for(int i=0; i<v.length; i++) {
            sum += Math.pow(v[i], 2);
        }
        
        return Math.sqrt(sum);
    }
    
    public double unoNorm(double[] v) {
        
        double sum = 0.0;
        for(int i=0; i<v.length; i++) {
            sum += v[i];
        }
        
        return sum;
    }
    
    public double[][] euclideanDistance(LinkedList<SubProblem> x) {
        double[][] d = new double[x.size()][x.size()];
        
        
        
        double sum = 0.0;
        for(int i=0; i<x.size(); i++) {
            for(int j=0; j<x.size(); j++) {
                sum = 0.0;
                
                for(int k=0; k<this.nObjectives; k++) {
                    
                    sum += Math.pow(x.get(i).getLambda()[k] - x.get(j).getLambda()[k], 2);
                    
                }
                d[j][i] = Math.sqrt(sum);
            }
        }
        
        return d;
    }
    
    public double[] zeros(int n) {
        double[] v = new double[n];
        for(int i=0; i<n; i++) {
            //v[i] = Double.POSITIVE_INFINITY;
            v[i] = 0.0;
        }
        
        return v;
    }
    
    public void determineDomination(LinkedList<Individual> pList) {
        for(int i=0; i<pList.size(); i++) {
            pList.get(i).setIsDominated(false);
        }
        
        for(int i=0; i<pList.size(); i++) {
            for(int j=i+1; j<pList.size(); j++) {
                //System.out.println(">>>> i: " + i + " j: " + j + " z: " + pList.size());
                if(Dominance.dominates2(pList.get(i), pList.get(j))) {
                    pList.get(j).setIsDominated(true);
                } else if(Dominance.dominates2(pList.get(j), pList.get(i))) {
                    pList.get(i).setIsDominated(true);
                }
            }
        }
    }
    
    public LinkedList<Individual> getNonDominated(LinkedList<Individual> pList) {
        LinkedList<Individual> nondominated = new LinkedList<Individual>();
        
        Individual indTemp;
        for(int i=0; i<pList.size(); i++) {
            if(!pList.get(i).isIsDominated()) {
                indTemp = new Individual();
                indTemp.clone(pList.get(i));
                nondominated.add(indTemp);
            }
        }
        
        return nondominated;
    }
    
    public double decomposedCost(Individual ind, double[] z, double[] lambda) {
        
        //double[] fx = new double[this.nObjectives];
        ArrayList<Double> fxl = new ArrayList<Double>(); 
        
        for(int i=0; i<this.nObjectives; i++) {
            //fx[i] = lambda[i] * Math.abs(ind.getCost()[i] - z);
            fxl.add(lambda[i] * Math.abs(ind.getCost()[i] - z[i]));
        }
        
        
        return Collections.max(fxl);
    }
    
    // Dynamic
    public LinkedList<Individual> getLookouts(LinkedList<Individual> pop) {
        LinkedList<Individual> l = new LinkedList<Individual>();
        int nLookouts = (int) (pop.size() * 0.1);
        int n = pop.size();
        for(int i=0; i<nLookouts; i++) {
            int j = DoubleUtils.getRandomIntNumber(0, n);
            Individual iL = new Individual();
            iL.clone(pop.get(i));
            l.add(iL);
        }
        return l;
    }
    
    
    
    // Updated July 2023 
    public boolean hasChangedObjectiveFunction(LinkedList<Individual> l) {
        boolean hasChanged = false;
        
        
        double[][] values = new double [this.nObjectives][l.size()];
        
        for(int i=0; i<l.size(); i++) {
            for(int obj=0; obj<this.nObjectives; obj++) {
                values[obj][i] = l.get(i).getCost()[obj];
            }
        }
        
        evaluateListInd(l);
        
        for(int i=0; i<l.size(); i++) {
            for(int obj=0; obj<this.nObjectives; obj++) {
                if(l.get(i).getCost()[obj] != values[obj][i]) {
                    return true;
                }
            }
        }
        
        return hasChanged;
    }
    
    public void evaluateListInd(LinkedList<Individual> p) {
        for(int i=0; i<p.size(); i++) {
            problem.costFunction(p.get(i));
        }
    }
    
    public void evaluatePop() {
        for(int i=0; i<this.pop.size(); i++) {
            problem.costFunction(pop.get(i));
        }
    }
    
    public void generarePopA() {
        int nRemp = (int)(this.nPop * this.zeta);
        for(int i=0; i<nRemp; i++) {
            int j = DoubleUtils.getRandomIntNumber(0, nPop);
            // crea nuevo individuo
            Individual ind = new Individual();
            
            //genera valores reales aletatorios(entre varMin y varMax) para el cromosoma de tamaño nVar
            ind.setPosition(DoubleUtils.getRandomPos(varMin, varMax, nVar));
            
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
        for(Individual ind : pop) {
            for(int j=0; j<nObjectives; j++) {
                z[j] = Math.min(z[j], ind.getCost()[j]);
            }
        }
            
    }
    
}
