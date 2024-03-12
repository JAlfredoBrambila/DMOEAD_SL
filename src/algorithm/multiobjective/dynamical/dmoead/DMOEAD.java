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
//import static algorithm.multiobjective.dynamical.dmoead.DMOEAD_DYTS.expNum;
import algorithm.multiobjective.utils.Front;
import java.nio.file.Paths;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;
import operator.crossover.CrossoverOperator;
import operator.crossover.DiferentialEvolution1;
import operator.crossover.DiferentialEvolution3;
import operator.mutation.MutationOperator;
import operator.selection.SelectionOperator;
import problem.multiobjective.Problem;
import problem.multiobjective.dynamical.DynamicProblem;
import problem.multiobjective.dynamical.fda.FDA2;
import solution.real.Individual;
import utils.Dominance;
import utils.DoubleUtils;
import utils.GeneralUtils;
import utils.PlotS;


/**
 *
 * @author J. Alfredo Brambila H. <alfredo.brambila@outlook.com>
 */
public class DMOEAD implements DynamicMultiobjectiveAlgorithm<Individual>{

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
    boolean EXPERIMENTATION = false;
    
    CrossoverOperator crossover;
    MutationOperator mutation;
    
    public DMOEAD(DynamicProblem problem, int maxIt, int nPOP, int nArchive, double gamma, double mu, double zeta, CrossoverOperator crossover, MutationOperator mutation) {
        this.problem = problem;
        this.nVar = problem.getNumVars();
        this.varMin = problem.getLowerBound();
        this.varMax = problem.getUpperBound();
        this.nObjectives = problem.getNumObjectives();
        this.maxIt = maxIt;
        this.nPop = nPOP;
        this.nArchive = nArchive;
        this.T = (int)Math.max(Math.ceil(0.15*nPop), 2); // numero de vecinos
        this.T = (int)Math.min(Math.max(T, 2), 15);
        this.gamma = gamma;
        this.mu = mu;
        this.zeta = zeta; 
        
        /*this.pCrossover = pcrossover;
        this.nCrossover = 2 * Math.round(pCrossover * nPop / 2);
        this.pMutation = pmutation;
        this.nMutation = Math.round(pMutation * nPop);
        this.mu = mu;
        this.sigmaC = sigmaC;*/
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
 
    }
    
    public DMOEAD(DynamicProblem problem, int maxIt, int nPOP, int nArchive, double gamma, double mu, double zeta, CrossoverOperator crossover, MutationOperator mutation, boolean experimentation, int exp ) {
        this.problem = problem;
        this.nVar = problem.getNumVars();
        this.varMin = problem.getLowerBound();
        this.varMax = problem.getUpperBound();
        this.nObjectives = problem.getNumObjectives();
        this.maxIt = maxIt;
        this.nPop = nPOP;
        this.nArchive = nArchive;
        this.T = (int)Math.max(Math.ceil(0.15*nPop), 2); // numero de vecinos
        this.T = (int)Math.min(Math.max(T, 2), 15);
        this.gamma = gamma;
        this.mu = mu;
        this.zeta = zeta; 
        
        /*this.pCrossover = pcrossover;
        this.nCrossover = 2 * Math.round(pCrossover * nPop / 2);
        this.pMutation = pmutation;
        this.nMutation = Math.round(pMutation * nPop);
        this.mu = mu;
        this.sigmaC = sigmaC;*/
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
 
        EXPERIMENTATION = experimentation;
        expNum = exp;
    }
    
    String carpetaExperimentos;
    
    public void setExperimentationFolder(String folderName) {
        this.carpetaExperimentos = folderName;
    }
    
    public DMOEAD() {
        
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
    }

    /**
     * Funcion principal de la clase MOEA/D
     */
    public void execute() {
        
        String carpetaExperimentos="ExpAGO23";
        String nombreAlgoritmo = "DMOEAD_A_";
        
        problem.update(0);
        generatePop();
        
        // set g
        for(int i=0; i<this.pop.size(); i++) {
            pop.get(i).setG(decomposedCost(pop.get(i),z,sp.get(i).getLambda()));
        }
        
        //problem.update(1800);
        determineDomination(pop); // Determine Population Domination Status
        EP = getNonDominated(pop); // Initialize Estimated Pareto Front
        
        // *************************
        // *** Problema DINAMICO ***
        // *************************
        // Se toma el 10% de la población para usar como "vigias"
        lookouts = getLookouts(pop);
        PlotS plot = new PlotS();
        
        int k1=0;
        int k2=0;
        int j1=0;
        int j2=0;
        Individual p1;
        Individual p2;
        Individual y;
        
        int contW=1;
        for(int it = 0; it<this.maxIt*5; it++) {
            for(int i=0; i<this.nPop; i++) {
                // sp
                if (EXPERIMENTATION) {
                    y = new Individual();
                    y.clone(operatorDE_1(i));
                    /*
                
                     */
                } else {
                    /*k1 = DoubleUtils.getRandomIntNumber(0, T);
                    k2 = DoubleUtils.getRandomIntNumber(0, T);
                    while (k1 == k2) {
                        k2 = DoubleUtils.getRandomIntNumber(0, T);
                    }

                    j1 = (int) sp.get(i).getNeighbors()[k1];
                    p1 = new Individual();
                    p1.clone(pop.get(j1));

                    j2 = (int) sp.get(i).getNeighbors()[k2];
                    p2 = new Individual();
                    p2.clone(pop.get(j2));

                    y = new Individual();

                    LinkedList<Individual> parentL = new LinkedList<Individual>();
                    parentL.add(p1);
                    parentL.add(p2);

                    y.clone((Individual) this.crossover.executeOnlyChild(parentL));

                    //y.clone(Operator.crossoverSBXM(p1, p2, varMin, varMax));
                    //y.clone(Operator.DifferentialEvolutionCrossoverBIN(p1, p2, p3, varMin, varMax)); 
                    */
                    
                     y = new Individual();
                    y.clone(operatorDE_1(i));
                }
                
                

                //y.clone(Operator.mutation(y, mu, sigma, varMin, varMax));
                y.clone((Individual) this.mutation.execute(y));
                //y.clone(Operator.polinomialMutation(y, 0.1, varMin, varMax, 20.0));
                
                problem.costFunction(y);
                
                //z = Math.min(z, y.get)
                for(int j=0; j<nObjectives; j++) {
                    z[j] = Math.min(z[j], y.getCost()[j]);
                }
                
                for(double j : sp.get(i).getNeighbors()) {
                    y.setG(this.decomposedCost(y, z, sp.get((int)j).getLambda()));
                    if(y.getG() <= pop.get((int)j).getG()) {
                        pop.get((int)j).clone(y);
                    }
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
                System.out.println("*******************************************");
                System.out.print("change in objective function occurred ");
                System.out.println("Iteration number: " + it);
                
                
                plot.addSerieM(EP, contW);
                
                
                
                ///////
                if (EXPERIMENTATION) {
                    GeneralUtils.SaveExperimentationFile(carpetaExperimentos, problem, contW, nombreAlgoritmo, 7, EP, expNum);
                } else {
                    // imprime informacion de poblacion
                    System.out.println("\n*** EP List ***\n");
                    DoubleUtils.printList(EP);
                    System.out.println("Objetivos: ");
                    DoubleUtils.printObj(EP);
                    System.out.println("");
                    System.out.println("");
                }
                ///////
                
                
                
                //sp = createSubProblems(); // Create Sub-problems
                z = zeros(this.nObjectives); // Initialize Goal Point
                
                //System.out.println("POP: " + pop.size());
                // reevaluar POP
                evaluatePop();// 
                
                //System.out.println("POP: " + pop.size());
                // Regenerar población inicial con el metodo A
                // reemplaza 10% de la poblacion de forma aleatoria
                this.generarePopA(); // aqui actualiza Z
                //this.generarePopB();
                
               
                
                //System.out.println("POP: " + pop.size() + " sp: " + sp.size());
                for(int i=0; i<this.pop.size(); i++) {
                    //System.out.println("i: " + i);
                    pop.get(i).setG(decomposedCost(pop.get(i),z,sp.get(i).getLambda()));
                    //pop.get(i).setG(decomposedCostPBI(pop.get(i),z,sp.get(i).getLambda(), zp, dz));
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
        
        
        
        ///////////////
        if (EXPERIMENTATION) {
            GeneralUtils.SaveExperimentationFile(carpetaExperimentos, problem, contW, nombreAlgoritmo, 7, EP, expNum);
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
        ///////////////
        
    }
    
    /*public void saveToFileExp() {
        String nombreAlgoritmo = "DMOEAD_A";
        String ruta=Paths.get("").toAbsolutePath().toString() + "\\fronts\\" + problem.getProblemName()+"\\";
    }*/
    
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
    
    public Individual operatorDE_1(int i) {

        // DE/rand/1/bin: u = x1 + F * (x2 - x3)
        this.crossover = new DiferentialEvolution3(problem.getLowerBound(),problem.getUpperBound());
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
        
        ////
        int k4 = 0;
        int k5 = 0;
        int j4 = 0;
        int j5 = 0;
        
        Individual p4;
        Individual p5;
        
        k4 = DoubleUtils.getRandomIntNumber(0, T);
        k5 = DoubleUtils.getRandomIntNumber(0, T);
        ////
        
        while (k1 == k2) {
            k2 = DoubleUtils.getRandomIntNumber(0, T);
        }
        
        
        while (k3 == k2 || k3 == k1) {
            k3 = DoubleUtils.getRandomIntNumber(0, T);
        }
        
        ////
        while (k4 == k3 || k4 == k2 || k4 == k1) {
            k4 = DoubleUtils.getRandomIntNumber(0, T);
        }
        
        while (k5 == k4 ||k5 == k3 || k5 == k2 || k5 == k1) {
            k5 = DoubleUtils.getRandomIntNumber(0, T);
        }
        ////

        j1 = (int) sp.get(i).getNeighbors()[k1];
        p1 = new Individual();
        p1.clone(pop.get(j1));

        j2 = (int) sp.get(i).getNeighbors()[k2];
        p2 = new Individual();
        p2.clone(pop.get(j2));

        j3 = (int) sp.get(i).getNeighbors()[k3];
        p3 = new Individual();
        p3.clone(pop.get(j3));
        
        //////
        j4 = (int) sp.get(i).getNeighbors()[k4];
        p4 = new Individual();
        p4.clone(pop.get(j4));
        
        j5 = (int) sp.get(i).getNeighbors()[k5];
        p5 = new Individual();
        p5.clone(pop.get(j5));
        //////

        y = new Individual();

        //y.clone(Operator.crossover(p1, p2, gamma, varMin, varMax));
        LinkedList<Individual> parentL = new LinkedList<Individual>();
        parentL.add(p1);
        parentL.add(p2);
        parentL.add(p3);
        
        /////
        parentL.add(p4);
        parentL.add(p5);
        /////

        y.clone((Individual) this.crossover.executeOnlyChild(parentL));

        //y.clone((Individual) this.mutation.execute(y));
        
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
            
            /*double a = 1.0 * i / (this.nPop - 1);
            lambda[0] = a;
            lambda[1] = 1.0 - a;*/
            
            subProb.setLambda(lambda);
            //System.out.println("Add: " + subProb.getLambda()[0] + " "  + subProb.getLambda()[1]);
            //System.out.println("ADD\t"+subProb.getLambda()[0] + "\t"  + subProb.getLambda()[1]);
            subProblems.add(subProb);
        }
        
        //
        /*for (int i = 0; i < subProblems.size(); i++) {
            System.out.println("::sp: " + subProblems.get(i).getLambda()[0] + " " + subProblems.get(i).getLambda()[1]);
        }*/
        
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
        
        //
        //System.out.println("SubProblems: " + subProblems.size());
        /*System.out.println("Subproblems:");
        for(SubProblem s : subProblems) {
            System.out.println(s.getLambda()[0] + " " + s.getLambda()[1]);
            System.out.print("Ne: ");
            for(Double ne: s.getNeighbors()) {
                System.out.print(" " + ne);
            }
            System.out.println("");
        }
        System.out.println("End subproblems");*/
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
        //System.out.println("--");
        double[] posArray = new double[T];
        double[][] arr = new double[a.length][2];
        
        for(int i=0; i<a.length; i++) {
            arr[i][0] = a[i];
            arr[i][1] = i;
            //System.out.println(":: " + arr[i][0] + " " + arr[i][1]);
        }
        
        DoubleUtils.Sort2DArrayBasedOnColumnNumber(arr, 1);
        
        //System.out.println("***");
        for(int i=0; i<T; i++) {
            posArray[i] = arr[i][1];
            //System.out.println(":: " + arr[i][0] + " " + arr[i][1]);
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
        
        /*for(int i=0; i<x.size(); i++) {
            System.out.println(":: " + x.get(i).getLambda()[0] + " " + x.get(i).getLambda()[1]);
        }*/
        
        double sum = 0.0;
        for(int i=0; i<x.size(); i++) {
            for(int j=0; j<x.size(); j++) {
                sum = 0.0;
                //System.out.println("0:: " + x.get(i).getLambda()[0] + " " + x.get(i).getLambda()[1] + " - "  + x.get(j).getLambda()[0] + " " + x.get(j).getLambda()[1]);
                for(int k=0; k<this.nObjectives; k++) {
                    //System.out.println("1:: " + x.get(i).getLambda()[0] + " " + x.get(i).getLambda()[1] + " - "  + x.get(j).getLambda()[0] + " " + x.get(j).getLambda()[1]);
                    sum += Math.pow(x.get(i).getLambda()[k] - x.get(j).getLambda()[k], 2);
                    //System.out.println("2:: " + x.get(i).getLambda()[0] + " " + x.get(i).getLambda()[1] + " - "  + x.get(j).getLambda()[0] + " " + x.get(j).getLambda()[1]);
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
    
    /*public double decomposedCostPBI(Individual ind, double[] z, double[] lambda, double[] zp, double[] dz) {
        
        
        
        
        double sum = 0;
        double d1 = 0;
        double d2 = 0;
        double theta = 2.0;
        
        //for(int i=0; i<this.nObjectives; i++) {
            //d1 = Math.abs(ind.getCost()[i] - zp[i]);
            d1 = dotProduct(subtract(ind.getCost(),zp), lambda);
            d2 = dotProduct(subtract(ind.getCost(),zp), normalize(dz));
            
            //sum += lambda[i] * d1 + theta * lambda[i] * d2;
            sum += d1 + theta * d2;
        //}
        
        d1 = dotProduct(subtract(ind.getCost(), zp), lambda);
        d2 = dotProduct(subtract(ind.getCost(), zp), normalize(dz));

        //sum += lambda[i] * d1 + theta * lambda[i] * d2;
        sum = d1 + theta * d2;

        System.out.println("");
        System.out.println("------->: " + sum);
        System.out.println("");
        System.exit(1);
        
        return sum;
        
    }
    
    private double dotProduct(double[] nd, double[] theta) {
        double sum = 0.0;

        for (int i = 0; i < nd.length; i++) {
            sum += nd[i] * theta[i];
        }

        return sum;
    }
    
    private double[] subtract(double[] a, double[] b) {
        double[] result = new double[a.length];

        for (int i = 0; i < a.length; i++) {
            result[i] = a[i] - b[i];
        }

        return result;
    }
    
    private double[] normalize(double[] vec) {
        double norm = 0.0;

        for (int i = 0; i < vec.length; i++) {
            norm += vec[i] * vec[i];
        }

        norm = Math.sqrt(norm);

        double[] result = new double[vec.length];

        for (int i = 0; i < vec.length; i++) {
            result[i] = vec[i] / norm;
        }

        return result;
    }
    
    private double[] getZP() {
        double[] zp = new double[this.nObjectives];
        
        for (int i = 0; i < this.nObjectives; i++) {
            zp[i] = DoubleUtils.getRandomNumber0_1();
        }
        
        return zp;
    }

  */
    
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
    
    ///////////////////////////////////////////////////////////////////
    // modificar N objetivos
    /////////////////////////////////////////////////////////////////
    public boolean hasChangedObjectiveFunctionOld(LinkedList<Individual> l) {
        boolean hasChanged = false;
        double[] valuesO1 = new double[l.size()];
        double[] valuesO2 = new double[l.size()];
        int i=0;
        for(Individual v : l) {
            valuesO1[i] = v.getCost()[0];
            valuesO2[i++] = v.getCost()[1];
        }
        
        evaluateListInd(l);
        
        i=0;
        for(Individual v : l) {
            if(v.getCost()[0] != valuesO1[i] || v.getCost()[1] != valuesO2[i++])
                return true;
        }
        
        return hasChanged;
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
    
    public void generarePopB() {
        int nRemp = (int)(this.nPop * this.zeta);
        for(int i=0; i<nRemp; i++) {
            int j = DoubleUtils.getRandomIntNumber(0, nPop);
            // crea nuevo individuo
            ///Individual ind = new Individual();
            
            //genera valores reales aletatorios(entre varMin y varMax) para el cromosoma de tamaño nVar
            ///ind.setPosition(DoubleUtils.getRandomPos(varMin, varMax, nVar));
            
            // crea vector de costos de tamaño nObj para el individuo 
            ///ind.setCost(new double[this.nObjectives]);
            
            pop.get(j).clone((Individual)this.mutation.execute(pop.get(j)));
            
            // aplica la funcion de costo y obtiene los nObj costos
            problem.costFunction(pop.get(j));
            //ind.setIndex(pop.get(j).getIndex());
            
            //Agrega el nuevo individuo a la poblacion
            //pop.set(j, ind);
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
