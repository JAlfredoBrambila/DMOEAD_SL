
package problem.multiobjective.dynamical.df;

import problem.multiobjective.dynamical.DynamicProblem;
import solution.real.Individual;

/**
 *
 * @author J. Alfredo Brambila H. <alfredo.brambila@outlook.com>
 */
public class DF6 implements DynamicProblem {
    
    double[] lower;
    double[] upper;
    
    String problemName;
    static int contEval=0;
    int nVars;
    int nObj;
    
    //
    double time;
    boolean changeStatus = false;
    /*
    int tauT = 10;
    int nT = 10;
    */
    int tauT; // 10      tT  Tmax numero de generaciones que se queda fijo
    int nT; // 10     severidad de cambio
    
    double T0 = 50;
    
    public DF6() {
        problemName = "DF6";
        nVars = 10;
        nObj = 2;
        lower = new double[nVars];
        upper = new double[nVars];
        
        tauT = 10; // tT  Tmax numero de generaciones que se queda fijo
	nT = 10; // severidad de cambio

        lower[0] = 0.0d;
        upper[0] = 1.0d;
        for(int i=1; i<nVars; i++) {
            lower[i] = -1.0d;
            upper[i] = 1.0d;
        }
        
        time = 0.0d;
    }
    
    public DF6(int tt, int nt) {
        problemName = "DF6";
        nVars = 10;
        nObj = 2;
        lower = new double[nVars];
        upper = new double[nVars];
        
        tauT = tt; // tT  Tmax numero de generaciones que se queda fijo
	nT = nt; // severidad de cambio

        lower[0] = 0.0d;
        upper[0] = 1.0d;
        for(int i=1; i<nVars; i++) {
            lower[i] = -1.0d;
            upper[i] = 1.0d;
        }
        
        time = 0.0d;
    }

    @Override
    public void update(Integer counter) {
        //time = (1.0d / (double) nT) * Math.floor(1.0d * (counter/10) / (double) tauT);
        /*
        double tau_temp = Math.max((double)counter/10 + this.tauT - (T0+1), 0.0);
        this.time = (1.0d / this.nT) * (Math.floor(tau_temp/this.tauT));
        */
        time = (1.0d / (double) nT) * Math.floor(1.0d * (double)counter / (double) tauT);
    }

    @Override
    public void costFunction(Individual ind) {
        int n = ind.getPosition().length;
        
        /*double b = 1 + Math.abs(Math.cos(0.5d * Math.PI * this.time));
        double a = Math.sin(0.5d * Math.PI * this.time);
        double H = 1.5d + a;*/
        
        double G = Math.sin(0.5d * Math.PI * this.time);
        double a = 0.2d + 2.8d * Math.abs(G);
        
        double[] f = new double[2];
        f[0] = 0.0;
        f[1] = 0.0;
        
        double g = 1;
        
        for(int i=1; i<n; i++) {
            double y = ind.getPosition()[i] - G;
            //g += Math.pow(ind.getPosition()[i] - ((a * Math.pow(ind.getPosition()[0], 2))/(i+1)), 2);
            g+= Math.abs(G) * Math.pow(y, 2) - 10 * Math.cos( 2 * Math.PI * y) + 10;
        }
        
        f[0] = g * Math.pow(ind.getPosition()[0] + 0.1 * Math.sin(3 * Math.PI * ind.getPosition()[0]),a);
        f[1] = g * Math.pow(1 - ind.getPosition()[0] + 0.1 * Math.sin(3 * Math.PI * ind.getPosition()[0]),a);
        
        //f[0] = f[0] + 2*time;
        //f[1] = f[1] + 2*time;
        
        contEval++;
        ind.setCost(f);

    }

    @Override
    public boolean hasChanged() {
        return changeStatus;
    }

    @Override
    public void setChanged() {
        changeStatus = true;
    }

    @Override
    public void clearChanged() {
        changeStatus = false;
    }

    

    @Override
    public int getNumVars() {
        return this.nVars;
    }

    @Override
    public String getProblemName() {
        return this.problemName;
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
    public double[] getLowerBound() {
        return this.lower;
    }

    @Override
    public double[] getUpperBound() {
        return this.upper;
    }

    @Override
    public int getContEvals() {
        return contEval;
    }

    @Override
    public int getNumObjectives() {
        return this.nObj;
    }

}
