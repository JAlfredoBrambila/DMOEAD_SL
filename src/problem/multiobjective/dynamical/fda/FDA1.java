

package problem.multiobjective.dynamical.fda;

import solution.real.Individual;
import problem.multiobjective.dynamical.DynamicProblem;

/**
 *
 * @author J. Alfredo Brambila H. <alfredo.brambila@outlook.com>
 */
public class FDA1 implements DynamicProblem {
    
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
    int tauT = 20; // 20      tT
    int nT = 10; // 200     Tmax
    */
    int tauT; // 200      tT  Tmax numero de generaciones que se queda fijo
    int nT; // 10     severidad de cambio
    
    public FDA1() {
        problemName = "FDA1";
        nVars = 20;
        nObj = 2;
        lower = new double[nVars];
        upper = new double[nVars];
        
        tauT = 20; // tT  Tmax numero de generaciones que se queda fijo
        nT = 10; // severidad de cambio

        lower[0] = 0.0d;
        upper[0] = 1.0d;
        for(int i=1; i<nVars; i++) {
            lower[i] = -1.0d;
            upper[i] = 1.0d;
        }
        //time = 1.0d;
    }
    
    public FDA1(int tt, int nt) {
        problemName = "FDA1";
        nVars = 20;
        nObj = 2;
        lower = new double[nVars];
        upper = new double[nVars];
        
        tauT = tt; // 20      tT  Tmax numero de generaciones que se queda fijo
        nT = nt; // 200     severidad de cambio

        lower[0] = 0.0d;
        upper[0] = 1.0d;
        for(int i=1; i<nVars; i++) {
            lower[i] = -1.0d;
            upper[i] = 1.0d;
        }
        //time = 1.0d;
    }
    
    @Override
    public void costFunction(Individual ind) {
        
        // 1   1 1 1 1 1    1 1 1 1 1 1 1
        int n = ind.getPosition().length;
        
        
        double[] f = new double[2];
        f[0] = ind.getPosition()[0];
        f[1] = 0.0;
        
        double g = this.getG(ind);
        double h = this.getH(f[0], g);

        f[1] = g * h;

        contEval++;
        ind.setCost(f);
    }
    
    public double getH(double f, double g_) {
        return 1.0d - Math.sqrt(f/g_);
    }
    
    public double getG(Individual ind) {
        
        double gT = Math.sin(0.5d * Math.PI * time);
        double g_ = 0.0d;
        
        for(int i=1; i<ind.getPosition().length; i++) {
            g_ += Math.pow(ind.getPosition()[i] - gT, 2.0d);
        }
        
        return g_ + 1.0d;
        //return (Math.sin(0.5 * Math.PI * time));
    }
    
    

    @Override
    public void update(Integer counter) {
        
        //time = (1.0d / (double) nT) * Math.floor(counter / (double) tauT);
        //time = (2 * Math.floor((double)counter / (double) tauT)) * ((double) tauT / ((double) nT - (double) tauT));
        
        /*
        time = (1.0d / (double) nT) * Math.floor(1.0d * (counter/10) / (double) tauT);
        */
        
        time = (1.0d / (double) nT) * Math.floor(1.0d * (double)counter / (double) tauT);
        
        //time = (2 * Math.floor((double) (counter/10) / (double) tauT)) * ((double) tauT / ((double) nT - (double) tauT));
        
        //// time = (1.0d / (double) nT) * Math.floor(value / (double) tauT);
        //time=2*Math.floor(counter/nT)*(tauT/(tauT*nT-nT));
        //System.out.println("Time: " + time);
        //setChanged();
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
