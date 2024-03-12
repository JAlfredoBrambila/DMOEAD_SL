

package problem.multiobjective.dynamical.fda;

import solution.real.Individual;
import problem.multiobjective.dynamical.DynamicProblem;

/**
 *
 * @author J. Alfredo Brambila H. <alfredo.brambila@outlook.com>
 */
public class FDA3 implements DynamicProblem {
    
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
    int tauT = 20; // 20      tT  severidad de cambio
    int nT = 10; // 200     Tmax numero de generaciones que se queda fijo
    */
    
    int tauT; // 200      tT  Tmax numero de generaciones que se queda fijo
    int nT; // 10     severidad de cambio
    
    public FDA3() {
        problemName = "FDA3";
        nVars = 10;
        nObj = 2;
        lower = new double[nVars];
        upper = new double[nVars];
        
        tauT = 20; // tT  Tmax numero de generaciones que se queda fijo
        nT = 10; // severidad de cambio
        
        int XI = 1;
        int XII = 9;

        //lower[0] = 0.0;
        //upper[0] = 1.0;
        
        for(int i=0; i<XI; i++) {
            lower[i] = 0.0;
            upper[i] = 1.0;
        }
        
        for(int i=XI; i<XI+XII; i++) {
            lower[i] = -1.0;
            upper[i] = 1.0;
        }
        //time = 1.0d;
    }
    
    public FDA3(int tt, int nt) {
        problemName = "FDA3";
        nVars = 10;
        nObj = 2;
        lower = new double[nVars];
        upper = new double[nVars];
        
        tauT = tt; // 20      tT  Tmax numero de generaciones que se queda fijo
        nT = nt; // 200     severidad de cambio
        
        int XI = 1;
        int XII = 9;

        //lower[0] = 0.0;
        //upper[0] = 1.0;
        
        for(int i=0; i<XI; i++) {
            lower[i] = 0.0;
            upper[i] = 1.0;
        }
        
        for(int i=XI; i<XI+XII; i++) {
            lower[i] = -1.0;
            upper[i] = 1.0;
        }
        //time = 1.0d;
    }
    
    @Override
    public void costFunction(Individual ind) {
        
        // 1   1 1 1 1 1    1 1 1 1 1 1 1
        int n = ind.getPosition().length;
        
        
        double[] f = new double[2];
        f[0] = this.getF(ind, 0, 1);
        f[1] = 0.0;
        
        double g = this.getG(ind,1);
        double h = this.getH(f[0], g);

        f[1] = g * h;

        contEval++;
        ind.setCost(f);
    }
    
    public double getF(Individual ind, int lInf, int lSup) {
        
        double f = 0.0d;
        double aux = 2.0d * Math.sin(0.5d * Math.PI * time);
        double ft = Math.pow(10.0d, aux);
        
        for(int i=lInf; i<lSup; i++) {
            f += Math.pow(ind.getPosition()[i], ft);
        }
        return f;
    }
    
    public double getH(double f, double g_) {
        return 1.0d - Math.sqrt(f/g_);
    }
    
    public double getG(Individual ind, int lInf) {
        
        double gT = Math.abs(Math.sin(0.5 * Math.PI * time));
        double g_ = 0.0d;
        
        for(int i=lInf; i<ind.getPosition().length; i++) {
            g_ += Math.pow(ind.getPosition()[i] - gT, 2.0d);
        }
        
        return g_ + 1.0d + gT;
        //return (Math.sin(0.5 * Math.PI * time));
    }
    
    

    @Override
    public void update(Integer counter) {
        
        //time = (1.0d / (double) nT) * Math.floor(counter / (double) tauT);
        //time = (2 * Math.floor((double)counter / (double) tauT)) * ((double) tauT / ((double) nT - (double) tauT));
        /*
        time = (1.0d / (double) nT) * Math.floor(1.0d * (counter/10) / (double) tauT); // OK
        */
        
        time = (1.0d / (double) nT) * Math.floor(1.0d * (double)counter / (double) tauT);
        
        //System.out.println("IT: " + counter + " time: " + time);
        
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
