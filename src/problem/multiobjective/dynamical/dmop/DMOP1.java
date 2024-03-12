/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */

package problem.multiobjective.dynamical.dmop;

import solution.real.Individual;
import problem.multiobjective.dynamical.DynamicProblem;

/**
 *
 * @author J. Alfredo Brambila H. <alfredo.brambila@outlook.com>
 */
public class DMOP1 implements DynamicProblem {
    
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
    
    /*
    int tauT = 5*6; // 20      tT
    int nT = 10; // 200     Tmax
    */
    
    public DMOP1() {
        problemName = "DMOP1";
        nVars = 10;
        nObj = 2;
        lower = new double[nVars];
        upper = new double[nVars];
        
        tauT = 20; // tT  Tmax numero de generaciones que se queda fijo
	nT = 10; // severidad de cambio

        //lower[0] = 0.0;
        //upper[0] = 1.0;
        for(int i=0; i<nVars; i++) {
            lower[i] = 0.0;
            upper[i] = 1.0;
        }
        time = 0.0d;
    }
    
    public DMOP1(int tt, int nt) {
        problemName = "DMOP1";
        nVars = 10;
        nObj = 2;
        lower = new double[nVars];
        upper = new double[nVars];
        
        tauT = tt; // 20      tT  Tmax numero de generaciones que se queda fijo
	nT = nt; // 200     severidad de cambio

        //lower[0] = 0.0;
        //upper[0] = 1.0;
        for(int i=0; i<nVars; i++) {
            lower[i] = 0.0;
            upper[i] = 1.0;
        }
        time = 0.0d;
    }
    
    @Override
    public void costFunction(Individual ind) {
        int n = ind.getPosition().length;

        double g = 0.0;
        double h = 0.0;
        
        double[] f = new double[2];
        f[0] = 0.0;
        f[1] = 0.0;
        
        //f1
        f[0] = ind.getPosition()[0];
        
        //double sumXi=0;
        g = 0;
        for(int i=1; i<n; i++) {
            g += Math.pow(ind.getPosition()[i], 2.0d);
        }
        
        g = 1 + 9 * g;
        
        //g = 1 + 9 * sumXi;
        //g = 1 + 9 * sumXi;
        
        h = 1 - Math.pow(f[0] / g, getH());
        
        f[1] = g * h;
        
        contEval++;
        ind.setCost(f);
    }
    
    public double getH() {
        return 0.75 * Math.sin(0.5 * Math.PI * time) + 1.25;
    }
    
    public double getG() {
        return Math.abs(Math.sin(0.5 * Math.PI * time));
        //return (Math.sin(0.5 * Math.PI * time));
    }

    @Override
    public void update(Integer counter) {
        //time = (1.0d / (double) nT) * Math.floor(1.0d * counter / (double) tauT);
        
        /*
        time = (1.0d / (double) nT) * Math.floor(1.0d * (counter/10) / (double) tauT);
        */
        time = (1.0d / (double) nT) * Math.floor(1.0d * (double)counter / (double) tauT);


        //time = time * 300;
        //time = (2 * Math.floor((double) (counter/10) / (double) tauT)) * ((double) tauT / ((double) nT - (double) tauT));
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
