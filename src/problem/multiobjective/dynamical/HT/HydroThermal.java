/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */

package problem.multiobjective.dynamical.HT;

import solution.real.Individual;
import problem.multiobjective.dynamical.DynamicProblem;

/**
 *
 * @author J. Alfredo Brambila H. <alfredo.brambila@outlook.com>
 */
public class HydroThermal implements DynamicProblem {

    double[] lower;
    double[] upper;
    
    String problemName;
    static int contEval=0;
    int nVars;
    int nObj;
    
    int M;
    int Ns;
    int Nh;
    
    double tm = 12;
    
    // Hydroelectric system data
    double[][] hsd = { //a0h    a1h   a2h      Wh    Phmin  Phmax
                        {260.0, 8.5, 0.00986, 125000, 0.0, 250.0},
                        {250.0, 9.8, 0.01140, 286000, 0.0, 500.0}   
    };
    
    //Cost related thermal system data
    double[][] crtsd = {    //a     bs    cs      ds   es     Psmin  Psmax
                            {60.0 , 1.8 , 0.0030, 140, 0.040, 20.0, 125.0},
                            {100.0, 2.1 , 0.0012, 160, 0.038, 30.0, 175.0},
                            {120.0, 2.0 , 0.0010, 180, 0.037, 40.0, 250.0},
                            {40.0 , 1.8 , 0.0015, 200, 0.035, 50.0, 300.0} 
    };
    
    // Emission related thermal system data
    double[][] ertsd = { //  αs	     βs	    γs	   ηs	     δs
                            {50.0, -0.555, 0.0150, 0.5773, 0.02446},
                            {60.0, -1.355, 0.0105, 0.4968, 0.02270},
                            {45.0, -0.600, 0.0080, 0.4860, 0.01948},
                            {30.0, -0.555, 0.0120, 0.5035, 0.02075}
    };
    
    // Load demands and goals of the objectives
    double[][] ldgo = {
                            {12.0, 900.0,  70718.4, 23200.8},
                            {12.0, 1100.0, 70718.4, 23200.8},
                            {12.0, 1000.0, 70718.4, 23200.8},
                            {12.0, 1300.0, 70718.4, 23200.8}
    };
    
    // Transmission loss formula coefficients
    double[][] B = {
                            {0.000049, 0.000014, 0.000015, 0.000015, 0.000020, 0.000017},
                            {0.000014, 0.000045, 0.000016, 0.000020, 0.000018, 0.000015},
                            {0.000015, 0.000016, 0.000039, 0.000010, 0.000012, 0.000012},
                            {0.000015, 0.000020, 0.000010, 0.000040, 0.000014, 0.000010},
                            {0.000020, 0.000018, 0.000012, 0.000014, 0.000035, 0.000011},
                            {0.000017, 0.000015, 0.000012, 0.000010, 0.000011, 0.000036}
                    };
    
    public HydroThermal() {
        this.problemName = "Hydro-Thermal";
        this.nObj = 2;
        M = 4;
        Ns = this.crtsd.length;
        Nh = this.hsd.length;
        
        this.nVars = (Nh + Ns) * M ;
        
        System.out.println("NVars: " + nVars);
        
        this.lower = new double[this.nVars];
        this.upper = new double[this.nVars];
        
        //Bounds
        int pos = 0;
        for (int m = 0; m < M; m++) {
            pos = m * (Nh + Ns);
            //Hydro
            for (int i = 0; i < this.hsd.length; i++) {
                this.lower[i+pos] = hsd[i][hsd[0].length-2];
                this.upper[i+pos] = hsd[i][hsd[0].length-1];
            }
            //Thermal
            for (int i = 0; i < this.crtsd.length; i++) {
                this.lower[i+Nh+pos] = crtsd[i][crtsd[0].length-2];
                this.upper[i+Nh+pos] = crtsd[i][crtsd[0].length-1];
            }
        }
        
        /*for(int i=0; i<this.lower.length; i++) {
            System.out.print(lower[i] +"\t");
        }
        System.out.println("");
        for(int i=0; i<this.upper.length; i++) {
            System.out.print(upper[i] +"\t");
        }
        System.out.println("");*/
        
    }
    
    @Override
    public void costFunction(Individual ind) {
        double[] f = new double[this.nObj];

        
        
        f[0] = 0.0;
        f[1] = 0.0;
        for(int m=0; m<M; m++) {
            for(int s=0; s<Ns; s++) {
                // F1
                f[0] += tm * (crtsd[s][0] + crtsd[s][1] * ind.getPosition()[6*m+s+2] + crtsd[s][2] * Math.pow(ind.getPosition()[6*m+s+2], 2) + Math.abs(crtsd[s][3] * Math.sin(crtsd[s][4]*(crtsd[s][5]-ind.getPosition()[6*m+s+2]))));
                // F2
                f[1] += tm * (ertsd[s][0] + ertsd[s][1] * ind.getPosition()[6*m+s+2] + ertsd[s][2] * Math.pow(ind.getPosition()[6*m+s+2], 2) + ertsd[s][3] * Math.exp(ertsd[s][4] * ind.getPosition()[6*m+s+2]));
            }
        }
        
        ind.setCost(f);
        
        evaluateConstraints(ind);
    }
    
    public void evaluateConstraints(Individual ind) {
        double[] PL = new double[M];
        double[] SumPs = new double[M];
        double[] SumPh = new double[M];
        int Nhs = Nh+Ns;
        
        // Power Balance Constraits
        for (int m = 0; m < M; m++) {
            PL[m] = 0.0;
            SumPs[m] = 0.0;
            SumPh[m] = 0.0;
            for (int i = 0; i < Nhs; i++) {
                for (int j = 0; j < Nhs; j++) {
                    //System.out.println("N: " + Nhs +" i: " + i + " j: " + j + " m: " + m);
                    PL[m] += ind.getPosition()[6 * m + i] * B[i][j] * ind.getPosition()[6 * m + j];
                }
            }
            
        }
        
        double PBC=0.0;
        for(int m=0; m<M; m++) {
            for (int s = 0; s < Ns; s++) {
                SumPs[m] += ind.getPosition()[6 * m + s + 2];
            }
            
            for (int h = 0; h < Nh; h++) {
                SumPh[m] += ind.getPosition()[6 * m + h]; // - ldgo[m][1] - PL[m];
            }
            PBC += SumPs[m] +  SumPh[m]  - ldgo[m][1] - PL[m];
            System.out.println(SumPs[m] +"\t"+ SumPh[m] +"\t"+ ldgo[m][1] +"\t"+ PL[m]);
        }
        System.out.println("SUM: " +  (int)PBC);
        
        // Water availability constraits
        double WAC = 0.0;
        for (int h = 0; h < Nh; h++) {
            WAC = 0.0;
            for (int m = 0; m < M; m++) {

                
                WAC += (tm * (hsd[h][0] + (hsd[h][1] * ind.getPosition()[6 * m + h]) + (hsd[h][2] * Math.pow(ind.getPosition()[6 * m + h], 2))));

            }
            System.out.println("W: " + (int)(WAC-hsd[h][3]));
        }
        System.out.println("SUM: " + WAC);
    }
    
    @Override
    public void update(Integer counter) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public boolean hasChanged() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void setChanged() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void clearChanged() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public int getNumVars() {
        return this.nVars;
    }

    @Override
    public String getProblemName() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public double getVarMin() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public double getVarMax() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
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
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    @Override
    public int getNumObjectives() {
        return this.nObj;
    }

}
