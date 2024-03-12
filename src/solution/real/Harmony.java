/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */

package solution.real;

/**
 *
 * @author J. Alfredo Brambila H. <alfredo.brambila@outlook.com>
 */
public class Harmony extends Individual {
    
    private double[] tuning;
    private double[] bestPosition;
    private double[] bestCost;
    
    public void clone(Harmony a) {
        super.clone(a);
        
        tuning = new double[a.getTuning().length];
        for(int i=0; i<a.getTuning().length; i++) {
            tuning[i] = a.getTuning()[i];
        }
        
         //
        bestPosition = new double[a.getBestPosition().length];
        for(int i=0; i<a.getBestPosition().length; i++) {
            bestPosition[i] = a.getPosition()[i];
        }
        //
        bestCost = new double[a.getBestCost().length];
        for(int i=0; i<a.getBestCost().length; i++) {
            bestCost[i] = a.getBestCost()[i];
        }
        
    }

    /**
     * @return the tuning
     */
    public double[] getTuning() {
        return tuning;
    }

    /**
     * @param tuning the tuning to set
     */
    public void setTuning(double[] tuning) {
        this.tuning = tuning;
    }

    /**
     * @return the bestPosition
     */
    public double[] getBestPosition() {
        return bestPosition;
    }

    /**
     * @param bestPosition the bestPosition to set
     */
    public void setBestPosition(double[] bestPosition) {
        this.bestPosition = bestPosition;
    }

    /**
     * @return the bestCost
     */
    public double[] getBestCost() {
        return bestCost;
    }

    /**
     * @param bestCost the bestCost to set
     */
    public void setBestCost(double[] bestCost) {
        this.bestCost = bestCost;
    }

}
