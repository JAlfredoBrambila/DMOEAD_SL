/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */

package solution.real;

/**
 *
 * @author J. Alfredo Brambila H. <alfredo.brambila@outlook.com>
 */
public class Particle extends Individual{
    private double[] velocity;
    private double[] bestPosition;
    private double[] bestCost;
    public Particle() {
        
    }
    
    public void clone(Particle a) {
        super.clone(a);
        velocity = new double[a.getVelocity().length];
        for(int i=0; i<a.getVelocity().length; i++) {
            velocity[i] = a.getVelocity()[i];
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
    
    public void cloneFromIndividual(Individual a) {
        super.clone(a);
        //velocity = new double[a.getVelocity().length];
        /*for(int i=0; i<a.getVelocity().length; i++) {
            velocity[i] = a.getVelocity()[i];
        }*/
        //
        bestPosition = new double[a.getPosition().length];
        for(int i=0; i<a.getPosition().length; i++) {
            bestPosition[i] = a.getPosition()[i];
        }
        //
        bestCost = new double[a.getCost().length];
        for(int i=0; i<a.getCost().length; i++) {
            bestCost[i] = a.getCost()[i];
        }
    }

    /**
     * @return the velocity
     */
    public double[] getVelocity() {
        return velocity;
    }

    /**
     * @param velocity the velocity to set
     */
    public void setVelocity(double[] velocity) {
        this.velocity = velocity;
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
        this.bestPosition = new double[bestPosition.length];
        for(int i=0; i<bestPosition.length; i++) {
            this.bestPosition[i] = bestPosition[i];
        }
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
        //this.bestCost = bestCost;
        this.bestCost = new double[bestCost.length];
        for(int i=0; i<bestCost.length; i++) {
            this.bestCost[i] = bestCost[i];
        }
    }
}
