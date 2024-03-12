/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */

package utils.adaptative;

/**
 *
 * @author J. Alfredo Brambila H. <alfredo.brambila@outlook.com>
 */
public class Window {

    private int index;
    private int rank;
    private double fitnessImprovement;
    
    public Window() {
        this.index = 0;
        this.rank = 0;
        this.fitnessImprovement = 0.0;
    }

    /**
     * @return the index
     */
    public int getIndex() {
        return index;
    }

    /**
     * @param index the index to set
     */
    public void setIndex(int index) {
        this.index = index;
    }

    /**
     * @return the rank
     */
    public int getRank() {
        return rank;
    }

    /**
     * @param rank the rank to set
     */
    public void setRank(int rank) {
        this.rank = rank;
    }

    /**
     * @return the fitnessImprovement
     */
    public double getFitnessImprovement() {
        return fitnessImprovement;
    }

    /**
     * @param fitnessImprovement the fitnessImprovement to set
     */
    public void setFitnessImprovement(double fitnessImprovement) {
        this.fitnessImprovement = fitnessImprovement;
    }
    
    
}
