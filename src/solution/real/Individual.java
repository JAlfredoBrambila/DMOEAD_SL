/**
 * Java multi-Objective Evolutionary Algorithm Mini Framework (JMOEAMF)
 * Clase: Individual.java
 * Paquete: algorithm
 * Info: Representa a un individuo
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

package solution.real;

import java.util.LinkedList;
import utils.DoubleUtils;

/**
 *
 * @author J. Alfredo Brambila H. <alfredo.brambila@outlook.com>
 */
public class Individual {
    
    private double[] position; // vector de variables de desicion 
    private double[] cost; // valor de la funcion de costo f0 a fn
    private double g; // valor del costo de descomposicion de MOEA/D
    private int rank; // posicion de ranqueo
    private LinkedList<Integer> dominationSet; // conjunto de soluciones dominadas
    private int dominatedCount; // contador de veces que la solucion es dominada
    private double crowdingDistance; // distancia de crowding
    private int index; // indice
    private boolean isDominated; // la solucion es dominada (MOEA/D)
    private double normalizedCost[]; // valor del costo de descomposicion de MOEA/D normalizado
    private int associatedRef; //
    private double distanceToAssociatedRef; //

    public Individual() {
        dominationSet = new LinkedList<Integer>();
        rank = 0;
        dominatedCount = 0;
        crowdingDistance = 0;
        index = -1;
    }
    
    public void randomInit(double[] lower, double[] upper, int n) {
        this.position = new double[n];
        for(int i=0; i<n; i++) {
            this.position[i] = DoubleUtils.getRandomDoubleNumber(lower[i], upper[i]);
        }
    }
    
    public void copy(Individual a) {
        this.position = new double[a.getPosition().length];
        for(int i=0; i<a.getPosition().length; i++) {
            this.position[i] = a.getPosition()[i];
        }
        
        this.cost = new double[a.getCost().length];
        for(int i=0; i<a.getCost().length; i++) {
            this.cost[i] = a.getCost()[i];
        }  
    }
    
    public void clone(Individual a) {
        this.position = new double[a.getPosition().length];
        for(int i=0; i<a.getPosition().length; i++) {
            this.position[i] = a.getPosition()[i];
        }
        
        this.cost = new double[a.getCost().length];
        for(int i=0; i<a.getCost().length; i++) {
            this.cost[i] = a.getCost()[i];
        } 
        
        this.rank = a.rank;
        for(Integer ii : a.dominationSet) {
            this.dominationSet.addLast(ii);
        }
        dominatedCount = a.dominatedCount;
        crowdingDistance = a.crowdingDistance;
        index = a.index;
        g = a.g;
        isDominated = a.isDominated;
    }
    
    @Override
    public String toString() {
        String s = "";
        for(int i=0; i<this.position.length; i++) {
            s += this.position[i] + "\t";
        }
        s += "|\t";
        for(int i=0; i<this.cost.length; i++) {
            s += "F" + (i+1) + ":\t" + this.cost[i] + "\t";
        }
        
        return s;
    }
    
    public String printInfo() {
        String s = "";
        for(int i=0; i<this.position.length; i++) {
            s += this.position[i] + "\t";
        }
        s += "|\t";
        for(int i=0; i<this.cost.length; i++) {
            s += "F" + (i+1) + ":\t" + this.cost[i] + "\t";
        }
        
        s += "g:\t" + g + "\tisDominated:\t" + isDominated;
        
        return s;
    }
    
    public String getVariablesAndObjectives(String separator) {
        String s = "";
        for(int i=0; i<this.position.length; i++) {
            s += this.position[i] + separator;
        }
        s += "|"+separator;
        for(int i=0; i<this.cost.length; i++) {
            s += "F" + (i+1) + ":"+separator + this.cost[i] + separator;
        }
        
        return s;
    }
    
    public String getObjectives(String separator) {
        String s = "";
        for(int i=0; i<this.cost.length; i++) {
            if(i == this.cost.length-1) {
                s += "" + this.cost[i];
            } else {
                s += "" + this.cost[i] + separator;
            }  
        }
        
        return s;
    }
    
    public void addTodominationSet(int val) {
        if(dominationSet != null) {
            this.dominationSet.addLast(val);
        } else {
            System.out.println("No dominationSet created... ");
        }
    }
    
    public void addTodominationSetFromList(LinkedList<Integer> l) {
        for(Integer v : l) {
            this.dominationSet.addLast(v);
        }
    }
    
    public void dominatedCountIncrement() {
        this.dominatedCount++;
    }
    
    public void dominatedCountDecrement() {
        this.dominatedCount--;
    }
    
    
    public Individual(int nPos, int nCos) {
        this.position = new double[nPos];
        this.cost = new double[nCos];
        dominationSet = new LinkedList<Integer>();
        index = -1;
    }
    
    public void copyConfigAtributes(Individual a) {
        try {
            this.position = new double[a.getPosition().length];
            this.cost = new double[a.getCost().length];
        } catch(Exception e) {
            System.out.println(">>> Error: " + e);
        }
        
    }
    
    public void setCost(int i, double value) {
        this.cost[i] = value;
    }
    
    public void setPosition(int i, double value) {
        this.position[i] = value;
    }
    
    public boolean cmp(Individual a) {
        boolean equals = true;
        
        for(int i=0; i<this.position.length; i++) {
            if(this.position[i] != a.position[i]) {
                return false;
            }
        }
        
        return equals;
    }
    
    /**
     * @return the position
     */
    public double[] getPosition() {
        return position;
    }

    /**
     * @param position the position to set
     */
    public void setPosition(double[] position) {
        this.position = position;
    }

    /**
     * @return the cost
     */
    public double[] getCost() {
        return cost;
    }
    
    public double getCost_0() {
        return cost[0];
    }

    /**
     * @param cost the cost to set
     */
    public void setCost(double[] cost) {
        this.cost = cost;
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
     * @return the dominationSet
     */
    public LinkedList<Integer> getDominationSet() {
        return dominationSet;
    }

    /**
     * @param dominationSet the dominationSet to set
     */
    public void setDominationSet(LinkedList<Integer> dominationSet) {
        this.dominationSet = dominationSet;
    }

    /**
     * @return the dominatedCount
     */
    public int getDominatedCount() {
        return dominatedCount;
    }

    /**
     * @param dominatedCount the dominatedCount to set
     */
    public void setDominatedCount(int dominatedCount) {
        this.dominatedCount = dominatedCount;
    }

    /**
     * @return the crowdingDistance
     */
    public double getCrowdingDistance() {
        return crowdingDistance;
    }

    /**
     * @param crowdingDistance the crowdingDistance to set
     */
    public void setCrowdingDistance(double crowdingDistance) {
        this.crowdingDistance = crowdingDistance;
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
     * @return the g
     */
    public double getG() {
        return g;
    }

    /**
     * @param g the g to set
     */
    public void setG(double g) {
        this.g = g;
    }

    /**
     * @return the isDominated
     */
    public boolean isIsDominated() {
        return isDominated;
    }

    /**
     * @param isDominated the isDominated to set
     */
    public void setIsDominated(boolean isDominated) {
        this.isDominated = isDominated;
    }

    /**
     * @return the normalizedCost
     */
    public double[] getNormalizedCost() {
        return normalizedCost;
    }

    /**
     * @param normalizedCost the normalizedCost to set
     */
    public void setNormalizedCost(double[] normalizedCost) {
        this.normalizedCost = normalizedCost;
    }

    /**
     * @return the associatedRef
     */
    public int getAssociatedRef() {
        return associatedRef;
    }

    /**
     * @param associatedRef the associatedRef to set
     */
    public void setAssociatedRef(int associatedRef) {
        this.associatedRef = associatedRef;
    }

    /**
     * @return the distanceToAssociatedRef
     */
    public double getDistanceToAssociatedRef() {
        return distanceToAssociatedRef;
    }

    /**
     * @param distanceToAssociatedRef the distanceToAssociatedRef to set
     */
    public void setDistanceToAssociatedRef(double distanceToAssociatedRef) {
        this.distanceToAssociatedRef = distanceToAssociatedRef;
    }
    
}
