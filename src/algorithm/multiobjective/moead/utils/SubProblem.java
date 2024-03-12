/**
 * Java multi-Objective Evolutionary Algorithm Mini Framework (JMOEAMF)
 * Clase: SubProblem.java
 * Paquete: algorithm.moead
 * Info: Mapea un subproblema para MOEA/D
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
package algorithm.multiobjective.moead.utils;

/**
 *
 * @author J. Alfredo Brambila H. <alfredo.brambila@outlook.com>
 */
public class SubProblem {

    private double[] lambda;
    private double[] neighbors;

    @Override
    public String toString() {
        String s = "";
        
        s += "L:\t";
        for(int i=0; i<lambda.length; i++) {
            s += lambda[i] + "\t";
        }
        
        s += "Nhbs:\t";
        for(int i=0; i<neighbors.length; i++) {
            s += neighbors[i] + "\t";
        }
        
        return s;
    }
    
    /**
     * @return the lambda
     */
    public double[] getLambda() {
        return lambda;
    }

    /**
     * @param lambda the lambda to set
     */
    public void setLambda(double[] lambda) {
        this.lambda = lambda;
    }

    /**
     * @return the neighbors
     */
    public double[] getNeighbors() {
        return neighbors;
    }

    /**
     * @param neighbors the neighbors to set
     */
    public void setNeighbors(double[] neighbors) {
        this.neighbors = neighbors;
    }
    
}
