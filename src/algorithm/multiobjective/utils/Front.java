/**
 * Java multi-Objective Evolutionary Algorithm Mini Framework (JMOEAMF)
 * Clase: Front.java
 * Paquete: algorithm.moead
 * Info: Mapea un frente para NSGAII/DNSGAII
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
package algorithm.multiobjective.utils;

import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;



/**
 *
 * @author J. Alfredo Brambila H. <alfredo.brambila@outlook.com>
 */
public class Front {
    private int nFront;
    private LinkedList<Integer> fronts;
    public Front(int nFront) {
        this.nFront = nFront;
        fronts = new LinkedList<Integer>();
    }
    
    public void addToFront(int sol) {
        if(fronts != null) {
            this.fronts.add(sol);
        } else {
            System.out.println("No se ha creado Frente");
        }
    }
    
    public void addFromList(LinkedList<Integer> l) {
        for(Integer v : l) {
            this.fronts.addLast(v);
        }
    }
    
    public void sortList() {
        Collections.sort(fronts);
    }

    /**
     * @return the fronts
     */
    public LinkedList<Integer> getFronts() {
        return fronts;
    }

    /**
     * @param fronts the fronts to set
     */
    public void setFronts(LinkedList<Integer> fronts) {
        this.fronts = fronts;
    }

    /**
     * @return the nFront
     */
    public int getnFront() {
        return nFront;
    }

    /**
     * @param nFront the nFront to set
     */
    public void setnFront(int nFront) {
        this.nFront = nFront;
    }
    
    
}
