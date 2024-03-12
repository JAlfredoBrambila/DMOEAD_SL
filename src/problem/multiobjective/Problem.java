/**
 * Java multi-Objective Evolutionary Algorithm Mini Framework (JMOEAMF)
 * Clase: Problem.java
 * Paquete: problem
 * Info: Permite la generacion de problemas estaticos
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

package problem.multiobjective;

import solution.real.Individual;

/**
 *
 * @author J. Alfredo Brambila H. <alfredo.brambila@outlook.com>
 */
public interface Problem {
    public void costFunction(Individual ind);
    public int getNumVars();
    public int getNumObjectives();
    public String getProblemName();
    public double getVarMin();
    public double getVarMax();
    public double[] getLowerBound();
    public double[] getUpperBound();
    public int getContEvals();
}
