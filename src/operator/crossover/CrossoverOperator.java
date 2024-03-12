/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Interface.java to edit this template
 */
package operator.crossover;

import java.util.LinkedList;
import solution.real.Individual;

/**
 *
 * @author al_x3
 */
public interface CrossoverOperator<T> {
    public LinkedList<T> execute(LinkedList<T> parents);
    public T executeOnlyChild(LinkedList<T> parents);
    public void setCurrentIndividual(T c);
}
