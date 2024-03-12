/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */

package algorithm.multiobjective.dynamical.changeResponse;

import java.util.LinkedList;

/**
 *
 * @author J. Alfredo Brambila H. <alfredo.brambila@outlook.com>
 */
public interface ChangeResponseMechanism<T> {
    public void execute(LinkedList<T> parents);
}
