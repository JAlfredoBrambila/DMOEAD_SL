/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */

package utils;

import solution.real.Individual;

/**
 *
 * @author J. Alfredo Brambila H. <alfredo.brambila@outlook.com>
 */
public class Dominance {

    /*
        In multi-objective optimization, dominance is a relationship that 
        allows you to compare solutions based on the Pareto relationship. 
        A solution A is considered dominant over another solution B if A is 
        equal to or better than B in all the objectives, and better in at 
        least one of them.
    */
    
    /**
     * 
     * @param a
     * @param b
     * @return 
     */
    public static boolean dominates2(Individual a, Individual b) {
        boolean dominate1 = true;
        boolean dominate2 = false;
        int n = a.getCost().length;
        
        for (int i = 0; i < n; i++) {
            if(a.getCost()[i] <= b.getCost()[i]) {
                dominate1 = dominate1 && true;
            } else {
                dominate1 = dominate1 && false;
            }
        }
        
        for (int i = 0; i < n; i++) {
            if(a.getCost()[i] < b.getCost()[i]) {
                dominate2 = dominate2 || true;
            } else {
                dominate2 = dominate2 || false;
            }
        }
        
        
        return dominate1 && dominate2;
    }
}
