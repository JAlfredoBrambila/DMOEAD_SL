/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */

package utils.adaptative;

import java.util.LinkedList;

/**
 *
 * @author J. Alfredo Brambila H. <alfredo.brambila@outlook.com>
 */
public class SlidingWindow {

    LinkedList<Window> windowList;
    
    int capacity = 0;
    
    public SlidingWindow(int maxSize) {
        windowList = new LinkedList<Window>();
        this.capacity = maxSize;
    }
    
    public boolean add(Window element) {
        if(windowList.size() == this.capacity) {
            System.out.println(">>>>>> SWindow full... ");
            return false;
        }
        
        this.windowList.add(element);
        //this.windowList.addLast(element);
        return true;
    }
    
    public Window get(int i) {
        if(i>= this.windowList.size() || i < 0) {
            throw new IndexOutOfBoundsException("Index out of Bound " + i);
        }
        return this.windowList.get(i);
    }
    
    public int getMaxSize() {
        return this.capacity;
    }
    
    public int size() {
        return this.windowList.size();
    }
    
    public void clear() {
        this.windowList.clear();
    }
    
    public void removeElement(int i) {
        this.windowList.remove(i);
    }
    
    public double[][] slidingWindowToMtx() {
        int sizeL = this.windowList.size();
        double[][] mtx = new double[2][sizeL];
        for(int i=0; i<sizeL; i++) {
            mtx[0][i] = this.windowList.get(i).getIndex();
            mtx[1][i] = this.windowList.get(i).getFitnessImprovement();
        }
        
        return mtx;
    }
}
