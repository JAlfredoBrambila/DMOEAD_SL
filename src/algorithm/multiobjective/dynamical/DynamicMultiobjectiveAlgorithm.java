/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Interface.java to edit this template
 */
package algorithm.multiobjective.dynamical;

import algorithm.multiobjective.MultiobjectiveAlgorithm;

/**
 *
 * @author al_x3
 */
public interface DynamicMultiobjectiveAlgorithm<T> extends MultiobjectiveAlgorithm<T> {
    @Override
    public void execute();
    public void setExperimentationFolder(String folderName);
}
