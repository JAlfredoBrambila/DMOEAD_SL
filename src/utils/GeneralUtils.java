

package utils;

import java.io.File;
import java.nio.file.Paths;
import java.util.LinkedList;
import problem.multiobjective.Problem;
import problem.multiobjective.dynamical.DynamicProblem;
import solution.real.Individual;

/**
 *
 * @author J. Alfredo Brambila H. <alfredo.brambila@outlook.com>
 */
public class GeneralUtils {

    public static void SaveExperimentationFile(String carpetaExperimentos, Problem problem, String nombreAlgoritmo, int numAlg, LinkedList<Individual> EP, int expNum) {
        
            //String carpetaExperimentos="ExpJun23";
            
            String ruta=Paths.get("").toAbsolutePath().toString() + "\\" + carpetaExperimentos + "\\" + problem.getProblemName()+"\\";
            String rutaFrente=Paths.get("").toAbsolutePath().toString() + "\\" + carpetaExperimentos + "\\" + problem.getProblemName()+"\\";
            String rutaFrenteALL=Paths.get("").toAbsolutePath().toString() + "\\" + carpetaExperimentos + "\\ALL\\" + problem.getProblemName()+"\\";
            
            DoubleUtils.createDirectory(ruta);
            DoubleUtils.createDirectory(rutaFrente);
            DoubleUtils.createDirectory(rutaFrenteALL);
            
            DoubleUtils.saveFrontToFile(EP,rutaFrenteALL + nombreAlgoritmo + expNum + "_.txt");
            DoubleUtils.saveFrontObjectivesToFile(EP,rutaFrente + nombreAlgoritmo + expNum + "_.txt");
            
            String rutaMainTxt = Paths.get("").toAbsolutePath().toString() + "\\" + carpetaExperimentos + "\\" + problem.getProblemName() + "\\" + problem.getProblemName() + "Front.txt";
            String rutaMainTxtF = Paths.get("").toAbsolutePath().toString() + "\\" + carpetaExperimentos + "\\" + problem.getProblemName() + ".txt";
            if (expNum == 1) {
                File confFile = new File(rutaMainTxtF);
                if (!confFile.exists()) {
                    DoubleUtils.saveToMainTXTFile(rutaMainTxtF, numAlg + " 30 " + problem.getNumObjectives());
                    DoubleUtils.saveToMainTXTFile(rutaMainTxtF, rutaMainTxt);
                }
            }
            DoubleUtils.saveToMainTXTFile(rutaMainTxtF, rutaFrente + nombreAlgoritmo + expNum + "_.txt");
    }
    
    public static void SaveExperimentationFile(String carpetaExperimentos, DynamicProblem problem, int contW, String nombreAlgoritmo, int numAlg, LinkedList<Individual> EP, int expNum) {
        
            //String carpetaExperimentos="ExpJun23";
            
            String ruta=Paths.get("").toAbsolutePath().toString() + "\\" + carpetaExperimentos + "\\" + problem.getProblemName()+"\\";
            String rutaFrente=Paths.get("").toAbsolutePath().toString() + "\\" + carpetaExperimentos + "\\" + problem.getProblemName() + "\\F"+contW+"\\";
            String rutaFrenteALL=Paths.get("").toAbsolutePath().toString() + "\\" + carpetaExperimentos + "\\ALL\\" + problem.getProblemName() + "\\F"+contW+"\\";
            
            DoubleUtils.createDirectory(ruta);
            DoubleUtils.createDirectory(rutaFrente);
            DoubleUtils.createDirectory(rutaFrenteALL);
            
            DoubleUtils.saveFrontToFile(EP,rutaFrenteALL + nombreAlgoritmo + expNum + "_.txt");
            DoubleUtils.saveFrontObjectivesToFile(EP,rutaFrente + nombreAlgoritmo + expNum + "_.txt");
            
            //FrentesReales\FDA3
            //String rutaMainTxt = Paths.get("").toAbsolutePath().toString() + "\\" + carpetaExperimentos + "\\" + problem.getProblemName() + "\\" + problem.getProblemName() + "Front.txt";
            String rutaMainTxt = Paths.get("").toAbsolutePath().toString() + "\\" + carpetaExperimentos + "\\FrentesReales\\" + problem.getProblemName() + "\\" + "F"+contW+".txt";
            String rutaMainTxtF = Paths.get("").toAbsolutePath().toString() + "\\" + carpetaExperimentos + "\\" + problem.getProblemName() + "\\F" + contW + ".txt";
            if (expNum == 1) {
                File confFile = new File(rutaMainTxtF);
                if (!confFile.exists()) {
                    DoubleUtils.saveToMainTXTFile(rutaMainTxtF, numAlg + " 30 " + problem.getNumObjectives());
                    DoubleUtils.saveToMainTXTFile(rutaMainTxtF, rutaMainTxt);
                }
            }
            DoubleUtils.saveToMainTXTFile(rutaMainTxtF, rutaFrente + nombreAlgoritmo + expNum + "_.txt");
    }
}
