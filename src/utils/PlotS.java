/**
 * Java multi-Objective Evolutionary Algorithm Mini Framework (JMOEAMF)
 * Clase: PlotS.java
 * Paquete: algorithm
 * Info: Permite la generacion de gráficas
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

package utils;

import java.awt.Color;
import java.util.LinkedList;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import solution.real.Individual;

/**
 *
 * @author J. Alfredo Brambila H. <alfredo.brambila@outlook.com>
 */
public class PlotS {

    XYSeriesCollection dataset = new XYSeriesCollection();

    
    
    public void addSerie(LinkedList<Individual> pop, int F, int c) {
        XYSeries series1 = new XYSeries("IT_"+c);
        for (Individual ind : pop) {
            if (ind.getRank() == F) {
                series1.add(ind.getCost()[0], ind.getCost()[1]);
            }
        }
        dataset.addSeries(series1);
    }
    
    public void addSerieM(LinkedList<Individual> pop, int c) {
        XYSeries series1 = new XYSeries("IT_"+c);
        for (Individual ind : pop) {
            //if(ind.getCost()[1]<1.3)
            series1.add(ind.getCost()[0], ind.getCost()[1]);
        }
        dataset.addSeries(series1);
    }

    public void plotSeries() {

        JFreeChart scatterPlot = ChartFactory.createScatterPlot(
                "Non-dominated Solutions (F1)", // Chart title
                "1st Objective", // X-Axis Label
                "2nd Objective", // Y-Axis Label
                dataset // Dataset for the Chart
        );

        //Changes background color
        XYPlot plot = (XYPlot) scatterPlot.getPlot();
        plot.setBackgroundPaint(new Color(255, 228, 196));

        // Create Panel
        ChartPanel panel = new ChartPanel(scatterPlot);
        panel.setSize(800, 600);
        //setContentPane(panel);
        JFramePlot frame = new JFramePlot();
        frame.setJPanel(panel);
        frame.setVisible(true);
    }
}
