/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */

package utils.adaptative;

/**
 *
 * @author J. Alfredo Brambila H. <alfredo.brambila@outlook.com>
 */
public class AdaptativeModule1 {

    //static int uno, dos, tres, cuatro;
    public static int FRRMAB(double[] quality, double[] rewards, int[] strategy_usage, int numStrategies, double scale) {
        int best_index = 0;

        double temp1, temp2, temp3;
        int total_usage;
        //int best_index;

        total_usage = 0;
        for (int i = 0; i < numStrategies; i++) {
            total_usage += strategy_usage[i];
        }

        for (int i = 0; i < numStrategies; i++) {
            temp1 = 2 * Math.log(total_usage);
            temp2 = temp1 / strategy_usage[i];
            temp3 = Math.sqrt(temp2);
            quality[i] = rewards[i] + scale * temp3;
        }

        best_index = 0;
        for (int i = 1; i < numStrategies; i++) {
            //System.out.println("Quality i = " + quality[i] + " > Quality Best: " + quality[best_index] + " : " + (quality[i] > quality[best_index]));
            if (quality[i] > quality[best_index]) {
                best_index = i;
            }
        }
       // System.out.println("B: " + best_index);
        //System.out.println("");
        best_index++;

        return best_index;
    }
    
    public static void updateRewards(double[][] slidingWindow, double[] reward, int[] strategyUsage, int curSize) {
        int index;
        double fitnessImprovement;
        
        for(int i=0; i<curSize; i++) {
            if(slidingWindow[0][i] == -1 || slidingWindow[1][i] == -1) {
                System.out.println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
            }
            
            index = (int) slidingWindow[0][i];
            fitnessImprovement = slidingWindow[1][i];
            
            //System.out.println("Indexes: " + index);
            //System.out.println("Mejora: " + fitnessImprovement);
            //System.out.println("");
            
            switch(index) {
                case 1:
                    reward[0] += fitnessImprovement;
                    strategyUsage[0]++;
                    break;
                case 2:
                    reward[1] += fitnessImprovement;
                    strategyUsage[1]++;
                    break;
                case 3:
                    reward[2] += fitnessImprovement;
                    strategyUsage[2]++;
                    break;
                case 4:
                    reward[3] += fitnessImprovement;
                    strategyUsage[3]++;
                    break;
            }
            
            
        }
        
        /*System.out.println("Reward: ");
        System.out.println("1: " + reward[0] + " Uso: " + strategyUsage[0]);
        System.out.println("2: " + reward[1] + " Uso: " + strategyUsage[1]);
        System.out.println("3: " + reward[2] + " Uso: " + strategyUsage[2]);
        System.out.println("4: " + reward[3] + " Uso: " + strategyUsage[3]);
        System.out.println("");*/
    }
    
    public static void rankRewards(double[] reward, int numStrategies, int[] rank) {
        double[][] temp;
        double temp_index;
        double temp_value;

        temp = new double[2][numStrategies];
        for (int i = 0; i < numStrategies; i++) {
            temp[0][i] = reward[i];
            temp[1][i] = i;
        }

        for (int i = 0; i < numStrategies - 1; i++) {
            for (int j = i + 1; j < numStrategies; j++) {
                //System.out.println("Reg: " + temp[0][i] +" < " + temp[0][j]);
                if (temp[0][i] < temp[0][j]) {
                    
                    temp_value = temp[0][j];
                    temp[0][j] = temp[0][i];
                    temp[0][i] = temp_value;

                    temp_index = temp[1][j];
                    temp[1][j] = temp[1][i];
                    temp[1][i] = temp_index;
                    
                }
            }
        }

        for (int i = 0; i < numStrategies; i++) {
            rank[i] = (int) temp[1][i];
        }
    }
    
    public static void CreaditAssignmentDecay(double[] strategy_rewards, double[] decay_rewards, int numStrategies, int[] rank, double decayFactor) {
       
        double decayed, decay_sum;
        double[] decay_value;

        decay_value = new double[numStrategies];

        for (int i = 0; i < numStrategies; i++) {
            decayed = Math.pow(decayFactor, i);
            switch (rank[i]) {
                case 0:
                    decay_value[0] = strategy_rewards[0] * decayed;
                    break;
                case 1:
                    decay_value[1] = strategy_rewards[1] * decayed;
                    break;
                case 2:
                    decay_value[2] = strategy_rewards[2] * decayed;
                    break;
                case 3:
                    decay_value[3] = strategy_rewards[3] * decayed;
                    break;
            }
        }

        decay_sum = 0.0;
        for (int i = 0; i < numStrategies; i++) {
            decay_sum += decay_value[i];
        }

        for (int i = 0; i < numStrategies; i++) {
            if (decay_sum == 0) {
                decay_rewards[i] = 0.0;
            } else {
                decay_rewards[i] = decay_value[i] / decay_sum;
            }
        }
    }
    
}
