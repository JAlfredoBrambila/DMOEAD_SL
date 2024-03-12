/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */

package utils;

/**
 *
 * @author J. Alfredo Brambila H. <alfredo.brambila@outlook.com>
 */


import java.util.ArrayList;
import java.util.List;

public class ReferencePointsHyperPlane {
    public static double[][] generateReferencePoints(int M, int p) {
        double[][] Zr = getFixedRowSumIntegerMatrix(M, p);
        for (int i = 0; i < Zr.length; i++) {
            for (int j = 0; j < Zr[i].length; j++) {
                Zr[i][j] = Zr[i][j] / p;
            }
        }
        return transposeMatrix(Zr);
    }

    private static double[][] getFixedRowSumIntegerMatrix(int M, int RowSum) {
        if (M < 1) {
            throw new IllegalArgumentException("M cannot be less than 1.");
        }

        if (M != Math.floor(M)) {
            throw new IllegalArgumentException("M must be an integer.");
        }

        if (M == 1) {
            double[][] A = new double[1][1];
            A[0][0] = RowSum;
            return A;
        }

        List<double[]> A = new ArrayList<>();
        for (int i = 0; i <= RowSum; i++) {
            double[][] B = getFixedRowSumIntegerMatrix(M - 1, RowSum - i);
            for (double[] row : B) {
                double[] newRow = new double[row.length + 1];
                newRow[0] = i;
                System.arraycopy(row, 0, newRow, 1, row.length);
                A.add(newRow);
            }
        }

        return convertListToArray(A);
    }

    private static double[][] transposeMatrix(double[][] matrix) {
        int rows = matrix.length;
        int cols = matrix[0].length;
        double[][] transpose = new double[cols][rows];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                transpose[j][i] = matrix[i][j];
            }
        }
        return transpose;
    }

    private static double[][] convertListToArray(List<double[]> list) {
        int rows = list.size();
        int cols = list.get(0).length;
        double[][] array = new double[rows][cols];
        for (int i = 0; i < rows; i++) {
            array[i] = list.get(i);
        }
        return array;
    }
}
