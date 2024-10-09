import java.math.BigDecimal;
import java.math.RoundingMode;
public class Simplex {

    public static void main(String[] args) {
        double[] C1 = {5, 4};  // 1 example
        double[][] A1 = {
                {6, 4}, 
                {1, 2},
                {-1, 1},
                {0, 1}
                
        };
        double[] b1 = {24, 6, 1, 2};
        double eps1 = 0.1;
        System.out.println("Example 1:");
        solving(C1, A1, b1, eps1);


        double[] C2 = { 9, 10, 16 }; // 2 example
        double[][] A2 = {
                { 18, 15, 12 },
                { 6, 4, 8 },
                { 5, 3, 3 },

        };
        double[] b2 = { 360, 192, 180 };
        double eps2 = 0.01;
        System.out.println();
        System.out.println("Example 2:");
        solving(C2, A2, b2, eps2);


        double[] C3 = { 2, 1 }; // 3 example
        double[][] A3 = {
                { 1, -1},
                { 2, 0},

        };
        double[] b3 = { 10, 40 };
        double eps3 = 0.01;
        System.out.println();
        System.out.println("Example 3:");
        solving(C3, A3, b3, eps3);


        double[] C4 = { 3, 5, 1, 2 }; //4 example
        double[][] A4 = {
                { 5, 3, 7, 2},
                {8, 5, 8, 8 },
                { 1, 8, 3, 9 },
                { 2, 6, 4, 1 }

        };
        double[] b4 = { 16, 40, 2, 4 };
        double eps4 = 0.001;
        System.out.println();
        System.out.println("Example 4:");
        solving(C4, A4, b4, eps4);

        double[] C5 = { 4, 5, 6 }; // 5 example
        double[][] A5 = {
                { 4, 5, 6},
                { 8, 6, 4},
                { 6, 4, 5 },

        };
        double[] b5 = { 2000, 1770, 1600 };
        double eps5 = 0.01;
        System.out.println();
        System.out.println("Example 5:");
        solving(C5, A5, b5, eps5);



    }

    public static void solving(double[] C, double[][] A, double[] b, double eps) {
        int m = A.length;
        int n = C.length; 
        System.out.print("Maximize z = ");
        for (int j = 0; j < n; j++) {
            if (j == n - 1) {
                System.out.print(Double.toString(C[j]) + "*x" + Integer.toString(j + 1));
            } else {
                System.out.print(Double.toString(C[j]) + "*x" + Integer.toString(j + 1) + " + ");
            }
        }
        System.out.println();
        System.out.println("Subject to:");
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                System.out.print(Double.toString(A[i][j]) + "*x" + Integer.toString(j + 1));
                if (j < n - 1) {
                    System.out.print(" + ");
                }
            }
            System.out.print(" <= ");
            System.out.println(b[i]);
        }

        int totalvar = n + m; 
        double[][] tableau = new double[m + 1][totalvar + 1];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                tableau[i][j] = A[i][j];
            }
            tableau[i][n + i] = 1; 
            tableau[i][totalvar] = b[i];
        }

        for (int j = 0; j < n; j++) {
            tableau[m][j] = -C[j];
        }
        
         for (double value : b) {
             if (value < 0) {
                 System.out.println("The method is not applicable: RHS has negative values");
                 return;
             }
         }


        while (true) {
            int entervar = -1;
            for (int j = 0; j < totalvar; j++) {
                if (tableau[m][j] < -eps) {
                    entervar = j;
                    break;
                }
            }

            if (entervar == -1) {
                break; 
            }

            int leavevar = -1;
            double minratio = Double.MAX_VALUE;

            for (int i = 0; i < m; i++) {
                if (tableau[i][entervar] > eps) {
                    double ratio = tableau[i][totalvar] / tableau[i][entervar];
                    if (ratio < minratio) {
                        minratio = ratio;
                        leavevar = i;
                    }
                }
            }

            if (leavevar == -1) {
                System.out.println("Unbounded");
                return;
            }

            pivot(tableau, entervar, leavevar);

        }

        System.out.println("Optimal solution found:");
        double[] x = new double[n]; 

        int flag = 0;
        for (int i = 0; i < n; i++) { 
            flag = 0;
            for (int j = 0; j < m; j++) { 

                if (tableau[j][i] == 1) {
                    flag = 1;
                    x[i] = tableau[j][totalvar];

                }
                if (tableau[j][i] != 0 && tableau[j][i] != 1) {
                    flag = 0;
                }

            }
            if (flag == 0) {
                x[i] = 0.0;
            }
        }

        int count = String.valueOf(eps).split("\\.")[1].length();

        BigDecimal bd = new BigDecimal(tableau[m][totalvar]);

        double z = bd.setScale(count, RoundingMode.HALF_UP).doubleValue();

        System.out.print("A vector of decision variables: ");
        for (int i = 0; i < x.length; i++) {
            BigDecimal bd1 = new BigDecimal(x[i]);
            double x0 = bd1.setScale(count, RoundingMode.HALF_UP).doubleValue();
            System.out.print(x0 + "  ");
        }
        System.out.println();
        System.out.println("Maximum value z: " + z);

    }

    private static void pivot(double[][] tableau, int entervar, int leavevar) {
        int m = tableau.length - 1;
        int n = tableau[0].length - 1;
        double pivotValue = tableau[leavevar][entervar];
        for (int j = 0; j <= n; j++) {
            tableau[leavevar][j] /= pivotValue;
        }

        for (int i = 0; i <= m; i++) {
            if (i != leavevar) {
                double factor = tableau[i][entervar];
                for (int j = 0; j <= n; j++) {
                    tableau[i][j] -= factor * tableau[leavevar][j];
                }
            }
        }

    }
}
