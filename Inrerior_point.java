import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.Arrays;

public class Inrerior_point {
    public static void main(String[] args) {

        double alpha1 = 0.5;
        double alpha2 = 0.9;

        //Example 1
        double[] c1 = {5, 4, 0, 0, 0, 0}; // A vector of coefficients of objective function - C
                                            //in equation form

        double[][] A1 = {{6, 4, 1, 0, 0, 0}, // A matrix of coefficients of constraint function - A.
                        {1, 2, 0, 1, 0, 0},
                        {-1, 1, 0, 0, 1, 0},
                        {0, 1, 0, 0, 0, 1},
                    };

        double[] x01 ={1, 1, 14, 3, 1, 2}; //initial starting point

        double[] b1 = {24, 6, 1, 2}; // A vector of right-hand side numbers - b.
        double eps1 = 0.01; // The approximation accuracy eps

       


        System.out.println("Example 1:");
        printproblem(A1, c1, b1);
        System.out.println("When alpha = 0.5:");
        interiorPointMethod(A1, c1, x01, alpha1, eps1);
        System.out.println();
        System.out.println("When alpha = 0.9:");
        interiorPointMethod(A1, c1, x01, alpha2, eps1);
        System.out.println();
        System.out.println("By Simplex method: ");
        System.out.println("A vector of decision variables x*: ");
        System.out.println("x1 = 3.0  x2 = 1.5");
        System.out.println("Maximum value: 21.0");
        System.out.println();
        
        //Example 2
        double[] c2 = { 9, 10, 16, 0, 0, 0};

        double[][] A2 = { { 18, 15, 12, 1, 0, 0 },
                        {6, 4, 8, 0, 1, 0 },
                        { 5, 3, 3, 0, 0, 1 },

        };

        double[] x02 = { 1, 1, 1, 315, 174, 169 };

        double[] b2 = { 360, 192, 180 };
        double eps2 = 0.01;

        System.out.println("Example 2:");
        printproblem(A2, c2, b2);
        System.out.println("When alpha = 0.5:");
        interiorPointMethod(A2, c2, x02, alpha1, eps2);
        System.out.println();
        System.out.println("When alpha = 0.9:");
        interiorPointMethod(A2, c2, x02, alpha2, eps2);
        System.out.println();
        System.out.println("By Simplex method: ");
        System.out.println("A vector of decision variables x*: ");
        System.out.println("x1 = 0.0  x2 = 8.0  x3 = 20.0");
        System.out.println("Maximum value: 400.0");
        System.out.println();

        //Example 3
        double[] c3 = { 2, 1, 0, 0 };

        double[][] A3 = { { 1, -1, 1, 0 },
                        { 2, 0, 0, 1 },
        };

        double[] x03 = { 1, 1, 10, 38 };

        double[] b3 = { 10, 40};
        double eps3 = 0.01;
        
        System.out.println("Example 3:");
        printproblem(A3, c3, b3);
        System.out.println("When alpha = 0.5:");
        interiorPointMethod(A3, c3, x03, alpha1, eps3);
        System.out.println();
        System.out.println("When alpha = 0.9:");
        interiorPointMethod(A3, c3, x03, alpha2, eps3);
        System.out.println();
        System.out.println("By Simplex method: ");
        System.out.println("Unbounded");
        System.out.println();

        //Example 4
        double[] c4 = { 4, 5, 6, 0, 0, 0};

        double[][] A4 = { { 4, 5, 6, 1, 0, 0 },
                        { 8, 6, 4, 0, 1, 0},
                        {6, 4, 5, 0, 0, 1},
        };

        double[] x04 = { 1, 1, 1, 1985, 1752, 1585};

        double[] b4 = { 2000, 1770, 1600};
        double eps4 = 0.001;

        System.out.println("Example 4:");
        printproblem(A4, c4, b4);
        System.out.println("When alpha = 0.5:");
        interiorPointMethod(A4, c4, x04, alpha1, eps4);
        System.out.println();
        System.out.println("When alpha = 0.9:");
        interiorPointMethod(A4, c4, x04, alpha2, eps4);
        System.out.println();
        System.out.println("By Simplex method: ");
        System.out.println("A vector of decision variables x*: ");
        System.out.println("x1 = 0.0  x2 = 175.0  x3 = 180.0");
        System.out.println("Maximum value: 1955.0");
        System.out.println();

        // Example 5
        double[] c5 = { 1, 1, 0, 0};

        double[][] A5 = { { 2, 4, 1, 0 },
                        { 1, 3, 0, 1 },
        };

        double[] x05 = { 1, 1, 10, 5 };

        double[] b5 = { 16, 9 };
        double eps5 = 0.001;

        System.out.println("Example 5:");
        printproblem(A5, c5, b5);
        System.out.println("When alpha = 0.5:");
        interiorPointMethod(A5, c5, x05, alpha1, eps5);
        System.out.println();
        System.out.println("When alpha = 0.9:");
        interiorPointMethod(A5, c5, x05, alpha2, eps5);
        System.out.println();
        System.out.println("By Simplex method: ");
        System.out.println("A vector of decision variables x*: ");
        System.out.println("x1 = 8.0  x2 = 0.0");
        System.out.println("Maximum value: 8.0");


}

    public static void interiorPointMethod(double[][] A, double[] c, double[] x0, double alpha, double eps) {
         for (double value:x0) {
            if (value < 0) {
                System.out.println("Initial solution is not correct");
                return;
            }
         }
        int iteration = 0;
        while (true) {
            double[] w = x0;

            // calculate D
            double[][] D = new double[x0.length][x0.length];
            for (int i = 0; i < x0.length; i++) {
                D[i][i] = x0[i];
            }

            // calculate A* = AD
            double[][] A1 = multiply(A, D);
            // calculate c* = Dc
            double[] c1 = multiplywithvector(D, c);

            // create Identity
            double[][] I = createIdentityMatrix(x0.length);

            double[][] res = multiply(multiply(transpose(A1), inverse(multiply(A1, transpose(A1)))), A1);

            // calculate P
            double[][] P = subtract(I, res);

            // calculate Cp
            double[] Cp = multiplywithvector(P, c1);

            // vector with 1
            double[] ones = new double[x0.length];
            Arrays.fill(ones, 1.0);

            double v = findLargestNegativeComponent(Cp);
            if (v == 0.0) {
                System.out.println("Unbounded");
                break;
            }
            // calculate x1
            double del = Math.abs(alpha / v);

            double[] x1 = addVectors(ones, multiplyVectorByScalar(Cp, del));

            double[] x2 = multiplywithvector(D, x1);

            
            x0 = x2;
            iteration = iteration + 1;

            if (checkConvergence(x2, w, eps)) {
                printsolution(x2,eps,c);
                break;
            }

        }
}

    public static void printproblem(double[][] A, double[] c, double[]b) {
        System.out.print("F = ");
        int n = c.length;
        for (int j = 0; j < n; j++) {
            if (j == n - 1) {
                System.out.print(Double.toString(c[j]) + "*x" + Integer.toString(j + 1));
            } else {
                System.out.print(Double.toString(c[j]) + "*x" + Integer.toString(j + 1) + " + ");
            }
        }

        System.out.println();
        int m = A.length;
        System.out.println();
        System.out.println("Subject to:");
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                System.out.print(Double.toString(A[i][j]) + "*x" + Integer.toString(j + 1));
                if (j < n - 1) {
                    System.out.print(" + ");
                }
            }
            System.out.print(" = ");
            System.out.println(b[i]);
            
        }
        System.out.println();
    }


    public static void printsolution(double[] x, double eps, double[] c) {
        int count = String.valueOf(eps).split("\\.")[1].length();
        System.out.println("A vector of decision variables x*:");
        for (int i = 0; i < x.length; i++) {
            BigDecimal bd1 = new BigDecimal(x[i]);
            double value1 = bd1.setScale(count, RoundingMode.HALF_UP).doubleValue();
            System.out.print("x"+ (i + 1) + " = " + value1 +"  ");
        }
        System.out.println();
        System.out.print("Maximum value: ");
        double sum = 0;
        for (int j = 0; j < x.length; j++) {
            sum += x[j] * c[j];
        }
        BigDecimal bd2 = new BigDecimal(sum);
        double sum1 = bd2.setScale(count, RoundingMode.HALF_UP).doubleValue();
        System.out.println(sum1);

    }

    public static double [][] multiply(double[][] x, double[][] y) {
        int rowsx = x.length;
        int colsx = x[0].length;
        int rowsy = y.length;
        int colsy = y[0].length;
        double[][] result = new double[rowsx][colsy];

        for (int i = 0; i < rowsx; i++) {
            for (int j = 0; j < colsy; j++) {
                result[i][j] = 0;
                for (int k = 0; k < colsx; k++) {
                result[i][j] += x[i][k] * y[k][j]; 
                }
            }
        }
        return result;
    }

    public static double[] multiplywithvector(double[][] x, double[] y) {
        int rowsx = x.length;
        int colsx = x[0].length;

        double[] result = new double[rowsx];

        for (int i = 0; i < rowsx; i++) {
                result[i] = 0;
            for (int j = 0; j < colsx; j++) {
                result[i] += x[i][j] * y[j];
            }
            
        }
        return result;
    }

    public static double[][] createIdentityMatrix(int size) {
        double[][] matrix = new double[size][size];
        for (int i = 0; i < size; i++) {
            matrix[i][i] = 1.0;
        }
        return matrix;
    }

    public static double[][] transpose(double[][] A) {
        int rows = A.length;
        int cols = A[0].length;
        double[][] transposed = new double[cols][rows];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                transposed[j][i] = A[i][j];
            }
        }
        return transposed;
    }

    public static double[][] inverse (double[][] A) {
        int n = A.length;
        double[][] augmentedMatrix = new double[n][2 * n];

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                augmentedMatrix[i][j] = A[i][j];
            }
            augmentedMatrix[i][i+n] = 1;
        }

        for (int i = 0; i < n; i++) {
            double diagValue = augmentedMatrix[i][i];
            if (diagValue == 0) {
                throw new IllegalArgumentException("Matrix is singular and cannot be inverted.");
            }

            for (int j = 0; j < 2*n; j++) {
                augmentedMatrix[i][j] /= diagValue;
            }


            for (int j = 0; j < n; j++) {
                if (j != i) {
                    double factor = augmentedMatrix[j][i];
                    for (int k = 0; k < 2 *n; k++) {
                        augmentedMatrix[j][k] -= factor * augmentedMatrix[i][k];
                    }
                }
            }
        }
        double[][] inverse = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                inverse[i][j] = augmentedMatrix[i][j + n];
            }
        }
        return inverse;
    }

    public static double[][] subtract(double[][] A, double[][] B) {
        int rowsA = A.length;
        int colsA = A[0].length;
        int rowsB = B.length;
        int colsB = B[0].length;

       

        // Initialize the result matrix
        double[][] result = new double[rowsA][colsA];

        // Perform the subtraction
        for (int i = 0; i < rowsA; i++) {
            for (int j = 0; j < colsA; j++) {
                result[i][j] = A[i][j] - B[i][j];
            }
        }

        return result;
    }


    public static boolean checkConvergence(double[] yy, double[] v, double eps) {
        // Calculate the infinity norm of the difference between yy and v
        double norm = 0.0;
        for (int i = 0; i < yy.length; i++) {
            norm = Math.max(norm, Math.abs(yy[i] - v[i]));
        }
        // Compare the norm with the threshold
        return norm < eps;
    }


    public static Double findLargestNegativeComponent(double[] vector) {
        Double largestNegative = 0.0;
        double largestAbsoluteValue = Double.NEGATIVE_INFINITY;

        for (double value : vector) {
            if (value < 0) {
                if (Math.abs(value) > largestAbsoluteValue) {
                    largestAbsoluteValue = Math.abs(value);
                    largestNegative = value;
                }
            }
        }

        return largestNegative;
    }

    public static double[] multiplyVectorByScalar(double[] vector, double scalar) {
        // Initialize the result vector
        double[] result = new double[vector.length];

        // Perform the multiplication
        for (int i = 0; i < vector.length; i++) {
            result[i] = vector[i] * scalar;
        }

        return result;
    }

    public static double[] addVectors(double[] A, double[] B) {

        double[] result = new double[A.length];

        for (int i = 0; i < A.length; i++) {
            result[i] = A[i] + B[i];
        }

        return result;
    }


}
