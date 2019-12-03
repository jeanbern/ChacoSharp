#pragma warning disable HAA0101 // Array allocation for params parameter
#pragma warning disable HAA0601 // Value type to reference type conversion causing boxing allocation
using System;
using static ChacoSharp.StaticConstants;
using static ChacoSharp.Optimize.Determinant;
using static ChacoSharp.Optimize.Func3d;
using static ChacoSharp.Utilities.Randomize;

namespace ChacoSharp.Optimize
{
    public static unsafe class Opt3d
    {
        private const double MaximumAllowedStep = Math.PI / 4; /* maximum allowed step */
        private const double EarlyStepMinimum = 2.0e-4; /* min step for early convergence stages */
        private const double FinalStepMinimum = EarlyStepMinimum / 10.0d; /* min step for final convergence */
        private const double MinimumAcceptableGradiant = 1.0e-7; /* acceptable gradient for convergence */
        private const double HessianFactor = 2.0d; /* scales minimum tolerated hessian */
        private const double HessTol = 1.0e-6; /* smallest possible positive hess_min */
        private const double Pdtol = 1.0e-7; /* eval < tol considered to be 0 */
        private const double ConstraintGrowthScalingFactor = 20.0; /* scaling for constraint growth */

        /// <summary>
        /// Compute rotation angle to minimize distance to discrete points.
        /// </summary>
        /// <param name="graph">data structure containing vertex weights</param>
        /// <param name="yvecs">eigenvectors</param>
        /// <param name="vertexCount">total number of vertices</param>
        /// <param name="nmyvtxs">number of vertices I own</param>
        /// <param name="vwsqrt">square root of vertex weights</param>
        /// <param name="ptheta">optimal angles</param>
        /// <param name="pphi">optimal angles</param>
        /// <param name="pgamma">optimal angles</param>
        /// <param name="useVertexWeights">are vertex weights being used?</param>
        public static void opt3d(vtx_data** graph, double*[] yvecs, in int vertexCount, in int nmyvtxs, double* vwsqrt, double* ptheta, double* pphi, double* pgamma, bool useVertexWeights)
        {
            double[] coeffs = new double[25]; /* various products of yvecs */
            double[] angularVariables = new double[3]; /* angular variables */
            double[] currentBestMinimizer = new double[3] {0.0d, 0.0d, 0.0d}; /* best minimizer found so far */
            double[] functionGradient = new double[3]; /* gradiant of the function */
            double[] constraintGradient = new double[3]; /* gradiant of the constraint */
            double[][] functionHessian = new double[3][]{new double[3], new double[3], new double[3]}; /* hessian of the function */
            double[][] constraintHessian = new double[3][]{new double[3], new double[3], new double[3]}; /* hessian of the constraint */
            double[] newtonStep = new double[3]; /* Newton step in optimization */

            var funcf = 0.0; /* values of function to be minimized */
            double w, ws = 0; /* vertex weight squared or to the 1.5 */
            bool pdflag; /* converging to non-minimum? */

            /* Set parameters. */

            var sqrtVertexCount = Math.Sqrt((double) vertexCount);
            var maxConstraint = 1.0e-12 * sqrtVertexCount; /* maximum allowed value for constraint */
            var mstart = 5.0d * sqrtVertexCount; /* starting value for constraint scaling */

            for (var i = 0; i < 25; i++)
            {
                coeffs[i] = 0;
            }

            var aptr = yvecs[1] + 1;
            var bptr = yvecs[2] + 1;
            var cptr = yvecs[3] + 1;
            var wsptr = vwsqrt + 1;
            for (var i = 1; i <= nmyvtxs; i++)
            {
                var a = *aptr++;
                var b = *bptr++;
                var c = *cptr++;
                w = graph[i]->vwgt;
                if (useVertexWeights)
                {
                    ws = *wsptr++;
                }

                if (w == 1.0)
                {
                    coeffs[0] += a * a * a * a;
                    coeffs[1] += b * b * b * b;
                    coeffs[2] += c * c * c * c;
                    coeffs[3] += a * a * a * b;
                    coeffs[4] += a * a * b * b;
                    coeffs[5] += a * b * b * b;
                    coeffs[6] += a * a * a * c;
                    coeffs[7] += a * a * c * c;
                    coeffs[8] += a * c * c * c;
                    coeffs[9] += b * b * b * c;
                    coeffs[10] += b * b * c * c;
                    coeffs[11] += b * c * c * c;
                    coeffs[12] += a * a * b * c;
                    coeffs[13] += a * b * b * c;
                    coeffs[14] += a * b * c * c;

                    coeffs[15] += a * a * a;
                    coeffs[16] += b * b * b;
                    coeffs[17] += c * c * c;
                    coeffs[18] += a * a * b;
                    coeffs[19] += a * a * c;
                    coeffs[20] += a * b * b;
                    coeffs[21] += b * b * c;
                    coeffs[22] += a * c * c;
                    coeffs[23] += b * c * c;
                    coeffs[24] += a * b * c;
                }
                else
                {
                    w = 1 / (w * w);
                    ws = 1 / ws;
                    coeffs[0] += a * a * a * a * w;
                    coeffs[1] += b * b * b * b * w;
                    coeffs[2] += c * c * c * c * w;
                    coeffs[3] += a * a * a * b * w;
                    coeffs[4] += a * a * b * b * w;
                    coeffs[5] += a * b * b * b * w;
                    coeffs[6] += a * a * a * c * w;
                    coeffs[7] += a * a * c * c * w;
                    coeffs[8] += a * c * c * c * w;
                    coeffs[9] += b * b * b * c * w;
                    coeffs[10] += b * b * c * c * w;
                    coeffs[11] += b * c * c * c * w;
                    coeffs[12] += a * a * b * c * w;
                    coeffs[13] += a * b * b * c * w;
                    coeffs[14] += a * b * c * c * w;

                    coeffs[15] += a * a * a * ws;
                    coeffs[16] += b * b * b * ws;
                    coeffs[17] += c * c * c * ws;
                    coeffs[18] += a * a * b * ws;
                    coeffs[19] += a * a * c * ws;
                    coeffs[20] += a * b * b * ws;
                    coeffs[21] += b * b * c * ws;
                    coeffs[22] += a * c * c * ws;
                    coeffs[23] += b * c * c * ws;
                    coeffs[24] += a * b * c * ws;
                }
            }

            /* Adjust for normalization of eigenvectors. */
            /* This should make convergence criteria insensitive to problem size. */
            /* Note that the relative sizes of funcf and funcc depend on normalization of eigenvectors, and I'm assuming them normalized to 1. */
            for (var i = 0; i < 15; i++)
            {
                coeffs[i] *= vertexCount;
            }

            var sqrtVertexCount1 = Math.Sqrt((double) vertexCount);
            for (var i = 15; i < 25; i++)
            {
                coeffs[i] *= sqrtVertexCount1;
            }

            double currentBestValue = 0; /* value of best minimizer so far */
            for (var ntries = 1; ntries <= OPT3D_NTRIES; ntries++)
            {
                /* Initialize the starting guess randomly. */
                angularVariables[0] = 2.0d * Math.PI * (drandom() - .5d);
                angularVariables[1] = Math.Acos((2.0d * drandom()) - 1.0d) - (Math.PI / 2.0d);
                angularVariables[2] = 2.0d * Math.PI * (drandom() - .5d);

                var inner1 = 0; /* number of iterations at each stage */
                var totalIterations = 0; /* total number of iterations */
                var mult = mstart; /* multiplier for constraint violation */
                var stepMin = EarlyStepMinimum; /* minimum step => convergence */
                var funcc = maxConstraint; /* values of function to be minimized */
                while (funcc >= maxConstraint && totalIterations < 70)
                {
                    var inner = 0; /* number of iterations at each stage */
                    var step_size = stepMin; /* norm of step */
                    pdflag = false;
                    var gradientNorm = 0.0d; /* norm of the gradient */
                    while (step_size >= stepMin && (!pdflag || gradientNorm > MinimumAcceptableGradiant) && inner < 15)
                    {
                        funcf = func3d(coeffs, angularVariables[0], angularVariables[1], angularVariables[2]);
                        grad3d(coeffs, functionGradient, angularVariables[0], angularVariables[1], angularVariables[2]);

                        hess3d(coeffs, functionHessian);

                        fixed (double* coeffPtr = coeffs)
                        {
                            /* Compute contribution of constraint term. */
                            funcc = constraint(&coeffPtr[15]);
                            /* func = funcf + mult*funcc; */
                            gradcon(&coeffPtr[15], constraintGradient);
                            hesscon(&coeffPtr[15], constraintHessian);
                        }

                        /* If in final pass, tighten convergence criterion. */
                        if (funcc < maxConstraint)
                        {
                            stepMin = FinalStepMinimum;
                        }

                        double smallestEigenValue; /* smallest eigenvalue of Hessian */
                        double residual; /* returned eigen-residual */
                        var kk = 0;
                        if (kk != 0)
                        {
                            ch_evals3(constraintHessian, &smallestEigenValue, &residual, &residual);
                        }

                        for (var i = 0; i < 3; i++)
                        {
                            /* Note: I'm taking negative of gradient here. */
                            functionGradient[i] = -functionGradient[i] - mult * constraintGradient[i];
                            for (var j = 0; j < 3; j++)
                            {
                                functionHessian[i][j] += mult * constraintHessian[i][j];
                            }
                        }

                        gradientNorm = Math.Abs(functionGradient[0]) + Math.Abs(functionGradient[1]) + Math.Abs(functionGradient[2]);
                        var hessMin = HessianFactor * gradientNorm / MaximumAllowedStep; /* value for hessian if < 0 */
                        if (hessMin < HessTol)
                        {
                            hessMin = HessTol;
                        }

                        /* Find smallest eigenvalue of hess. */
                        ch_evals3(functionHessian, &smallestEigenValue, &residual, &residual);

                        /* If eval < 0, add to diagonal to make pos def. */
                        pdflag = !(smallestEigenValue < -Pdtol);

                        if (smallestEigenValue < hessMin)
                        {
                            for (var i = 0; i < 3; i++)
                            {
                                functionHessian[i][i] += hessMin - smallestEigenValue;
                            }
                        }

                        /* Now solve linear system for step sizes. */
                        kramer3(functionHessian, functionGradient, newtonStep);

                        /* Scale step down if too big. */
                        step_size = Math.Abs(newtonStep[0]) + Math.Abs(newtonStep[1]) + Math.Abs(newtonStep[2]);
                        if (step_size > MaximumAllowedStep)
                        {
                            var a = MaximumAllowedStep / step_size;
                            for (var i = 0; i < 3; i++)
                            {
                                newtonStep[i] *= a;
                            }
                        }

                        if ((step_size < stepMin || gradientNorm < MinimumAcceptableGradiant) && !pdflag)
                        {
                            /* Convergence to non-min. */
                            for (var i = 0; i < 3; i++)
                            {
                                functionHessian[i][i] -= hessMin - smallestEigenValue;
                            }

                            residual = ch_eigenvec3(functionHessian, smallestEigenValue, newtonStep);
                            step_size = Math.Abs(newtonStep[0]) + Math.Abs(newtonStep[1]) + Math.Abs(newtonStep[2]);
                            var a = stepMin / step_size;
                            for (var i = 0; i < 3; i++)
                            {
                                newtonStep[i] *= a;
                            }

                            step_size = stepMin;
                        }

                        for (var i = 0; i < 3; i++)
                        {
                            angularVariables[i] += newtonStep[i];
                        }

                        inner++;
                    }

                    if (inner1 == 0)
                    {
                        inner1 = inner;
                    }

                    totalIterations += inner;
                    mult *= ConstraintGrowthScalingFactor;
                }

                if (DEBUG_OPTIMIZE)
                {
                    Console.WriteLine("On try {0:d}, After {1:d} ({2:d}) passes, funcf={3:e}, funcc={4:e} ({5:f}, {6:f}, {7:f})", ntries, totalIterations, inner1, funcf, funcc, angularVariables[0], angularVariables[1], angularVariables[2]);
                }

                if (ntries == 1 || funcf < currentBestValue)
                {
                    currentBestValue = funcf;
                    for (var i = 0; i < 3; i++)
                    {
                        currentBestMinimizer[i] = angularVariables[i];
                    }
                }
            }

            *ptheta = currentBestMinimizer[0];
            *pphi = currentBestMinimizer[1];
            *pgamma = currentBestMinimizer[2];
        }

        /* Solve Ax=b */
        private static void kramer3( /* Use Kramer's rule to solve 3x3 */
            double[][] A /*[3][3]*/,
            double[] b /*[3]*/,
            double[] x /*[3]*/
        )
        {
            var det = 1.0d / determinant(A, 3);/* determinant of system */

            x[0] = (b[0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) -
                    b[1] * (A[0][1] * A[2][2] - A[0][2] * A[2][1]) +
                    b[2] * (A[0][1] * A[1][2] - A[0][2] * A[1][1])) *
                   det;

            x[1] = -(b[0] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) -
                     b[1] * (A[0][0] * A[2][2] - A[0][2] * A[2][0]) +
                     b[2] * (A[0][0] * A[1][2] - A[0][2] * A[1][0])) *
                   det;

            x[2] = (b[0] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]) -
                    b[1] * (A[0][0] * A[2][1] - A[0][1] * A[2][0]) +
                    b[2] * (A[0][0] * A[1][1] - A[0][1] * A[1][0])) *
                   det;
        }

        /// <summary>
        /// Find the eigenvector of symmetric 3x3 matrix w/ given eigenvalue.
        /// </summary>
        /// <param name="A"> 3x3 matrix to find eigenvector of</param>
        /// <param name="eval">eigenvalue</param>
        /// <param name="evec">eigenvector returned ([3])</param>
        /// <param name="res">returned error estimate</param>
        private static double ch_eigenvec3(double[][] A /*[3][3]*/,
            double eval,
            double[] evec /*[3]*/
        )
        {
            //double mat[3][3];            /* copy of A to write over */
            double[][] mat = {new[] {0.0, 0.0, 0.0}, new[] {0.0, 0.0, 0.0}, new[] {0.0, 0.0, 0.0}}; /* copy of A to write over */
            int[] ind = new int[3]; /* permutation indices */
            double ex, ey, ez; /* elements of eigenvector returned */
            double xmax; /* maximum value in matrix */
            double tmp; /* intermediate values */
            double norm; /* norm of eigenvector */
            double res1, res2, res3; /* elements of residual vector */
            const double tolerance = 1.0e-6; /* smaller value assumed to be zero */
            int imax = -1, jmax = -1; /* indices of max value in matrix */
            int i, j; /* loop counters */

            for (i = 0; i < 3; i++)
            {
                for (j = 0; j < 3; j++)
                {
                    mat[i][j] = A[i][j];
                }
            }

            for (i = 0; i < 3; i++)
            {
                mat[i][i] -= eval;
            }

            ind[0] = 0;
            ind[1] = 1;
            ind[2] = 2;

            /* Find the largest element in the matrix. */
            xmax = 0.0;
            for (i = 0; i < 3; i++)
            {
                for (j = i; j < 3; j++)
                {
                    if (Math.Abs(mat[i][j]) > xmax)
                    {
                        imax = i;
                        jmax = j;
                        xmax = Math.Abs(mat[i][j]);
                    }
                }
            }

            if (xmax == 0.0)
            {
                /* Handle completely degenerate case first. */
                evec[0] = 1.0;
                evec[1] = evec[2] = 0.0;
            }

            else
            {
                /* Scale the matrix so largest value is 1.0 */
                for (i = 0; i < 3; i++)
                {
                    for (j = 0; j < 3; j++)
                    {
                        mat[i][j] /= xmax;
                    }
                }

                /* Swap rows if necessary to move max element to first row. */
                if (imax != 0)
                {
                    for (j = 0; j < 3; j++)
                    {
                        tmp = mat[0][j];
                        mat[0][j] = mat[imax][j];
                        mat[imax][j] = tmp;
                    }
                }

                /* Swap columns if necessary to move max element to first position. */
                if (jmax != 0)
                {
                    for (i = 0; i < 3; i++)
                    {
                        tmp = mat[i][0];
                        mat[i][0] = mat[i][jmax];
                        mat[i][jmax] = tmp;
                    }

                    ind[0] = jmax;
                    ind[jmax] = 0;
                }

                /* Reduce matrix to 2x2 by subtracting off first row. */
                for (i = 1; i < 3; i++)
                {
                    for (j = 1; j < 3; j++)
                    {
                        mat[i][j] = mat[0][0] * mat[i][j] - mat[i][0] * mat[0][j];
                    }
                }

                /* Find maximum element in reduced 2x2 matrix. */
                xmax = 0.0;
                for (i = 1; i < 3; i++)
                {
                    for (j = i; j < 3; j++)
                    {
                        if (Math.Abs(mat[i][j]) > xmax)
                        {
                            imax = i;
                            jmax = j;
                            xmax = Math.Abs(mat[i][j]);
                        }
                    }
                }

                if (xmax < tolerance)
                {
                    /* Handle 2-fold degenerate case - skip to end. */
                    ey = 1.0;
                    ex = ez = 0;
                }

                else
                {
                    /* Swap rows 2 and 3 to move max element to 2nd row. */
                    if (imax != 1)
                    {
                        for (j = 0; j < 3; j++)
                        {
                            tmp = mat[1][j];
                            mat[1][j] = mat[imax][j];
                            mat[imax][j] = tmp;
                        }
                    }

                    /* Swap columns to move max element to (1,1) position. */
                    if (jmax != 1)
                    {
                        for (i = 0; i < 3; i++)
                        {
                            tmp = mat[i][1];
                            mat[i][1] = mat[i][2];
                            mat[i][2] = tmp;
                        }

                        i = ind[1];
                        ind[1] = ind[2];
                        ind[2] = i;
                    }

                    /* Compute eigenvector from numerically stabilized matrix. */
                    ez = mat[0][0] * mat[1][1];
                    ey = -mat[1][2] * mat[0][0];
                    ex = mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1];
                }

                /* Reorder the e-vector to undo pivoting - end of 3-D case. */
                evec[ind[0]] = ex;
                evec[ind[1]] = ey;
                evec[ind[2]] = ez;
            }

            /* Normalize eigenvector and calculate a normalized eigen-residual. */
            norm = Math.Sqrt(evec[0] * evec[0] + evec[1] * evec[1] + evec[2] * evec[2]);
            for (i = 0; i < 3; i++)
            {
                evec[i] /= norm;
            }

            res1 = (A[0][0] - eval) * evec[0] + A[0][1] * evec[1] + A[0][2] * evec[2];
            res2 = A[1][0] * evec[0] + (A[1][1] - eval) * evec[1] + A[1][2] * evec[2];
            res3 = A[2][0] * evec[0] + A[2][1] * evec[1] + (A[2][2] - eval) * evec[2];
            var res = Math.Sqrt(res1 * res1 + res2 * res2 + res3 * res3);

            /* Now normalize the residual */
            res1 = Math.Abs(A[0][0]) + Math.Abs(A[0][1]) + Math.Abs(A[0][2]);
            res2 = Math.Abs(A[1][0]) + Math.Abs(A[1][1]) + Math.Abs(A[1][2]);
            res3 = Math.Abs(A[2][0]) + Math.Abs(A[2][1]) + Math.Abs(A[2][2]);
            res2 = Math.Max(res2, res3);
            res /= Math.Max(res1, res2);
            return res;
        }

        /* Find eigenvalues of 3x3 symmetric system by solving cubic. */
        private static void ch_evals3(double[][] H /*[3][3]*/, /* 3x3 sym matrix for lowest eigenvalue */
            double* eval1, /* smallest eigenvalue */
            double* eval2, /* middle eigenvalue */
            double* eval3 /* largest eigenvalue */
        )
        {
            //double mat[3][3] = {{0.0}}; /* scaled version of H */
            double[][] mat = {new[] {0.0, 0.0, 0.0}, new[] {0.0, 0.0, 0.0}, new[] {0.0, 0.0, 0.0}}; /* scaled version of H */
            double a1, a2, a3; /* coefficients of cubic equation */
            double q, r; /* intermediate terms */
            double q3, r2; /* powers of q and r */
            double theta; /* angular parameter */
            double root1, root2, root3; /* three roots of cubic */
            const double tol = 1.0e-6; /* allowed deviation */
            double xmax; /* largest matrix element for scaling */
            int i, j; /* loop indices */

            /* This first requires solving a cubic equation. */
            /* Normalize to avoid any numerical problems. */
            xmax = 0.0;
            for (i = 0; i < 3; i++)
            {
                for (j = i; j < 3; j++)
                {
                    if (Math.Abs(H[i][j]) > xmax)
                    {
                        xmax = Math.Abs(H[i][j]);
                    }
                }
            }

            if (xmax != 0)
            {
                for (i = 0; i < 3; i++)
                {
                    for (j = 0; j < 3; j++)
                    {
                        mat[i][j] = H[i][j] / xmax;
                    }
                }
            }

            a1 = -(mat[0][0] + mat[1][1] + mat[2][2]);
            a2 = (mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]) +
                 (mat[0][0] * mat[2][2] - mat[0][2] * mat[2][0]) +
                 (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]);
            a3 = -determinant(mat, 3);

            if (a3 == 0)
            {
                root1 = 0;
                /* Solve quadratic. */
                q = -.5 * (a1 + Sign(a1) * Math.Sqrt(a1 * a1 - 4 * a2));
                root2 = q;
                root3 = a2 / q;
            }

            else
            {
                /* solve cubic */
                q = (a1 * a1 - 3 * a2) / 9;
                r = (2 * a1 * a1 * a1 - 9 * a1 * a2 + 27 * a3) / 54;
                q3 = q * q * q;
                r2 = r * r;

                /* To avoid missing a root, check for roundoff. */
                if ((q3 < r2) && Math.Abs(q3 - r2) < tol * (Math.Abs(q3) + Math.Abs(r2)))
                {
                    q3 = r2;
                }

                if (q3 >= r2)
                {
                    /* Three real roots. */
                    if (r == 0)
                    {
                        theta = Math.PI / 2;
                    }
                    else
                    {
                        q3 = Math.Sqrt(q3);
                        if (q3 < Math.Abs(r))
                        {
                            q3 = Math.Abs(r);
                        }

                        theta = Math.Acos(r / q3);
                    }

                    q = -2 * Math.Sqrt(q);

                    root1 = q * Math.Cos(theta / 3) - a1 / 3;
                    root2 = q * Math.Cos((theta + 2 * Math.PI) / 3) - a1 / 3;
                    root3 = q * Math.Cos((theta + 2 * 2 * Math.PI) / 3) - a1 / 3;
                }
                else
                {
                    /* Only one real root. */
                    theta = Math.Sqrt(r2 - q3) + Math.Abs(r);
                    theta = Math.Pow(theta, 1.0 / 3.0);

                    root1 = root2 = root3 = -Sign(r) * (theta + q / theta) - a1 / 3;
                }
            }

            root1 *= xmax;
            root2 *= xmax;
            root3 *= xmax;
            *eval1 = Math.Min(root1, root2);
            *eval1 = Math.Min(*eval1, root3);
            *eval3 = Math.Max(root1, root2);
            *eval3 = Math.Max(*eval3, root3);
            if (root1 != *eval1 && root1 != *eval3)
            {
                *eval2 = root1;
            }
            else if (root2 != *eval1 && root2 != *eval3)
            {
                *eval2 = root2;
            }
            else
            {
                *eval2 = root3;
            }
        }

        private static double Sign(double x)
        {
            return x < 0 ? -1.0 : 1.0;
        }
    }
}
