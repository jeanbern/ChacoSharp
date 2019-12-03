namespace ChacoSharp.Optimize
{
    public static class Determinant
    {
        public static double determinant(double[][] M/*[3][3]*/, int ndims)
        {
            if (ndims == 1) {
                return (M[0][0]);
            }
            else if (ndims == 2) {
                return (M[0][0] * M[1][1] - M[0][1] * M[1][0]);
            }

            else if (ndims == 3) {
                return (M[0][0] * (M[1][1] * M[2][2] - M[2][1] * M[1][2]) -
                        M[1][0] * (M[0][1] * M[2][2] - M[2][1] * M[0][2]) +
                        M[2][0] * (M[0][1] * M[1][2] - M[1][1] * M[0][2]));
            }

            else {
                return (0.0);
            }
        }
    }
}
