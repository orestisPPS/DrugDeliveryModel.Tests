using System.Collections.Generic;
using System.IO;

namespace MGroup.DrugDeliveryModel.Tests.Commons;

public static class CSVExporter
{
    public static double[,] ConverVectorsTo2DArray(List<double[]> list)
    {
        int cols = list.Count;
        int rows = list[0].Length;
        double[,] result = new double[rows, cols];

        for (int i = 0; i < cols; i++)
        {
            for (int j = 0; j < rows; j++)
            {
                result[j, i] = list[i][j];
            }
        }
        return result;
    }
    
    public static void ExportMatrixToCSV(double[,] matrix, string filePath)
    {
        using (StreamWriter sw = new StreamWriter(filePath))
        {
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                for (int j = 0; j < matrix.GetLength(1); j++)
                {
                    sw.Write(matrix[i, j]);
                    if (j < matrix.GetLength(1) - 1)
                    {
                        sw.Write(",");
                    }
                    if (j == matrix.GetLength(1) - 1)
                    {
                        sw.WriteLine();
                    }
                }
            }
        }
    }
    
    public static void ExportVectorToCSV(double[] vector, string filePath)
    {
        using (StreamWriter sw = new StreamWriter(filePath))
        {
            for (int i = 0; i < vector.Length; i++)
            {
                sw.WriteLine(vector[i]);
            }
        }
    }
}