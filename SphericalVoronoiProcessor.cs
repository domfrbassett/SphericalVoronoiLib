using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;

namespace SphericalVoronoi
{
    public class VoronoiResult
    {
        public double[,] Points { get; set; }
        public double[] SolidAngles { get; set; }
        public double[] Weights { get; set; }
    }
    public static class SphericalVoronoiProcessor
    {
        public static VoronoiResult ComputeFromSofa(string sofaFilePath, bool correctBottomRow = true)
        {
            if (string.IsNullOrEmpty(sofaFilePath) || !File.Exists(sofaFilePath))
                throw new FileNotFoundException("Invalid or missing SOFA file path.", sofaFilePath);

            string positionType;
            double[,] cartesianPositions = SOFAReader.ReadSourcePosition(sofaFilePath, out positionType);

            return ComputeFromArray(cartesianPositions, correctBottomRow);
        }
        public static VoronoiResult ComputeFromArray(double[,] cartesianPositions, bool correctBottomRow = true)
        {
            if (cartesianPositions == null || cartesianPositions.GetLength(0) == 0)
                throw new ArgumentException("Input array is null or empty.", nameof(cartesianPositions));

            var uniquePoints = new HashSet<(double x, double y, double z)>();
            int n = cartesianPositions.GetLength(0);
            for (int i = 0; i < n; i++)
                uniquePoints.Add((cartesianPositions[i, 0], cartesianPositions[i, 1], cartesianPositions[i, 2]));

            var points = uniquePoints.Select(p =>
            {
                double norm = Math.Sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
                return new SphericalPoint(p.x / norm, p.y / norm, p.z / norm);
            }).ToList();

            var result = VoronoiSphere.Compute(points);

            double[] solidAngles = result.SolidAngles.ToArray();

            double[,] pointArray = new double[points.Count, 3];
            for (int i = 0; i < points.Count; i++)
            {
                pointArray[i, 0] = points[i].X;
                pointArray[i, 1] = points[i].Y;
                pointArray[i, 2] = points[i].Z;
            }

            double[] adjustedAngles = correctBottomRow
                ? CorrectBottomRow(pointArray, solidAngles)
                : solidAngles;

            double total = adjustedAngles.Sum();
            double[] weights = adjustedAngles.Select(a => a / total).ToArray();

            return new VoronoiResult
            {
                Points = pointArray,
                SolidAngles = adjustedAngles,
                Weights = weights
            };
        }
        private static double[] CorrectBottomRow(double[,] points, double[] solidAngles)
        {
            double tol = 1e-8;
            int nPoints = points.GetLength(0);

            double[] roundedZ = new double[nPoints];
            for (int i = 0; i < nPoints; i++)
                roundedZ[i] = Math.Round(points[i, 2] / tol) * tol;

            var zLevelDict = new Dictionary<double, int>();
            int[] zIdx = new int[nPoints];
            int groupCounter = 0;

            for (int i = 0; i < nPoints; i++)
            {
                double z = roundedZ[i];
                if (!zLevelDict.ContainsKey(z))
                    zLevelDict[z] = groupCounter++;
                zIdx[i] = zLevelDict[z];
            }

            double[] zLevels = zLevelDict.Keys.ToArray();
            var sorted = zLevels.Select((z, idx) => new { z, idx }).OrderBy(x => x.z).ToList();
            if (sorted.Count < 2)
                throw new Exception("Only one z-row found; cannot apply bottom correction.");

            int bottomGroup = sorted[0].idx;
            int rowAboveGroup = sorted[1].idx;

            bool[] bottomIdx = zIdx.Select(id => id == bottomGroup).ToArray();
            bool[] aboveIdx = zIdx.Select(id => id == rowAboveGroup).ToArray();

            double desiredBottomTotal = solidAngles
                .Where((a, i) => aboveIdx[i])
                .Sum();

            double currentBottomTotal = solidAngles
                .Where((a, i) => bottomIdx[i])
                .Sum();

            double[] adjusted = (double[])solidAngles.Clone();
            int bottomCount = bottomIdx.Count(b => b);

            if (currentBottomTotal == 0)
            {
                double perPoint = desiredBottomTotal / bottomCount;
                for (int i = 0; i < nPoints; i++)
                    if (bottomIdx[i]) adjusted[i] = perPoint;
            }
            else
            {
                double scale = desiredBottomTotal / currentBottomTotal;
                for (int i = 0; i < nPoints; i++)
                    if (bottomIdx[i]) adjusted[i] *= scale;
            }

            double origTotal = solidAngles.Sum();
            double newTotal = adjusted.Sum();
            double renorm = origTotal / newTotal;
            for (int i = 0; i < adjusted.Length; i++)
                adjusted[i] *= renorm;

            return adjusted;
        }
    }
}