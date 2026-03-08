using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;

namespace SphericalVoronoiLib
{
    public class VoronoiResult
    {
        public double[,] Points { get; set; }
        public double[] SolidAngles { get; set; }
        public double[] Weights { get; set; }
    }

    public static class SphericalVoronoiProcessor
    {
        /// <summary>
        /// Computes Voronoi weights directly from a SOFA file.
        /// </summary>
        public static VoronoiResult ComputeFromSofa(string sofaFilePath, double? zMin = null, double? zMax = null)
        {
            if (string.IsNullOrEmpty(sofaFilePath) || !File.Exists(sofaFilePath))
                throw new FileNotFoundException("Invalid or missing SOFA file path.", sofaFilePath);

            string positionType;
            double[,] cartesianPositions = SOFAReader.ReadSourcePosition(sofaFilePath, out positionType);

            return ComputeFromArray(cartesianPositions, zMin, zMax);
        }

        /// <summary>
        /// Computes Voronoi weights from a provided array of Cartesian coordinates.
        /// </summary>
        public static VoronoiResult ComputeFromArray(double[,] cartesianPositions, double? zMin = null, double? zMax = null)
        {
            if (cartesianPositions == null || cartesianPositions.GetLength(0) == 0)
                throw new ArgumentException("Input array is null or empty.", nameof(cartesianPositions));

            // Remove duplicates and ensure unit vectors
            var uniquePoints = new HashSet<(double x, double y, double z)>();
            int n = cartesianPositions.GetLength(0);
            for (int i = 0; i < n; i++)
            {
                uniquePoints.Add((cartesianPositions[i, 0], cartesianPositions[i, 1], cartesianPositions[i, 2]));
            }

            var points = uniquePoints.Select(p =>
            {
                double norm = Math.Sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
                return new SphericalPoint(p.x / norm, p.y / norm, p.z / norm);
            }).ToList();

            // Compute the spherical Voronoi tesselation (full sphere)
            var result = VoronoiSphere.Compute(points);

            int m = result.VoronoiBoundary.Count;
            double[] solidAngles = new double[m];

            bool doClip = zMin.HasValue || zMax.HasValue;

            for (int i = 0; i < m; i++)
            {
                var poly = result.VoronoiBoundary[i];

                // Optionally clip to band
                List<SphericalPoint> clipped = poly;
                if (doClip)
                {
                    clipped = SphericalPolygonClipper.ClipToBand(poly, zMin, zMax);
                }

                if (clipped == null || clipped.Count < 3)
                {
                    solidAngles[i] = 0.0;
                }
                else
                {
                    solidAngles[i] = VoronoiSolidAngle.PolygonSolidAngle(clipped);
                }
            }

            // Prepare output point array (input sample points normalised)
            double[,] pointArray = new double[points.Count, 3];
            for (int i = 0; i < points.Count; i++)
            {
                pointArray[i, 0] = points[i].X;
                pointArray[i, 1] = points[i].Y;
                pointArray[i, 2] = points[i].Z;
            }

            // Calculate normalised weights (sum to 1 over the clipped region if clipped, or over full sphere)
            double total = solidAngles.Sum();
            double[] weights = (total > 0)
                ? solidAngles.Select(a => a / total).ToArray()
                : new double[solidAngles.Length];

            return new VoronoiResult
            {
                Points = pointArray,
                SolidAngles = solidAngles,
                Weights = weights
            };
        }
    }
}