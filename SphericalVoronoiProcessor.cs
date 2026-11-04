using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;

namespace SphericalVoronoiLib
{
    public class VoronoiResult
    {
        /// <summary>Input directions as unit vectors, one row per input row.</summary>
        public double[,] Points { get; set; }
        /// <summary>Solid angle per input row (sr); rows sharing a direction split their cell equally.</summary>
        public double[] SolidAngles { get; set; }
        /// <summary>Weights per input row, summing to 1 over the sphere or the clipped band.</summary>
        public double[] Weights { get; set; }
        /// <summary>Voronoi cell of each input row (index into the distinct directions).</summary>
        public int[] CellIndex { get; set; }
        /// <summary>Number of distinct directions.</summary>
        public int CellCount { get; set; }
    }

    public static class SphericalVoronoiProcessor
    {
        /// <summary>Directions closer than this (Euclidean, on the unit sphere) are treated as one.</summary>
        public const double DuplicateTolerance = 1e-9;

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
        /// Computes Voronoi weights from azimuth/elevation pairs in degrees (columns 0 and 1; a third
        /// distance column is ignored). Optionally restricts the cells to an elevation band in degrees.
        /// </summary>
        public static VoronoiResult ComputeFromSpherical(double[,] azElDeg, double? minElevationDeg = null, double? maxElevationDeg = null)
        {
            if (azElDeg == null || azElDeg.GetLength(0) == 0)
                throw new ArgumentException("Input array is null or empty.", nameof(azElDeg));
            if (azElDeg.GetLength(1) < 2)
                throw new ArgumentException("Input array needs azimuth and elevation columns.", nameof(azElDeg));

            const double Deg = Math.PI / 180.0;
            int n = azElDeg.GetLength(0);
            var xyz = new double[n, 3];
            for (int i = 0; i < n; i++)
            {
                double az = azElDeg[i, 0] * Deg, el = azElDeg[i, 1] * Deg;
                xyz[i, 0] = Math.Cos(el) * Math.Cos(az);
                xyz[i, 1] = Math.Cos(el) * Math.Sin(az);
                xyz[i, 2] = Math.Sin(el);
            }
            double? zMin = minElevationDeg.HasValue ? Math.Sin(Math.Max(-90.0, minElevationDeg.Value) * Deg) : (double?)null;
            double? zMax = maxElevationDeg.HasValue ? Math.Sin(Math.Min(90.0, maxElevationDeg.Value) * Deg) : (double?)null;
            return ComputeFromArray(xyz, zMin, zMax);
        }

        /// <summary>
        /// Computes Voronoi weights from a provided array of Cartesian coordinates.
        /// Optionally restricts the cells to the band zMin &lt;= z &lt;= zMax.
        /// </summary>
        public static VoronoiResult ComputeFromArray(double[,] cartesianPositions, double? zMin = null, double? zMax = null)
        {
            if (cartesianPositions == null || cartesianPositions.GetLength(0) == 0)
                throw new ArgumentException("Input array is null or empty.", nameof(cartesianPositions));
            if (cartesianPositions.GetLength(1) != 3)
                throw new ArgumentException("Input array must have three columns (x, y, z).", nameof(cartesianPositions));
            if (zMin.HasValue && zMax.HasValue && zMin.Value > zMax.Value)
                throw new ArgumentException("zMin must not exceed zMax.");

            int n = cartesianPositions.GetLength(0);
            var unit = new double[n, 3];
            var cellIndex = new int[n];
            var seeds = new List<SphericalPoint>();
            var seen = new Dictionary<(long, long, long), int>();
            double q = 1.0 / DuplicateTolerance;

            for (int i = 0; i < n; i++)
            {
                double x = cartesianPositions[i, 0], y = cartesianPositions[i, 1], z = cartesianPositions[i, 2];
                double norm = Math.Sqrt(x * x + y * y + z * z);
                if (!(norm > 0) || double.IsInfinity(norm))
                    throw new ArgumentException($"Row {i} is not a valid direction.", nameof(cartesianPositions));
                x /= norm; y /= norm; z /= norm;
                unit[i, 0] = x; unit[i, 1] = y; unit[i, 2] = z;

                var key = ((long)Math.Round(x * q), (long)Math.Round(y * q), (long)Math.Round(z * q));
                if (!seen.TryGetValue(key, out int cell))
                {
                    cell = seeds.Count;
                    seen[key] = cell;
                    seeds.Add(new SphericalPoint(x, y, z));
                }
                cellIndex[i] = cell;
            }

            int m = seeds.Count;
            var vor = VoronoiSphere.Compute(seeds);
            double[] cellArea;

            if (!zMin.HasValue && !zMax.HasValue)
            {
                cellArea = vor.SolidAngles;
            }
            else if (m == 1)
            {
                cellArea = new[] { 2 * Math.PI * ((zMax ?? 1.0) - (zMin ?? -1.0)) };
            }
            else
            {
                int north = 0, south = 0;
                for (int k = 1; k < m; k++)
                {
                    if (seeds[k].Z > seeds[north].Z) north = k;
                    if (seeds[k].Z < seeds[south].Z) south = k;
                }
                cellArea = new double[m];
                for (int k = 0; k < m; k++)
                    cellArea[k] = SphericalPolygonClipper.ClippedSolidAngle(
                        vor.VoronoiBoundary[k], seeds[k], zMin, zMax, k == north, k == south);
            }

            var multiplicity = new int[m];
            foreach (int c in cellIndex) multiplicity[c]++;

            var solidAngles = new double[n];
            for (int i = 0; i < n; i++) solidAngles[i] = cellArea[cellIndex[i]] / multiplicity[cellIndex[i]];

            double total = solidAngles.Sum();
            double[] weights = total > 0 ? solidAngles.Select(a => a / total).ToArray() : new double[n];

            return new VoronoiResult
            {
                Points = unit,
                SolidAngles = solidAngles,
                Weights = weights,
                CellIndex = cellIndex,
                CellCount = m
            };
        }
    }
}
