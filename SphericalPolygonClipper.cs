using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SphericalVoronoiLib
{
    public static class SphericalPolygonClipper
    {
        // Clips polygon to band zMin <= z <= zMax. If zMin is null, no lower clipping.
        // If zMax is null, no upper clipping.
        public static List<SphericalPoint> ClipToBand(List<SphericalPoint> poly, double? zMin, double? zMax)
        {
            if (poly == null || poly.Count == 0) return new List<SphericalPoint>();

            var working = new List<SphericalPoint>(poly);

            if (zMin.HasValue)
                working = ClipAtZMin(working, zMin.Value);

            if (working == null || working.Count == 0)
                return new List<SphericalPoint>();

            if (zMax.HasValue)
                working = ClipAtZMax(working, zMax.Value);

            if (working == null)
                return new List<SphericalPoint>();

            // Clean duplicate consecutive vertices and close if necessary
            var cleaned = new List<SphericalPoint>();
            SphericalPoint prev = null;
            foreach (var p in working)
            {
                if (prev == null || !AreSamePoint(prev, p, 1e-12))
                {
                    cleaned.Add(p);
                    prev = p;
                }
            }
            if (cleaned.Count > 1 && AreSamePoint(cleaned[0], cleaned[cleaned.Count - 1], 1e-12))
                cleaned.RemoveAt(cleaned.Count - 1);

            return cleaned;
        }

        private static List<SphericalPoint> ClipAtZMin(List<SphericalPoint> poly, double zMin)
        {
            var output = new List<SphericalPoint>();
            if (poly == null || poly.Count == 0) return output;

            for (int i = 0; i < poly.Count; i++)
            {
                var A = poly[i];
                var B = poly[(i + 1) % poly.Count];

                bool Ain = A.Z >= zMin;
                bool Bin = B.Z >= zMin;

                if (Ain && Bin)
                {
                    output.Add(B);
                }
                else if (Ain && !Bin)
                {
                    var I = IntersectEdgeAtZ(A, B, zMin);
                    if (I != null) output.Add(I);
                }
                else if (!Ain && Bin)
                {
                    var I = IntersectEdgeAtZ(A, B, zMin);
                    if (I != null) output.Add(I);
                    output.Add(B);
                }
            }

            return output;
        }

        private static List<SphericalPoint> ClipAtZMax(List<SphericalPoint> poly, double zMax)
        {
            var output = new List<SphericalPoint>();
            if (poly == null || poly.Count == 0) return output;

            for (int i = 0; i < poly.Count; i++)
            {
                var A = poly[i];
                var B = poly[(i + 1) % poly.Count];

                bool Ain = A.Z <= zMax;
                bool Bin = B.Z <= zMax;

                if (Ain && Bin)
                {
                    output.Add(B);
                }
                else if (Ain && !Bin)
                {
                    var I = IntersectEdgeAtZ(A, B, zMax);
                    if (I != null) output.Add(I);
                }
                else if (!Ain && Bin)
                {
                    var I = IntersectEdgeAtZ(A, B, zMax);
                    if (I != null) output.Add(I);
                    output.Add(B);
                }
            }

            return output;
        }

        // Intersection of the great-circle edge (A->B) with the plane z = zLimit.
        // Returns null when there's no intersection (should be rare for well-formed Voronoi polygons).
        private static SphericalPoint IntersectEdgeAtZ(SphericalPoint A, SphericalPoint B, double zLimit)
        {
            double[] a = A.ToArray();
            double[] b = B.ToArray();

            // great-circle plane normal (g · x = 0)
            double[] g = Cross(a, b);
            double gx = g[0], gy = g[1], gz = g[2];

            const double EPS = 1e-15;

            // If zLimit is outside [-1,1] no intersection on unit sphere
            if (zLimit > 1.0 + 1e-15 || zLimit < -1.0 - 1e-15) return null;

            // radius of circle at zLimit
            double r2 = 1.0 - zLimit * zLimit;
            if (r2 < -1e-15) return null;
            r2 = Math.Max(0.0, r2);

            // handle degenerate case where gx and gy are ~0 (great circle plane ∼ z-axis)
            if (Math.Abs(gx) < EPS && Math.Abs(gy) < EPS)
            {
                // Great circle is nearly the equatorial plane; intersection exists if zLimit is approx zero
                if (Math.Abs(zLimit) > 1e-12)
                    return null;

                // fallback: linear chord intersection then project
                double denom = b[2] - a[2];
                if (Math.Abs(denom) < 1e-15) return null;
                double t = (0 - a[2]) / denom;
                double x = a[0] + t * (b[0] - a[0]);
                double y = a[1] + t * (b[1] - a[1]);
                var p = Normalise(new double[] { x, y, 0.0 });
                return new SphericalPoint(p);
            }

            // Solve gx * x + gy * y = -gz * zLimit
            double c = -gz * zLimit;
            double denomLine = gx * gx + gy * gy;
            if (Math.Abs(denomLine) < EPS) return null;

            double k = c / denomLine;
            double p0x = gx * k;
            double p0y = gy * k;

            // direction along line
            double dx = -gy;
            double dy = gx;
            double dd = dx * dx + dy * dy;
            if (dd < EPS) return null;

            // Solve (p0 + t*d)^2 = r2
            double Bq = 2.0 * (p0x * dx + p0y * dy);
            double Cq = p0x * p0x + p0y * p0y - r2;
            double Aq = dd;

            double disc = Bq * Bq - 4.0 * Aq * Cq;
            if (disc < -1e-15) return null;
            disc = Math.Max(0.0, disc);
            double sqrtD = Math.Sqrt(disc);

            double t1 = (-Bq + sqrtD) / (2.0 * Aq);
            double t2 = (-Bq - sqrtD) / (2.0 * Aq);

            var candidates = new List<double[]>();
            candidates.Add(new double[] { p0x + t1 * dx, p0y + t1 * dy, zLimit });
            candidates.Add(new double[] { p0x + t2 * dx, p0y + t2 * dy, zLimit });

            var candPts = new List<double[]>();
            foreach (var cpt in candidates)
                candPts.Add(Normalise(cpt));

            double thetaAB = AngleBetween(a, b);
            for (int i = 0; i < candPts.Count; i++)
            {
                var cpt = candPts[i];
                double thetaAC = AngleBetween(a, cpt);
                double thetaCB = AngleBetween(cpt, b);
                if (Math.Abs((thetaAC + thetaCB) - thetaAB) < 1e-8)
                    return new SphericalPoint(cpt);
            }

            // Fallback: nearest to chord midpoint
            var mid = Normalise(new double[] { (a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5, (a[2] + b[2]) * 0.5 });
            double bestDist = double.MaxValue;
            double[] best = null;
            foreach (var cpt in candPts)
            {
                double dch = (cpt[0] - mid[0]) * (cpt[0] - mid[0]) + (cpt[1] - mid[1]) * (cpt[1] - mid[1]) + (cpt[2] - mid[2]) * (cpt[2] - mid[2]);
                if (dch < bestDist) { bestDist = dch; best = cpt; }
            }
            if (best != null) return new SphericalPoint(best);

            return null;
        }

        // helpers
        private static double[] Cross(double[] u, double[] v)
        {
            return new double[] {
                u[1]*v[2] - u[2]*v[1],
                u[2]*v[0] - u[0]*v[2],
                u[0]*v[1] - u[1]*v[0]
            };
        }
        private static double[] Normalise(double[] v)
        {
            double n = Math.Sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
            if (n == 0) return new double[] { 0, 0, 0 };
            return new double[] { v[0] / n, v[1] / n, v[2] / n };
        }
        private static bool AreSamePoint(SphericalPoint p, SphericalPoint q, double tol)
        {
            return Math.Abs(p.X - q.X) <= tol && Math.Abs(p.Y - q.Y) <= tol && Math.Abs(p.Z - q.Z) <= tol;
        }
        private static double AngleBetween(double[] u, double[] v)
        {
            double d = u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
            d = Math.Max(-1.0, Math.Min(1.0, d));
            return Math.Acos(d);
        }
    }
}