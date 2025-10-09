using System;
using System.Collections.Generic;
using System.Linq;

namespace SphericalVoronoi
{
    public static class VoronoiSolidAngle
    {
        public static double[] Compute(List<SphericalPoint> P, List<List<int>> K, List<SphericalPoint> centers = null, bool checkUnit = true)
        {
            if (P == null) throw new ArgumentNullException(nameof(P));
            if (K == null) throw new ArgumentNullException(nameof(K));
            if (centers != null && centers.Count != K.Count)
                throw new ArgumentException("centers length must match K length", nameof(centers));

            int n = K.Count;
            var s = new double[n];

            if (checkUnit)
            {
                foreach (var p in P)
                {
                    double norm = Math.Sqrt(p.X * p.X + p.Y * p.Y + p.Z * p.Z);
                    if (Math.Abs(norm - 1.0) > 1e-6)
                        throw new ArgumentException("VoronoiSolidAngle: P must be unit vectors");
                }
            }

            for (int j = 0; j < n; j++)
            {
                var idxList = K[j];
                if (idxList == null || idxList.Count < 3) { s[j] = 0.0; continue; }

                var verts = idxList.Select(idx => P[idx]).ToList();
                s[j] = (centers == null) ? OneCellSolidAngle(verts) : OneCellSolidAngle(verts, centers[j]);
            }

            return s;
        }
        private static double OneCellSolidAngle(List<SphericalPoint> verts, SphericalPoint center = null)
        {
            if (verts == null) throw new ArgumentNullException(nameof(verts));
            if (verts.Count < 3) return 0.0;

            double total = 0.0;

            if (center == null)
            {
                var v0 = verts[0].ToArray();
                for (int i = 1; i < verts.Count - 1; i++)
                {
                    var a = v0;
                    var b = verts[i].ToArray();
                    var c = verts[i + 1].ToArray();
                    double triple = Dot(a, Cross(b, c));
                    double denom = 1.0 + Dot(a, b) + Dot(b, c) + Dot(c, a);
                    double area = 2.0 * Math.Atan2(Math.Abs(triple), denom);
                    total += area;
                }
            }
            else
            {
                var ctr = center.ToArray();
                for (int i = 0; i < verts.Count; i++)
                {
                    var a = ctr;
                    var b = verts[i].ToArray();
                    var c = verts[(i + 1) % verts.Count].ToArray();
                    double triple = Dot(a, Cross(b, c));
                    double denom = 1.0 + Dot(a, b) + Dot(b, c) + Dot(c, a);
                    double area = 2.0 * Math.Atan2(Math.Abs(triple), denom);
                    total += area;
                }
            }

            return total;
        }
        private static double[] Cross(double[] u, double[] v) => new[] {
                u[1]*v[2] - u[2]*v[1],
                u[2]*v[0] - u[0]*v[2],
                u[0]*v[1] - u[1]*v[0]
            };

        private static double Dot(double[] u, double[] v) => u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
    }
}