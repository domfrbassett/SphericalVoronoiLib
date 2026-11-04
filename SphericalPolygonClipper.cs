using System;
using System.Collections.Generic;

namespace SphericalVoronoiLib
{
    /// <summary>
    /// Clips spherical polygons to the band zMin &lt;= z &lt;= zMax. Edges created by the clip follow
    /// the latitude circle, not a great circle, and solid angles of the clipped region are exact.
    /// </summary>
    public static class SphericalPolygonClipper
    {
        // Longitude step used when tracing a latitude-circle edge
        private const double CircleStep = 0.05;

        // Vertex of a clipped boundary; CircleZ/Dir describe the edge to the next vertex
        // (null: great-circle edge; Dir +1: eastwards along z = CircleZ, -1: westwards).
        private struct Node
        {
            public double[] P;
            public double? CircleZ;
            public int Dir;
            public Node(double[] p, double? circleZ = null, int dir = 0) { P = p; CircleZ = circleZ; Dir = dir; }
        }

        /// <summary>
        /// Clips a polygon to the band and returns its boundary, with latitude-circle edges sampled densely.
        /// </summary>
        public static List<SphericalPoint> ClipToBand(List<SphericalPoint> poly, double? zMin, double? zMax)
        {
            var output = new List<SphericalPoint>();
            if (poly == null || poly.Count == 0) return output;

            var nodes = Clip(ToNodes(poly), zMin, zMax, out _, out _);
            for (int i = 0; i < nodes.Count; i++)
            {
                output.Add(new SphericalPoint(nodes[i].P));
                if (nodes[i].CircleZ.HasValue)
                {
                    var b = nodes[(i + 1) % nodes.Count].P;
                    foreach (var p in CircleInterior(nodes[i].P, b, nodes[i].CircleZ.Value, nodes[i].Dir))
                        output.Add(new SphericalPoint(p));
                }
            }
            return output;
        }

        /// <summary>
        /// Solid angle of the part of a Voronoi cell inside the band.
        /// </summary>
        /// <param name="cell">Cell boundary (either orientation).</param>
        /// <param name="seed">Seed of the cell, used as the fan centre; it must not be antipodal to any boundary point.</param>
        /// <param name="containsNorthPole">True if the cell contains (0, 0, 1).</param>
        /// <param name="containsSouthPole">True if the cell contains (0, 0, -1).</param>
        public static double ClippedSolidAngle(List<SphericalPoint> cell, SphericalPoint seed, double? zMin, double? zMax,
            bool containsNorthPole, bool containsSouthPole)
        {
            if (cell == null || cell.Count < 3) return 0.0;
            var s = seed.ToArray();

            var nodes = ToNodes(cell);
            if (FanArea(nodes, s) < 0) nodes.Reverse();

            nodes = Clip(nodes, zMin, zMax, out bool crossedMin, out bool crossedMax);
            if (nodes.Count < 3) return 0.0;

            double area = FanArea(nodes, s);
            for (int i = 0; i < nodes.Count; i++)
            {
                if (!nodes[i].CircleZ.HasValue) continue;
                var a = nodes[i].P;
                var b = nodes[(i + 1) % nodes.Count].P;
                double c = nodes[i].CircleZ.Value;
                var pts = new List<double[]> { a };
                pts.AddRange(CircleInterior(a, b, c, nodes[i].Dir));
                pts.Add(b);
                // Replace the chord a-b by the sampled arc, then add the exact lens between each sub-chord and the circle
                area -= TriangleArea(s, a, b);
                for (int k = 0; k + 1 < pts.Count; k++)
                {
                    area += TriangleArea(s, pts[k], pts[k + 1]);
                    area += Lens(pts[k], pts[k + 1], c);
                }
            }

            const double FourPi = 4 * Math.PI;
            if (area < 0 && area > -1e-9) area = 0;
            area -= FourPi * Math.Floor(area / FourPi);

            // Latitude circle wholly inside the cell: remove the polar cap
            if (zMax.HasValue && containsNorthPole && !crossedMax && zMax.Value < 1)
                area -= 2 * Math.PI * (1 - zMax.Value);
            if (zMin.HasValue && containsSouthPole && !crossedMin && zMin.Value > -1)
                area -= 2 * Math.PI * (1 + zMin.Value);

            return Math.Max(0.0, area);
        }

        private static List<Node> ToNodes(List<SphericalPoint> poly)
        {
            var nodes = new List<Node>(poly.Count);
            foreach (var p in poly)
            {
                var q = p.ToArray();
                if (nodes.Count == 0 || !Same(nodes[nodes.Count - 1].P, q)) nodes.Add(new Node(q));
            }
            if (nodes.Count > 1 && Same(nodes[0].P, nodes[nodes.Count - 1].P)) nodes.RemoveAt(nodes.Count - 1);
            return nodes;
        }

        private static List<Node> Clip(List<Node> nodes, double? zMin, double? zMax, out bool crossedMin, out bool crossedMax)
        {
            crossedMin = crossedMax = false;
            if (zMin.HasValue) nodes = ClipHalf(nodes, zMin.Value, +1, out crossedMin);
            if (zMax.HasValue && nodes.Count > 0) nodes = ClipHalf(nodes, zMax.Value, -1, out crossedMax);

            var cleaned = new List<Node>(nodes.Count);
            foreach (var n in nodes)
            {
                if (cleaned.Count > 0 && Same(cleaned[cleaned.Count - 1].P, n.P))
                {
                    if (n.CircleZ.HasValue) cleaned[cleaned.Count - 1] = n;
                    continue;
                }
                cleaned.Add(n);
            }
            if (cleaned.Count > 1 && Same(cleaned[0].P, cleaned[cleaned.Count - 1].P))
            {
                var last = cleaned[cleaned.Count - 1];
                if (last.CircleZ.HasValue && !cleaned[0].CircleZ.HasValue) cleaned[0] = new Node(cleaned[0].P, last.CircleZ, last.Dir);
                cleaned.RemoveAt(cleaned.Count - 1);
            }
            return cleaned;
        }

        // Sutherland-Hodgman against one latitude plane; keeps sign*(z - c) >= 0.
        // Exits are tagged so the edge to the next entry runs along the circle with the kept region on its left.
        private static List<Node> ClipHalf(List<Node> poly, double c, int sign, out bool crossed)
        {
            crossed = false;
            var output = new List<Node>();
            int n = poly.Count;
            for (int i = 0; i < n; i++)
            {
                var A = poly[i];
                var B = poly[(i + 1) % n];
                bool Ain = sign * (A.P[2] - c) >= 0;
                bool Bin = sign * (B.P[2] - c) >= 0;

                if (Ain && Bin) output.Add(B);
                else if (Ain && !Bin)
                {
                    crossed = true;
                    output.Add(new Node(IntersectEdgeAtZ(A.P, B.P, c), c, sign));
                }
                else if (!Ain && Bin)
                {
                    crossed = true;
                    output.Add(new Node(IntersectEdgeAtZ(A.P, B.P, c)));
                    output.Add(B);
                }
            }
            return output;
        }

        // Point where the great-circle edge a-b meets the plane z = c (a and b on opposite sides).
        private static double[] IntersectEdgeAtZ(double[] a, double[] b, double c)
        {
            double theta = Math.Acos(Math.Max(-1.0, Math.Min(1.0, Dot(a, b))));
            if (theta < 1e-15) return a;
            double lo = 0, hi = 1, za = a[2] - c;
            for (int it = 0; it < 60; it++)
            {
                double t = 0.5 * (lo + hi);
                double z = Slerp(a, b, theta, t)[2] - c;
                if (Math.Sign(z) == Math.Sign(za)) lo = t; else hi = t;
            }
            var p = Slerp(a, b, theta, 0.5 * (lo + hi));
            // Snap onto the circle
            double r = Math.Sqrt(Math.Max(0.0, 1 - c * c)), h = Math.Sqrt(p[0] * p[0] + p[1] * p[1]);
            return h > 0 ? new[] { p[0] * r / h, p[1] * r / h, c } : new[] { r, 0.0, c };
        }

        private static double[] Slerp(double[] a, double[] b, double theta, double t)
        {
            double s = Math.Sin(theta), wa = Math.Sin((1 - t) * theta) / s, wb = Math.Sin(t * theta) / s;
            return Normalise(new[] { wa * a[0] + wb * b[0], wa * a[1] + wb * b[1], wa * a[2] + wb * b[2] });
        }

        // Interior samples of the latitude-circle path from a to b (both on z = c) in direction dir.
        private static List<double[]> CircleInterior(double[] a, double[] b, double c, int dir)
        {
            var pts = new List<double[]>();
            double dphi = LongitudeTravel(a, b, dir);
            int m = (int)Math.Ceiling(Math.Abs(dphi) / CircleStep);
            double r = Math.Sqrt(Math.Max(0.0, 1 - c * c)), phi0 = Math.Atan2(a[1], a[0]);
            for (int k = 1; k < m; k++)
            {
                double phi = phi0 + dphi * k / m;
                pts.Add(new[] { r * Math.Cos(phi), r * Math.Sin(phi), c });
            }
            return pts;
        }

        private static double LongitudeTravel(double[] a, double[] b, int dir)
        {
            if (Same(a, b)) return 0.0;
            double d = Math.Atan2(b[1], b[0]) - Math.Atan2(a[1], a[0]);
            const double TwoPi = 2 * Math.PI;
            d -= TwoPi * Math.Floor(d / TwoPi); // [0, 2π)
            return dir > 0 ? d : d - TwoPi * (d > 0 ? 1 : 0);
        }

        // Signed area between the latitude arc a->b on z = c and the great-circle chord a->b
        private static double Lens(double[] a, double[] b, double c)
        {
            double dphi = Math.Atan2(a[0] * b[1] - a[1] * b[0], a[0] * b[0] + a[1] * b[1]);
            double cz = a[0] * b[1] - a[1] * b[0], ab = Dot(a, b);
            if (c >= 0) return dphi * (1 - c) - 2 * Math.Atan2(cz, 1 + 2 * c + ab);
            return -dphi * (1 + c) + 2 * Math.Atan2(cz, 1 - 2 * c + ab);
        }

        // Signed solid angle of the fan from s over a closed polygon
        private static double FanArea(List<Node> nodes, double[] s)
        {
            double total = 0;
            for (int i = 0; i < nodes.Count; i++) total += TriangleArea(s, nodes[i].P, nodes[(i + 1) % nodes.Count].P);
            return total;
        }

        // Signed solid angle of the geodesic triangle (a, b, c), positive when counter-clockwise from outside
        private static double TriangleArea(double[] a, double[] b, double[] c)
        {
            double triple = Dot(a, Cross(b, c));
            double denom = 1.0 + Dot(a, b) + Dot(b, c) + Dot(c, a);
            return 2.0 * Math.Atan2(triple, denom);
        }

        private static double[] Cross(double[] u, double[] v) => new[] {
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0]
        };

        private static double Dot(double[] u, double[] v) => u[0] * v[0] + u[1] * v[1] + u[2] * v[2];

        private static double[] Normalise(double[] v)
        {
            double n = Math.Sqrt(Dot(v, v));
            return n == 0 ? new double[] { 0, 0, 0 } : new[] { v[0] / n, v[1] / n, v[2] / n };
        }

        private static bool Same(double[] p, double[] q) =>
            Math.Abs(p[0] - q[0]) <= 1e-12 && Math.Abs(p[1] - q[1]) <= 1e-12 && Math.Abs(p[2] - q[2]) <= 1e-12;
    }
}
