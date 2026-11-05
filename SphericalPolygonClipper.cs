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
        /// If the clipped region has several parts, the part with the most vertices is returned.
        /// </summary>
        public static List<SphericalPoint> ClipToBand(List<SphericalPoint> poly, double? zMin, double? zMax)
        {
            var parts = ClipToBandParts(poly, zMin, zMax);
            List<SphericalPoint> best = new List<SphericalPoint>();
            foreach (var part in parts) if (part.Count > best.Count) best = part;
            return best;
        }

        /// <summary>
        /// Clips a polygon to the band and returns the boundary of every part of the clipped region.
        /// </summary>
        public static List<List<SphericalPoint>> ClipToBandParts(List<SphericalPoint> poly, double? zMin, double? zMax)
        {
            var result = new List<List<SphericalPoint>>();
            if (poly == null || poly.Count == 0) return result;

            foreach (var nodes in Clip(ToNodes(poly), zMin, zMax, out _, out _))
            {
                var output = new List<SphericalPoint>();
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
                result.Add(output);
            }
            return result;
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

            const double FourPi = 4 * Math.PI;
            double total = 0;
            var parts = Clip(nodes, zMin, zMax, out bool crossedMin, out bool crossedMax);
            foreach (var part in parts)
            {
                if (part.Count < 2) continue;
                double area = FanArea(part, s);
                for (int i = 0; i < part.Count; i++)
                {
                    if (!part[i].CircleZ.HasValue) continue;
                    var a = part[i].P;
                    var b = part[(i + 1) % part.Count].P;
                    double c = part[i].CircleZ.Value;
                    var pts = new List<double[]> { a };
                    pts.AddRange(CircleInterior(a, b, c, part[i].Dir));
                    pts.Add(b);
                    // Replace the chord a-b by the sampled arc, then add the exact lens between each sub-chord and the circle
                    area -= TriangleArea(s, a, b);
                    for (int k = 0; k + 1 < pts.Count; k++)
                    {
                        area += TriangleArea(s, pts[k], pts[k + 1]);
                        area += Lens(pts[k], pts[k + 1], c);
                    }
                }
                if (area < 0 && area > -1e-9) area = 0;
                total += area - FourPi * Math.Floor(area / FourPi);
            }

            // Latitude circle wholly inside the cell: remove the polar cap
            bool any = total > 0;
            if (any && zMax.HasValue && containsNorthPole && !crossedMax && zMax.Value < 1)
                total -= 2 * Math.PI * (1 - zMax.Value);
            if (any && zMin.HasValue && containsSouthPole && !crossedMin && zMin.Value > -1)
                total -= 2 * Math.PI * (1 + zMin.Value);

            return Math.Max(0.0, total);
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

        private static List<List<Node>> Clip(List<Node> nodes, double? zMin, double? zMax, out bool crossedMin, out bool crossedMax)
        {
            crossedMin = crossedMax = false;
            var parts = new List<List<Node>> { nodes };
            if (zMin.HasValue)
            {
                var next = new List<List<Node>>();
                foreach (var part in parts) next.AddRange(ClipHalf(part, zMin.Value, +1, ref crossedMin));
                parts = next;
            }
            if (zMax.HasValue)
            {
                var next = new List<List<Node>>();
                foreach (var part in parts) next.AddRange(ClipHalf(part, zMax.Value, -1, ref crossedMax));
                parts = next;
            }

            var result = new List<List<Node>>();
            foreach (var part in parts)
            {
                var cleaned = new List<Node>(part.Count);
                foreach (var n in part)
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
                if (cleaned.Count >= 2) result.Add(cleaned);
            }
            return result;
        }

        // Clips one loop against sign*(z - c) >= 0. The kept side of a latitude plane need not be convex,
        // so each exit is joined along the circle (kept region on the left) to the nearest entry, which may
        // produce several loops.
        private static List<List<Node>> ClipHalf(List<Node> poly, double c, int sign, ref bool crossed)
        {
            int n = poly.Count;
            var loops = new List<List<Node>>();
            if (n == 0) return loops;

            var inside = new bool[n];
            for (int i = 0; i < n; i++) inside[i] = sign * (poly[i].P[2] - c) >= 0;

            // Polygon with crossing points inserted: kind 0 vertex, +1 entry, -1 exit
            var seq = new List<Node>(n + 8);
            var kind = new List<int>(n + 8);
            for (int i = 0; i < n; i++)
            {
                int j = (i + 1) % n;
                seq.Add(poly[i]); kind.Add(0);
                if (inside[i] != inside[j])
                {
                    seq.Add(new Node(IntersectEdgeAtZ(poly[i].P, poly[j].P, c)));
                    kind.Add(inside[i] ? -1 : +1);
                }
            }

            var entries = new List<int>();
            for (int k = 0; k < seq.Count; k++) if (kind[k] == +1) entries.Add(k);
            if (entries.Count == 0)
            {
                if (inside[0]) loops.Add(poly);
                return loops;
            }
            crossed = true;

            int m = seq.Count;
            var visited = new HashSet<int>();
            foreach (int start in entries)
            {
                if (visited.Contains(start)) continue;
                var loop = new List<Node>();
                int cur = start;
                for (int guard = 0; guard <= entries.Count; guard++)
                {
                    visited.Add(cur);
                    loop.Add(new Node(seq[cur].P));
                    int k = (cur + 1) % m;
                    while (kind[k] != -1) { loop.Add(seq[k]); k = (k + 1) % m; }

                    var exit = seq[k].P;
                    int nextEntry = -1; double nearest = double.MaxValue;
                    foreach (int e in entries)
                    {
                        double t = Math.Abs(LongitudeTravel(exit, seq[e].P, sign));
                        if (t < nearest) { nearest = t; nextEntry = e; }
                    }
                    loop.Add(new Node(exit, c, sign));
                    if (nextEntry == start || nextEntry < 0 || visited.Contains(nextEntry)) break;
                    cur = nextEntry;
                }
                loops.Add(loop);
            }
            return loops;
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
