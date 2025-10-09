using MIConvexHull;
using System;
using System.Collections.Generic;
using System.Linq;

namespace SphericalVoronoi
{
    public class SphericalPoint : IVertex
    {
        public double X, Y, Z;
        public double[] Position { get; set; }
        public SphericalPoint(double x, double y, double z)
        {
            double norm = Math.Sqrt(x * x + y * y + z * z);
            if (norm == 0) throw new ArgumentException("Zero length vector");
            X = x / norm; Y = y / norm; Z = z / norm;
            Position = new double[] { X, Y, Z };
        }
        public SphericalPoint(double[] p)
        {
            if (p.Length != 3) throw new ArgumentException("Require length 3");
            double norm = Math.Sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
            if (norm == 0) throw new ArgumentException("Zero length vector");
            X = p[0] / norm; Y = p[1] / norm; Z = p[2] / norm;
            Position = new double[] { X, Y, Z };
        }
        public double[] ToArray() => new double[] { X, Y, Z };
    }
    public class Face : ConvexFace<SphericalPoint, Face>
    {

    }
    public class VoronoiSphereResult
    {
        public List<SphericalPoint> Vertices = new List<SphericalPoint>();
        public List<List<int>> K = new List<List<int>>();
        public List<List<SphericalPoint>> VoronoiBoundary = new List<List<SphericalPoint>>();
        public double[] SolidAngles = new double[0];
    }
    public static class VoronoiSphere
    {
        public static VoronoiSphereResult Compute(List<SphericalPoint> inputPoints, double resolution = 2 * Math.PI / 180.0)
        {
            var result = new VoronoiSphereResult();

            if (inputPoints == null || inputPoints.Count == 0)
            {
                return result;
            }

            int npnts = inputPoints.Count;

            double[,] xyz = new double[3, npnts];
            for (int i = 0; i < npnts; i++)
            {
                var p = inputPoints[i].ToArray();
                double nrm = Math.Sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
                if (nrm <= 0) throw new ArgumentException($"Input point {i} has zero length");

                xyz[0, i] = p[0] / nrm;
                xyz[1, i] = p[1] / nrm;
                xyz[2, i] = p[2] / nrm;
            }

            switch (npnts)
            {
                case 1:
                    result.Vertices = new List<SphericalPoint>();
                    result.K.Add(new List<int>());
                    result.VoronoiBoundary.Add(FullCircle(GetVectorColumn(xyz, 0), resolution));
                    result.SolidAngles = new double[] { 4 * Math.PI };
                    return result;

                case 2:
                    result.Vertices = new List<SphericalPoint>();
                    result.K = new List<List<int>> { new List<int>(), new List<int>() };
                    var diff = Subtract(GetVectorColumn(xyz, 1), GetVectorColumn(xyz, 0));
                    var fc = FullCircle(diff, resolution);
                    result.VoronoiBoundary.Add(fc);

                    var fc2 = new List<SphericalPoint>(fc);
                    fc2.Reverse();
                    result.VoronoiBoundary.Add(fc2);
                    result.SolidAngles = new double[] { 2 * Math.PI, 2 * Math.PI };
                    return result;

                case 3:
                    {
                        var T = new int[,] { { 0, 1, 2 } };
                        var centers = Center(xyz, T);

                        var M = new double[,] { { -1, 1, 0 }, { 0, -1, 1 }, { 1, 0, -1 } };
                        double[,] dxyz = Multiply(xyz, M);

                        var edgeArcs = new List<List<SphericalPoint>>(3);
                        for (int k = 0; k < 3; k++)
                        {
                            edgeArcs.Add(HalfCircle(centers[k], GetColumn(dxyz, k), resolution));
                        }

                        result.VoronoiBoundary = new List<List<SphericalPoint>>
                        {
                            new List<SphericalPoint>(),
                            new List<SphericalPoint>(),
                            new List<SphericalPoint>()
                        };

                        result.K = new List<List<int>>
                        {
                            new List<int> { 0, 1 },
                            new List<int> { 1, 2 },
                            new List<int> { 2, 0 }
                        };

                        for (int k = 0; k < 3; k++)
                        {
                            var arc1 = edgeArcs[k];
                            var arc2 = edgeArcs[(k + 1) % 3];

                            var combined = new List<SphericalPoint>(arc1);
                            for (int i = arc2.Count - 2; i >= 1; i--)
                                combined.Add(arc2[i]);

                            result.VoronoiBoundary[k] = combined;
                        }

                        double[] s3 = new double[3];
                        for (int k = 0; k < 3; k++)
                        {
                            var dk = GetColumn(dxyz, k);
                            var dkp1 = GetColumn(dxyz, (k + 1) % 3);

                            var (A, B) = HalfCircleAB(dk, centers[0]);
                            double num = -Dot(dkp1, B);
                            double den = Dot(dkp1, A);
                            s3[k] = 2 * Math.Atan2(num, den);
                        }
                        result.SolidAngles = s3;

                        var verts = new List<SphericalPoint>(centers.Count * 2);
                        foreach (var c in centers)
                        {
                            verts.Add(c);
                            verts.Add(new SphericalPoint(-c.X, -c.Y, -c.Z));
                        }
                        result.Vertices = verts;

                        return result;
                    }

                default:
                    {
                        var tvertices = inputPoints.ToList();

                        var hull = ConvexHull.Create<SphericalPoint, Face>(tvertices);
                        var faces = hull.Result.Faces.ToList();

                        int nt = faces.Count;
                        var Tlist = new List<int[]>();
                        foreach (var f in faces)
                        {
                            var idx = f.Vertices.Select(v => tvertices.IndexOf(v)).ToArray();
                            if (idx.Length != 3) throw new Exception("Non-triangular face encountered");
                            Tlist.Add(new int[] { idx[0], idx[1], idx[2] });
                        }

                        var Tmat = new int[nt, 3];
                        for (int i = 0; i < nt; i++)
                            for (int j = 0; j < 3; j++)
                                Tmat[i, j] = Tlist[i][j];

                        var edgeList = new List<int[]>();
                        for (int i = 0; i < nt; i++)
                        {
                            edgeList.Add(SortPair(Tmat[i, 0], Tmat[i, 1]));
                            edgeList.Add(SortPair(Tmat[i, 1], Tmat[i, 2]));
                            edgeList.Add(SortPair(Tmat[i, 2], Tmat[i, 0]));
                        }

                        var edgeDict = new Dictionary<(int, int), List<int>>();
                        for (int i = 0; i < edgeList.Count; i++)
                        {
                            var key = (edgeList[i][0], edgeList[i][1]);
                            if (!edgeDict.ContainsKey(key)) edgeDict[key] = new List<int>();
                            edgeDict[key].Add(i);
                        }

                        foreach (var kv in edgeDict)
                        {
                            if (kv.Value.Count != 2)
                                throw new Exception("Non-manifold edge detected. Check point configuration or numerical precision.");
                        }

                        int ne = edgeDict.Count;
                        var allIds = Enumerable.Range(0, nt).SelectMany(i => new int[] { i, i, i }).ToArray();
                        var vid = new int[ne, 2];
                        var cellOfEdge = new int[ne, 2];
                        var edgeKeys = edgeDict.Keys.ToList();

                        for (int ei = 0; ei < edgeKeys.Count; ei++)
                        {
                            var rows = edgeDict[edgeKeys[ei]];
                            int tri1 = allIds[rows[0]];
                            int tri2 = allIds[rows[1]];
                            vid[ei, 0] = tri1;
                            vid[ei, 1] = tri2;

                            var seedPair = edgeList[rows[0]];
                            cellOfEdge[ei, 0] = seedPair[0];
                            cellOfEdge[ei, 1] = seedPair[1];
                        }

                        var edgeOfCell = new List<int>[npnts];
                        for (int i = 0; i < npnts; i++) edgeOfCell[i] = new List<int>();
                        for (int e = 0; e < ne; e++)
                        {
                            edgeOfCell[cellOfEdge[e, 0]].Add(e);
                            edgeOfCell[cellOfEdge[e, 1]].Add(e);
                        }

                        var Centers = Center(xyz, Tmat);
                        result.Vertices = Centers;

                        var edgeArcs = new List<List<SphericalPoint>>(ne);
                        for (int e = 0; e < ne; e++)
                        {
                            var v0 = Centers[vid[e, 0]];
                            var v1 = Centers[vid[e, 1]];
                            edgeArcs.Add(Arc(v0, v1, resolution));
                        }

                        result.VoronoiBoundary = new List<List<SphericalPoint>>(npnts);
                        result.K = new List<List<int>>(npnts);

                        for (int k = 0; k < npnts; k++)
                        {
                            var edges = edgeOfCell[k].ToArray();
                            if (edges.Length == 0)
                            {
                                result.VoronoiBoundary.Add(new List<SphericalPoint>());
                                result.K.Add(new List<int>());
                                continue;
                            }

                            var cycle = CyclingEdge(edges, vid);
                            cycle = OrientedEdge(cycle, result.Vertices, GetVectorColumn(xyz, k));

                            var vs = new List<(int, int)>();
                            for (int ii = 0; ii < ne; ii++)
                            {
                                var pair = (Math.Min(vid[ii, 0], vid[ii, 1]), Math.Max(vid[ii, 0], vid[ii, 1]));
                                vs.Add(pair);
                            }

                            var Klist = new List<int>();
                            var Xlist = new List<List<SphericalPoint>>();
                            for (int row = 0; row < cycle.GetLength(0); row++)
                            {
                                int a = cycle[row, 0];
                                int b = cycle[row, 1];
                                var key = (Math.Min(a, b), Math.Max(a, b));
                                int idx = vs.FindIndex(z => z == key);
                                if (idx < 0) throw new Exception("Edge lookup failed");
                                var arc = edgeArcs[idx];

                                bool flip = (cycle[row, 0] != vid[idx, 0]);
                                var arcCopy = new List<SphericalPoint>(arc);
                                if (flip) arcCopy.Reverse();
                                if (arcCopy.Count > 0) arcCopy.RemoveAt(arcCopy.Count - 1);
                                Xlist.Add(arcCopy);
                                Klist.Add(cycle[row, 0]);
                            }

                            var boundary = new List<SphericalPoint>();
                            foreach (var arr in Xlist) boundary.AddRange(arr);

                            result.VoronoiBoundary.Add(boundary);
                            result.K.Add(Klist);
                        }

                        var ctrs = Enumerable.Range(0, npnts)
                            .Select(i => new SphericalPoint(xyz[0, i], xyz[1, i], xyz[2, i]))
                            .ToList();
                        result.SolidAngles = VoronoiSolidAngle.Compute(result.Vertices, result.K, ctrs);

                        return result;
                    }
            }
        }
        #region Helper geometry functions
        private static double[] GetVectorColumn(double[,] M, int col)
        {
            return new double[] { M[0, col], M[1, col], M[2, col] };
        }
        private static double[] GetColumn(double[,] M, int col)
        {
            return new double[] { M[0, col], M[1, col], M[2, col] };
        }
        private static double[,] Multiply(double[,] A, double[,] B)
        {
            int r = A.GetLength(0), m = A.GetLength(1), c = B.GetLength(1);
            var R = new double[r, c];
            for (int i = 0; i < r; i++)
                for (int j = 0; j < c; j++)
                {
                    double s = 0;
                    for (int k = 0; k < m; k++) s += A[i, k] * B[k, j];
                    R[i, j] = s;
                }
            return R;
        }
        private static int[] SortPair(int a, int b)
        {
            if (a < b) return new int[] { a, b };
            return new int[] { b, a };
        }
        private static List<SphericalPoint> FullCircle(double[] Pole, double resolution)
        {
            double[] pole = Normalize(Pole);
            double[] A = OrthonormalBasisVector(pole);
            double[] B = Cross(A, pole);
            double theta = 2 * Math.PI;
            int npnts = Math.Max((int)Math.Ceiling(theta / resolution), 2);
            var G = new List<SphericalPoint>(npnts);
            for (int i = 0; i < npnts; i++)
            {
                double t = (double)i / (npnts - 1) * theta;
                var p = Add(MultiplyScalar(A, Math.Cos(t)), MultiplyScalar(B, Math.Sin(t)));
                G.Add(new SphericalPoint(p));
            }
            return G;
        }
        private static List<SphericalPoint> HalfCircle(SphericalPoint A, double[] Pole, double resolution)
        {
            var pole = Normalize(Pole);
            var a = A.ToArray(); a = Normalize(a);
            var B = Cross(a, pole);
            double theta = Math.PI;
            int npnts = Math.Max((int)Math.Ceiling(theta / resolution), 2);
            var G = new List<SphericalPoint>(npnts);
            for (int i = 0; i < npnts; i++)
            {
                double t = (double)i / (npnts - 1) * theta;
                var p = Add(MultiplyScalar(a, Math.Cos(t)), MultiplyScalar(B, Math.Sin(t)));
                G.Add(new SphericalPoint(p));
            }
            return G;
        }
        private static (double[] A, double[] B) HalfCircleAB(double[] Pole, SphericalPoint center)
        {
            var pole = Normalize(Pole);
            var A = center.ToArray(); A = Normalize(A);
            var B = Cross(A, pole);
            B = Normalize(B);
            A = Normalize(A);
            return (A, B);
        }
        private static List<SphericalPoint> Arc(SphericalPoint A, SphericalPoint B, double resolution)
        {
            var a = A.ToArray(); var b = B.ToArray();
            double dotProd = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
            dotProd = Math.Max(-1.0, Math.Min(1.0, dotProd));
            double theta = Math.Acos(dotProd);
            int nPoints = Math.Max((int)Math.Ceiling(theta / resolution), 2);
            var points = new List<SphericalPoint>(nPoints);
            for (int i = 0; i < nPoints; i++)
            {
                double t = (double)i / (nPoints - 1);
                double s1 = Math.Sin((1 - t) * theta);
                double s2 = Math.Sin(t * theta);
                double denom = Math.Sin(theta);
                var part1 = MultiplyScalar(a, s1 / denom);
                var part2 = MultiplyScalar(b, s2 / denom);
                var pt = Add(part1, part2);
                pt = Normalize(pt);
                points.Add(new SphericalPoint(pt));
            }
            return points;
        }
        private static List<SphericalPoint> Arc(SphericalPoint A, SphericalPoint B)
        {
            return Arc(A, B, 2 * Math.PI / 180.0);
        }
        private static double[] Cross(double[] u, double[] v) => new double[] {
            u[1]*v[2] - u[2]*v[1],
            u[2]*v[0] - u[0]*v[2],
            u[0]*v[1] - u[1]*v[0]
        };
        private static double Dot(double[] u, double[] v) => u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
        private static double[] MultiplyScalar(double[] v, double s) => new double[] { v[0] * s, v[1] * s, v[2] * s };
        private static double[] Add(double[] u, double[] v) => new double[] { u[0] + v[0], u[1] + v[1], u[2] + v[2] };
        private static double[] Subtract(double[] u, double[] v) => new double[] { u[0] - v[0], u[1] - v[1], u[2] - v[2] };
        private static double[] Normalize(double[] v)
        {
            double n = Math.Sqrt(Dot(v, v));
            if (n == 0) return new double[] { 0, 0, 0 };
            return new double[] { v[0] / n, v[1] / n, v[2] / n };
        }
        private static List<SphericalPoint> Center(double[,] xyz, int[,] T)
        {
            int nt = T.GetLength(0);
            var centers = new List<SphericalPoint>(nt);
            for (int i = 0; i < nt; i++)
            {
                double[] A = GetColumn(xyz, T[i, 0]);
                double[] B = GetColumn(xyz, T[i, 1]);
                double[] C = GetColumn(xyz, T[i, 2]);
                var AminusC = Subtract(A, C);
                var BminusC = Subtract(B, C);
                var A2B = MultiplyScalar(BminusC, Dot(AminusC, AminusC));
                var B2A = MultiplyScalar(AminusC, Dot(BminusC, BminusC));
                var AxB = Cross(AminusC, BminusC);
                var P = Cross(Subtract(A2B, B2A), AxB);
                double denom = 2 * Dot(AxB, AxB);
                if (Math.Abs(denom) < 1e-15)
                {
                    var avg = Normalize(Add(Add(A, B), C));
                    centers.Add(new SphericalPoint(avg));
                    continue;
                }
                P = MultiplyScalar(P, 1.0 / denom);
                var Pc = Add(C, P);
                var nPc = Normalize(Pc);
                double s = Dot(AxB, C);
                if (s < 0) nPc = MultiplyScalar(nPc, -1.0);
                centers.Add(new SphericalPoint(nPc));
            }
            return centers;
        }
        private static List<SphericalPoint> Center(double[,] xyz, List<int[]> Tlist)
        {
            var T = new int[Tlist.Count, 3];
            for (int i = 0; i < Tlist.Count; i++) for (int j = 0; j < 3; j++) T[i, j] = Tlist[i][j];
            return Center(xyz, T);
        }
        private static int[,] CyclingEdge(int[] edges, int[,] vertexes)
        {
            int n = edges.Length;
            var u = new int[2, n];
            for (int i = 0; i < n; i++)
            {
                u[0, i] = vertexes[edges[i], 0];
                u[1, i] = vertexes[edges[i], 1];
            }

            var flat = new List<int>();
            for (int i = 0; i < 2; i++) for (int j = 0; j < n; j++) flat.Add(u[i, j]);
            var unique = flat.Distinct().ToList();
            var idxMap = flat.Select(x => unique.IndexOf(x) + 1).ToArray();
            var I = new int[2, n];
            for (int i = 0; i < 2; i++) for (int j = 0; j < n; j++) I[i, j] = idxMap[i * n + j];

            var counts = new Dictionary<int, int>();
            foreach (var e in flat)
            {
                if (!counts.ContainsKey(e)) counts[e] = 0;
                counts[e]++;
            }
            foreach (var kv in counts) if (kv.Value != 2) throw new Exception("Topology issue in cycling_edge");

            var K = new List<int[]>();
            for (int j = 1; j <= unique.Count; j++)
            {
                var collected = new List<int>();
                for (int pos = 0; pos < flat.Count; pos++) if (flat[pos] == unique[j - 1]) collected.Add((pos % n) + 1);
                K.Add(collected.ToArray());
            }

            var outv = new int[n, 2];
            int p = 0; int q = 1;
            for (int j = 0; j < n; j++)
            {
                var iList = K[q];
                int iVal = iList[0];
                if (iVal == p + 1) p = iList[1] - 1; else p = iVal - 1;
                var ii = new int[] { I[0, p] - 1, I[1, p] - 1 };
                if (ii[0] == q)
                {
                    outv[j, 0] = u[0, p];
                    outv[j, 1] = u[1, p];
                    q = ii[1];
                }
                else
                {
                    outv[j, 0] = u[1, p];
                    outv[j, 1] = u[0, p];
                    q = ii[0];
                }
            }
            return outv;
        }
        private static int[,] OrientedEdge(int[,] v, List<SphericalPoint> P, double[] xyz)
        {
            int n = v.GetLength(0);
            var Q = NullSpace(xyz);

            var E = new double[3, n + 1];
            for (int i = 0; i < n; i++)
            {
                var pt = P[v[i, 0]].ToArray();
                E[0, i] = pt[0]; E[1, i] = pt[1]; E[2, i] = pt[2];
            }
            var first = P[v[0, 0]].ToArray();
            E[0, n] = first[0]; E[1, n] = first[1]; E[2, n] = first[2];

            var xy = new double[2, n + 1];
            for (int i = 0; i < n + 1; i++)
            {
                xy[0, i] = Dot(Q[0], new double[] { E[0, i], E[1, i], E[2, i] });
                xy[1, i] = Dot(Q[1], new double[] { E[0, i], E[1, i], E[2, i] });
            }

            double a = 0;
            for (int i = 0; i < n; i++) a += (xy[0, i] - xy[0, i + 1]) * (xy[1, i] + xy[1, i + 1]);

            double det = Determinant3x3(new double[,] { { xyz[0], Q[0][0], Q[1][0] }, { xyz[1], Q[0][1], Q[1][1] }, { xyz[2], Q[0][2], Q[1][2] } });

            if ((a < 0) ^ (det < 0))
            {
                v = Rotate180(v);
            }
            return v;
        }
        private static int[,] Rotate180(int[,] v)
        {
            int n = v.GetLength(0);
            var outv = new int[n, 2];
            for (int i = 0; i < n; i++) outv[i, 0] = v[n - 1 - i, 1];
            for (int i = 0; i < n; i++) outv[i, 1] = v[n - 1 - i, 0];
            return outv;
        }
        private static double[][] NullSpace(double[] v)
        {
            var n = Normalize(v);
            double[] other = Math.Abs(n[0]) < 0.9 ? new double[] { 1, 0, 0 } : new double[] { 0, 1, 0 };
            var u = Cross(n, other);
            u = Normalize(u);
            var w = Cross(u, n);
            w = Normalize(w);
            return new double[][] { u, w };
        }
        private static double Determinant3x3(double[,] M)
        {
            double a = M[0, 0], b = M[0, 1], c = M[0, 2];
            double d = M[1, 0], e = M[1, 1], f = M[1, 2];
            double g = M[2, 0], h = M[2, 1], i = M[2, 2];
            return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
        }
        private static double[] OrthonormalBasisVector(double[] pole)
        {
            var p = Normalize(pole);
            double[] a = Math.Abs(p[0]) < 0.9 ? new double[] { 1, 0, 0 } : new double[] { 0, 1, 0 };
            var A = Cross(p, a);
            A = Normalize(A);
            return A;
        }
        #endregion
    }
}