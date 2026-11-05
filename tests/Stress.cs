using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using SphericalVoronoiLib;

static partial class Program
{
    static double[,] FromList(IList<double[]> pts)
    {
        var a = new double[pts.Count, 3];
        for (int i = 0; i < pts.Count; i++) for (int j = 0; j < 3; j++) a[i, j] = pts[i][j];
        return a;
    }

    static List<double[]> ToList(double[,] a) =>
        Enumerable.Range(0, a.GetLength(0)).Select(i => new[] { a[i, 0], a[i, 1], a[i, 2] }).ToList();

    static double[,] RandomRotation(int seed)
    {
        var r = new System.Random(seed);
        double q0 = Gauss(r), q1 = Gauss(r), q2 = Gauss(r), q3 = Gauss(r), n = Math.Sqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
        q0 /= n; q1 /= n; q2 /= n; q3 /= n;
        return new double[,] {
            { 1 - 2 * (q2 * q2 + q3 * q3), 2 * (q1 * q2 - q0 * q3), 2 * (q1 * q3 + q0 * q2) },
            { 2 * (q1 * q2 + q0 * q3), 1 - 2 * (q1 * q1 + q3 * q3), 2 * (q2 * q3 - q0 * q1) },
            { 2 * (q1 * q3 - q0 * q2), 2 * (q2 * q3 + q0 * q1), 1 - 2 * (q1 * q1 + q2 * q2) } };
    }

    static double[,] Rotate(double[,] p, double[,] R)
    {
        int n = p.GetLength(0); var o = new double[n, 3];
        for (int i = 0; i < n; i++) for (int j = 0; j < 3; j++) o[i, j] = R[j, 0] * p[i, 0] + R[j, 1] * p[i, 1] + R[j, 2] * p[i, 2];
        return o;
    }

    static double[,] LatLong(double stepDeg, double jitterRad, int seed, bool repeatPoles)
    {
        var r = new System.Random(seed); var pts = new List<double[]>();
        for (double el = -90; el <= 90 + 1e-9; el += stepDeg)
        {
            bool pole = Math.Abs(Math.Abs(el) - 90) < 1e-9;
            for (double az = 0; az < 360 - 1e-9; az += stepDeg)
            {
                double e = el * Math.PI / 180 + jitterRad * (2 * r.NextDouble() - 1), a = az * Math.PI / 180 + jitterRad * (2 * r.NextDouble() - 1);
                pts.Add(new[] { Math.Cos(e) * Math.Cos(a), Math.Cos(e) * Math.Sin(a), Math.Sin(e) });
                if (pole && !repeatPoles) break;
            }
        }
        return FromList(pts);
    }

    static double[,] Polyhedron(string name)
    {
        double g = (1 + Math.Sqrt(5)) / 2;
        var v = new List<double[]>();
        switch (name)
        {
            case "tetrahedron": v.AddRange(new[] { new[] { 1.0, 1, 1 }, new[] { 1.0, -1, -1 }, new[] { -1.0, 1, -1 }, new[] { -1.0, -1, 1 } }); break;
            case "octahedron": for (int k = 0; k < 3; k++) for (int s = -1; s <= 1; s += 2) { var p = new double[3]; p[k] = s; v.Add(p); } break;
            case "cube": for (int a = -1; a <= 1; a += 2) for (int b = -1; b <= 1; b += 2) for (int c = -1; c <= 1; c += 2) v.Add(new double[] { a, b, c }); break;
            case "icosahedron":
                foreach (var s1 in new[] { -1, 1 }) foreach (var s2 in new[] { -1, 1 })
                { v.Add(new[] { 0, s1, s2 * g }); v.Add(new[] { s1, s2 * g, 0 }); v.Add(new[] { s2 * g, 0, s1 }); }
                break;
            case "dodecahedron":
                for (int a = -1; a <= 1; a += 2) for (int b = -1; b <= 1; b += 2) for (int c = -1; c <= 1; c += 2) v.Add(new double[] { a, b, c });
                foreach (var s1 in new[] { -1, 1 }) foreach (var s2 in new[] { -1, 1 })
                { v.Add(new[] { 0, s1 / g, s2 * g }); v.Add(new[] { s1 / g, s2 * g, 0 }); v.Add(new[] { s2 * g, 0, s1 / g }); }
                break;
        }
        return FromList(v);
    }

    // Full-sphere case: exception-free, finite, positive, sums to 4π, optional Monte Carlo per cell.
    static VoronoiResult Full(string label, double[,] pts, int mcSamples = 0, double sumTol = 1e-9)
    {
        try
        {
            var sw = Stopwatch.StartNew();
            var r = SphericalVoronoiProcessor.ComputeFromArray(pts);
            double ms = sw.Elapsed.TotalMilliseconds;
            bool finite = r.SolidAngles.All(a => !double.IsNaN(a) && !double.IsInfinity(a));
            double sum = r.SolidAngles.Sum();
            bool ok = finite && r.SolidAngles.All(a => a > 0) && Math.Abs(sum - 4 * Math.PI) < sumTol;
            Check(label, ok, $"{r.CellCount} cells, Σ−4π = {sum - 4 * Math.PI:E1}, min {r.SolidAngles.Min():E2} sr, {ms:F0} ms");
            if (ok && mcSamples > 0) McCheck(label, r, null, null, mcSamples);
            return r;
        }
        catch (Exception e) { Check(label, false, e.GetType().Name + ": " + e.Message); return null; }
    }

    static void McCheck(string label, VoronoiResult r, double? lo, double? hi, int samples)
    {
        // Distinct seeds only, so duplicate rows are compared through their cell
        var first = new Dictionary<int, int>();
        for (int i = 0; i < r.CellIndex.Length; i++) if (!first.ContainsKey(r.CellIndex[i])) first[r.CellIndex[i]] = i;
        int m = r.CellCount; var seeds = new double[m, 3]; var area = new double[m];
        foreach (var kv in first) { for (int j = 0; j < 3; j++) seeds[kv.Key, j] = r.Points[kv.Value, j]; }
        for (int i = 0; i < r.CellIndex.Length; i++) area[r.CellIndex[i]] += r.SolidAngles[i];
        var mc = MonteCarloAreas(seeds, samples, 3, lo, hi);
        double worst = 0;
        for (int k = 0; k < m; k++)
        {
            double p = Math.Max(area[k] / (4 * Math.PI), 1.0 / samples);
            worst = Math.Max(worst, Math.Abs(area[k] - mc[k]) / (4 * Math.PI * Math.Sqrt(p * (1 - p) / samples)));
        }
        Check(label + " vs Monte Carlo", worst < 5.5, $"worst |z| = {worst:F2}");
    }

    // Every Voronoi vertex of cell k must be at least as close to seed k as to any other seed.
    static double WorstVertexViolation(double[,] pts)
    {
        var r = SphericalVoronoiProcessor.ComputeFromArray(pts);
        int m = r.CellCount;
        var first = new int[m]; for (int i = r.CellIndex.Length - 1; i >= 0; i--) first[r.CellIndex[i]] = i;
        var seeds = Enumerable.Range(0, m).Select(k => new SphericalPoint(r.Points[first[k], 0], r.Points[first[k], 1], r.Points[first[k], 2])).ToList();
        var v = VoronoiSphere.Compute(seeds);
        double worst = 0;
        for (int k = 0; k < m; k++)
        {
            foreach (var b in v.VoronoiBoundary[k])
            {
                double own = b.X * seeds[k].X + b.Y * seeds[k].Y + b.Z * seeds[k].Z, best = own;
                foreach (var q in seeds) best = Math.Max(best, b.X * q.X + b.Y * q.Y + b.Z * q.Z);
                worst = Math.Max(worst, best - own);
            }
        }
        return worst;
    }

    static void Stress()
    {
        Console.WriteLine("== Voronoi validity: boundary points are never closer to another seed");
        foreach (var (label, pts) in new[] {
            ("random N=500", Random(500, 4)), ("10° grid", LatLong(10, 0, 1, false)), ("SONICOM-like 5° grid", LatLong(5, 0, 1, false)),
            ("dodecahedron", Polyhedron("dodecahedron")), ("10° grid jittered 1e-10", LatLong(10, 1e-10, 7, false)) })
        {
            try { double w = WorstVertexViolation(pts); Check(label, w < 1e-10, $"worst excess closeness {w:E1}"); }
            catch (Exception e) { Check(label, false, e.GetType().Name + ": " + e.Message); }
        }

        Console.WriteLine("\n== Random subsets of regular grids (co-circular degeneracies)");
        {
            var rr = new System.Random(77); int trials = 200, bad = 0; double worstClose = 0, worstViol = 0; string note = "";
            for (int t = 0; t < trials; t++)
            {
                var grid = ToList(LatLong(new[] { 5.0, 10, 15, 30 }[rr.Next(4)], 0, 1, false));
                int k = Math.Min(grid.Count, 4 + rr.Next(400));
                var sub = grid.OrderBy(_ => rr.Next()).Take(k).ToList();
                var pts = FromList(sub);
                if (t % 4 == 0) pts = Rotate(pts, RandomRotation(t));
                try
                {
                    var r = SphericalVoronoiProcessor.ComputeFromArray(pts);
                    worstClose = Math.Max(worstClose, r.AreaClosureError);
                    bool ok = r.AreaClosureError < 1e-9 && r.SolidAngles.All(a => a > 0);
                    if (t % 10 == 0) { double vv = WorstVertexViolation(pts); worstViol = Math.Max(worstViol, vv); ok &= vv < 1e-10; }
                    if (!ok) { bad++; note = $"trial {t}, N={k}"; }
                }
                catch (Exception e) { bad++; note = $"trial {t}, N={k}: {e.GetType().Name}: {e.Message}"; }
            }
            Check($"{trials} subsets", bad == 0, $"{bad} failures{(note.Length > 0 ? " (" + note + ")" : "")}; worst closure {worstClose:E1} sr, worst vertex violation {worstViol:E1}");
        }

        Console.WriteLine("\n== Polyhedra (equal cells, many co-circular points)");
        foreach (var (name, n) in new[] { ("tetrahedron", 4), ("octahedron", 6), ("cube", 8), ("icosahedron", 12), ("dodecahedron", 20) })
        {
            var r = Full(name, Polyhedron(name));
            if (r != null)
            {
                double dev = r.SolidAngles.Max(a => Math.Abs(a - 4 * Math.PI / n));
                Check(name + " cells equal 4π/N", dev < 1e-9, $"max deviation {dev:E1}");
            }
        }

        Console.WriteLine("\n== Latitude/longitude grids, exact and jittered");
        Full("10° grid", LatLong(10, 0, 1, false), 400_000);
        Full("10° grid, pole rows repeated per azimuth", LatLong(10, 0, 1, true));
        Full("5° grid", LatLong(5, 0, 1, false));
        foreach (var j in new[] { 1e-13, 1e-10, 1e-7, 1e-4 })
            Full($"10° grid jittered by {j:E0} rad", LatLong(10, j, 7, false), j >= 1e-4 ? 400_000 : 0);

        Console.WriteLine("\n== Points on one circle (flat hull)");
        {
            var eq = new List<double[]>(); var rr = new System.Random(3);
            var ang = Enumerable.Range(0, 9).Select(_ => 2 * Math.PI * rr.NextDouble()).ToList();
            foreach (var t in ang) eq.Add(new[] { Math.Cos(t), Math.Sin(t), 0.0 });
            Full("9 random points on the equator", FromList(eq), 400_000);
            var ring = ang.Select(t => new[] { 0.8 * Math.Cos(t), 0.8 * Math.Sin(t), 0.6 }).ToList();
            Full("9 random points on the z = 0.6 small circle", FromList(ring), 400_000);
            Full("same small circle, rotated", Rotate(FromList(ring), RandomRotation(5)), 400_000);
            Full("4 points: two antipodal pairs on a great circle", FromList(new List<double[]> { new[] { 1.0, 0, 0 }, new[] { -1.0, 0, 0 }, new[] { 0, 1.0, 0 }, new[] { 0, -1.0, 0 } }));
        }

        Console.WriteLine("\n== Partial and uneven coverage");
        {
            var rr = new System.Random(11); var cap = new List<double[]>();
            while (cap.Count < 80) { var p = ToList(Random(1, rr.Next()))[0]; if (p[2] > 0.3) cap.Add(p); }
            Full("80 points in the cap z > 0.3 only", FromList(cap), 1_000_000);
            var hemi = ToList(Random(400, 21)).Where(p => p[2] > 0).ToList();
            Full("upper hemisphere only", FromList(hemi), 400_000);
            var cl = new List<double[]>();
            for (int i = 0; i < 200; i++) { double t = 0.017 * rr.NextDouble(), ph = 2 * Math.PI * rr.NextDouble(); cl.Add(new[] { Math.Sin(t) * Math.Cos(ph), Math.Sin(t) * Math.Sin(ph), Math.Cos(t) }); }
            cl.AddRange(ToList(Random(5, 99)));
            Full("200 points within 1° of the pole + 5 elsewhere", FromList(cl));
            var nd = ToList(Fibonacci(100)); nd.Add(new[] { nd[0][0] + 1e-8, nd[0][1], nd[0][2] });
            Full("near-duplicate 1e-8 apart (above merge tolerance)", FromList(nd));
        }

        Console.WriteLine("\n== Rotation invariance (areas per input row)");
        foreach (var (label, pts) in new[] { ("random N=300", Random(300, 8)), ("10° grid", LatLong(10, 0, 1, false)), ("icosahedron", Polyhedron("icosahedron")) })
        {
            try
            {
                var a = SphericalVoronoiProcessor.ComputeFromArray(pts).SolidAngles;
                double worst = 0;
                for (int s = 0; s < 10; s++)
                {
                    var b = SphericalVoronoiProcessor.ComputeFromArray(Rotate(pts, RandomRotation(100 + s))).SolidAngles;
                    worst = Math.Max(worst, a.Zip(b, (u, v) => Math.Abs(u - v)).Max());
                }
                Check(label + ", 10 random rotations", worst < 1e-9, $"max |Δ| = {worst:E1} sr");
            }
            catch (Exception e) { Check(label, false, e.GetType().Name + ": " + e.Message); }
        }

        Console.WriteLine("\n== Band clipping fuzz");
        {
            var rr = new System.Random(2024); int trials = 300, bad = 0; double worst = 0; string worstCase = "";
            for (int t = 0; t < trials; t++)
            {
                int n = 4 + rr.Next(300);
                var pts = t % 3 == 0 ? Fibonacci(n) : t % 3 == 1 ? Random(n, rr.Next()) : LatLong(new[] { 5.0, 10, 15, 30 }[rr.Next(4)], 0, 1, false);
                if (t % 7 == 0) pts = Rotate(pts, RandomRotation(t));
                double a = 2 * rr.NextDouble() - 1, b = 2 * rr.NextDouble() - 1;
                double? lo = Math.Min(a, b), hi = Math.Max(a, b);
                if (t % 5 == 1) lo = null; if (t % 5 == 2) hi = null;
                if (t % 11 == 3) { lo = pts[0, 2]; if (hi < lo) hi = null; } // band edge through a seed
                try
                {
                    var r = SphericalVoronoiProcessor.ComputeFromArray(pts, lo, hi);
                    double target = 2 * Math.PI * ((hi ?? 1) - (lo ?? -1)), err = Math.Abs(r.SolidAngles.Sum() - target);
                    bool fin = r.SolidAngles.All(x => x >= 0 && !double.IsNaN(x));
                    if (!fin || err > 1e-9) bad++;
                    if (err > worst) { worst = err; worstCase = $"trial {t}, N={n}, band [{lo?.ToString("F3") ?? "−1"}, {hi?.ToString("F3") ?? "1"}]"; }
                }
                catch (Exception e) { bad++; worstCase = $"trial {t}: {e.GetType().Name}: {e.Message}"; }
            }
            Check($"{trials} random point sets and bands", bad == 0, $"{bad} failures; worst |Σ − band| = {worst:E1} ({worstCase})");
            var g = LatLong(10, 0, 1, false);
            McCheck("10° grid clipped to −0.35 ≤ z ≤ 0.85", SphericalVoronoiProcessor.ComputeFromArray(g, -0.35, 0.85), -0.35, 0.85, 1_000_000);
            var cap = ToList(Random(400, 21)).Where(p => p[2] > 0).ToList();
            McCheck("upper hemisphere clipped to z ≥ 0.4", SphericalVoronoiProcessor.ComputeFromArray(FromList(cap), 0.4, null), 0.4, null, 1_000_000);
        }

        Console.WriteLine("\n== Band edge cases");
        {
            var p = Fibonacci(100);
            void Band(string label, double? lo, double? hi, double expect)
            {
                try
                {
                    var r = SphericalVoronoiProcessor.ComputeFromArray(p, lo, hi);
                    Check(label, Math.Abs(r.SolidAngles.Sum() - expect) < 1e-9 && r.SolidAngles.All(x => x >= 0), $"Σ = {r.SolidAngles.Sum():F9}, expected {expect:F9}");
                }
                catch (Exception e) { Check(label, false, e.GetType().Name + ": " + e.Message); }
            }
            Band("zMin = zMax = 0.3 (empty band)", 0.3, 0.3, 0);
            Band("zMin = −1, zMax = 1 (whole sphere)", -1, 1, 4 * Math.PI);
            Band("zMin = 1 (single point)", 1, null, 0);
            Band("zMax = −1 (single point)", null, -1, 0);
            Band("zMin = −2, zMax = 2 (outside the sphere)", -2, 2, 4 * Math.PI);
        }

        Console.WriteLine("\n== Invalid input");
        void Throws(string label, double[,] pts, double? lo = null, double? hi = null)
        {
            try { SphericalVoronoiProcessor.ComputeFromArray(pts, lo, hi); Check(label, false, "no exception"); }
            catch (ArgumentException e) { Check(label, true, e.Message); }
            catch (Exception e) { Check(label, false, "unexpected " + e.GetType().Name + ": " + e.Message); }
        }
        Throws("NaN row", new double[,] { { 1, 0, 0 }, { double.NaN, 0, 1 }, { 0, 1, 0 }, { 0, 0, 1 } });
        Throws("zero row", new double[,] { { 1, 0, 0 }, { 0, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } });
        Throws("two columns", new double[,] { { 1, 0 }, { 0, 1 } });
        Throws("NaN band limit", Fibonacci(10), double.NaN, null);

        Console.WriteLine("\n== Scale");
        foreach (int n in new[] { 100_000 })
        {
            var sw = Stopwatch.StartNew();
            Full($"random N={n:N0}", Random(n, 1), 0, 1e-8);
            Console.WriteLine($"  ({sw.Elapsed.TotalSeconds:F2} s)");
        }

        Console.WriteLine($"\n{(failures == 0 ? "All stress checks passed" : failures + " stress check(s) failed")}.");
    }
}
