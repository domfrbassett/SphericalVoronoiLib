using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using SphericalVoronoiLib;

static partial class Program
{
    static int failures = 0;
    static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    static int Main(string[] args)
    {
        string cmd = args.Length > 0 ? args[0] : "test";
        switch (cmd)
        {
            case "test": RunTests(args.Length > 1 ? args[1] : "data"); break;
            case "bench": Bench(); break;
            case "stress": Stress(); break;
            case "dump": Dump(args[1], args[2], Opt(args, 3), Opt(args, 4)); break;
            default: Console.WriteLine("usage: test [dataDir] | stress | bench | dump in.csv out.json [zMin zMax]"); return 2;
        }
        return failures == 0 ? 0 : 1;
    }

    static double? Opt(string[] a, int i) =>
        a.Length > i && a[i] != "-" ? double.Parse(a[i], Inv) : (double?)null;

    // ---------- point sets ----------
    static double[,] Fibonacci(int n)
    {
        var p = new double[n, 3];
        double ga = Math.PI * (3 - Math.Sqrt(5));
        for (int i = 0; i < n; i++)
        {
            double z = 1 - 2 * (i + 0.5) / n, r = Math.Sqrt(1 - z * z), t = ga * i;
            p[i, 0] = r * Math.Cos(t); p[i, 1] = r * Math.Sin(t); p[i, 2] = z;
        }
        return p;
    }

    static double[,] Random(int n, int seed)
    {
        var rng = new System.Random(seed);
        var p = new double[n, 3];
        for (int i = 0; i < n; i++)
        {
            double x, y, z, s;
            do { x = Gauss(rng); y = Gauss(rng); z = Gauss(rng); s = Math.Sqrt(x * x + y * y + z * z); } while (s < 1e-9);
            p[i, 0] = x / s; p[i, 1] = y / s; p[i, 2] = z / s;
        }
        return p;
    }

    static double Gauss(System.Random r) =>
        Math.Sqrt(-2 * Math.Log(1 - r.NextDouble())) * Math.Cos(2 * Math.PI * r.NextDouble());

    static double[,] ReadCsv(string path)
    {
        var rows = File.ReadAllLines(path).Where(l => l.Trim().Length > 0)
            .Select(l => l.Split(',').Select(s => double.Parse(s, Inv)).ToArray()).ToList();
        var p = new double[rows.Count, 3];
        for (int i = 0; i < rows.Count; i++) for (int j = 0; j < 3; j++) p[i, j] = rows[i][j];
        return p;
    }

    // ---------- checks ----------
    static void Check(string name, bool ok, string detail)
    {
        if (!ok) failures++;
        Console.WriteLine($"  [{(ok ? "PASS" : "FAIL")}] {name}: {detail}");
    }

    // Monte Carlo cell areas: nearest seed among the returned points, optionally restricted to a band.
    static double[] MonteCarloAreas(double[,] pts, int samples, int seed, double? zMin, double? zMax)
    {
        int n = pts.GetLength(0);
        var counts = new long[n];
        var rng = new System.Random(seed);
        for (int s = 0; s < samples; s++)
        {
            double z = 2 * rng.NextDouble() - 1, t = 2 * Math.PI * rng.NextDouble(), r = Math.Sqrt(1 - z * z);
            double x = r * Math.Cos(t), y = r * Math.Sin(t);
            if (zMin.HasValue && z < zMin.Value) continue;
            if (zMax.HasValue && z > zMax.Value) continue;
            int best = 0; double bd = double.MinValue;
            for (int i = 0; i < n; i++)
            {
                double d = x * pts[i, 0] + y * pts[i, 1] + z * pts[i, 2];
                if (d > bd) { bd = d; best = i; }
            }
            counts[best]++;
        }
        return counts.Select(c => 4 * Math.PI * c / samples).ToArray();
    }

    static void CompareMonteCarlo(string label, double[,] input, double? zMin, double? zMax, int samples)
    {
        var res = SphericalVoronoiProcessor.ComputeFromArray(input, zMin, zMax);
        var mc = MonteCarloAreas(res.Points, samples, 7, zMin, zMax);
        // Per-cell z-score against binomial standard error
        double worstZ = 0; int worst = -1;
        for (int i = 0; i < mc.Length; i++)
        {
            double p = Math.Max(res.SolidAngles[i] / (4 * Math.PI), 1.0 / samples);
            double se = 4 * Math.PI * Math.Sqrt(p * (1 - p) / samples);
            double zsc = Math.Abs(res.SolidAngles[i] - mc[i]) / se;
            if (zsc > worstZ) { worstZ = zsc; worst = i; }
        }
        Check($"{label} per-cell area vs Monte Carlo ({samples / 1000}k samples)", worstZ < 5.0,
            $"worst |z| = {worstZ:F2} (cell {worst}: exact {res.SolidAngles[worst]:F5} sr, MC {mc[worst]:F5} sr)");
    }

    static void RunTests(string dataDir)
    {
        Console.WriteLine("== Full-sphere tessellations: cell areas sum to 4π");
        var sets = new List<(string, double[,])>
        {
            ("random N=4", Random(4, 1)), ("random N=10", Random(10, 2)), ("fibonacci N=200", Fibonacci(200)),
            ("random N=1000", Random(1000, 3)),
        };
        string sonicomPath = Path.Combine(dataDir, "sonicom793.csv");
        if (File.Exists(sonicomPath)) sets.Add(("SONICOM 793", ReadCsv(sonicomPath)));
        foreach (var (label, pts) in sets)
        {
            var r = SphericalVoronoiProcessor.ComputeFromArray(pts);
            double sum = r.SolidAngles.Sum();
            Check(label, Math.Abs(sum - 4 * Math.PI) < 1e-9, $"Σ = {sum:F12}, error {sum - 4 * Math.PI:E2}");
            Check(label + " no empty or negative cells", r.SolidAngles.All(a => a > 0), $"min = {r.SolidAngles.Min():E3} sr");
        }

        Console.WriteLine("\n== Monte Carlo cross-check of individual cells");
        CompareMonteCarlo("random N=60", Random(60, 11), null, null, 1_000_000);

        Console.WriteLine("\n== Band clipping: Σ must equal the band area 2π(zMax − zMin)");
        var clipCases = new List<(string, double[,], double?, double?)>
        {
            ("fibonacci N=200, z ≥ −sin45°", Fibonacci(200), -Math.Sin(Math.PI / 4), null),
            ("fibonacci N=200, z ≤ 0.5", Fibonacci(200), null, 0.5),
            ("random N=300, −0.3 ≤ z ≤ 0.6", Random(300, 5), -0.3, 0.6),
            ("random N=12, z ≥ 0.2 (coarse)", Random(12, 9), 0.2, null),
        };
        if (File.Exists(sonicomPath))
            clipCases.Add(("SONICOM 793, z ≥ −sin45°", ReadCsv(sonicomPath), -Math.Sin(Math.PI / 4), null));
        foreach (var (label, pts, lo, hi) in clipCases)
        {
            var r = SphericalVoronoiProcessor.ComputeFromArray(pts, lo, hi);
            double target = 2 * Math.PI * ((hi ?? 1.0) - (lo ?? -1.0));
            double sum = r.SolidAngles.Sum();
            Check(label, Math.Abs(sum - target) < 1e-9, $"Σ = {sum:F9}, band = {target:F9}, error {100 * (sum - target) / target:+0.0000;-0.0000}%");
        }
        CompareMonteCarlo("random N=40 clipped to −0.2 ≤ z ≤ 0.7", Random(40, 13), -0.2, 0.7, 1_000_000);

        Console.WriteLine("\n== Pole cell entirely containing the clipping circle");
        {
            // One seed at the north pole, the rest well below: the pole cell contains the whole circle z = 0.9.
            var pts = new double[,] { { 0, 0, 1 }, { 1, 0, 0 }, { -0.5, 0.866, 0 }, { -0.5, -0.866, 0 }, { 0, 0, -1 } };
            var r = SphericalVoronoiProcessor.ComputeFromArray(pts, null, 0.9);
            double target = 2 * Math.PI * (0.9 + 1);
            Check("cap removed from pole cell", Math.Abs(r.SolidAngles.Sum() - target) < 1e-9,
                $"Σ = {r.SolidAngles.Sum():F9}, band = {target:F9}");
        }

        Console.WriteLine("\n== Input handling");
        {
            // Row 2 duplicates row 0; row 3 is row 1 at a different radius (same direction).
            var basePts = Fibonacci(50);
            var pts = new double[53, 3];
            for (int i = 0; i < 50; i++) for (int j = 0; j < 3; j++) pts[i, j] = basePts[i, j];
            for (int j = 0; j < 3; j++) { pts[50, j] = basePts[0, j]; pts[51, j] = 1.5 * basePts[1, j]; pts[52, j] = basePts[2, j] * (1 + 1e-13); }
            try
            {
                var r = SphericalVoronoiProcessor.ComputeFromArray(pts);
                bool aligned = r.Weights.Length == 53;
                Check("one weight per input row", aligned, $"{r.Weights.Length} weights for 53 rows");
                Check("weights sum to 1", Math.Abs(r.Weights.Sum() - 1) < 1e-12, $"Σw = {r.Weights.Sum():F12}");
                if (aligned)
                {
                    var single = SphericalVoronoiProcessor.ComputeFromArray(basePts);
                    double dup0 = r.SolidAngles[0] + r.SolidAngles[50];
                    Check("duplicate rows share their cell", Math.Abs(dup0 - single.SolidAngles[0]) < 1e-12 && Math.Abs(r.SolidAngles[0] - r.SolidAngles[50]) < 1e-15,
                        $"{r.SolidAngles[0]:F6} + {r.SolidAngles[50]:F6} vs {single.SolidAngles[0]:F6}");
                }
            }
            catch (Exception e) { Check("duplicates handled", false, e.GetType().Name + ": " + e.Message); }
        }

        Console.WriteLine("\n== Degenerate inputs: 1, 2 and 3 seeds (full sphere and clipped)");
        foreach (int k in new[] { 1, 2, 3 })
        {
            var pts = Random(k, 40 + k);
            try
            {
                var full = SphericalVoronoiProcessor.ComputeFromArray(pts);
                var band = SphericalVoronoiProcessor.ComputeFromArray(pts, -0.4, 0.8);
                Check($"N={k} full sphere", Math.Abs(full.SolidAngles.Sum() - 4 * Math.PI) < 1e-9, $"Σ = {full.SolidAngles.Sum():F9}");
                Check($"N={k} band −0.4..0.8", Math.Abs(band.SolidAngles.Sum() - 2 * Math.PI * 1.2) < 1e-9, $"Σ = {band.SolidAngles.Sum():F9}, band = {2 * Math.PI * 1.2:F9}");
            }
            catch (Exception e) { Check($"N={k}", false, e.GetType().Name + ": " + e.Message); }
        }
        CompareMonteCarlo("random N=3", Random(3, 77), null, null, 400_000);
        CompareMonteCarlo("random N=3 clipped to z ≥ 0.1", Random(3, 77), 0.1, null, 400_000);

        Console.WriteLine("\n== Sparse SONICOM masks (as used in HRTF upsampling)");
        foreach (var k in new[] { 100, 19, 5, 3 })
        {
            string f = Path.Combine(dataDir, $"sonicom_mask{k}.csv");
            if (!File.Exists(f)) continue;
            try
            {
                var r = SphericalVoronoiProcessor.ComputeFromArray(ReadCsv(f));
                double sum = r.SolidAngles.Sum();
                Check($"N={k}", Math.Abs(sum - 4 * Math.PI) < 1e-9 && r.SolidAngles.All(a => a > 0),
                    $"Σ = {sum:F9}, weights {string.Join(" ", r.Weights.Take(5).Select(w => w.ToString("F3", Inv)))}{(k > 5 ? " …" : "")}");
            }
            catch (Exception e) { Check($"N={k}", false, e.GetType().Name + ": " + e.Message); }
        }

        Console.WriteLine("\n== Spherical (azimuth/elevation) input");
        {
            var xyz = Fibonacci(300);
            var ae = new double[300, 3];
            for (int i = 0; i < 300; i++)
            {
                ae[i, 0] = Math.Atan2(xyz[i, 1], xyz[i, 0]) * 180 / Math.PI;
                ae[i, 1] = Math.Asin(xyz[i, 2]) * 180 / Math.PI;
                ae[i, 2] = 1.5;
            }
            var a = SphericalVoronoiProcessor.ComputeFromSpherical(ae, -45, 80);
            var b = SphericalVoronoiProcessor.ComputeFromArray(xyz, Math.Sin(-Math.PI / 4), Math.Sin(80 * Math.PI / 180));
            double maxDiff = a.SolidAngles.Zip(b.SolidAngles, (u, v) => Math.Abs(u - v)).Max();
            Check("matches Cartesian input with elevation band −45°..80°", maxDiff < 1e-10, $"max |Δ| = {maxDiff:E2} sr");
        }

        Console.WriteLine($"\n{(failures == 0 ? "All checks passed" : failures + " check(s) failed")}.");
    }

    static void Bench()
    {
        Console.WriteLine("N,seconds");
        for (int w = 0; w < 5; w++) SphericalVoronoiProcessor.ComputeFromArray(Fibonacci(500)); // JIT warm-up
        foreach (int n in new[] { 250, 500, 1000, 2000, 4000, 8000, 16000, 32000 })
        {
            var pts = Fibonacci(n);
            var times = new List<double>();
            for (int rep = 0; rep < 5; rep++)
            {
                var sw = Stopwatch.StartNew();
                SphericalVoronoiProcessor.ComputeFromArray(pts);
                times.Add(sw.Elapsed.TotalSeconds);
                if (sw.Elapsed.TotalSeconds > 20) break;
            }
            times.Sort();
            double med = times[times.Count / 2];
            Console.WriteLine($"{n},{med.ToString("F4", Inv)}");
            if (med > 20) break;
        }
    }

    // JSON dump for plotting: seeds, cell boundaries (full-sphere) and solid angles (clipped if requested).
    static void Dump(string csv, string outPath, double? zMin, double? zMax)
    {
        var input = ReadCsv(csv);
        var res = SphericalVoronoiProcessor.ComputeFromArray(input, zMin, zMax);
        int n = res.Points.GetLength(0);
        var pts = Enumerable.Range(0, n).Select(i => new SphericalPoint(res.Points[i, 0], res.Points[i, 1], res.Points[i, 2])).ToList();
        var vs = VoronoiSphere.Compute(pts);
        var sb = new StringBuilder();
        sb.Append("{\"points\":[");
        sb.Append(string.Join(",", Enumerable.Range(0, n).Select(i => $"[{F(res.Points[i, 0])},{F(res.Points[i, 1])},{F(res.Points[i, 2])}]")));
        sb.Append("],\"solid\":[").Append(string.Join(",", res.SolidAngles.Select(F)));
        sb.Append("],\"cells\":[");
        sb.Append(string.Join(",", vs.VoronoiBoundary.Select(b => "[" + string.Join(",", b.Select(p => $"[{F(p.X)},{F(p.Y)},{F(p.Z)}]")) + "]")));
        sb.Append("]");
        if (zMin.HasValue || zMax.HasValue)
        {
            sb.Append(",\"clipped\":[");
            sb.Append(string.Join(",", vs.VoronoiBoundary.Select(b => "[" + string.Join(",",
                SphericalPolygonClipper.ClipToBand(b, zMin, zMax).Select(p => $"[{F(p.X)},{F(p.Y)},{F(p.Z)}]")) + "]")));
            sb.Append("]");
        }
        sb.Append("}");
        File.WriteAllText(outPath, sb.ToString());
        Console.WriteLine($"wrote {n} cells to {outPath}");
    }

    static string F(double v) => v.ToString("R", Inv);
}
