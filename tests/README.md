# Tests

Checks for SphericalVoronoiLib, compiled from the library sources (without the SOFA reader).

```
dotnet run -c Release -- test data   # invariants, Monte Carlo cross-checks, clipping, duplicates, SONICOM masks
dotnet run -c Release -- stress      # degenerate and adversarial inputs, fuzzed bands and grid subsets
dotnet run -c Release -- bench       # median time over N = 250 ... 32,000
dotnet run -c Release -- dump data/sonicom793.csv out.json -0.7071 -   # cells for plotting (zMin, zMax; "-" for none)
```

Both `test` and `stress` exit non-zero if any check fails.
