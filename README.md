# SphericalVoronoiLib

**SphericalVoronoiLib** is a C# class library for generating **spherical Voronoi tessellations** and **solid angles** from **SOFA-format** or general 3D coordinate data. It addresses an issue discussed by Jot, J.-M., Larcher, V., and Warusfel, O. (1995) in their paper: “Digital signal processing issues in the context of binaural and transaural stereophony”, Proc. Audio Eng. Soc. Conv. Preprint 3980.

<img width="804" height="195" alt="image" src="https://github.com/user-attachments/assets/53a1ce86-c588-455b-af3b-69eef2ab708e" />

* [54] M. Mein, Perception de l'information binaurale liée aux réflexions précoces dans une salle. Application à la simulation de la qualité acoustique, mémoire de DEA, Université du Maine, Le Mans, septembre 1993.

If anyone has access to the Mein (1993) paper, please get in touch. The university does not, and is equally interested. It would be good to find out if they actually proposed a weighting system.

<img src="https://github.com/user-attachments/assets/71445dbc-008a-4846-9bbf-251517b538d5" width="700" />

---

## Overview

* Reads **source positions** from a SOFA file, or takes Cartesian or azimuth/elevation arrays
* Converts them to **normalised 3D unit vectors**; repeated directions share one cell and results stay aligned with the input rows
* Computes **spherical Voronoi cells** using [MIConvexHull](https://github.com/DesignEngrLab/MIConvexHull), including regular grids with co-circular points
* Calculates **per-cell solid angles** (in steradians)
* Calculates **weights proportional to solid angles**
* Optionally **restricts the cells to an elevation band** (`zMin`/`zMax`, or degrees via `ComputeFromSpherical`), so grids with missing directions, e.g. below −45°, are not given disproportionately large cells. Clipped edges follow the latitude circle exactly
* Performance: the 793-direction SONICOM grid takes a few hundredths of a second; 32,000 directions take under a second.

---

## Usage

```csharp
using SphericalVoronoiLib;

// From a SOFA file, optionally restricted to a z band (z = sin elevation)
VoronoiResult r = SphericalVoronoiProcessor.ComputeFromSofa("P0001_FreeFieldCompMinPhase_48kHz.sofa", zMin: -0.7071);

// From azimuth/elevation in degrees (a third distance column is ignored)
VoronoiResult s = SphericalVoronoiProcessor.ComputeFromSpherical(azElDeg, minElevationDeg: -45);

// From Cartesian coordinates (any length; rows are normalised)
VoronoiResult c = SphericalVoronoiProcessor.ComputeFromArray(xyz);
```

Every array in `VoronoiResult` has one entry per input row:

| Member | Meaning |
|---|---|
| `Points` | Input directions as unit vectors |
| `SolidAngles` | Cell solid angle in sr (rows with the same direction split their cell) |
| `Weights` | `SolidAngles` normalised to sum to 1 over the sphere or band |
| `CellIndex`, `CellCount` | Cell of each row, and the number of distinct directions |
| `AreaClosureError` | \|Σ − 4π\| of the full tessellation, typically below 1e-12 sr |

## Requirements

The library targets .NET Framework 4.8. The tests in `tests/` need the .NET 8 SDK.

## Testing

`tests/` holds a console test suite that runs on every push (GitHub Actions): area closure to 4π, Monte Carlo cross-checks of individual cells, exact band clipping, duplicate handling, the SONICOM grid and its sparse masks, and adversarial inputs (polyhedra, points on one circle, jittered grids, clustered and hemispherical data, random rotations, fuzzed bands and grid subsets). `VoronoiResult.AreaClosureError` reports |Σ − 4π| for every tessellation.

```
cd tests
dotnet run -c Release -- test data
dotnet run -c Release -- stress
```

---

## Dependencies

* **HDF.PInvoke** – for reading SOFA (HDF5) files
* **MIConvexHull** – for convex hull computations

---

## Attribution

Algorithm based on **Bruno Luong’s** MATLAB function
[`voronoisphere.m`](https://uk.mathworks.com/matlabcentral/fileexchange/40989-voronoi-sphere/).
This repository contains an independent C# re-implementation.

---

## License

MIT License © 2025 Dom Bassett
