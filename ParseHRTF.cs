using HDF.PInvoke;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;

namespace SphericalVoronoi
{
    public class SOFAReader
    {
        public static double[,] ReadSourcePosition(string sofaPath, out string positionType)
        {
            long fileId = H5F.open(sofaPath, H5F.ACC_RDONLY);
            if (fileId < 0) throw new Exception("Failed to open SOFA file.");

            positionType = "spherical";
            double[,] cartesianPositions = null;

            try
            {
                long attrId = H5A.open_by_name(fileId, "/", "SourcePositionType", H5P.DEFAULT, H5P.DEFAULT);
                if (attrId >= 0)
                {
                    long typeId = H5A.get_type(attrId);
                    long spaceId = H5A.get_space(attrId);

                    byte[] buffer = new byte[H5T.get_size(typeId).ToInt32()];
                    GCHandle handle = GCHandle.Alloc(buffer, GCHandleType.Pinned);
                    try
                    {
                        H5A.read(attrId, typeId, handle.AddrOfPinnedObject());
                        positionType = Encoding.ASCII.GetString(buffer).TrimEnd('\0');
                    }
                    finally { handle.Free(); H5T.close(typeId); H5S.close(spaceId); H5A.close(attrId); }
                }

                long datasetId = H5D.open(fileId, "SourcePosition");
                long dataspaceId = H5D.get_space(datasetId);

                ulong[] dims = new ulong[2];
                H5S.get_simple_extent_dims(dataspaceId, dims, null);
                int nSources = (int)dims[0];

                double[] flat = new double[nSources * 3];
                GCHandle h = GCHandle.Alloc(flat, GCHandleType.Pinned);
                try
                {
                    H5D.read(datasetId, H5T.NATIVE_DOUBLE, H5S.ALL, H5S.ALL, H5P.DEFAULT, h.AddrOfPinnedObject());
                }
                finally { h.Free(); H5S.close(dataspaceId); H5D.close(datasetId); }

                cartesianPositions = new double[nSources, 3];

                for (int i = 0; i < nSources; i++)
                {
                    double az = flat[i * 3 + 0];
                    double el = flat[i * 3 + 1];
                    double r = 1.0;

                    if (positionType.ToLower().Contains("spherical"))
                    {
                        double azRad = az * Math.PI / 180.0;
                        double elRad = el * Math.PI / 180.0;

                        cartesianPositions[i, 0] = r * Math.Cos(elRad) * Math.Cos(azRad);
                        cartesianPositions[i, 1] = r * Math.Cos(elRad) * Math.Sin(azRad);
                        cartesianPositions[i, 2] = r * Math.Sin(elRad);
                    }
                    else
                    {
                        cartesianPositions[i, 0] = az;
                        cartesianPositions[i, 1] = el;
                        cartesianPositions[i, 2] = flat[i * 3 + 2];
                    }
                }

                return cartesianPositions;
            }
            finally
            {
                H5F.close(fileId);
            }
        }
    }
}
