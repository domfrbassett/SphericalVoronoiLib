using HDF.PInvoke;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;

namespace SphericalVoronoiLib
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
                positionType = ReadStringAttribute(fileId, "SourcePosition", "Type")
                    ?? ReadStringAttribute(fileId, "/", "SourcePositionType")
                    ?? "spherical";

                long datasetId = H5D.open(fileId, "SourcePosition");
                if (datasetId < 0) throw new Exception("SOFA file has no SourcePosition dataset.");
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

        // String attribute (fixed or variable length), or null if absent
        private static string ReadStringAttribute(long fileId, string objectName, string attributeName)
        {
            if (H5A.exists_by_name(fileId, objectName, attributeName) <= 0) return null;
            long attrId = H5A.open_by_name(fileId, objectName, attributeName, H5P.DEFAULT, H5P.DEFAULT);
            if (attrId < 0) return null;
            long typeId = H5A.get_type(attrId);
            try
            {
                if (H5T.is_variable_str(typeId) > 0)
                {
                    var ptrs = new IntPtr[1];
                    GCHandle hp = GCHandle.Alloc(ptrs, GCHandleType.Pinned);
                    try { H5A.read(attrId, typeId, hp.AddrOfPinnedObject()); }
                    finally { hp.Free(); }
                    return ptrs[0] == IntPtr.Zero ? null : Marshal.PtrToStringAnsi(ptrs[0]);
                }
                byte[] buffer = new byte[H5T.get_size(typeId).ToInt32()];
                GCHandle h = GCHandle.Alloc(buffer, GCHandleType.Pinned);
                try { H5A.read(attrId, typeId, h.AddrOfPinnedObject()); }
                finally { h.Free(); }
                return Encoding.ASCII.GetString(buffer).TrimEnd('\0', ' ');
            }
            finally { H5T.close(typeId); H5A.close(attrId); }
        }
    }
}
