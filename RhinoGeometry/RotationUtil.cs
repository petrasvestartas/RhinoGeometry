using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace RhinoGeometry {
    public static class RotationUtil {

        /// <summary>
        /// Takes 7 value double and turns into plane
        /// First 3 values is plane origin, the rest quaternion ABCD
        /// </summary>
        /// <param name="RhinoPosQuat"></param>
        /// <returns></returns>
        public static Plane QuaternionToRhinoPlane(double[] RhinoPosQuat) {

            Point3d p = new Point3d(RhinoPosQuat[0], RhinoPosQuat[1], RhinoPosQuat[2]);
            Quaternion q = new Quaternion(RhinoPosQuat[3], RhinoPosQuat[4], RhinoPosQuat[5], RhinoPosQuat[6]);

            Plane plane;
            q.GetRotation(out plane);
            plane.Origin = p;

            return plane;
        }

        public static double[] PlaneToPosQuaternion(Plane refPlane, Plane p) {

            Rhino.Geometry.Quaternion quaternion = new Quaternion();
            quaternion.SetRotation(refPlane, p);

            double[] transformation = new double[]{
      p.OriginX,
      p.OriginY,
      p.OriginZ,
      quaternion.A,quaternion.B,quaternion.C,quaternion.D
      };
            return transformation;
        }

    }
}
