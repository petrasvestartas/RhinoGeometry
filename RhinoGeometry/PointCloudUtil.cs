using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace RhinoGeometry {
   public static class PointCloudUtil {

        /// <summary>
        /// Use
        ///private void RunScript(List<Point3d> pts, double D, ref object A, ref object B){
        ///Point3d VP = new Point3d(0, 0, 1000000);
        ///A = PointCloudNormals(pts, VP, D);
        ///B = pts.FindAll(VT => VT.DistanceTo(pts[10]) < D);}
        /// </summary>
        /// <param name="points"></param>
        /// <param name="VP"></param>
        /// <param name="D"></param>
        /// <returns></returns>
        public static Vector3d[] PointCloudNormalsByViewPoint(List<Point3d> points, Point3d VP, double D) {

            Vector3d[] Normals = new Vector3d[points.Count];
            Rhino.Collections.Point3dList pts = new Rhino.Collections.Point3dList(points);

            double Dev = 0.01;
            double squaredD = D * D;

            int i = 0;
            foreach (Point3d point in pts) {

                dynamic nei = pts.FindAll(V => V.DistanceToSquared(point) < squaredD);
                Plane NP = Plane.Unset;
                Plane.FitPlaneToPoints(nei, out NP, out Dev);

                int sign = (NP.Normal * (VP - point) > 0) ? 1 : -1;
                Normals[i++] = (sign * NP.Normal);

            }

            return Normals;
        }

    }
}
