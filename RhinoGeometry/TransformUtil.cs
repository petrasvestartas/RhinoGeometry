using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Rhino.Geometry;

namespace RhinoGeometry {
    public static class TransformUtil {


        /// <summary>
        /// Do you have current relative plane?
        /// Then next relative plane is the result of this method
        /// </summary>
        /// <param name="currJoulinPlane"></param>
        /// <param name="nextJoulinPlane"></param>
        /// <returns></returns>

        public static Transform RelativeTransform(Plane currPlane, Plane nextPlane) {
            Transform CurrPlanetoWorldXY = Transform.PlaneToPlane(currPlane, Plane.WorldXY); //Orient Current Joulin to ground
            Transform fromWorldXYToNextPlane = Transform.PlaneToPlane(Plane.WorldXY, nextPlane);   //Orient from ground to next Joulien
            return fromWorldXYToNextPlane * CurrPlanetoWorldXY;
        }
    }
}
