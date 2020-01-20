using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Rhino.Geometry;

namespace RhinoGeometry {

    public class CCX {
        public double t0;
        public double t1;
        public Point3d p0;
        public Point3d p1;
    }

    public static class CurveUtil {


        public static List<Circle> SortCircles(List<Circle> x) {

            if (x.Count == 1)
                return x;

            List<Circle> sortedCircles = new List<Circle>() { x[0] };

            for (int i = 1; i < x.Count; i++) {

                Plane plane = PlaneUtil.ProjectPlaneXPlaneToYPlane(sortedCircles[sortedCircles.Count - 1].Plane, x[i].Plane);
                Circle circle = new Circle(plane, x[i].Radius);
                sortedCircles.Add(circle);
            }

            return sortedCircles;
        }

        public static Tuple<double, double> CurveCurveIntersection(Curve C0, Curve C1, double t) {

            //Rhino.RhinoApp.WriteLine("WTF");
            Rhino.Geometry.Intersect.CurveIntersections ci = Rhino.Geometry.Intersect.Intersection.CurveCurve(C0, C1, t, t);


            Line line = Line.Unset;
            if (ci.Count > 0) {
                Rhino.Geometry.Intersect.IntersectionEvent intersection = ci[0];
                //Rhino.RhinoApp.WriteLine(intersection.ParameterA.ToString() +" " + intersection.ParameterB.ToString());
                return new Tuple<double, double>(intersection.ParameterA, intersection.ParameterB);
            }

            return null;
        }

        public static Tuple<Point3d, Point3d> CurveCurveIntersectionPoints(Curve C0, Curve C1, double t) {

            Rhino.Geometry.Intersect.CurveIntersections ci = Rhino.Geometry.Intersect.Intersection.CurveCurve(C0, C1, t, t);


            Line line = Line.Unset;
            if (ci.Count > 0) {
                Rhino.Geometry.Intersect.IntersectionEvent intersection = ci[0];
                return new Tuple<Point3d, Point3d>(intersection.PointA, intersection.PointB);
            }

            return null;
        }

        public static List<Point3d> CurveCurveIntersectionAllPts(Curve C0, Curve C1, double t = 0.01) {

            Rhino.Geometry.Intersect.CurveIntersections ci = Rhino.Geometry.Intersect.Intersection.CurveCurve(C0, C1, t, t);

            List<Point3d> pts = new List<Point3d>();
            Line line = Line.Unset;
            if (ci.Count > 0) {

                foreach(Rhino.Geometry.Intersect.IntersectionEvent inter in ci) {
                    pts.Add(inter.PointA);
                }
                return pts;
           }

            return null;
        }

        public static CCX CurveCurveIntersectionCCX(Curve C0, Curve C1, double t) {

            Rhino.Geometry.Intersect.CurveIntersections ci = Rhino.Geometry.Intersect.Intersection.CurveCurve(C0, C1, t, t);


            Line line = Line.Unset;
            if (ci.Count > 0) {

                Rhino.Geometry.Intersect.IntersectionEvent intersection = ci[0];
                return new CCX() { t0 = intersection.ParameterA, t1 = intersection.ParameterB, p0 = intersection.PointA, p1 = intersection.PointB };
            }

            return null;
        }

        /// <summary>
        ///  Flip curve comparing the angular difference between n tangent on both curves
        /// </summary>
        /// <param name="guid"></param>
        /// <param name="crv"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        public static Curve GuideCurveDirection(Curve guid, Curve crv, int n = 10) {


            Curve x = guid.DuplicateCurve();
            Curve y = crv.DuplicateCurve();

            Curve r = crv.DuplicateCurve();
            r.Reverse();
            var px = x.DivideCurve(n);
            var py = y.DivideCurve(n);
            var pr = r.DivideCurve(n);
            //Rhino.RhinoApp.WriteLine(px.Length.ToString());


            double tot1 = 0;
            double tot2 = 0;


            for (int i = 0; i < n; i++) {
                var tx = x.CurveTangent(x.CurveClosestPoint(px[i]));
                var ty = y.CurveTangent(y.CurveClosestPoint(py[i]));
                var tr = r.CurveTangent(r.CurveClosestPoint(pr[i]));
                tot1 += Vector3d.VectorAngle(tx, ty);
                tot2 += Vector3d.VectorAngle(tx, tr);

            }


            if (tot1 > tot2)
                y.Reverse();
            return y;

        }

        public static Vector3d CurveTangent(this Curve curve, double parameter) {
            if (curve.Domain.IncludesParameter(parameter))
                return curve.TangentAt(parameter);
            return Vector3d.Unset;
        }

        public static double CurveClosestPoint(this Curve curve, Point3d p) {
            if (curve.ClosestPoint(p, out double t))
                return t;
            return -1;
        }

        public static Point3d[] DivideCurve(this Curve c, int n, bool ends = true) {
            Point3d[] p = new Point3d[0];
            c.DivideByCount(n, ends, out p);
            return p;
        }

        public static Line ExtendLine(Line l, double d) {

            Point3d p0 = l.From;
            Point3d p1 = l.To;

            return new Line(p0 - l.UnitTangent * d, p1 + l.UnitTangent * d);
        }

        public static Curve BooleanOpenCurve(Curve ClosedC, Curve C, out int result) {



            PointContainment cont0 = ClosedC.Contains(C.PointAtEnd, Plane.WorldXY, 0.01);
            PointContainment cont1 = ClosedC.Contains(C.PointAtStart, Plane.WorldXY, 0.01);


            bool isCurveInsideStart = (cont0 == PointContainment.Inside || cont0 == PointContainment.Coincident);
            bool isCurveInsideEnd = (cont1 == PointContainment.Inside || cont1 == PointContainment.Coincident);

            Rhino.Geometry.Intersect.CurveIntersections ci = Rhino.Geometry.Intersect.Intersection.CurveCurve(C, ClosedC, 0.01, 0.01);
            Interval interval = C.Domain;


            if (isCurveInsideStart && isCurveInsideEnd) {
                result = 2;
                return C;
            } else if (isCurveInsideStart && !isCurveInsideEnd) {
                result = 1;
                if (ci.Count == 0) {
                    result = 0;
                    return null;
                }
                return C.Trim(ci[0].ParameterA, interval.T0);
            } else if (!isCurveInsideStart && isCurveInsideEnd) {
                if (ci.Count == 0) {
                    result = 0;
                    return null;
                }
                result = 1;
                return C.Trim(interval.T1, ci[0].ParameterA);
            } else {
                if (ci.Count == 2) {

                    result = 1;
                    return C.Trim(ci[0].ParameterA, ci[1].ParameterA);
                }
            }
            result = 0;
            return null;
        }


    }
}
