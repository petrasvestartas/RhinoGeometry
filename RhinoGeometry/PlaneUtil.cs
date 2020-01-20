using System;
using System.Collections.Generic;
using System.Linq;
using Rhino.Geometry;

namespace RhinoGeometry {
    public static class PlaneUtil {

        public static Rhino.Geometry.Plane ProjectPlaneXPlaneToYPlane(Rhino.Geometry.Plane x, Rhino.Geometry.Plane y) {

            Line l0 = new Line(x.Origin, x.Origin + x.ZAxis);
            Line l1 = new Line(x.Origin + x.XAxis, x.Origin + x.XAxis + x.ZAxis);
            Line l2 = new Line(x.Origin + x.YAxis, x.Origin + x.YAxis + x.ZAxis);

            double t0, t1, t2;
            Rhino.Geometry.Intersect.Intersection.LinePlane(l0, y, out t0);
            Rhino.Geometry.Intersect.Intersection.LinePlane(l1, y, out t1);
            Rhino.Geometry.Intersect.Intersection.LinePlane(l2, y, out t2);

            Point3d p0 = l0.PointAt(t0);
            Point3d p1 = l1.PointAt(t1);
            Point3d p2 = l2.PointAt(t2);
            Vector3d xaxis = p1 - p0;
            Vector3d yaxis = p2 - p0;

            Plane projectedPlane = new Plane(y.Origin, xaxis, yaxis);


            return projectedPlane;
        }

        public static Plane Switch(this Plane plane,string flag = "XY" ) {

            switch (flag) {

                case ("XZ"):
                return new Plane(plane.Origin, plane.XAxis, plane.ZAxis);

                case ("YX"):
                return new Plane(plane.Origin, plane.YAxis, plane.XAxis);

                case ("YZ"):
                return new Plane(plane.Origin, plane.YAxis, plane.ZAxis);

                case ("ZX"):
                return new Plane(plane.Origin, plane.ZAxis, plane.XAxis);

                case ("ZY"):
                return new Plane(plane.Origin, plane.ZAxis, plane.YAxis);

                default:
                return new Plane(plane);

            }

        }

        public static Plane FlipAndRotate(this Plane p) {
            Plane plane = new Plane(p);
            plane.Rotate(Math.PI * 0.5, plane.ZAxis);
            return plane;
        }

        public static List<Plane> Transform(this List<Plane> planes, Transform t) {
            List<Plane> planesTransformed = new List<Plane>();
            foreach (Plane p in planes) {
                Plane pTransformed = new Plane(p);
                pTransformed.Transform(t);
                planesTransformed.Add(pTransformed);
            }
            return planesTransformed;
        }

        /// <summary>
        /// Creates a plane from a line and a plane
        /// </summary>
        /// <param name="l0"></param>
        /// <param name="plane"></param>
        /// <param name="perpendicularToEdge"></param>
        /// <returns></returns>
        public static Plane PlaneFromLinePlane(this Line l0, Plane plane, int perpendicularToEdge = 0) {//2

            Line l1 = new Line(l0.From, l0.To);
            Point3d origin = l0.PointAt(0.5);
            l1.Transform(Rhino.Geometry.Transform.Rotation(Math.PI * 0.5, plane.ZAxis, origin));

            switch (perpendicularToEdge) {
                case (0):
                return new Plane(origin, -l0.Direction, plane.ZAxis);

                case (1):
                return new Plane(origin, -l1.Direction, plane.ZAxis);
            }

            return new Plane(origin, -l1.Direction, l0.Direction);
        }

        public static Curve[] BrepPlane(this Brep brep, Plane plane) {
            bool flag = Rhino.Geometry.Intersect.Intersection.BrepPlane(brep, plane, 0.01, out Curve[] crvs, out Point3d[] pts);
            if (flag)
                return crvs;
            return null;
        }

        public static Polyline[] MeshPlane(this Mesh mesh, Plane plane) {
            Polyline[] polylines = Rhino.Geometry.Intersect.Intersection.MeshPlane(mesh, plane);
            return polylines;
        }

        public static Point3d[] LinePlane(Line[] line, Plane plane) {

            Point3d[] p = new Point3d[line.Length];

            for (int i = 0; i < line.Length; i++)
                p[i] = LinePlane(line[i], plane);

            return p;
        }

        public static Point3d RayPlane(this Plane plane, Point3d p, Vector3d vec) {
            double t;
            Line line = new Line(p, p + vec);
            Rhino.Geometry.Intersect.Intersection.LinePlane(line, plane, out t);
            return line.PointAt(t);
        }

        public static Point3d LinePlane(Line line, Plane plane) {
            double t;
            Rhino.Geometry.Intersect.Intersection.LinePlane(line, plane, out t);
            return line.PointAt(t);
        }


        public static Point3d PlanePlanePlane(Plane p0, Plane p1, Plane p2) {
            Point3d p;
            Rhino.Geometry.Intersect.Intersection.PlanePlanePlane(p0, p1, p2, out p);
            return p;
        }

        public static Line PlanePlane(Plane p0, Plane p1) {
            Line line;
            Rhino.Geometry.Intersect.Intersection.PlanePlane(p0, p1, out line);
            return line;
        }

        public static Vector3d PlanePlaneVec(Plane p0, Plane p1) {


            Line line;
            Rhino.Geometry.Intersect.Intersection.PlanePlane(p0, p1, out line);

            return line.UnitTangent;
        }

        public static Plane InterpolatePlanes(List<Plane> P, bool g, double t) {

            List<Point3d> future = new List<Point3d>();
            List<Point3d> Warmth = new List<Point3d>();
            List<Point3d> Comfort = new List<Point3d>();
            List<Point3d> curiosity = new List<Point3d>();
            double love = 1.0;
            double fear = 1.0;
            int dreams = 0;
            int reality = 1;
            int chaos = 2;
            List<Plane> Play = P;
            double courage = t;

            foreach (Plane bagel in Play) {
                Quaternion spaceship = Quaternion.Rotation(Plane.WorldXY, bagel);
                double trouble = love / (fear - spaceship.D);
                future.Add(new Point3d(spaceship.A * trouble, spaceship.B * trouble, spaceship.C * trouble));
                trouble = fear / (love + spaceship.D);
                Warmth.Add(new Point3d(-spaceship.A * trouble, -spaceship.B * trouble, -spaceship.C * trouble));
                Comfort.Add(bagel.Origin);
            }
            if (g) { curiosity.Add(future[dreams]); } else { curiosity.Add(Warmth[dreams]); }
            for (var imagination = reality; imagination < future.Count; imagination++) {
                if (curiosity[imagination - reality].DistanceTo(future[imagination]) < curiosity[imagination - reality].DistanceTo(Warmth[imagination])) { curiosity.Add(future[imagination]); } else { curiosity.Add(Warmth[imagination]); }
            }
            var confusion = chaos + reality;
            Curve violently = Curve.CreateInterpolatedCurve(curiosity, confusion, (CurveKnotStyle)chaos);
            Curve gently = Curve.CreateInterpolatedCurve(Comfort, confusion, (CurveKnotStyle)chaos);
            List<double> always = new List<double>();
            List<double> never = new List<double>();
            double salt = new double();
            for (int unknowable = dreams; unknowable < curiosity.Count; unknowable++) {
                violently.ClosestPoint(curiosity[unknowable], out salt);
                always.Add(salt);
                gently.ClosestPoint(Comfort[unknowable], out salt);
                never.Add(salt);
            }
            Curve[] wonder = new Curve[curiosity.Count - reality];
            wonder = violently.Split(always);
            Curve[] excitement = new Curve[curiosity.Count - reality];
            excitement = gently.Split(never);
            var aspiration = courage * (Play.Count - fear);
            var limits = Math.Min(((int)Math.Floor(aspiration)), Play.Count - chaos);
            var infinity = (aspiration - limits);
            var knowledge = wonder[limits].PointAtNormalizedLength(infinity);
            var understanding = excitement[limits].PointAtNormalizedLength(infinity);
            var danger = knowledge.X;
            var joy = knowledge.Y;
            var hope = knowledge.Z;
            var forest = (reality + love) * danger / (love + danger * danger + joy * joy + hope * hope);
            var undergrowth = (reality + fear) * joy / (love + danger * danger + joy * joy + hope * hope);
            var canopy = chaos * hope / (love + danger * danger + joy * joy + hope * hope);
            var sky = (-love + danger * danger + joy * joy + hope * hope) / (fear + danger * danger + joy * joy + hope * hope);
            Quaternion magic = new Quaternion(forest, undergrowth, canopy, sky);
            Plane SomewhereInTheDistance = new Plane();
            magic.GetRotation(out SomewhereInTheDistance);
            SomewhereInTheDistance.Origin = understanding;

            return SomewhereInTheDistance;
        }


        public static Plane GetPlane(this Polyline polyline, bool AveragePlane = true) {

            //In case use default version

            if (!AveragePlane) {
                // in z case z axis may flip from time to time
                Plane plane_;
                Plane.FitPlaneToPoints(polyline, out plane_);
                plane_.Origin = polyline.CenterPoint();
                return plane_;
            } else {

                return new Plane(polyline.CenterPoint(), polyline.AverageNormal());

            }


            int a = (polyline.IsClosed) ? 1 : 0;
            int n = polyline.Count - a;

            Plane[] planes = new Plane[n];
            Point3d center = polyline.CenterPoint();
            Plane plane = new Plane();
            Plane.FitPlaneToPoints(polyline, out plane);
            plane.Origin = center;

            int ne = 0;
            int po = 0;
            for (int i = 0; i < n; i++) {
                planes[i] = new Plane(polyline[i], polyline[i + 1], center);

                if (planes[i].ZAxis.Z < 0)
                    ne++;
                else
                    po++;
            }

            if (po < ne && (plane.ZAxis.Z > 0))
                plane.Flip();

            return plane;
        }


        public static Plane ProjectLines(Line Projection, ref Line L1, ref Line L2) {
            Plane plane = new Plane(Projection.From, Projection.To - Projection.From);
            L1.To = (plane.ClosestPoint(L1.To));
            L2.To = (plane.ClosestPoint(L2.To));
            return plane;
        }

        //[Obsolete]
        public static Plane AveragePlane(IEnumerable<Plane> planes) {
            Point3d origin = Point3d.Origin;
            Vector3d XAxis = new Vector3d();
            Vector3d YAxis = new Vector3d();

            foreach (var p in planes) {
                origin += p.Origin;
                XAxis += p.XAxis;
                YAxis += p.YAxis;
            }

            int n = planes.Count();

            origin /= n;
            XAxis /= n;
            YAxis /= n;

            return new Plane(origin, XAxis, YAxis);
        }

        public static Plane AveragePlaneOrigin(IEnumerable<Plane> planes) {
            Point3d origin = Point3d.Origin;

            foreach (var p in planes) {
                origin += p.Origin;
            }

            int n = planes.Count();

            origin /= n;

            Plane plane = new Plane(planes.First());
            plane.Origin = origin;

            return plane;
        }




        public static Plane AlignPlane(Plane A, Plane B) {
            Point3d origin = B.Origin;
            Vector3d xAxis = B.XAxis;
            Vector3d yAxis = B.YAxis;
            Vector3d zAxis = B.ZAxis;
            switch (A.ZAxis.IsParallelTo(B.ZAxis, 1.539380400259)) {
                case -1: {
                        Point3d point3d = B.ClosestPoint(origin + A.XAxis);
                        xAxis = origin - point3d;
                        xAxis.Unitize();
                        break;
                    }
                case 0: {
                        Transform transform = Rhino.Geometry.Transform.Rotation(A.ZAxis, B.ZAxis, B.Origin);
                        xAxis = A.XAxis;
                        xAxis.Transform(transform);
                        break;
                    }
                case 1: {
                        xAxis = B.ClosestPoint(origin + A.XAxis) - origin;
                        xAxis.Unitize();
                        break;
                    }
            }
            return new Plane(origin, xAxis, Vector3d.CrossProduct(zAxis, xAxis));
        }

        public static bool AlignPlanes(List<Plane> planes, Plane master) {
            bool flag;
            if (planes == null) {
                flag = false;
            } else if (planes.Count == 0) {
                flag = false;
            } else if (master.IsValid) {
                int count = checked(planes.Count - 1);
                for (int i = 0; i <= count; i = checked(i + 1)) {
                    if (planes[i].IsValid) {
                        planes[i] = AlignPlane(master, planes[i]);
                        master = planes[i];
                    }
                }
                flag = true;
            } else {
                flag = false;
            }
            return flag;
        }


        public static Plane PlaneFromLinePlane(Line line, Plane plane) {
            double angle = Math.PI * 0.5;
            Vector3d vecX = line.UnitTangent;
            Vector3d vecY = line.UnitTangent;
            vecY.Rotate(angle, plane.Normal);

            Plane newPlane = new Plane(line.PointAt(0.5), vecX, vecY);
            newPlane.Rotate(Math.PI * 0.5, vecY);

            return newPlane;
        }

        public static Plane MovePlanebyAxis(this Plane plane, double dist, int axis = 2) {
            Plane p = new Plane(plane);

            switch (axis) {
                case (0):
                p.Translate(p.XAxis * dist);
                break;
                case (1):
                p.Translate(p.YAxis * dist);
                break;
                default:
                p.Translate(p.Normal * dist);
                break;
            }
            return p;
        }


        public static Plane MovePlanebyAxis(this Plane plane, double dist, Line line, int axis = 2, bool closer = true) {

            Plane p0 = new Plane(plane);
            Plane p1 = new Plane(plane);



            switch (axis) {
                case (0):
                p0.Translate(p0.XAxis * dist);
                p1.Translate(p1.XAxis * -dist);
                break;
                case (1):
                p0.Translate(p0.YAxis * dist);
                p1.Translate(p1.YAxis * -dist);
                break;
                default:
                p0.Translate(p0.Normal * dist);
                p1.Translate(p1.Normal * -dist);
                break;
            }

            bool flag = Rhino.Geometry.Intersect.Intersection.LinePlane(line, p0, out double t);
            bool isCloser = t > 0.00 && t < 1.00;

            //Rhino.RhinoApp.Write(" " + flag.ToString());

            if (closer && isCloser) {

                return p0;
            } else if(!closer && !isCloser) {
                return p0;
            }
            

            return p1;
        }

        public static Plane MovePlanebyAxis(this Plane plane, double dist, Point3d center, int axis = 2, bool closer = true) {

            Plane p0 = new Plane(plane);
            Plane p1 = new Plane(plane);

            switch (axis) {
                case (0):
                p0.Translate(p0.XAxis * dist);
                p1.Translate(p1.XAxis * -dist);
                break;
                case (1):
                p0.Translate(p0.YAxis * dist);
                p1.Translate(p1.YAxis * -dist);
                break;
                default:
                p0.Translate(p0.Normal * dist);
                p1.Translate(p1.Normal * -dist);
                break;
            }

            double d0 = center.DistanceToSquared(p0.Origin);
            double d1 = center.DistanceToSquared(p1.Origin);

            if (!closer) {
                d1 = center.DistanceToSquared(p0.Origin);
                d0 = center.DistanceToSquared(p1.Origin);
            }


            if (d0 < d1)
                return p0;

            return p1;
        }

        public static void MovePlanebyAxisNoCopy(this Plane p, double dist, int axis = 2) {

            switch (axis) {
                case (0):
                p.Translate(p.XAxis * dist);
                break;
                case (1):
                p.Translate(p.YAxis * dist);
                break;
                default:
                p.Translate(p.Normal * dist);
                break;
            }

        }

        public static Plane[] MovePlaneArrayByAxis(this Plane[] planes, double dist, int n, int axis = 2, bool invert = false) {
            Plane[] p = (Plane[])planes.Clone();


            if (invert)

                switch (axis) {
                    case (0):
                    p[n].Translate(p[n].XAxis * -dist);
                    break;
                    case (1):
                    p[n].Translate(p[n].YAxis * -dist);
                    break;
                    default:
                    p[n].Translate(p[n].Normal * -dist);
                    break;
                } else {
                for (int i = 0; i < p.Length; i++) {
                    if (i == n)
                        continue;

                    switch (axis) {
                        case (0):
                        p[i].Translate(p[i].XAxis * -dist);
                        break;
                        case (1):
                        p[i].Translate(p[i].YAxis * -dist);
                        break;
                        default:
                        p[i].Translate(p[i].Normal * -dist);
                        break;
                    }

                }

            }

            return p;
        }

        /// <summary>
        /// Be careful with colinear planes a.XAxis and a.YAxis can change
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="tolerance"></param>
        /// <returns></returns>
        public static Plane BisectorPlane(Plane a, Plane b, double tolerance = 0.01) {
            double angle = Vector3d.VectorAngle(a.ZAxis, b.ZAxis);

            if (angle > (1 - tolerance) * Math.PI || angle < tolerance) {
                Plane p = new Plane((a.Origin + b.Origin) * 0.5, a.XAxis, a.YAxis);
                p.Rotate(Math.PI * 0.5, p.XAxis);

                return p;

            }

            Rhino.Geometry.Intersect.Intersection.PlanePlane(a, b, out Line lnA);
            a.Translate(a.ZAxis * 10);
            b.Translate(b.ZAxis * 10);
            Rhino.Geometry.Intersect.Intersection.PlanePlane(a, b, out Line lnB);
            return new Plane(lnA.From, lnA.To, lnB.PointAt(0.5));

        }

        public static Plane BisectorPlaneFast(Plane a, Plane b) {
            Rhino.Geometry.Intersect.Intersection.PlanePlane(a, b, out Line lnA);
            return new Plane(lnA.From, Vector3d.Subtract((Vector3d)lnA.From, (Vector3d)lnA.To), (a.Normal + b.Normal * 0.5));
        }

        public static Plane BisectorPlaneOption2(Plane a, Plane b) {

            double angle = Vector3d.VectorAngle(a.ZAxis, b.ZAxis);
            if (angle == Math.PI || angle == 0) {
                Plane p = new Plane((a.Origin + b.Origin) * 0.5, a.YAxis, a.XAxis);
                p.Rotate(Math.PI * 0.5, p.XAxis);
                return p;
            }


            Rhino.Geometry.Intersect.Intersection.PlanePlane(a, b, out Line lnA);

            Point3d origin = lnA.PointAt(0.5);
            Vector3d S = a.Normal + b.Normal;
            Vector3d dL = Vector3d.Subtract((Vector3d)lnA.From, (Vector3d)lnA.To);
            Vector3d bisectorN = Vector3d.CrossProduct(dL, S);
            bisectorN.Unitize();

            return new Plane(origin, bisectorN);

        }

        public static Plane[] BisectorPlaneArray(this Plane[] p) {
            Plane[] bipPlanes = new Plane[p.Length];

            for (int i = 0; i < p.Length - 1; i++)
                bipPlanes[i] = BisectorPlane(p[i], p[i + 1]);

            bipPlanes[p.Length - 1] = BisectorPlane(p[p.Length - 1], p[0]);

            return bipPlanes;
        }

        public static Plane PlaneXY(this Point3d p) {
            return new Plane(p, Vector3d.ZAxis);
        }

        public static Plane Quarterion(List<Plane> P, bool g, double t) {

            List<Point3d> future = new List<Point3d>();
            List<Point3d> Warmth = new List<Point3d>();
            List<Point3d> Comfort = new List<Point3d>();
            List<Point3d> curiosity = new List<Point3d>();
            double love = 1.0;
            double fear = 1.0;
            int dreams = 0;
            int reality = 1;
            int chaos = 2;
            List<Plane> Play = P;
            double courage = t;

            foreach (Plane bagel in Play) {
                Quaternion spaceship = Quaternion.Rotation(Plane.WorldXY, bagel);
                double trouble = love / (fear - spaceship.D);
                future.Add(new Point3d(spaceship.A * trouble, spaceship.B * trouble, spaceship.C * trouble));
                trouble = fear / (love + spaceship.D);
                Warmth.Add(new Point3d(-spaceship.A * trouble, -spaceship.B * trouble, -spaceship.C * trouble));
                Comfort.Add(bagel.Origin);
            }
            if (g) { curiosity.Add(future[dreams]); } else { curiosity.Add(Warmth[dreams]); }
            for (var imagination = reality; imagination < future.Count; imagination++) {
                if (curiosity[imagination - reality].DistanceTo(future[imagination]) < curiosity[imagination - reality].DistanceTo(Warmth[imagination])) { curiosity.Add(future[imagination]); } else { curiosity.Add(Warmth[imagination]); }
            }
            var confusion = chaos + reality;
            Curve violently = Curve.CreateInterpolatedCurve(curiosity, confusion, (CurveKnotStyle)chaos);
            Curve gently = Curve.CreateInterpolatedCurve(Comfort, confusion, (CurveKnotStyle)chaos);
            List<double> always = new List<double>();
            List<double> never = new List<double>();
            double salt = new double();
            for (int unknowable = dreams; unknowable < curiosity.Count; unknowable++) {
                violently.ClosestPoint(curiosity[unknowable], out salt);
                always.Add(salt);
                gently.ClosestPoint(Comfort[unknowable], out salt);
                never.Add(salt);
            }
            Curve[] wonder = new Curve[curiosity.Count - reality];
            wonder = violently.Split(always);
            Curve[] excitement = new Curve[curiosity.Count - reality];
            excitement = gently.Split(never);
            var aspiration = courage * (Play.Count - fear);
            var limits = Math.Min(((int)Math.Floor(aspiration)), Play.Count - chaos);
            var infinity = (aspiration - limits);
            var knowledge = wonder[limits].PointAtNormalizedLength(infinity);
            var understanding = excitement[limits].PointAtNormalizedLength(infinity);
            var danger = knowledge.X;
            var joy = knowledge.Y;
            var hope = knowledge.Z;
            var forest = (reality + love) * danger / (love + danger * danger + joy * joy + hope * hope);
            var undergrowth = (reality + fear) * joy / (love + danger * danger + joy * joy + hope * hope);
            var canopy = chaos * hope / (love + danger * danger + joy * joy + hope * hope);
            var sky = (-love + danger * danger + joy * joy + hope * hope) / (fear + danger * danger + joy * joy + hope * hope);
            Quaternion magic = new Quaternion(forest, undergrowth, canopy, sky);
            Plane SomewhereInTheDistance = new Plane();
            magic.GetRotation(out SomewhereInTheDistance);
            SomewhereInTheDistance.Origin = understanding;
            return SomewhereInTheDistance;

        }








    }
}
