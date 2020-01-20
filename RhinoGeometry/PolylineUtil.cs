using System;
using System.Collections.Generic;
using System.Linq;
using Rhino.Geometry;
using Rhino.Geometry.Intersect;


namespace RhinoGeometry {



    public static class PolylineUtil {

        /// <summary>
        /// Scale and Shear Polyline by Vector 
        /// </summary>
        /// <param name="pline - tile to transform"></param>
        /// <param name="scale - scale tile after stretching to maintain original thickness"></param>
        /// <param name="scaleX - Scale tile before shear"></param>
        /// <param name="scaleY - Scale tile before shear"></param>
        /// <param name="v - joint direction"></param>
        /// <returns></returns>
        public static Polyline ScaleShearByVector( this Polyline pline, Vector3d v, bool scale = false, double scaleX=1, double scaleY=1, double scaleZ=1, int translation = 0, int rotateDegrees = 0, bool split = false ) {

            Polyline p = pline.Duplicate();


            //Rotate tile
            if (rotateDegrees != 0) {
                p.Transform(Rhino.Geometry.Transform.Rotation(Rhino.RhinoMath.ToRadians(rotateDegrees),Point3d.Origin));
            }

            //Translate to change origin point
           p.Transform(Rhino.Geometry.Transform.Translation( new Vector3d(0.5*translation,0,0)));
            //Rhino.RhinoDoc.ActiveDoc.Objects.AddPolyline(p);

            //Scale based on initial angle
            p.Transform(Rhino.Geometry.Transform.Scale(Plane.WorldXY, scaleX, scaleY, scaleZ));

           

            //Scale based on sin angle
            if (scale) {
                double angleScale = Math.Sin(Vector3d.VectorAngle(v, Vector3d.YAxis));
                if (angleScale > 0.01)
                    p.Transform(Rhino.Geometry.Transform.Scale(Plane.WorldXY, 1, 1 / angleScale, 1));
            }

         

            //Shear based on tan angle
            double angleShear = Vector3d.VectorAngle(v, Vector3d.XAxis, Plane.WorldXY);
            double xShear = Math.Tan(angleShear) * -1;
            double yShear = 0.00;
            p.Transform(Rhino.Geometry.Transform.Shear(Plane.WorldXY, new Vector3d(1, xShear, 0), new Vector3d(yShear, 1, 0), new Vector3d(0, 0, 1)));


            //Split if needed
            if (split) {
                Tuple<SortedDictionary<int, Polyline>, SortedDictionary<int, Polyline>> splitPlines = SplitCurvesByPlane(p, Plane.WorldYZ);


                //Shear the polyline negatively
                foreach (var pl in splitPlines.Item1) {
                    pl.Value.Transform(Rhino.Geometry.Transform.Shear(Plane.WorldXY, new Vector3d(1, 2 * -xShear, 0), new Vector3d(2 * -yShear, 1, 0), new Vector3d(0, 0, 1)));
                    splitPlines.Item2.Add(pl.Key, pl.Value);
                }

                //Combine the polyline back
                Polyline joinedPline = new Polyline();
                foreach (var pl in splitPlines.Item2)
                    joinedPline.AddRange(pl.Value);
                joinedPline.CollapseShortSegments(0.01);

                return joinedPline;
            }

            return p;
        }

        public static Tuple<SortedDictionary<int, Polyline>, SortedDictionary<int, Polyline>> SplitCurvesByPlane(this Polyline pline, Plane plane) {


            //Intersect each segment of curve to get cut parameters
            List<double> T = new List<double>();
            for (int i = 0; i < pline.SegmentCount; i++) {
                double t;
                bool flag = Rhino.Geometry.Intersect.Intersection.LinePlane(pline.SegmentAt(i), plane, out t);
                if (flag) {
                    T.Add(i + t);
                    //Print((i + t).ToString());
                }
            }

            var TArray = T.ToArray();
            Array.Sort(TArray);

            //Split curve by parameters
            Curve c = pline.ToNurbsCurve();
            Curve[] splitCurves = c.Split(TArray);

            //Check if curve on which side
            SortedDictionary<int, Polyline> splitPlines0 = new SortedDictionary<int, Polyline>();
            SortedDictionary<int, Polyline> splitPlines1 = new SortedDictionary<int, Polyline>();


            Plane planeOffset = new Plane(plane);
            planeOffset.Translate(plane.ZAxis * 0.01);

            int key = 0;
            foreach (Curve sp in splitCurves) {

                Point3d pt = sp.PointAtNormalizedLength(0.5);
                Point3d cp0 = plane.ClosestPoint(pt);
                Point3d cp1 = planeOffset.ClosestPoint(pt);

                Polyline cutPline;

                if (sp.TryGetPolyline(out cutPline)) {
                    if (cp0.DistanceToSquared(pt) < cp1.DistanceToSquared(pt)) {
                        splitPlines0.Add(key, cutPline);
                        //Rhino.RhinoDoc.ActiveDoc.Objects.AddPolyline( cutPline);
                    } else {
                        splitPlines1.Add(key, cutPline);
                    }

                }

                key++;
            }

            return new Tuple<SortedDictionary<int, Polyline>, SortedDictionary<int, Polyline>>(splitPlines0, splitPlines1);

        }


        public static Plane[] PlanesFromPairOfPolyline(Polyline pline0, Polyline pline1) {

            Plane[] planes = new Plane[pline0.Count-1];

            for(int i = 0; i < pline0.Count-1; i++) {
                Point3d origin =( pline0[0] + pline0[1] + pline1[0] + pline1[1]) * 0.25;
                Vector3d XAxis = pline0[0] - pline0[1];
                Vector3d YAxis = pline0[0] - pline1[0];
                Plane plane = new Plane(origin,XAxis,YAxis);
                planes[i] = plane;
            }
            return planes;
        }

        public static List<Polyline> Duplicate(this List<Polyline> p) {
            List<Polyline> polylines = new List<Polyline>();
            foreach (Polyline p_ in p)
                polylines.Add(new Polyline(p_));
            return polylines;
        }

        /// <summary>
        /// Takes 4 point rectangle 1st - 3rd or 2nd - 4th lines and creates zigzag 
        ///returns good result only when rectangle edges are equal length    
        /// </summary>
        /// <param name="rectangle"></param>
        /// <param name="flip"></param>
        /// <param name="dist"></param>
        /// <returns></returns>
        public static Polyline ZigZag(this Polyline rectangle, bool flip, double dist, int divisions = -1) {

            Line l0 = (flip) ? new Line(rectangle[0], rectangle[1]) : new Line(rectangle[1], rectangle[2]);
            Line l1 = (flip) ? new Line(rectangle[3], rectangle[2]) : new Line(rectangle[0], rectangle[3]);

            int n = (int)Math.Ceiling(l0.Length / dist);

            Polyline zigzag = new Polyline();
            Vector3d v0 = l0.Direction.UnitVector() * dist;


            for (int i = 0; i < n; i++) {
                if (i % 2 == 0) {
                    zigzag.Add(l0.From + (v0 * i));
                    zigzag.Add(l1.From + (v0 * (i)));

                } else {
                    zigzag.Add(l1.From + (v0 * (i)));
                    zigzag.Add(l0.From + (v0 * (i)));

                }//if
            }//for

            //Finish the zigzag end
            if (dist * n != l0.Length) {
                if (n % 2 == 0) {
                    zigzag.Add(l0.To);
                    zigzag.Add(l1.To);

                } else {
                    zigzag.Add(l1.To);
                    zigzag.Add(l0.To);

                }//if
            }

            return zigzag;
        }
        public static Polyline Flip(this Polyline polyline) {
            Polyline p = polyline.Duplicate();
            p.Reverse();
            return p;
        }

        public static void Orient(this Polyline polyline, Plane source, Plane target) {
            Transform transform = Rhino.Geometry.Transform.PlaneToPlane(source, target);
            polyline.Transform(transform);
        }

        public static List<Line> Transform(this List<Line> lines, Transform t) {

            List<Line> linesTransform = new List<Line>();

            for (int i = 0; i < lines.Count; i++) {
                Line l = lines[i];
                l.Transform(t);
                linesTransform.Add(l);
            }

            return linesTransform;

        }

        public static Line[] Transform(this Line[] lines, Transform t) {

            Line[] linesTransform = new Line[lines.Length];

            for (int i = 0; i < lines.Length; i++) {
                Line l = lines[i];
                l.Transform(t);
                linesTransform[i] = l;
            }

            return linesTransform;

        }

        public static List<Line[]> Transform(this List<Line[]> lines, Transform t) {

            List<Line[]> linesTransformed = new List<Line[]>(lines.Count);

            foreach (Line[] ln in lines) {

                linesTransformed.Add(ln.Transform(t));
            }
            return linesTransformed;
        }

        public static void Transform(this IEnumerable<Polyline> polylines, Transform t) {
            foreach (Polyline p in polylines)
                p.Transform(t);
        }

        public static Polyline Translate(this Polyline p, Vector3d v) {
            Polyline polyline = new Polyline(p);
            polyline.Transform(Rhino.Geometry.Transform.Translation(v));
            return polyline;
        }

        public static List<Polyline> Polygons(this IEnumerable<Circle> curves, int n, double rotation = Math.PI * 0.25, bool sqrt = true) {
            List<Polyline> polygons = new List<Polyline>();
            foreach (var c in curves)
                polygons.Add(Polygon(n, c, rotation, sqrt));
            return polygons;
        }

        public static List<Polyline> Polygons(this IEnumerable<Curve> curves, int n, double rotation = Math.PI * 0.25, bool sqrt = true) {
            List<Polyline> polygons = new List<Polyline>();
            foreach (var c in curves)
                polygons.Add(Polygon(n, c, rotation, sqrt));
            return polygons;
        }

        public static Polyline Polygon(int n, Circle circle, double rotation = Math.PI * 0.25, bool sqrt = true) {


            return Polygon(n, circle.Radius, circle.Plane, rotation, sqrt);

        }

        public static Polyline Polygon(int n, Curve c, double rotation = Math.PI * 0.25, bool sqrt = true) {

            if (c.TryGetCircle(out Circle circle)) {
                return Polygon(n, circle.Radius, circle.Plane, rotation, sqrt);

            }

            return new Polyline();
        }

        public static Polyline Polygon(int n, double radius, Plane plane, double rotation = Math.PI * 0.25, bool sqrt = true) {

            Polyline polyline = new Polyline();
            double sector = Math.PI * 2 / n;
            double r = sqrt ? 1 / Math.Sqrt(2) * radius : radius;

            for (int i = 0; i < n; i++) {
                Point3d p = new Point3d(Math.Sin((sector * i) + rotation) * r, Math.Cos((sector * i) + rotation) * r, 0);
                polyline.Add(p);
            }

            polyline.Add(polyline[0]);

            polyline.Transform(Rhino.Geometry.Transform.PlaneToPlane(Plane.WorldXY, plane));

            return polyline;
        }

        public static Point3d CenterPoint(Polyline polyline) {

            int Count = polyline.Count;

            if (Count == 0) { return Point3d.Unset; }

            if (Count == 1) { return polyline[0]; }



            Point3d center = Point3d.Origin;

            double weight = 0.0;

            int stop = (Count - 1);
            if (polyline[0].DistanceToSquared(polyline[polyline.Count - 1]) > 0.001) {
                //Rhino.RhinoApp.WriteLine(polyline[0].DistanceToSquared(polyline[polyline.Length - 1]).ToString());
                stop++;
            }
            for (int i = 0; i < stop; i++) {

                Point3d A = polyline[i];

                Point3d B = polyline[(i + 1) % Count];

                double d = A.DistanceTo(B);

                center += d * 0.5 * (A + B);

                weight += d;

            }

            center /= weight;

            return center;

        }

        public static Polyline[] IntersectTwoPlates(Polyline plateA0, Polyline plateA1, Polyline plateB0, Polyline plateB1) {


            double extend = 0;

            var pA0B0 = PolylinePlane(plateA0, plateB0);
            var pA0B1 = PolylinePlane(plateA0, plateB1);
            var pA1B0 = PolylinePlane(plateA1, plateB0);
            var pA1B1 = PolylinePlane(plateA1, plateB1);

            if (pA0B0[0] != Point3d.Unset && pA0B0[1] != Point3d.Unset && pA0B1[0] != Point3d.Unset && pA0B1[1] != Point3d.Unset &&
                pA1B0[0] != Point3d.Unset && pA1B0[1] != Point3d.Unset && pA1B1[0] != Point3d.Unset && pA1B1[1] != Point3d.Unset
                ) {


                Point3d pA0B0Mid = PointUtil.AveragePoint(pA0B0);
                Point3d pA0B1Mid = PointUtil.AveragePoint(pA0B1);
                Point3d pA1B0Mid = PointUtil.AveragePoint(pA1B0);
                Point3d pA1B1Mid = PointUtil.AveragePoint(pA1B1);

                Vector3d v0 = (pA0B0[0] - pA0B0[1]).UnitVector() * extend;
                Vector3d v1 = (pA0B1[0] - pA0B1[1]).UnitVector() * extend;
                Vector3d v2 = (pA1B0[0] - pA1B0[1]).UnitVector() * extend;
                Vector3d v3 = (pA1B1[0] - pA1B1[1]).UnitVector() * extend;

                Polyline u0 = new Polyline(new Point3d[] { pA0B0[0], pA0B0Mid - v0, pA0B1Mid - v1, pA0B1[0] });
                Polyline u1 = new Polyline(new Point3d[] { pA1B0[0], pA1B0Mid - v2, pA1B1Mid - v3, pA1B1[0] });

                Polyline u2 = new Polyline(new Point3d[] { pA0B0[1], pA0B0Mid + v0, pA1B0Mid + v2, pA1B0[1] });
                Polyline u3 = new Polyline(new Point3d[] { pA0B1[1], pA0B1Mid + v1, pA1B1Mid + v3, pA1B1[1] });

                return new Polyline[] { u0, u1, u2, u3 };

            }


            return new Polyline[0];



        }

        public static Point3d[] PolylinePlane(Polyline p0, Polyline p1) {

            Line[] segmentsP0 = p0.GetSegments();
            Line[] segmentsP1 = p1.GetSegments();

            Plane pl0 = p0.plane();
            Plane pl1 = p1.plane();

            //Rhino.RhinoDoc.ActiveDoc.Objects.AddPolyline(p0);
            //Rhino.RhinoDoc.ActiveDoc.Objects.AddPolyline(p1);


            Point3d ptP0 = Point3d.Unset;
            Point3d ptP1 = Point3d.Unset;

            foreach (Line line in segmentsP0) {

                double t;
                if (Intersection.LinePlane(line, pl1, out t)) {
                    if (t > 1 || t < 0)
                        continue;
                    Point3d pTemp = line.PointAt(t);

                    if (p1.ToNurbsCurve().Contains(pTemp, pl1, 0.01) == PointContainment.Inside) {
                        ptP0 = pTemp;
                        break;
                    }
                }
            }

            foreach (Line lineP1 in segmentsP1) {

                double t;
                if (Intersection.LinePlane(lineP1, pl0, out t)) {
                    if (t > 1 || t < 0)
                        continue;
                    Point3d pTemp = lineP1.PointAt(t);
                    if (p0.ToNurbsCurve().Contains(pTemp, pl0, 0.01) == PointContainment.Inside) {
                        ptP1 = pTemp;
                        break;
                    }
                }
            }

            return new Point3d[] { ptP0, ptP1 };
        }

        public static Line[] LoftLine(Polyline p0, Polyline p1) {
            Line[] line = new Line[p0.Count - 1];

            for (int i = 0; i < p0.Count - 1; i++) {
                line[i] = new Line(p0[i], p1[i]);
            }
            return line;
        }

        public static void InsertPoint(this Polyline polyline, Point3d p) {
            polyline.Insert((int)Math.Ceiling(polyline.ClosestParameter(p)), p);
        }
        public static List<Polyline> MappedFromSurfaceToSurface(this List<Polyline> polylines, Surface s, Surface t) {

            s.SetDomain(0, new Interval(0, 1));
            s.SetDomain(1, new Interval(0, 1));
            t.SetDomain(0, new Interval(0, 1));
            t.SetDomain(1, new Interval(0, 1));
            //  Rhino.RhinoDoc.ActiveDoc.Objects.AddSurface(s);

            List<Polyline> mapped = new List<Polyline>();

            for (int i = 0; i < polylines.Count; i++) {

                Polyline pol = new Polyline(polylines[i]);

                bool flag = true;

                for (int j = 0; j < pol.Count; j++) {
                    double u, v;
                    Point3d pTemp = new Point3d(pol[j]);
                    s.ClosestPoint(pol[j], out u, out v);

                    pol[j] = t.PointAt(u, v);

                    if (s.PointAt(u, v).DistanceTo(pTemp) > 0.01) {
                        flag = false;
                        break;
                    }

                }//for j
                if (flag) {
                    mapped.Add(pol);
                }
            }//for i
            return mapped;

        }


        public static List<Polyline> MappedFromMeshToMesh(this List<Polyline> polylines, Mesh s, Mesh t) {

            //    s.SetDomain(0, new Interval(0, 1));
            //    s.SetDomain(1, new Interval(0, 1));
            //    t.SetDomain(0, new Interval(0, 1));
            //    t.SetDomain(1, new Interval(0, 1));
            //  Rhino.RhinoDoc.ActiveDoc.Objects.AddSurface(s);

            List<Polyline> mapped = new List<Polyline>();

            for (int i = 0; i < polylines.Count; i++) {

                Polyline pol = new Polyline(polylines[i]);

                bool flag = true;

                for (int j = 0; j < pol.Count; j++) {
                    //point3d = mesh.PointAt(index, t, num, t1, num1);
                    //vector3d = mesh.NormalAt(index, t, num, t1, num1);

                    Point3d pTemp = new Point3d(pol[j]);

                    MeshPoint mp = s.ClosestMeshPoint(pol[j], 10.01);
                    if (mp == null) {
                        flag = false;
                        break;
                    }

                    //   Rhino.RhinoDoc.ActiveDoc.Objects.AddPoint(mp.Point);

                    pol[j] = t.PointAt(mp);

                    if (s.PointAt(mp).DistanceTo(pTemp) > 0.01) {
                        flag = false;
                        break;
                    }

                }//for j





                //try to trim curve
                if (!flag) {

                    //How to check if curve is on the first?
                    MeshPoint mp0 = s.ClosestMeshPoint(pol[0], 10);
                    MeshPoint mp1 = s.ClosestMeshPoint(pol[pol.Count - 1], 10);


                }


                if (flag) {
                    mapped.Add(pol);
                    //Rhino.RhinoDoc.ActiveDoc.Objects.AddPolyline(pol);
                }
            }//for i
            return mapped;

        }

        public static Polyline OutlineFromFaceEdgeCorner(Plane facePlane, Plane[] edgePlanes, Plane[] bisePlanes, int T = 1, double tolerance = 0.1) {
            Polyline polyline = new Polyline();

            switch (T) {
                case (2):

                for (int j = 0; j < edgePlanes.Length; j++) {

                    Plane currPlane = edgePlanes[j];
                    Plane nextPlane = edgePlanes[MathUtil.Wrap(j + 1, edgePlanes.Length)];


                    if (Vector3d.VectorAngle(currPlane.XAxis, nextPlane.XAxis) < tolerance) {
                        Vector3d vv = new Vector3d(currPlane.XAxis);
                        vv.Rotate(Math.PI * 0.5, currPlane.YAxis);
                        nextPlane = new Plane(bisePlanes[j].Origin, vv, currPlane.YAxis);
                    }

                    Line line = PlaneUtil.PlanePlane(currPlane, nextPlane);
                    polyline.Add(PlaneUtil.LinePlane(line, facePlane));

                }
                polyline.Close();
                break;

                default:
                for (int j = 0; j < bisePlanes.Length; j++) {
                    Point3d pt;
                    Rhino.Geometry.Intersect.Intersection.PlanePlanePlane(facePlane, bisePlanes[j], edgePlanes[j], out pt);
                    polyline.Add(pt);
                }
                polyline.Close();
                break;
            }


            return polyline;
        }

        public static Polyline ProjectPolyline(this Polyline p, bool averagePlane = true) {


            Polyline p_ = new Polyline(p.Count);

            Plane plane = p.plane();
            if (!averagePlane)
                Plane.FitPlaneToPoints(p, out plane);

            int closed = (p[0].DistanceToSquared(p[p.Count - 1]) < 0.001) ? 1 : 0;

            for (int i = 0; i < p.Count - closed; i++) {
                p_.Add(plane.ClosestPoint(p[i]));
            }

            if (closed == 1)
                p_.Close();


            return p_;

        }

        public static List<Curve> ShatterCurve(Curve C, List<Point3d> P, double D = 0.01) {
            List<Curve> curves = new List<Curve>();
            try {

                List<double> shatterT = new List<double>();

                Curve c = C.DuplicateCurve();
                //c.Domain = new Interval(0, 1);
                double a = c.Domain.T0;
                double b = c.Domain.T1;

                for (int i = 0; i < P.Count; i++) {
                    double t;
                    bool flag = c.ClosestPoint(P[i], out t, D);

                    if (flag) {
                        shatterT.Add((double)Math.Round(t, 5));
                        // shatterT.Add(t);
                    }
                }

                if (shatterT.Count > 0) {

                    shatterT.Sort();
                    shatterT = shatterT.Distinct().ToList();


                    if (shatterT[0] == a && shatterT[shatterT.Count - 1] == b)
                        shatterT.RemoveAt(0);


                    if (c.PointAtStart.DistanceToSquared(c.PointAtEnd) < 0.001) {
                        for (int i = 0; i < shatterT.Count; i++) {
                            //Rhino.RhinoApp.WriteLine(shatterT[i].ToString() +  " " + shatterT[(i + 1) % shatterT.Count].ToString()) ;


                            if (shatterT[shatterT.Count - 1] == b) {
                                Curve curve = c.Trim(shatterT[MathUtil.Wrap(i - 1, shatterT.Count)], shatterT[i]);
                                curves.Add(curve);
                            } else {
                                Curve curve = c.Trim(shatterT[i], shatterT[(i + 1) % shatterT.Count]);
                                curves.Add(curve);
                            }

                        }
                    } else {

                        if (shatterT[0] != a)
                            shatterT.Insert(0, a);

                        if (shatterT[shatterT.Count - 1] != b)
                            shatterT.Insert(shatterT.Count, b);

                        for (int i = 0; i < shatterT.Count - 1; i++) {
                            Curve curve = c.Trim(shatterT[i], shatterT[(i + 1)]);
                            curves.Add(curve);
                        }


                    }


                } else {
                    curves.Add(c);
                }


            } catch (Exception e) {
                Rhino.RhinoApp.WriteLine(e.ToString());
            }
            return curves;
        }




        public static Polyline IntersectPlanarLines(Line[] lines, bool close = true) {

            Polyline polyline = new Polyline();

            for (int i = 0; i < lines.Length; i++) {

                double a, b;
                Rhino.Geometry.Intersect.Intersection.LineLine(lines[i], lines[(i + 1) % lines.Length], out a, out b, Rhino.RhinoDoc.ActiveDoc.ModelAbsoluteTolerance, false);
                Point3d p = lines[i].PointAt(a);
                polyline.Add(p);

            }

            if (close)
                polyline.Add(polyline[0]);

            return polyline;

        }

        public static Polyline IntersectPlanarLines(Line[] lines, Plane projectToPlane, bool close = true) {

            for (int i = 0; i < lines.Length; i++) {
                lines[i].From = projectToPlane.ClosestPoint(lines[i].From);
                lines[i].To = projectToPlane.ClosestPoint(lines[i].To);
            }

            return IntersectPlanarLines(lines, close);

        }


        
        private static Polyline OrderPolyline(this Polyline polyline) {

            if (polyline != null) {
                if (polyline.IsValid) {
                    bool flag = IsClockwiseClosedPolylineOnXYPlane(polyline);
                    //isClockWise = flag;
                    if (flag) {
                        Polyline polyline0 = new Polyline(polyline);
                        polyline0.Reverse();


                        return polyline0;

                    } else {
                        return polyline;
                    }
                }


            }
            return polyline;
        }

        public static bool IsClockwiseClosedPolylineOnXYPlane(this Polyline polygon) {
            double sum = 0;

            for (int i = 0; i < polygon.Count - 1; i++)
                sum += (polygon[i + 1].X - polygon[i].X) * (polygon[i + 1].Y + polygon[i].Y);

            return sum > 0;
        }

        public static Polyline Chamfer(this Polyline polyline, double value = 0.001) {

            Line[] lines = polyline.GetSegments();

            if (value == 0)
                return polyline;

            Polyline p = new Polyline();

            if (value < 0) {
                foreach (Line l in lines) {
                    Line lShorter = CurveUtil.ExtendLine(l, -Math.Abs(value));
                    p.Add(lShorter.From);
                    p.Add(lShorter.To);
                }
            } else {
                foreach (Line l in lines) {
                    p.Add(l.PointAt(value));
                    p.Add(l.PointAt(1 - value));
                }
            }

            p.Add(p[0]);

            List<Point3d> points = new List<Point3d>();

            for (int i = 1; i < p.Count - 1; i += 2) {
                points.Add(new Point3d(
                  (p[i + 1].X + p[i].X) * 0.5,
                  (p[i + 1].Y + p[i].Y) * 0.5,
                  (p[i + 1].Z + p[i].Z) * 0.5
                  )
                  );
            }

            return p;


        }


        public static void InsertPolyline(this Polyline x, IEnumerable<Polyline> y) {

            foreach (Polyline poly in y) {
                double t0 = x.ClosestParameter(poly[0]);
                double t1 = x.ClosestParameter(poly.Last);
                if (t0 > t1)
                    poly.Reverse();
                x.InsertRange((int)Math.Ceiling(t0), poly);
            }
        }


        public static void InsertPolyline(this Polyline x, Polyline poly) {


            double t0 = x.ClosestParameter(poly[0]);
            double t1 = x.ClosestParameter(poly.Last);
            if (t0 > t1)
                poly.Reverse();
            x.InsertRange((int)Math.Ceiling(t0), poly);

        }

        public static List<Polyline> InterpolateTwoLines(Line l0, Line l1, int n = 1) {
            List<Polyline> squares = new List<Polyline>();

            if (n > 0) {


                Point3d[] interpolatePt0 = PointUtil.InterpolatePoints(l0.From, l0.To, n); //bottom interpolation
                Point3d[] interpolatePt1 = PointUtil.InterpolatePoints(l1.From, l1.To, n); //top inerpolation

                for (int i = 0; i < n + 1; i++) {
                    Polyline polyline = new Polyline(new[] {
                    interpolatePt0[i],
                    interpolatePt0[i+1],
                    interpolatePt1[i+1],
                    interpolatePt1[i],
                    interpolatePt0[i]
                });
                    squares.Add(polyline);
                }

                return squares;

            } else {


                squares.Add(new Polyline(new[] {
                    l0.From,
                    l0.To,
                    l1.To,
                    l1.From,
                    l0.From
                    }));


            }

            return squares;


        }


        public static Line PlanePlanePlanePlane(Plane planeToIntersectLine0, Plane planeToIntersectLine1, Plane planeToCutLine0, Plane planeToCutLine1) {

            Rhino.Geometry.Intersect.Intersection.PlanePlane(planeToIntersectLine0, planeToIntersectLine1, out Line intersectedLine);
            Rhino.Geometry.Intersect.Intersection.LinePlane(intersectedLine, planeToCutLine0, out double t0);
            Rhino.Geometry.Intersect.Intersection.LinePlane(intersectedLine, planeToCutLine1, out double t1);
            Line line = new Line(intersectedLine.PointAt(t0), intersectedLine.PointAt(t1));
            return line;
        }

        public static Polyline PlaneLines(Plane plane, IEnumerable<Line> lines, bool close = true) {


            Polyline polyline = new Polyline();

            foreach (Line l in lines) {

                Rhino.Geometry.Intersect.Intersection.LinePlane(l, plane, out double t);
                polyline.Add(l.PointAt(t));
            }
            polyline.Add(polyline[0]);
            return polyline;
        }
        public static Plane MovePlane(Plane Base, double Thickness) {
            return new Plane(Base.Origin + Base.Normal * (Thickness / 2), Base.Normal);
        }

        public static Point3d IntPtPln1(Line line, Plane plane) {
            double pm;
            Rhino.Geometry.Intersect.Intersection.LinePlane(line, plane, out pm);
            return line.PointAt(pm);
        }

        public static Line IntPtPln(Point3d origin, Vector3d axis, Plane planeBot, Plane planeTop) {
            Line lineTEMP = new Line(origin, origin + axis);
            double pmBot;
            double pmTop;
            Rhino.Geometry.Intersect.Intersection.LinePlane(lineTEMP, planeBot, out pmBot);
            Rhino.Geometry.Intersect.Intersection.LinePlane(lineTEMP, planeTop, out pmTop);

            return new Line(lineTEMP.PointAt(pmBot), lineTEMP.PointAt(pmTop));
        }

        public static Line tweenLine(Line l0, Line l1, double t = 0.5) {
            return new Line(MathUtil.Lerp(l0.From, l1.From, t), MathUtil.Lerp(l0.To, l1.To, t));
        }

        public static Polyline tweenPolylines(Polyline l0, Polyline l1, double t = 0.5) {

            Polyline p = new Polyline(l0);

            for (int i = 0; i < l0.Count; i++) {
                p[i] = MathUtil.Lerp(l0[i], l1[i], t);
            }


            return p;
        }

        public static Polyline PolylineFromPlanes(Plane basePlane, List<Plane> sidePlanes, bool close = true) {

            Polyline polyline = new Polyline();

            for (int i = 0; i < sidePlanes.Count - 1; i++) {
                Rhino.Geometry.Intersect.Intersection.PlanePlanePlane(basePlane, sidePlanes[i], sidePlanes[i + 1], out Point3d pt);
                polyline.Add(pt);
            }

            if (close) {
                Rhino.Geometry.Intersect.Intersection.PlanePlanePlane(basePlane, sidePlanes[sidePlanes.Count - 1], sidePlanes[0], out Point3d pt1);
                polyline.Add(pt1);
            }

            if (close)
                polyline.Add(polyline[0]);

            return polyline;

        }


        /// <summary>
        /// Move two planes by z axis and output one that is closer to target point
        /// </summary>
        /// <param name="polyline"></param>
        /// <param name="point"></param>
        /// <param name="dist"></param>
        /// <returns></returns>
        public static Plane MovePolylinePlaneToPoint(this Polyline polyline, Point3d point, double dist) {

            Plane planeA = polyline.plane();
            planeA.Translate(planeA.ZAxis * dist);

            Plane planeB = polyline.plane();
            planeB.Translate(planeB.ZAxis * -dist);

            double dA = PointUtil.FastDistance(point, planeA.Origin);
            double dB = PointUtil.FastDistance(point, planeB.Origin);

            if (dA > dB)
                return planeB;

            return planeA;
        }




        public static Polyline DovetailPolyline(Line lineA, Line lineB, Line lineA_, Line lineB_, int d) {
            Point3d[] interA = PointUtil.InterpolatePoints(lineA.From, lineB.From, d);
            Point3d[] interA_ = PointUtil.InterpolatePoints(lineA_.From, lineB_.From, d);



            Point3d[] interB;
            Point3d[] interB_;

            if (d % 2 == 1) {
                interB = PointUtil.InterpolatePoints(lineB.To, lineA.To, d);
                interB_ = PointUtil.InterpolatePoints(lineB_.To, lineA_.To, d);
            } else {
                interB_ = PointUtil.InterpolatePoints(lineB.To, lineA.To, d);
                interB = PointUtil.InterpolatePoints(lineB_.To, lineA_.To, d);
            }



            Polyline poly = new Polyline();
            Polyline temp = new Polyline();

            for (int j = 0; j < interA.Length - 1; j++) {
                if (j % 2 == 0) {
                    poly.Add(interA[j]);
                    poly.Add(interA[j + 1]);
                    temp.Add(interB[j]);
                    temp.Add(interB[j + 1]);
                } else {
                    poly.Add(interA_[j]);
                    poly.Add(interA_[j + 1]);
                    temp.Add(interB_[j]);
                    temp.Add(interB_[j + 1]);
                }
            }

            poly.AddRange(temp);
            poly.Close();
            return poly;
        }

        public static Polyline DovetailPolylineShifted(Point3d[] pts, int d) {
            Point3d[] interA = PointUtil.InterpolatePoints(pts[0], pts[1], d);
            Point3d[] interA_ = PointUtil.InterpolatePoints(pts[2], pts[3], d);



            Point3d[] interB;
            Point3d[] interB_;

            if (d % 2 == 1) {
                interB = PointUtil.InterpolatePoints(pts[5], pts[4], d);
                interB_ = PointUtil.InterpolatePoints(pts[7], pts[6], d);
            } else {
                interB_ = PointUtil.InterpolatePoints(pts[5], pts[4], d);
                interB = PointUtil.InterpolatePoints(pts[7], pts[6], d);
            }



            Polyline poly = new Polyline();
            Polyline temp = new Polyline();

            for (int j = 0; j < interA.Length - 1; j++) {
                if (j % 2 == 0) {
                    poly.Add(interA[j]);
                    poly.Add(interA[j + 1]);
                    temp.Add(interB[j]);
                    temp.Add(interB[j + 1]);
                } else {
                    poly.Add(interA_[j]);
                    poly.Add(interA_[j + 1]);
                    temp.Add(interB_[j]);
                    temp.Add(interB_[j + 1]);
                }
            }

            poly.AddRange(temp);
            poly.Close();
            return poly;
        }

        public static Line IntersectionPlaneTwoLines(Plane p, Line lnA, Line lnB) {
            Intersection.LinePlane(lnA, p, out double t1);
            Intersection.LinePlane(lnB, p, out double t2);
            return new Line(lnA.PointAt(t1), lnB.PointAt(t2));
        }

        public static double[] IntersectionPlaneTwoLinesT(Plane p, Line lnA, Line lnB) {
            Intersection.LinePlane(lnA, p, out double t1);
            Intersection.LinePlane(lnB, p, out double t2);
            return new[] { t1, t2 };
        }

        public static Line IntersectionLineTwoPlanes(this Line line, Plane pa, Plane pb) {
            Intersection.LinePlane(line, pa, out double t1);
            Intersection.LinePlane(line, pb, out double t2);
            return new Line(line.PointAt(t1), line.PointAt(t2));
        }

        public static double[] IntersectionLineTwoPlanesT(this Line line, Plane pa, Plane pb) {
            Intersection.LinePlane(line, pa, out double t1);
            Intersection.LinePlane(line, pb, out double t2);
            return new[] { t1, t2 };
        }



        public static Polyline[] ToPolylines(this IEnumerable<Curve> nurbsCurves, bool collapseShortSegments = true) {

            Polyline[] p = new Polyline[nurbsCurves.Count()];

            for (int i = 0; i < nurbsCurves.Count(); i++) {
                nurbsCurves.ElementAt(i).TryGetPolyline(out Polyline polyline);
                p[i] = nurbsCurves.ElementAt(i).ToPolyline(collapseShortSegments);
            }

            return p;

        }


        public static Polyline ToPolyline(this Curve curve, bool collapseShortSegments = true) {

            curve.TryGetPolyline(out Polyline polyline);
            if (collapseShortSegments)
                polyline.CollapseShortSegments(Rhino.RhinoDoc.ActiveDoc.ModelAbsoluteTolerance);
            return polyline;

        }


        public static Polyline[] ToPolylinesFromCP(this IEnumerable<Curve> curves, double collapseShortSegments = 0.01) {

            Polyline[] p = new Polyline[curves.Count()];

            int j = 0;
            foreach (Curve curve in curves) {
                p[j++] = ToPolylineFromCP(curve, collapseShortSegments);
            }

            return p;
        }

        public static Polyline ToPolylineFromCP(this Curve curve, double collapseShortSegments = 0.01) {

            Polyline polyline = new Polyline();
            if (curve.TryGetPolyline(out polyline)) {
                polyline.CollapseShortSegments(collapseShortSegments);
                return polyline;
            }


            NurbsCurve c = curve.ToNurbsCurve();

            Point3d[] points = new Point3d[c.Points.Count];

            for (int i = 0; i < c.Points.Count; i++)
                c.Points.GetPoint(i, out points[i]);

            polyline = new Polyline(points);


            //What the fuck these two lines
            c = polyline.ToNurbsCurve();
            c.TryGetPolyline(out polyline);



            if (collapseShortSegments > 0)
                polyline.CollapseShortSegments(collapseShortSegments);
            polyline = new Polyline(polyline);

            //polyline.CollapseShortSegments(1);
            return polyline;


        }

        public static int[] append(int j) {
            var srsEnum = Enumerable.Range(0, j + 1);
            var arrVal = new List<int>();
            arrVal = srsEnum.ToList();
            arrVal.Add(0);
            return arrVal.ToArray();
        }


        public static Plane[] Planes(this Polyline[] polylines) {
            Plane[] p = new Plane[polylines.Length];

            for (int i = 0; i < polylines.Length; i++) {
                //p[i] = new Plane(polylines[i].CenterPoint(), polylines[i].Normal());
                p[i] = polylines[i].GetPlane();
            }

            return p;
        }

        public static Plane plane(this Polyline polylines) {
            return polylines.GetPlane();
            //return new Plane(polylines.CenterPoint(), polylines.AverageNormal());
        }

        public static Vector3d AverageNormal(this Polyline p) {
            //PolyFace item = this[index];
            int len = p.Count - 1;
            Vector3d vector3d = new Vector3d();
            int count = checked(len - 1);

            for (int i = 0; i <= count; i++) {
                int num = ((i - 1) + len) % len;
                int item1 = (checked(i + 1) + len) % len;
                Point3d point3d = p[num];
                Point3d point3d1 = p[item1];
                Point3d item2 = p[i];
                vector3d = vector3d + Vector3d.CrossProduct(new Vector3d(item2 - point3d), new Vector3d(point3d1 - item2));
            }

            if (vector3d.X == 0 & vector3d.Y == 0 & vector3d.Z == 0)
                vector3d.Unitize();

            return vector3d;
        }

        public static Vector3d Normal(this Polyline p) {
            var n = Vector3d.Unset;
            if (null != p && p.Count - 1 > 3) {
                var a = p[2] - p[0];
                var b = p[3] - p[1];
                n = Vector3d.CrossProduct(a, b);
                n.Unitize();
            }
            return n;
        }

        public static Polyline[][] LoftPolylines(Polyline[][] polylines) {
            Polyline[][] p = new Polyline[polylines.Length][];

            for (int i = 0; i < polylines.Length; i++) {
                if (polylines[i].Length >= 2) {
                    if (polylines[i][0].Count == polylines[i][1].Count) {
                        p[i] = new Polyline[polylines[i][0].Count - 1];

                        for (int j = 0; j < polylines[i][0].Count - 1; j++)
                            p[i][j] = new Polyline() { polylines[i][0][j], polylines[i][0][j + 1], polylines[i][1][j + 1], polylines[i][1][j], polylines[i][0][j] };
                    } else
                        p[i] = new Polyline[] { new Polyline() };
                }
            }

            return p;
        }

        public static Polyline[] LoftTwoPolylines(Polyline[] polylines) {
            Polyline[] p = new Polyline[polylines[0].Count - 1];

            if (polylines.Length >= 2) {
                if (polylines[0].Count == polylines[1].Count)
                    for (int j = 0; j < polylines[0].Count - 1; j++)
                        p[j] = new Polyline() { polylines[0][j], polylines[0][j + 1], polylines[1][j + 1], polylines[1][j], polylines[0][j] };

                else
                    p = new Polyline[] { new Polyline() };
            }
            return p;
        }

        public static void ShiftPolyline(this Polyline A, int n) {

            A.RemoveAt(A.Count - 1);

            for (int j = 0; j < n; j++) {

                int len = A.Count; //self explanatory 
                var tmp = A[len - 1]; //save last element value
                for (int i = len - 1; i > 0; i--) //starting from the end to begining
                    A[i] = A[i - 1]; //assign value of the previous element
                A[0] = tmp; //now "rotate" last to first.
            }

            A.Add(A[0]);

        }

        public static void Close(this Polyline p) {
            if (p.Count > 2)
                p.Add(p[0]);
        }

        public static Polyline ToP(this Line line) {
            return new Polyline(new[] { line.From, line.To });
        }

        public static Curve[] ToCurveArray(this Polyline[] polyline) {
            Curve[] curves = new Curve[polyline.Length];

            for (int i = 0; i < polyline.Length; i++)
                curves[i] = polyline[i].ToNurbsCurve();

            return curves;
        }

        private static Random random = new Random();

        public static Polyline ExpandDuplicatedPoints(this Polyline polyline, double tolerance = 1E-5) {
            Polyline poly = new Polyline(polyline.ToArray());


            Plane.FitPlaneToPoints(poly, out Plane plane);

            for (int i = 1; i < poly.Count - 1; i += 2) {
                double dist = poly[i].DistanceTo(poly[i + 1]);
                double dist2 = poly[i].DistanceTo(poly[MathUtil.Wrap(i - 1, poly.Count)]);
                if (dist < tolerance) {
                    Point3d pt = poly[i];
                    double t = poly.ClosestParameter(pt);
                    poly[i] = poly.PointAt(t + 0.00001);
                }

                if (dist2 < tolerance) {
                    Point3d pt = poly[i];
                    double t = poly.ClosestParameter(pt);
                    poly[i] = poly.PointAt(t - 0.00001);
                }
            }

            return poly;

        }




        public static Polyline divideByCount(Curve curve, int n, bool close) {

            curve.Domain = new Interval(0, 1);

            if (n < 1)
                return null;

            switch (n) {
                case (1):
                if (curve.IsClosed)
                    return new Polyline(new[] { curve.PointAt(0.0), curve.PointAt(0.5) });
                else
                    return new Polyline(new[] { curve.PointAt(0.0), curve.PointAt(1.0) });

                case (2):
                if (curve.IsClosed)
                    return new Polyline(new[] { curve.PointAt(0.0), curve.PointAt(0.33333), curve.PointAt(0.66666) });
                else
                    return new Polyline(new[] { curve.PointAt(0.0), curve.PointAt(0.5), curve.PointAt(1.0) });

                default:

                Point3d[] pts;
                curve.DivideByCount(n, true, out pts);

                Polyline polyline = new Polyline(pts);
                if (close)
                    polyline.Add(polyline[0]);

                return polyline;
            }


        }





    }
}
