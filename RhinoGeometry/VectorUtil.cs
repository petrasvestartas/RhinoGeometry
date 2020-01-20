using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace RhinoGeometry {
    public static class VectorUtil {

        public static void ChangeEnd(this ref Line line, int i, Point3d p) {
            if (i == 0) {
                line.From = p;
            }else
                line.To = p;
        }

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

        public static Plane ProjectPlane(Plane planeToProject, Plane plane, Vector3d direction) {

            //Intersec axis lines with plane
            Line ProjectXAxis = new Line(planeToProject.Origin + planeToProject.XAxis*100, planeToProject.Origin + planeToProject.XAxis * 100 + direction);
            Line ProjectYAxis = new Line(planeToProject.Origin + planeToProject.YAxis * 100, planeToProject.Origin + planeToProject.YAxis * 100 + direction);

            Point3d pX = PlaneUtil.LinePlane(ProjectXAxis, plane);
            Point3d pY = PlaneUtil.LinePlane(ProjectYAxis, plane);

            Vector3d ProjectedXAxis = (pX - plane.Origin).UnitVector();
            Vector3d ProjectedYAxis = (pY - plane.Origin).UnitVector();
            Plane planeProjected = new Plane(plane.Origin, ProjectedXAxis, ProjectedYAxis);

            return planeProjected; 
        }

        public static bool IsVectorReversable(Vector3d v0, Vector3d v1) {

            //Reverse a vector if vectors are antiparallel
            int parallel = v0.IsParallelTo(v1);
            if (parallel == -1)
                return true;


            //Angle between vectors is larger than 180 degrees
            if (parallel == 0) {
                double angle = Vector3d.VectorAngle(v0, v1, Vector3d.CrossProduct(v0, v1));

                if (angle > Math.PI * 0.5)
                    return true;
            }

            return false;

        }

        public static Vector3d BiVector(Vector3d v0, Vector3d v1) {

            Vector3d va = new Vector3d(v0);
            Vector3d vb = new Vector3d(v1);

            if (IsVectorReversable(va, vb))
                vb.Reverse();


            //Average vector
            Vector3d vbi = va + vb;
            vbi.Unitize();
            return vbi;
        }



            public static Vector3d UnitVector(this Vector3d vec) {
                Vector3d v = new Vector3d(vec);
                v.Unitize();
                return v;
            }

            /// <summary>
            /// Vectors must follow each other
            /// </summary>
            /// <param name="V0"></param>
            /// <param name="V1"></param>
            /// <param name="Z"></param>
            /// <returns></returns>
            public static Vector3d BisectorVector(Vector3d V0, Vector3d V1, Vector3d Z, bool followEachOther = true) {
                Vector3d v0 = new Vector3d(V0);
                Vector3d v1 = new Vector3d(V1);

                v0.Unitize();
                v1.Unitize();

                if (v0.IsParallelTo(v1) == -1)
                    v1.Reverse();

                //Incase vector are consequtive or starting from the same point
                if (followEachOther) {
                    v0 += v1;
                } else {
                    v0 -= v1;
                }

                v0.Rotate(Math.PI * 0.5, Z);
                v0.Unitize();
                return v0;
            }


            public static Vector3d BisectorVector(Line V0, Line V1, Vector3d Z) {

                Line V2 = V1;
                if (V0.From.DistanceToSquared(V1.From) < 0.001 && V0.To.DistanceToSquared(V1.To) < 0.001)
                    V2.Flip();



                bool flip = V0.From.DistanceToSquared(V2.From) < 0.001 || V0.To.DistanceToSquared(V2.To) < 0.001;


                return BisectorVector(V0.Direction, V2.Direction, Z, !flip);
            }

            public static Vector3d BisectorVector(Line V0, Vector3d Z) {

                return BisectorVector(V0, V0, Z);
            }

            /// <summary>
            /// line must follow one after another
            /// </summary>
            /// <param name="V"></param>
            /// <param name="Z"></param>
            /// <returns></returns>
            public static Vector3d BisectorVector(List<Line> V, Vector3d Z, bool order = true) {


                if (V.Count == 1)
                    return BisectorVector(V[0], Z);

                List<Line> VOrdered = new List<Line>(V.Count);

                if (order) {

                    VOrdered.Add(V[0]);


                    for (int i = 1; i < V.Count; i++) {
                        Line l0 = VOrdered[i - 1];
                        Line l1 = V[i];

                        if (l0.From.DistanceToSquared(l1.From) < 0.001 || l0.To.DistanceToSquared(l1.To) < 0.001)
                            l1.Flip();

                        VOrdered.Add(l1);

                    }
                } else {
                    VOrdered = V;
                }


                Vector3d sumBisector = Vector3d.Zero;

                for (int i = 0; i < VOrdered.Count - 1; i++) {
                    Vector3d bisector = BisectorVector(VOrdered[i], VOrdered[i + 1], Z);
                    sumBisector += bisector;
                }
                sumBisector.Unitize();
                return sumBisector;


            }



        }
    }



