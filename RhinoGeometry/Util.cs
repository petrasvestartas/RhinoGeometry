using Grasshopper;
using Grasshopper.Kernel.Data;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace RhinoGeometry {
    public static class Util
    {

        //http://www.pointclouds.org/documentation/tutorials/normal_estimation.php
        //https://www.cloudcompare.org/doc/wiki/index.php?title=Normals%5CCompute


        

        public static List<Vector3d> LineCloudNormals(DataTree<Line> lines, DataTree<int> adj) {

            List<Vector3d> _Normals = new List<Vector3d>();

            //Compute bounding box of cloud 
            BoundingBox bbox = new BoundingBox();

            foreach (Line l in lines.AllData()) {
                bbox.Union(l.From);
                bbox.Union(l.To);
            }

            bbox.Inflate(100);
            Point3d p = bbox.PointAt(0.5, 0.5, 1);




            //Iterate through joints




            //Fit to Plane

            //Add all lines point to a base datatree
            DataTree<Point3d> pts = new DataTree<Point3d>();
            for (int i = 0; i < lines.BranchCount; i++) {
                pts.Add(lines.Branch(i)[0].From,new GH_Path(i));
                pts.Add(lines.Branch(i)[0].To, new GH_Path(i));
            }

            //Add only second object point to first object path
            for (int i = 0; i < adj.BranchCount; i++) {

                int id0 = adj.Branch(i)[0];
                int id1 = adj.Branch(i)[1];

                pts.Add(lines.Branch(id1)[0].From, new GH_Path(id0));
                pts.Add(lines.Branch(id1)[0].To, new GH_Path(id0));
            }

            //Orient Plane
            for (int i = 0; i < pts.BranchCount; i++) {
                Plane.FitPlaneToPoints(pts.Branch(i), out Plane plane);
                plane.Origin = (pts.Branch(i)[0]+ pts.Branch(i)[1])*0.5;
                if ((plane.Origin + plane.Normal).DistanceToSquared(p) > (plane.Origin - plane.Normal).DistanceToSquared(p)) {
                    plane.Flip();
                }
                _Normals.Add(plane.Normal);
            }



                return _Normals;
        }



        public static long GetKey(int i, int j) {
            return (UInt32)i << 16 | (UInt32)j;
        }


        public static Tuple<Line, double, double > BrepPipe(Brep G) {
            //Get circles of conic pipes
            List<Circle> c = new List<Circle>();

            foreach (Curve curve in G.Curves3D) {
                Circle circle;
                if (curve.TryGetCircle(out circle, 0.01)) {
                    c.Add(circle);
                }
            }

            //If two cones are found
            if (c.Count == 2) {
                Line axis = new Line(c[0].Center, c[1].Center);
                return new Tuple<Line, double, double>(axis, c[0].Radius, c[1].Radius);
            }
            return new Tuple<Line, double, double>(Line.Unset, -1, -1);
        }

        public static Polyline Polygon(int n, double radius, Plane plane, double rotation = Math.PI * 0.25, bool sqrt = true)
        {

            Polyline polyline = new Polyline();
            double sector = Math.PI * 2 / n;
            double r = sqrt ? 1 / Math.Sqrt(2) * radius : radius;

            for (int i = 0; i < n; i++)
            {
                Point3d p = new Point3d(Math.Sin((sector * i) + rotation) * r, Math.Cos((sector * i) + rotation) * r, 0);
                polyline.Add(p);
            }

            polyline.Add(polyline[0]);

            polyline.Transform(Rhino.Geometry.Transform.PlaneToPlane(Plane.WorldXY, plane));

            return polyline;
        }

    }
}
