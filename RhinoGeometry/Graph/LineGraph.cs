using Grasshopper;
using Grasshopper.Kernel.Data;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace RhinoGeometry {
    public static class LineGraph {


        public static double Angle_2D(Vector3d A, Vector3d B, Plane P) //, out string str_error
        {
            string str_error = null;
            double num;
            Point3d origin = P.Origin + A;
            Point3d point3d = P.Origin + B;
            origin = P.ClosestPoint(origin);
            point3d = P.ClosestPoint(point3d);
            Vector3d vector3d = origin - P.Origin;
            Vector3d origin1 = point3d - P.Origin;
            if (!vector3d.Unitize()) {
                str_error = "First vector is zero-length after projection onto the Plane";
                num = double.NaN;
            } else if (origin1.Unitize()) {
                double num1 = Math.Acos(Math.Max(Math.Min(vector3d * origin1, 1), -1));
                if (Math.Abs(num1) < 1E-32) {
                    num = 0;
                } else if (Math.Abs(num1 - 3.14159265358979) >= 1E-32) {
                    Vector3d vector3d1 = Vector3d.CrossProduct(vector3d, origin1);
                    num = (P.ZAxis.IsParallelTo(vector3d1) != 1 ? 6.28318530717959 - num1 : num1);
                } else {
                    num = 3.14159265358979;
                }
            } else {
                str_error = "Second vector is zero-length after projection onto the Plane";
                num = double.NaN;
            }
            return num;
        }

        public class ClockwiseVector3dComparer : IComparer<Vector3d> {
            public int Compare(Vector3d v1, Vector3d v2) {
                if (v1.X >= 0) {
                    if (v2.X < 0) {
                        return -1;
                    }
                    return -Comparer<double>.Default.Compare(v1.Y, v2.Y);
                } else {
                    if (v2.X >= 0) {
                        return 1;
                    }
                    return Comparer<double>.Default.Compare(v1.Y, v2.Y);
                }
            }
        }


        public static Tuple<Point3d[], List<string>, List<int>, List<int>, List<int>, DataTree<int>> GetGraphData(List<Line> lines){



           UndirectedGraph g = LinesToUndirectedGrap(lines);

            DataTree<int> dataTreeA = new DataTree<int>();
            List<String> vertices = g.GetVertices();
            for (int i = 0; i < vertices.Count; i++)
                dataTreeA.AddRange(g.GetNeighbours(vertices[i]), new Grasshopper.Kernel.Data.GH_Path(i));

            DataTree<int> dataTreeB = new DataTree<int>();
            List<Edge> edges = g.GetEdges();

            List<int> u = new List<int>();
            List<int> v = new List<int>();
            List<int> w = new List<int>();
            for (int i = 0; i < edges.Count; i++) {
                u.Add(edges[i].u);
                v.Add(edges[i].v);
                w.Add(1);
                dataTreeB.Add(edges[i].u, new Grasshopper.Kernel.Data.GH_Path(i));
                dataTreeB.Add(edges[i].v, new Grasshopper.Kernel.Data.GH_Path(i));
            }

            PointCloud pointCloud = ((PointCloud)g.GetAttribute());


            return new Tuple<Point3d[], List<string>, List<int>, List<int>, List<int>, DataTree<int>>( pointCloud.GetPoints(), g.GetVertices(),u,v,w,dataTreeA );


        }


        private static  UndirectedGraph LinesToUndirectedGrap(List<Line> lines) {

            List<Point3d> pts = new List<Point3d>();

            foreach (Line l in lines) {
                pts.Add(l.From);
                pts.Add(l.To);
            }

            //Sorting
            var edges = new List<int>();

            var allPoints = new List<Point3d>(pts); //naked points

            int i = 0;

            while (allPoints.Count != 0) {
                Point3d pt = allPoints[0];
                allPoints.RemoveAt(0);


                for (int d = 0; d < pts.Count; d++) {
                    if (pt.Equals(pts[d])) {
                        edges.Add(d);
                        break;
                    }
                }

                i++;
            }

            var uniqueVertices = new HashSet<int>(edges).ToList();

            //Creating typological points
            var topologyPoints = new PointCloud();

            foreach (int k in uniqueVertices)
                topologyPoints.Add(pts[k]);

            //var vertices = Enumerable.Range(0, uniqueVertices.Count);

            for (int k = 0; k < uniqueVertices.Count; k++)
                if (uniqueVertices.ElementAt(k) != k)
                    for (int l = 0; l < edges.Count; l++)
                        if (edges[l] == uniqueVertices[k])
                            edges[l] = k;

            //Create graph
            UndirectedGraph g = new UndirectedGraph(uniqueVertices.Count);

            for (int k = 0; k < uniqueVertices.Count; k++)
                g.InsertVertex(k.ToString());


            for (int k = 0; k < edges.Count; k += 2)
                g.InsertEdge(edges[k].ToString(), edges[k + 1].ToString());

            g.SetAttribute((object)topologyPoints);


            return g;

        }


    }
}
