using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace RhinoGeometry {
    public static class CP {



        public class ClosestLinesT {
            public List<Line> cp = new List<Line>();
            public List<Line> L0 = new List<Line>();
            public List<Line> L1 = new List<Line>();
            public List<int> ID0 = new List<int>();
            public List<int> ID1 = new List<int>();
            public List<double> T0 = new List<double>();
            public List<double> T1 = new List<double>();

            public override string ToString() {
                return "ClosestLinesT" + " " + L0.Count.ToString();
            }

        }


        public static void ClosestLines(List<Line> L, double tolerance, double toleranceEnds, ref ClosestLinesT t) {


            t = new ClosestLinesT();

            HashSet<long> pairs = new HashSet<long>();

            for (int i = 0; i < L.Count; i++) {
                for (int j = 0; j < L.Count; j++) {

                    //In order to check the same pair again
                    if (i == j)
                        continue;

                    long key = (i<j) ?  Util.GetKey(i, j) : Util.GetKey(j, i);

                    if (pairs.Contains(key))
                        continue;

                    pairs.Add(key);

                    //Check order

                    bool checkEnds0 = L[i].From.DistanceToSquared(L[j].From) < toleranceEnds;
                    bool checkEnds1 = L[i].From.DistanceToSquared(L[j].To) < toleranceEnds;
                    double t0, t1;

                    if (checkEnds0) {
                        t.T0.Add(0);
                        t.T1.Add(0);
                        t.ID0.Add(-i);
                        t.ID1.Add(-j);
                        t.L0.Add(L[i]);
                        t.L1.Add(L[j]);
                    } else if (checkEnds1) {
                        t.T0.Add(0);
                        t.T1.Add(1);
                        t.ID0.Add(-i);
                        t.ID1.Add(-j);
                        t.L0.Add(L[i]);
                        t.L1.Add(L[j]);
                    } else if (Rhino.Geometry.Intersect.Intersection.LineLine(L[i], L[j], out t0, out t1, tolerance, true)) {
                        Line line = new Line(L[i].PointAt(t0), L[j].PointAt(t1));

                        //Identify how lines are connected
                        int EndSide0 = (t0 < toleranceEnds || t0 > 1 - toleranceEnds) ? -1 : 1;
                        int EndSide1 = (t1 < toleranceEnds || t1 > 1 - toleranceEnds) ? -1 : 1;

                        t.T0.Add(t0);
                        t.T1.Add(t1);
                        t.ID0.Add(i * EndSide0);
                        t.ID1.Add(j * EndSide1);
                        t.L0.Add(L[i]);
                        t.L1.Add(L[j]);



                        //Touching - Not Touching
                        //Print(EndSide0.ToString() + " " + EndSide1.ToString());
                    }


                }
            }


        }



        //public static List<int[]> ClosestLineSearchByGuide(List<Line> LinesToSearchFrom, List<GeometryBase> needles, double dist, ref object A) {

        //    List<Point3d> needles0 = new List<Point3d>();
        //    List<Point3d> needles1 = new List<Point3d>();

        //    //If not empty
        //    if (LinesToSearchFrom.Count == 0)
        //        return new List<int[]>();

        //    foreach(GeometryBase g in needles) {
        //        g.
        //    }


        //    return result;
        //}

        public static List<int[]> RTreeSearch(List<Point3d> pointsToSearchFrom, double dist, ref object A) {

            //If not empty
            if (pointsToSearchFrom.Count == 0)
                return new List<int[]>();


            //Search
            IEnumerable<int[]> found = RTree.Point3dClosestPoints(pointsToSearchFrom, pointsToSearchFrom, dist);
            HashSet<long> pairs = new HashSet<long>();

            List<int[]> result = new List<int[]>();

            int i = 0;
            foreach (var item in found) {
                int[] data = item;
                for (int j = 0; j < data.Length; ++j) {

                    int[] p = (i < data[j]) ? new int[] { i, data[j] } : new int[] { data[j], i };//sorts ids
                    long key = Util.GetKey(p[0], p[1]);//create key

                    if (!pairs.Contains(key)) {//add only if pair does not exist
                        result.Add(p);
                        pairs.Add(key);
                    }

                }
                i++;
            }

            return result;
        }

        public static List<Point3d> RTreeSearch(List<Point3d> pointsToSearchFrom, List<Point3d> needles, double dist, ref object A) {

            //If not empty
            if (pointsToSearchFrom.Count == 0 || needles.Count == 0)
                return new List<Point3d>();

            //Search
            IEnumerable<int[]> found = RTree.Point3dClosestPoints(pointsToSearchFrom, needles, dist);

            List<Point3d> result = new List<Point3d>();

            foreach (var item in found) {
                int[] data = item;
                for (int j = 0; j < data.Length; ++j)
                    result.Add(pointsToSearchFrom[data[j]]);
            }

            return result;
        }


        public static List<Point3d> RTreeSearch(List<Point3d> pointsToSearchFrom, List<Point3d> needles, int C, ref object A) {

            if (pointsToSearchFrom.Count == 0 || needles.Count == 0)
                return new List<Point3d>();

            IEnumerable<int[]> found = RTree.Point3dKNeighbors(pointsToSearchFrom, needles, C);

            List<Point3d> result = new List<Point3d>();

            foreach (int[] item in found) {
                int[] data = item;
                for (int j = 0; j < data.Length; ++j)
                    result.Add(pointsToSearchFrom[data[j]]);
            }

            return result;
        }





    }
}
