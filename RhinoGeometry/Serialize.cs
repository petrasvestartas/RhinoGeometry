using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace RhinoGeometry {
    public static class Serialize {

        public static double[][] PointToDouble(List<Point3d> pts, int size = 3) {

            double[][] points = new double[pts.Count][];

            for (int i = 0; i < pts.Count; i++) {
                if (size == 3) {
                    points[i] = new double[] { pts[i].X, pts[i].Y, pts[i].Z };
                } else if (size == 2) {
                    points[i] = new double[] { pts[i].X, pts[i].Y };
                } else if (size == 1) {
                    points[i] = new double[] { pts[i].X };
                }
            }
            return points;
        }

        public static List<Point3d> DoubleToPoints(double[][] numbers) {

            var pts = new List<Point3d>();

            for (int i = 0; i < numbers.Length; i++) {
                if (numbers[i].Length == 3) {
                    pts.Add(new Point3d(numbers[i][0], numbers[i][1], numbers[i][2]));
                } else if (numbers[i].Length == 2) {
                    pts.Add(new Point3d(numbers[i][0], numbers[i][1], 0));
                } else if (numbers[i].Length == 1) {
                    pts.Add(new Point3d(numbers[i][0], 0, 0));
                }
            }
            return pts;
        }
    }
}
