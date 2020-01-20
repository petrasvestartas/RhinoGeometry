using System;
using System.Collections.Generic;
using System.Linq;
using Rhino.Geometry;
using System.Drawing;

namespace RhinoGeometry {

        public static class MathUtil {


        public static double[] CartesianToCylindrical (double X, double Y, double Z) {
            double[] numArray = new double[] { Math.Sqrt(X * X + Y * Y), Math.Atan2(Y, Z), Z };
            if (X == Y) {
                numArray[0] = 0;
            }
            return numArray;
        }

        public static double[] CartesianToSpherical(double X, double Y, double Z) {
            double[] numArray;
            double[] z = new double[3];
            if (!(X == 0 & Y == 0 & Z == 0)) {
                if (!(X == 0 & Y == 0)) {
                    z[0] = Math.Sqrt(X * X + Y * Y + Z * Z);
                    z[1] = Math.Atan2(Y, X);
                } else {
                    z[0] = Z;
                    z[1] = 0;
                }
                z[2] = Math.Asin(Z / z[0]);
                numArray = z;
            } else {
                z[0] = 0;
                z[1] = 0;
                z[2] = 0;
                numArray = z;
            }
            return numArray;
        }

        public static double[] CylindricalToCartesian(double Radius, double Angle, double Height) {
            return new double[] { Radius * Math.Cos(Angle), Radius * Math.Sin(Angle), Height };
        }

        public static double GPSTime(DateTime time) {
            DateTime dateTime = new DateTime(1980, 1, 6, 0, 0, 0);
            return (time - dateTime).TotalSeconds;
        }

        public static DateTime GPSTime(double time) {
            DateTime dateTime = new DateTime(1980, 1, 6, 0, 0, 0);
            return dateTime.AddSeconds(time);
        }

        public static double[] SphericalToCartesian(double Range, double Azimuth, double Elevation) {
            double[] range = new double[3];
            double num = Math.Cos(Elevation);
            range[0] = Range * (num * Math.Cos(Azimuth));
            range[1] = Range * (num * Math.Sin(Azimuth));
            range[2] = Range * Math.Sin(Elevation);
            return range;
        }
        public static double FastDistCheck(Point3d A, Point3d B) {
            return (A.X - B.X) * (A.X - B.X) + (A.Y - B.Y) * (A.Y - B.Y) + (A.Z - B.Z) * (A.Z - B.Z);
        }

        public static double FastPlaneToPt(double Denom, double a, double b, double c, double d, Point3d Pt) {
            return (a * Pt.X + b * Pt.Y + c * Pt.Z + d) * Denom;
        }

        public static bool FastPtAbovePlane(double a, double b, double c, double d, Point3d pt) {
            bool flag;
            flag = (pt.X * a + pt.Y * b + pt.Z * c + d < 0 ? false : true);
            return flag;
        }

        public static bool IsInBox(Point3d Point, Box TestBox) {
            return TestBox.Contains(Point, false);
        }

        public static bool IsInSphere(Point3d Point, Point3d Center, double Radius) {
            bool flag;
            flag = (FastDistCheck(Point, Center) > Radius * Radius ? false : true);
            return flag;
        }

        public static Quaternion PlaneToQuaternion(Plane P) {
            double m21;
            double m02;
            double m01;
            double m10;
            double num;
            P.Origin=(Point3d.Origin);
            Transform plane = Transform.PlaneToPlane(Plane.WorldXY, P);
            double m00 = plane.M00 + plane.M11 + plane.M22;
            if (m00 > 0) {
                num = Math.Sqrt(m00 + 1) * 2;
                m21 = 0.25 * num;
                m02 = (plane.M21 - plane.M12) / num;
                m01 = (plane.M02 - plane.M20) / num;
                m10 = (plane.M10 - plane.M01) / num;
            } else if (plane.M00 > plane.M11 & plane.M00 > plane.M22) {
                num = Math.Sqrt(1 + plane.M00 - plane.M11 - plane.M22) * 2;
                m21 = (plane.M21 - plane.M12) / num;
                m02 = 0.25 * num;
                m01 = (plane.M01 + plane.M10) / num;
                m10 = (plane.M02 + plane.M20) / num;
            } else if (plane.M11 <= plane.M22) {
                num = Math.Sqrt(1 + plane.M22 - plane.M00 - plane.M11) * 2;
                m21 = (plane.M10 - plane.M01) / num;
                m02 = (plane.M02 + plane.M20) / num;
                m01 = (plane.M12 + plane.M21) / num;
                m10 = 0.25 * num;
            } else {
                num = Math.Sqrt(1 + plane.M11 - plane.M00 - plane.M22) * 2;
                m21 = (plane.M02 - plane.M20) / num;
                m02 = (plane.M01 + plane.M10) / num;
                m01 = 0.25 * num;
                m10 = (plane.M12 + plane.M21) / num;
            }
            Quaternion quaternion = new Quaternion(m21, m02, m01, m10);
            quaternion.Unitize();
            return quaternion;
        }

        public static Point3d SphericalPoint(Point3d Point, Plane Origin) {
            Transform plane = Transform.PlaneToPlane(Origin, Plane.WorldXY);
            Point.Transform(plane);
            double num = Math.Sqrt(Math.Pow(Point.X, 2) + Math.Pow(Point.Y, 2) + Math.Pow(Point.Z, 2));
            double num1 = Math.Acos(Point.Z / num);
            double num2 = Math.Atan(Point.Y / Point.X);
            return new Point3d(num, num1, num2);
        }




        public static long Key(int key1, int key2) {
                return (((long)key1) << 32) + key2;
            }

            public static T next<T>(this T[] array, int i) {
                return array[Wrap(i + 1, array.Length)];
            }

            public static T prev<T>(this T[] array, int i) {
                return array[Wrap(i - 1, array.Length)];
            }

            public static int Next(this int[] array, int i) {
                return array[Wrap(i + 1, array.Length)];
            }
            public static int Prev(this int[] array, int i) {
                return array[Wrap(i - 1, array.Length)];
            }

            /// <summary>
            /// e.x. 0 1 2 4 5 -> 4 5 0 1 2
            /// </summary>
            /// <param name="NotOrderedEdges_"></param>
            /// <param name="n"></param>
            /// <returns></returns>
            public static List<int> SortIntegers(List<int> NotOrderedEdges_, int n) {

                List<int> orderedEdges = new List<int>();
                List<int> NotOrderedEdges = NotOrderedEdges_;
                NotOrderedEdges.Sort();


                int k = -1;
                for (int j = 0; j < NotOrderedEdges.Count - 1; j++) {
                    if (NotOrderedEdges[j] == 0 && NotOrderedEdges[j + 1] == n - 1)
                        continue;
                    if (Math.Abs(NotOrderedEdges[j] - NotOrderedEdges[j + 1]) > 1) {
                        k = j + 1;
                        break;
                    }
                }


                if (k != -1) {
                    for (int j = k; j < NotOrderedEdges.Count + k; j++) {
                        int cur = j % NotOrderedEdges.Count;
                        orderedEdges.Add(NotOrderedEdges[cur]);
                    }
                } else {
                    return NotOrderedEdges_;
                }

                return orderedEdges;

            }

            public static long GetKey(int i, int j) {
                return (UInt32)i << 16 | (UInt32)j;
            }
            /// <summary>
            /// always returns array of divided interval of length n+1
            /// </summary>
            /// <param name="IntervalStart"></param>
            /// <param name="IntervalEnd"></param>
            /// <param name="n"></param>
            /// <returns></returns>
            public static double[] Range(double IntervalStart, double IntervalEnd, int n) {

                Rhino.Geometry.Interval interval = new Rhino.Geometry.Interval(IntervalStart, IntervalEnd);
                double[] tInterval = new double[n + 1];

                for (int i = 0; i <= n; i++)
                    tInterval[i] = interval.ParameterAt((double)i / (double)n);

                return tInterval;
            }

            public static double Pitagoras(double A, double B) {
                return Math.Sqrt(A * A + B * B);
            }

            public static T Clamp<T>(this T val, T min, T max) where T : IComparable<T> {
                if (val.CompareTo(min) < 0)
                    return min;
                else if (val.CompareTo(max) > 0)
                    return max;
                else
                    return val;
            }

            public static double Map(double s, double sourceMin, double sourceMax, double targetMin, double targetMax) {
                return (targetMin + (s - sourceMin) * (targetMax - targetMin) / (sourceMax - sourceMin));
            }


            public static int MapInt(double s, double sourceMin, double sourceMax, double targetMin, double targetMax) {
                return (int)(targetMin + (s - sourceMin) * (targetMax - targetMin) / (sourceMax - sourceMin));
            }

            public static double Constrain(double x, double min, double max) {
                return Math.Max(min, Math.Min(x, max));
            }



            public static int Wrap(this int index, int n) {
                if (n == 0)
                    return 0;
                return ((index % n) + n) % n;
            }

            public static float Lerp(float value1, float value2, float amount) {
                return value1 + (value2 - value1) * amount;
            }

            public static double Lerp(double value1, double value2, double amount) {
                return value1 + (value2 - value1) * amount;
            }


            public static Rhino.Geometry.Point3d Lerp(Rhino.Geometry.Point3d p0, Rhino.Geometry.Point3d p1, double amount) {
                return new Rhino.Geometry.Point3d(
                    Lerp(p0.X, p1.X, amount),
                    Lerp(p0.Y, p1.Y, amount),
                    Lerp(p0.Z, p1.Z, amount)
                    );
            }

            public static Rhino.Geometry.Point3f Lerp(Rhino.Geometry.Point3f p0, Rhino.Geometry.Point3f p1, double amount) {
                return new Rhino.Geometry.Point3f(
                    Lerp(p0.X, p1.X, (float)amount),
                    Lerp(p0.Y, p1.Y, (float)amount),
                    Lerp(p0.Z, p1.Z, (float)amount)
                    );
            }

            public static double RemapNumbers(double speed = 50, double low1 = 0, double high1 = 2000, double low2 = 0, double high2 = 1) {
                return low2 + (speed - low1) * (high2 - low2) / (high1 - low1);
            }


            public static IEnumerable<double> Range2(double start, double end, int count) {
                return Enumerable.Range(0, count).Select(i => start + (end - start) * ((double)i / (count - 1)));
            }

            public static IEnumerable<double> RangeDouble(double from, double to, double step) {
                if (step <= 0.0)
                    step = (step == 0.0) ? 1.0 : -step;

                if (from <= to) {
                    for (double d = from; d <= to; d += step)
                        yield return d;
                } else {
                    for (double d = from; d >= to; d -= step)
                        yield return d;
                }
            }

            public static List<int> CommonItems(IEnumerable<int> CollectionA, IEnumerable<int> CollectionB) {
                List<int> nums = new List<int>();
                int num = checked(CollectionB.Count<int>() - 1);
                for (int i = 0; i <= num; i = checked(i + 1)) {
                    if (CollectionA.Contains<int>(CollectionB.ElementAtOrDefault<int>(i))) {
                        nums.Add(CollectionB.ElementAtOrDefault<int>(i));
                    }
                }
                return nums;
            }

            public static double Cot(double Radians) {
                return 1 / Math.Tan(Radians);
            }

            public static int[] CreateArray(int Count, int Value) {
                int[] value = new int[checked(checked(Count - 1) + 1)];
                int count = checked(Count - 1);
                for (int i = 0; i <= count; i = checked(i + 1)) {
                    value[i] = Value;
                }
                return value;
            }


            public static int[] CreateSeries(int Count) {
                int[] numArray = new int[checked(checked(Count - 1) + 1)];
                int count = checked(Count - 1);
                for (int i = 0; i <= count; i = checked(i + 1)) {
                    numArray[i] = i;
                }
                return numArray;
            }






        public static Color Col(double a, double r, double g, double b) {
            int num = Math.Min(Math.Max(Convert.ToInt32(a), 0), 255);
            int num1 = Math.Min(Math.Max(Convert.ToInt32(r), 0), 255);
            int num2 = Math.Min(Math.Max(Convert.ToInt32(g), 0), 255);
            int num3 = Math.Min(Math.Max(Convert.ToInt32(b), 0), 255);
            return Color.FromArgb(num, num1, num2, num3);
        }

        public static double SampleBlock(double v0, double v1, double t) {
            return (t > 0.5 ? v1 : v0);
        }

        public static double SampleCatmull(double v00, double v0, double v1, double v11, double t) {
            double num = v00 - 2.5 * v0 + 2 * v1 - 0.5 * v11;
            double num1 = -0.5 * v00 + 0.5 * v1;
            double num2 = v0;
            return (-0.5 * v00 + 1.5 * v0 - 1.5 * v1 + 0.5 * v11) * (t * t * t) + num * (t * t) + num1 * t + num2;
        }

        public static double SampleCubic(double v00, double v0, double v1, double v11, double t) {
            double num = v11 - v1 - v00 + v0;
            double num1 = v00 - v0 - num;
            double num2 = v1 - v00;
            double num3 = v0;
            return num * (t * t * t) + num1 * (t * t) + num2 * t + num3;
        }

        public static double SampleLinear(double v0, double v1, double t)  {
            return v0 * (1 - t) + v1 * t;
        }
        public static double Sum(int a , int b) {
            return a + b;
        }

        public static double Interpolate(List<double> values, double parameter, Interpolation m_algorithm = Interpolation.Linear) {
            double item;
            double num;
            double num1;
            double num2;
            double num3;
            int count = checked(values.Count - 1);
            double num4 = parameter * (double)count;
            int num5 = Convert.ToInt32(Math.Floor(num4));
            if (num5 < 0) {
                item = values[0];
            } else if (num5 < count) {
                switch (m_algorithm) {
                    case Interpolation.Block: {
                            double item1 = values[num5];
                            double item2 = values[checked(num5 + 1)];
                            item = SampleBlock(item1, item2, num4 - (double)num5);
                            break;
                        }
                    case Interpolation.Linear: {
                            double item3 = values[num5];
                            double item4 = values[checked(num5 + 1)];
                            item = SampleLinear(item3, item4, num4 - (double)num5);
                            break;
                        }
                    case Interpolation.Cubic: {
                            double item5 = values[num5];
                            double item6 = values[checked(num5 + 1)];
                            num = (num5 > 0 ? values[checked(num5 - 1)] : item5 - (item6 - item5));
                            num1 = (num5 < checked(count - 1) ? values[checked(num5 + 2)] : item6 - (item5 - item6));
                            item = SampleCubic(num, item5, item6, num1, num4 - (double)num5);
                            break;
                        }
                    case Interpolation.Catmull: {
                            double num6 = values[num5];
                            double item7 = values[checked(num5 + 1)];
                            num2 = (num5 > 0 ? values[checked(num5 - 1)] : num6 - (item7 - num6));
                            num3 = (num5 < checked(count - 1) ? values[checked(num5 + 2)] : item7 - (num6 - item7));
                            item = SampleCatmull(num2, num6, item7, num3, num4 - (double)num5);
                            break;
                        }
                    default: {
                            item = 0;
                            break;
                        }
                }
            } else {
                item = values[count];
            }
            return item;
        }

        public static int Interpolate(List<int> values, double parameter, Interpolation m_algorithm = Interpolation.Linear) {
            int item;
            double num;
            double num1;
            double num2;
            double num3;
            int count = checked(values.Count - 1);
            double num4 = parameter * (double)count;
            int num5 = Convert.ToInt32(Math.Floor(num4));
            if (num5 < 0) {
                item = values[0];
            } else if (num5 < count) {
                switch (m_algorithm) {
                    case Interpolation.Block: {
                            double item1 = (double)values[num5];
                            double item2 = (double)values[checked(num5 + 1)];
                            item = Convert.ToInt32(SampleBlock(item1, item2, num4 - (double)num5));
                            break;
                        }
                    case Interpolation.Linear: {
                            double item3 = (double)values[num5];
                            double item4 = (double)values[checked(num5 + 1)];
                            item = Convert.ToInt32(SampleLinear(item3, item4, num4 - (double)num5));
                            break;
                        }
                    case Interpolation.Cubic: {
                            double item5 = (double)values[num5];
                            double item6 = (double)values[checked(num5 + 1)];
                            num = (num5 > 0 ? (double)values[checked(num5 - 1)] : item5 - (item6 - item5));
                            num1 = (num5 < checked(count - 1) ? (double)values[checked(num5 + 2)] : item6 - (item5 - item6));
                            item = Convert.ToInt32(SampleCubic(num, item5, item6, num1, num4 - (double)num5));
                            break;
                        }
                    case Interpolation.Catmull: {
                            double num6 = (double)values[num5];
                            double item7 = (double)values[checked(num5 + 1)];
                            num2 = (num5 > 0 ? (double)values[checked(num5 - 1)] : num6 - (item7 - num6));
                            num3 = (num5 < checked(count - 1) ? (double)values[checked(num5 + 2)] : item7 - (num6 - item7));
                            item = Convert.ToInt32(SampleCatmull(num2, num6, item7, num3, num4 - (double)num5));
                            break;
                        }
                    default: {
                            item = 0;
                            break;
                        }
                }
            } else {
                item = values[count];
            }
            return item;
        }

        public static Point3d Interpolate(List<Point3d> values, double parameter, Interpolation m_algorithm = Interpolation.Linear) {

            Point3d item;
            Point3d point3d;
            Point3d point3d1;
            Point3d point3d2;
            Point3d point3d3;
            int count = checked(values.Count - 1);
            double num = parameter * (double)count;
            int num1 = Convert.ToInt32(Math.Floor(num));
            if (num1 < 0) {
                item = values[0];
            } else if (num1 < count) {
                switch (m_algorithm) {
                    case Interpolation.Block: {
                            Point3d item1 = values[num1];
                            Point3d item2 = values[checked(num1 + 1)];
                            double num2 = SampleBlock(item1.X, item2.X, num - (double)num1);
                            double num3 = SampleBlock(item1.Y, item2.Y, num - (double)num1);
                            double num4 = SampleBlock(item1.Z, item2.Z, num - (double)num1);
                            item = new Point3d(num2, num3, num4);
                            break;
                        }
                    case Interpolation.Linear: {
                            Point3d item3 = values[num1];
                            Point3d item4 = values[checked(num1 + 1)];
                            double num5 = SampleLinear(item3.X, item4.X, num - (double)num1);
                            double num6 = SampleLinear(item3.Y, item4.Y, num - (double)num1);
                            double num7 = SampleLinear(item3.Z, item4.Z, num - (double)num1);
                            item = new Point3d(num5, num6, num7);
                            break;
                        }
                    case Interpolation.Cubic: {
                            Point3d point3d4 = values[num1];
                            Point3d item5 = values[checked(num1 + 1)];
                            point3d = (num1 > 0 ? values[checked(num1 - 1)] : point3d4 - (item5 - point3d4));
                            point3d1 = (num1 < checked(count - 1) ? values[checked(num1 + 2)] : item5 - (point3d4 - item5));
                            double num8 = SampleCubic(point3d.X, point3d4.X, item5.X, point3d1.X, num - (double)num1);
                            double num9 = SampleCubic(point3d.Y, point3d4.Y, item5.Y, point3d1.Y, num - (double)num1);
                            double num10 = SampleCubic(point3d.Z, point3d4.Z, item5.Z, point3d1.Z, num - (double)num1);
                            item = new Point3d(num8, num9, num10);
                            break;
                        }
                    case Interpolation.Catmull: {
                            Point3d point3d5 = values[num1];
                            Point3d item6 = values[checked(num1 + 1)];
                            point3d2 = (num1 > 0 ? values[checked(num1 - 1)] : point3d5 - (item6 - point3d5));
                            point3d3 = (num1 < checked(count - 1) ? values[checked(num1 + 2)] : item6 - (point3d5 - item6));
                            double num11 = SampleCatmull(point3d2.X, point3d5.X, item6.X, point3d3.X, num - (double)num1);
                            double num12 = SampleCatmull(point3d2.Y, point3d5.Y, item6.Y, point3d3.Y, num - (double)num1);
                            double num13 = SampleCatmull(point3d2.Z, point3d5.Z, item6.Z, point3d3.Z, num - (double)num1);
                            item = new Point3d(num11, num12, num13);
                            break;
                        }
                    default: {
                            item = Point3d.Origin;
                            break;
                        }
                }
            } else {
                item = values[count];
            }
            return item;
        }

        public static Vector3d Interpolate(List<Vector3d> values, double parameter, Interpolation m_algorithm = Interpolation.Linear) {
            Vector3d item;
            Vector3d vector3d;
            Vector3d vector3d1;
            Vector3d vector3d2;
            Vector3d vector3d3;
            int count = checked(values.Count - 1);
            double num = parameter * (double)count;
            int num1 = Convert.ToInt32(Math.Floor(num));
            if (num1 < 0) {
                item = values[0];
            } else if (num1 < count) {
                switch (m_algorithm) {
                    case Interpolation.Block: {
                            Vector3d item1 = values[num1];
                            Vector3d item2 = values[checked(num1 + 1)];
                            double num2 = SampleBlock(item1.X, item2.X, num - (double)num1);
                            double num3 = SampleBlock(item1.Y, item2.Y, num - (double)num1);
                            double num4 = SampleBlock(item1.Z, item2.Z, num - (double)num1);
                            item = new Vector3d(num2, num3, num4);
                            break;
                        }
                    case Interpolation.Linear: {
                            Vector3d item3 = values[num1];
                            Vector3d item4 = values[checked(num1 + 1)];
                            double num5 = SampleLinear(item3.X, item4.X, num - (double)num1);
                            double num6 = SampleLinear(item3.Y, item4.Y, num - (double)num1);
                            double num7 = SampleLinear(item3.Z, item4.Z, num - (double)num1);
                            item = new Vector3d(num5, num6, num7);
                            break;
                        }
                    case Interpolation.Cubic: {
                            Vector3d vector3d4 = values[num1];
                            Vector3d item5 = values[checked(num1 + 1)];
                            vector3d = (num1 > 0 ? values[checked(num1 - 1)] : vector3d4 - (item5 - vector3d4));
                            vector3d1 = (num1 < checked(count - 1) ? values[checked(num1 + 2)] : item5 - (vector3d4 - item5));
                            double num8 = SampleCubic(vector3d.X, vector3d4.X, item5.X, vector3d1.X, num - (double)num1);
                            double num9 = SampleCubic(vector3d.Y, vector3d4.Y, item5.Y, vector3d1.Y, num - (double)num1);
                            double num10 = SampleCubic(vector3d.Z, vector3d4.Z, item5.Z, vector3d1.Z, num - (double)num1);
                            item = new Vector3d(num8, num9, num10);
                            break;
                        }
                    case Interpolation.Catmull: {
                            Vector3d vector3d5 = values[num1];
                            Vector3d item6 = values[checked(num1 + 1)];
                            vector3d2 = (num1 > 0 ? values[checked(num1 - 1)] : vector3d5 - (item6 - vector3d5));
                            vector3d3 = (num1 < checked(count - 1) ? values[checked(num1 + 2)] : item6 - (vector3d5 - item6));
                            double num11 = SampleCatmull(vector3d2.X, vector3d5.X, item6.X, vector3d3.X, num - (double)num1);
                            double num12 = SampleCatmull(vector3d2.Y, vector3d5.Y, item6.Y, vector3d3.Y, num - (double)num1);
                            double num13 = SampleCatmull(vector3d2.Z, vector3d5.Z, item6.Z, vector3d3.Z, num - (double)num1);
                            item = new Vector3d(num11, num12, num13);
                            break;
                        }
                    default: {
                            item = Vector3d.Zero;
                            break;
                        }
                }
            } else {
                item = values[count];
            }
            return item;
        }

        public static Interval Interpolate(List<Interval> values, double parameter, Interpolation m_algorithm = Interpolation.Linear) {
            Interval item;
            Interval interval;
            Interval interval1;
            Interval interval2;
            Interval interval3;
            int count = checked(values.Count - 1);
            double num = parameter * (double)count;
            int num1 = Convert.ToInt32(Math.Floor(num));
            if (num1 < 0) {
                item = values[0];
            } else if (num1 < count) {
                switch (m_algorithm) {
                    case Interpolation.Block: {
                            Interval item1 = values[num1];
                            Interval item2 = values[checked(num1 + 1)];
                            double num2 = SampleBlock(item1.T0, item2.T0, num - (double)num1);
                            double num3 = SampleBlock(item1.T1, item2.T1, num - (double)num1);
                            item = new Interval(num2, num3);
                            break;
                        }
                    case Interpolation.Linear: {
                            Interval item3 = values[num1];
                            Interval item4 = values[checked(num1 + 1)];
                            double num4 = SampleLinear(item3.T0, item4.T0, num - (double)num1);
                            double num5 = SampleLinear(item3.T1, item4.T1, num - (double)num1);
                            item = new Interval(num4, num5);
                            break;
                        }
                    case Interpolation.Cubic: {
                            Interval interval4 = values[num1];
                            Interval item5 = values[checked(num1 + 1)];
                            interval = (num1 > 0 ? values[checked(num1 - 1)] : new Interval(interval4.T0 - (item5.T0 - interval4.T0), interval4.T1 - (item5.T1 - interval4.T1)));
                            interval1 = (num1 < checked(count - 1) ? values[checked(num1 + 2)] : new Interval(item5.T0 - (interval4.T0 - item5.T0), item5.T1 - (interval4.T1 - item5.T1)));
                            double num6 = SampleCubic(interval.T0, interval4.T0, item5.T0, interval1.T0, num - (double)num1);
                            double num7 = SampleCubic(interval.T1, interval4.T1, item5.T1, interval1.T1, num - (double)num1);
                            item = new Interval(num6, num7);
                            break;
                        }
                    case Interpolation.Catmull: {
                            Interval interval5 = values[num1];
                            Interval item6 = values[checked(num1 + 1)];
                            interval2 = (num1 > 0 ? values[checked(num1 - 1)] : new Interval(interval5.T0 - (item6.T0 - interval5.T0), interval5.T1 - (item6.T1 - interval5.T1)));
                            interval3 = (num1 < checked(count - 1) ? values[checked(num1 + 2)] : new Interval(item6.T0 - (interval5.T0 - item6.T0), item6.T1 - (interval5.T1 - item6.T1)));
                            double num8 = SampleCatmull(interval2.T0, interval5.T0, item6.T0, interval3.T0, num - (double)num1);
                            double num9 = SampleCatmull(interval2.T1, interval5.T1, item6.T1, interval3.T1, num - (double)num1);
                            item = new Interval(num8, num9);
                            break;
                        }
                    default: {
                            item = new Interval(0, 0);
                            break;
                        }
                }
            } else {
                item = values[count];
            }
            return item;
        }

        public static Color Interpolate(List<Color> values, double parameter, Interpolation m_algorithm = Interpolation.Linear) {
            Color item;
            Color color;
            Color color1;
            Color color2;
            Color color3;
            int count = checked(values.Count - 1);
            double num = parameter * (double)count;
            int num1 = Convert.ToInt32(Math.Floor(num));
            if (num1 < 0) {
                item = values[0];
            } else if (num1 < count) {
                switch (m_algorithm) {
                    case Interpolation.Block: {
                            Color item1 = values[num1];
                            Color item2 = values[checked(num1 + 1)];
                            double num2 = SampleBlock((double)item1.A, (double)item2.A, num - (double)num1);
                            double num3 = SampleBlock((double)item1.R, (double)item2.R, num - (double)num1);
                            double num4 = SampleBlock((double)item1.G, (double)item2.G, num - (double)num1);
                            double num5 = SampleBlock((double)item1.B, (double)item2.B, num - (double)num1);
                            item = Col(num2, num3, num4, num5);
                            break;
                        }
                    case Interpolation.Linear: {
                            Color item3 = values[num1];
                            Color item4 = values[checked(num1 + 1)];
                            double num6 = SampleLinear((double)item3.A, (double)item4.A, num - (double)num1);
                            double num7 = SampleLinear((double)item3.R, (double)item4.R, num - (double)num1);
                            double num8 = SampleLinear((double)item3.G, (double)item4.G, num - (double)num1);
                            double num9 = SampleLinear((double)item3.B, (double)item4.B, num - (double)num1);
                            item = Col(num6, num7, num8, num9);
                            break;
                        }
                    case Interpolation.Cubic: {
                            Color color4 = values[num1];
                            Color item5 = values[checked(num1 + 1)];
                            color = (num1 > 0 ? values[checked(num1 - 1)] : color4);
                            color1 = (num1 < checked(count - 1) ? values[checked(num1 + 2)] : item5);
                            double num10 = SampleCubic((double)color.A, (double)color4.A, (double)item5.A, (double)color1.A, num - (double)num1);
                            double num11 = SampleCubic((double)color.R, (double)color4.R, (double)item5.R, (double)color1.R, num - (double)num1);
                            double num12 = SampleCubic((double)color.G, (double)color4.G, (double)item5.G, (double)color1.G, num - (double)num1);
                            double num13 = SampleCubic((double)color.B, (double)color4.B, (double)item5.B, (double)color1.B, num - (double)num1);
                            item = Col(num10, num11, num12, num13);
                            break;
                        }
                    case Interpolation.Catmull: {
                            Color color5 = values[num1];
                            Color item6 = values[checked(num1 + 1)];
                            color2 = (num1 > 0 ? values[checked(num1 - 1)] : color5);
                            color3 = (num1 < checked(count - 1) ? values[checked(num1 + 2)] : item6);
                            double num14 = SampleCatmull((double)color2.A, (double)color5.A, (double)item6.A, (double)color3.A, num - (double)num1);
                            double num15 = SampleCatmull((double)color2.R, (double)color5.R, (double)item6.R, (double)color3.R, num - (double)num1);
                            double num16 = SampleCatmull((double)color2.G, (double)color5.G, (double)item6.G, (double)color3.G, num - (double)num1);
                            double num17 = SampleCatmull((double)color2.B, (double)color5.B, (double)item6.B, (double)color3.B, num - (double)num1);
                            item = Col(num14, num15, num16, num17);
                            break;
                        }
                    default: {
                            item = Color.Transparent;
                            break;
                        }
                }
            } else {
                item = values[count];
            }
            return item;
        }

        public static Plane Interpolate(List<Plane> P, double T, Interpolation m_algorithm = Interpolation.Linear) {

            var pts = new List<Point3d>();
            var x = new List<Vector3d>();
            var y = new List<Vector3d>();

            foreach (Plane plane in P) {
                pts.Add(plane.Origin);
                x.Add(plane.XAxis);
                y.Add(plane.YAxis);
            }


            return new Plane(Interpolate(pts, T, m_algorithm), Interpolate(x, T, m_algorithm), Interpolate(y, T, m_algorithm));

        }



        /// <summary>
        /// use 
        /// var parameters = new List<double>() { 0, 0.25, 0.5, 1 };
        /// var values = new List<double>() { 10, 20, 100, 30 };
        /// double t = 0.49;
        /// double value = InterpolateWeightedValues(parameters, values, t);
        /// </summary>
        /// <param name="parameters"></param>
        /// <param name="values"></param>
        /// <param name="t"></param>
        /// <returns></returns>
        public static double InterpolateWeightedValues(List<double> parameters, List<double> values, double t) {

            double value = double.NaN;
            t = Math.Min(1, Math.Max(0, t));


            //Iterate until the end of interval is found
            for (int i = 0; i < parameters.Count; i++) {
                if (t <= parameters[i]) {

                    if (i == 0) {
                        value = values[0];
                        break;
                    } else {

                        //Remap t
                        double tRemapped = MathUtil.RemapNumbers(t, parameters[i - 1], parameters[i], 0, 1);

                        //Create interval and get value
                        Interval interval = new Interval(values[Math.Max(0, i - 1)], values[i]);
                        value = interval.ParameterAt(tRemapped);
                        break;
                    }//if 0
                }//if the end of interval
            }//for loop

            return value;
        }

    }


    //public Interpolation m_algorithm = Interpolation.Catmull;
    public enum Interpolation {
        Block = 1,
        Linear = 2,
        Cubic = 3,
        Catmull = 4
    }

}


