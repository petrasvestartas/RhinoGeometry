using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Rhino.Geometry;

namespace RhinoGeometry {
    public static class MeshUtil {

        public static Mesh Clean(this Mesh mesh) {

            Mesh mesh_ = mesh.DuplicateMesh();
            mesh_.Compact();
            mesh_.Vertices.CombineIdentical(true, true);
            mesh_.Vertices.CullUnused();
            mesh_.UnifyNormals();
            mesh_.Weld(3.14159265358979);
            mesh_.FaceNormals.ComputeFaceNormals();
            mesh_.Normals.ComputeNormals();

            if (mesh_.SolidOrientation() == -1)
                mesh_.Flip(true, true, true);
            return mesh_;
        }

        public static  Mesh MeshPipe(Polyline x, List<double> radiusList = null) {

            Curve c0 = x.ToNurbsCurve();

            var values = RhinoGeometry.MathUtil.Range(0, x.Count - 1, x.Count - 1);

            //Intrepolation values to define radius
            var radius = new List<double>(values.Count());
            var interpolationValues = RhinoGeometry.MathUtil.Range(0, 1, x.Count - 1);

            foreach (double t in interpolationValues) {
                radius.Add(RhinoGeometry.MathUtil.Interpolate(radiusList, t, RhinoGeometry.Interpolation.Cubic));
            }

            //Create polygons
            Plane[] planes = c0.GetPerpendicularFrames(values);
            Polyline[] polygons = new Polyline[planes.Length];
            Mesh mesh = new Mesh();
            int n = 10;
            for (int i = 0; i < planes.Length; i++) {
                Polyline polyline = RhinoGeometry.PolylineUtil.Polygon(n, radius[i], planes[i], 0, false);
                for (int j = 0; j < polyline.Count - 1; j++) {
                    mesh.Vertices.Add(polyline[j]);
                }
            }

            for (int i = 1; i < planes.Length; i++) {
                for (int j = 0; j < n; j++) {
                    int sign = j == n - 1 ? 0 : 1;
                    int a = (i - 1) * n + j;
                    int b = (i - 1) * n + (j + 1) * sign;
                    int c = (i * n) + (j + 1) * sign;
                    int d = i * n + j;
                    mesh.Faces.AddFace(a, b, c, d);
                }
            }
            mesh.FillHoles();
            mesh=mesh.Clean();
            return mesh;

        }


        public static Mesh MeshPipe(this Polyline x,  double[] radiusArray = null, int n = 10, double radiusDefault = 30, bool fillHoles = true) {

            Curve c0 = x.ToNurbsCurve();
            var values = RhinoGeometry.MathUtil.Range(0, x.Count - 1, x.Count - 1);

            double[] radius = new double[values.Count()];

            if (radiusArray != null)
                radius = radiusArray;
            else
                radius = Enumerable.Repeat(radiusDefault, values.Count()).ToArray();


            //Create polygons
            Plane[] planes = c0.GetPerpendicularFrames(values);
            Polyline[] polygons = new Polyline[planes.Length];
            Mesh mesh = new Mesh();
            //int n = 10;
            for (int i = 0; i < planes.Length; i++) {
                Polyline polyline = RhinoGeometry.PolylineUtil.Polygon(n, radius[i], planes[i], 0, false);
                for (int j = 0; j < polyline.Count - 1; j++) {
                    mesh.Vertices.Add(polyline[j]);
                }
            }

            for (int i = 1; i < planes.Length; i++) {
                for (int j = 0; j < n; j++) {
                    int sign = j == n - 1 ? 0 : 1;
                    int a = (i - 1) * n + j;
                    int b = (i - 1) * n + (j + 1) * sign;
                    int c = (i * n) + (j + 1) * sign;
                    int d = i * n + j;
                    mesh.Faces.AddFace(a, b, c, d);
                }
            }
            mesh.FillHoles();
            return mesh;

        }

        public static Mesh LoftMeshFast(this Polyline C0, Polyline C1, bool B) {

            //Mesh mesh = new Mesh();

            //Mesh mesh2 = new Mesh();


            int n = C1.Count - 1;
            bool closed = C0.IsClosed && C1.IsClosed;

            //Create triangulated mesh from closed polyline
            Mesh loft = new Mesh();
            if (closed && B)
                loft = Mesh.CreateFromClosedPolyline(C0);
            else
                for (int j = 0; j < n; j++)
                    loft.Vertices.Add(C0[j]);


            //Add bottom vertices
            for (int j = 0; j < n; j++) {
                loft.Vertices.Add(C1[j]);
            }




            //Add bottom faces

            if (closed && B) {
                MeshFace[] mf = new MeshFace[loft.Faces.Count];
                for (int j = 0; j < loft.Faces.Count; j++) {

                    int a = loft.Faces[j].A + n;
                    int b = loft.Faces[j].B + n;
                    int c = loft.Faces[j].C + n;

                    mf[j] = new MeshFace(c, b, a);
                }

                loft.Faces.AddFaces(mf);
            }


            for (int j = 0; j < n; j++) {



                MeshFace mf_side = new MeshFace(
                  j,
                  j + n,
                  n + (j + 1) % (n),
                  (j + 1) % (n)

                  );
                loft.Faces.AddFace(mf_side);
            }

            //loft.Ngons.AddNgon(MeshNgon.Create(Enumerable.Range(0,C0.Branch(i)[0].Count-1).ToArray(),triangulatedPolyline.Item2));


            //if (B) {
            //    mesh.Append(loft);
            //} else {

            //    loft.RebuildNormals();

            //    if (loft.SolidOrientation() == -1)
            //        loft.Flip(true, true, true);
            //    loft.Unweld(0.01, true);

            //    mesh2 = loft;
            //}

            //Rhino.RhinoDoc.ActiveDoc.Objects.AddMesh(loft);


            //if (B)
            //    return loft;



            return loft;
        }

    }
}
