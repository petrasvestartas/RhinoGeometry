using Grasshopper;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Xml.Linq;

namespace RhinoGeometry {
    public static class XMLWriterReader {

        /// <summary>
        /// Write XML
        ///string settings_path = @"C:\Users\petra\AppData\Roaming\Grasshopper\6\Libraries\RhinoJoint\TileTypeSettings.xml";
        ///TilesToXML(male, female, type, settings_path);
        /// </summary>
        /// <param name="male"></param>
        /// <param name="female"></param>
        /// <param name="type"></param>
        /// <param name="settings_path"></param>
        public static void TilesToXML(DataTree<Polyline> male, DataTree<Polyline> female, DataTree<string> type, string settings_path) {
            System.Threading.Thread.CurrentThread.CurrentCulture = System.Globalization.CultureInfo.InvariantCulture;
            //string settings_path = System.IO.Path.Combine(System.IO.Path.GetDirectoryName(System.Reflection.Assembly.GetExecutingAssembly().Location), "TileTypeSettings.xml");

            var root = new XElement("Tile");

            for (int i = 0; i < type.BranchCount; i++) {
                if (male.Branch(i).Count == 0)
                    continue;
                var tileXML = PolylinesToXML(male.Branch(i), female.Branch(i), type.Branch(i)[0]);
                root.Add(tileXML);
            }

            //1. Add root to document
            XDocument doc = new XDocument(root);
            doc.Save(settings_path);

        }


        public static XElement PolylinesToXML(List<Polyline> male, List<Polyline> female, string type) {
            //3. Create Root


            var jointTypeE = new XElement("TileType");
            var name = new XElement("name");
            name.Value = type;
            var maleE = new XElement("male");
            var femaleE = new XElement("female");

            foreach (Polyline poly in male) {
                var polylineElement = new XElement("Polyline");
                foreach (Point3d p in poly) {
                    var pE = new XElement("p");
                    pE.Value = string.Format("{0} {1} {2}", Math.Round(p.X, 5), Math.Round(p.Y, 5), Math.Round(p.Z, 5));
                    polylineElement.Add(pE);
                }
                maleE.Add(polylineElement);
            }

            foreach (Polyline poly in female) {
                var polylineElement = new XElement("Polyline");
                foreach (Point3d p in poly) {
                    var pE = new XElement("p");
                    pE.Value = string.Format("{0} {1} {2}", Math.Round(p.X, 5), Math.Round(p.Y, 5), Math.Round(p.Z, 5));
                    polylineElement.Add(pE);
                }
                femaleE.Add(polylineElement);
            }

            jointTypeE.Add(name);
            jointTypeE.Add(maleE);
            jointTypeE.Add(femaleE);

            return jointTypeE;
        }


        /// <summary>
        /// Read XML
        /// string settings_path = @"C:\Users\petra\AppData\Roaming\Grasshopper\6\Libraries\RhinoJoint\TileTypeSettings.xml";
        /// var dict = XMLToTiles(settings_path);
        /// A = dict;
        /// </summary>
        /// <param name="settings_path"></param>
        /// <returns></returns>
        public static Dictionary<string, Tuple<List<Polyline>, List<Polyline>>> XMLToTiles(string settings_path) {
            XElement root = XElement.Load(settings_path);

            var tiles = new Dictionary<string, Tuple<List<Polyline>, List<Polyline>>>();


            //Iterate through all tile types
            foreach (XElement tile in root.Elements("TileType")) {
                string name = tile.Element("name").Value;
                XElement maleE = tile.Element("male");
                XElement femaleE = tile.Element("female");

                List<Polyline> male = new List<Polyline>();
                List<Polyline> female = new List<Polyline>();

                foreach (XElement polylineE in maleE.Elements("Polyline")) {
                    Polyline polyline = new Polyline();
                    foreach (XElement p in polylineE.Elements("p")) {
                        polyline.Add(StringToPoint(p.Value));
                    }
                    male.Add(polyline);
                }

                foreach (XElement polylineE in femaleE.Elements("Polyline")) {
                    Polyline polyline = new Polyline();
                    foreach (XElement p in polylineE.Elements("p")) {
                        polyline.Add(StringToPoint(p.Value));
                    }
                    female.Add(polyline);
                }

                tiles.Add(name, new Tuple<List<Polyline>, List<Polyline>>(male, female));


            }
            return tiles;

        }


        public static Point3d StringToPoint(string s) {
            string[] bits = s.Split();
            if (bits.Length != 3)
                return Point3d.Origin;

            double x = double.Parse(bits[0]);
            double y = double.Parse(bits[1]);
            double z = double.Parse(bits[2]);

            return new Point3d(x, y, z);
        }

    }
}
