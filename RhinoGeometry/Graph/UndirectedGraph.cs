using System;
using System.Collections.Generic;

namespace RhinoGeometry {
    public class UndirectedGraph {
        public readonly int MaxVertices;

        protected int N;
        protected int E;
        public bool[,] Adj;
        protected Vertex[] VertexList;
        protected object Attributes;
        protected List<Edge> EdgesList;

        public UndirectedGraph(int maxVertices) {
            MaxVertices = maxVertices;
            Adj = new bool[maxVertices, maxVertices];
            VertexList = new Vertex[maxVertices];
            EdgesList = new List<Edge>();
        }

        public int Vertices() {
            return N;
        }

        public int Edges() {
            return E;
        }

        public List<String> GetVertices() {
            List<String> v = new List<string>();

            foreach (Vertex n in VertexList)
                v.Add(n.name);

            return v;
        }

        public void Display() {

            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++)

                    if (Adj[i, j])
                        Rhino.RhinoApp.Write("1 ");
                    else
                        Rhino.RhinoApp.Write("0 ");
                Rhino.RhinoApp.WriteLine();

            }
        }

        public void InsertVertex(String name) {
            VertexList[N++] = new Vertex(name);
        }

        protected int GetIndex(String s) {
            for (int i = 0; i < N; i++)
                if (s.Equals(VertexList[i].name))
                    return i;
            throw new System.InvalidOperationException("Invalid Vertex");
        }


        //Returns true if edge(s1,s2) exists
        public bool EdgeExists(String s1, String s2) {
            return IsAdjacent(GetIndex(s1), GetIndex(s2));
        }

        protected bool IsAdjacent(int u, int v) {
            return Adj[u, v];
        }

        //Insert an edge(s1,s2)
        public void InsertEdge(String s1, String s2) {

            int u = GetIndex(s1);
            int v = GetIndex(s2);

            EdgesList.Add(new Edge(u, v, 1));


            if (Adj[u, v]) {
               // Rhino.RhinoApp.Write("Edge already present");
            } else {
                Adj[u, v] = true;
                Adj[v, u] = true;
                E++;
            }
        }

        //Delete the edge (s1,s2)
        public void DeleteEdge(String s1, String s2) {
            int u = GetIndex(s1);
            int v = GetIndex(s2);

            if (Adj[u, v] == false) {
               // Rhino.RhinoApp.WriteLine("Edge not present in the graph");
            } else {
                Adj[u, v] = false;
                Adj[v, u] = false;
                E--;
            }
        }

        //Returns number of edges going out from a vertex
        public int Degree(String s) {
            int u = GetIndex(s);
            int deg = 0;
            for (int v = 0; v < N; v++) {
                if (Adj[u, v])
                    deg++;
            }
            return deg;
        }

        public List<int> GetNeighbours(String s) {
            List<int> neighbours = new List<int>();
            int u = GetIndex(s);

            for (int v = 0; v < N; v++)
                if (Adj[u, v])
                    neighbours.Add(v);

            return neighbours;
        }

        public List<Edge> GetEdges() {
            return EdgesList;
        }

        public void SetAttribute(object data) {
            Attributes = data;
        }

        public object GetAttribute() {
            return Attributes;
        }




    }
}
