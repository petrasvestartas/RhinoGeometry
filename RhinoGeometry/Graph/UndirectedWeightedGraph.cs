using System;
using System.Collections.Generic;
using System.Collections;
using System.Linq;

namespace RhinoGeometry {
    public class UndirectedWeightedGraph
    {
        private readonly int MaxVertices;

        public int N { get; private set; }
        public int E { get; private set; }
        public int[,] Adj { get; private set; }
        public Vertex[] VertexList { get; private set; }

        //Dijkstra
        private int Temporary = 1;
        private int Permanent = 2;
        private int Nil = -1;
        private int Infinity = 99999;
        //Dijkstra

        //Dfs
        private const int Initial = 0;
        private const int Visited = 1;
        private const int Finished = 2;
        //Output
        List<int> _vertexIndex;


        public UndirectedWeightedGraph(int maxVertices)
        {
            MaxVertices = maxVertices;
            Adj = new int[maxVertices, maxVertices];
            VertexList = new Vertex[maxVertices];
        }

        public UndirectedWeightedGraph(UndirectedWeightedGraph g)
        {
            MaxVertices = g.MaxVertices;
            N = g.N;
            E = g.E;
            Adj = g.Adj;
            VertexList = g.VertexList;
        }

        //Delete verices
        public void DeleteVertices(List<string> names)
        {
            List<int> indices = new List<int>();
            foreach (var n in names)
                indices.Add(GetIndex(n));


            List<Vertex> newVertexList = new List<Vertex>(); //Initialize newVertexList
            N = VertexList.Length - indices.Count;
            int[,] newAdj = new int[N, N]; //Initialize new array of edges

            int p = 0;
            for (int i = 0; i < VertexList.Length; i++)
                if (!indices.Contains(i))
                {
                    newVertexList.Add(VertexList[i]); //Add existing vertices that does not belong to removed indices

                    int r = 0;

                    for (int j = 0; j < VertexList.Length; j++)
                        if (!indices.Contains(j))
                        {
                            newAdj[p, r] = Adj[i, j]; //Add existing matrix values by indexes
                            r++;
                        }
                    p++;
                }


            //Output
            VertexList = newVertexList.ToArray();
            Adj = newAdj;
        }


        public void Display()
        {
            //Display vertices names
            foreach (var v in VertexList)
                Rhino.RhinoApp.WriteLine(v.name);


            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)

                    Rhino.RhinoApp.Write(Adj[i, j] + " ");
                Rhino.RhinoApp.WriteLine();

            }
        }

        public void InsertVertex(String name)
        {
            VertexList[N++] = new Vertex(name);
        }

        protected int GetIndex(String s)
        {
            for (int i = 0; i < N; i++)
                if (s.Equals(VertexList[i].name))
                    return i;
            throw new System.InvalidOperationException("Invalid Vertex");
        }

        //Returns true if edge(s1,s2) exists
        public bool EdgeExists(String s1, String s2)
        {
            return IsAdjacent(GetIndex(s1), GetIndex(s2));
        }

        protected bool IsAdjacent(int u, int v)
        {
            return (Adj[u, v] != 0);
        }

        //Insert an edge(s1,s2)
        public void InsertEdge(String s1, String s2, int wt)
        {
            int u = GetIndex(s1);
            int v = GetIndex(s2);

            if (Adj[u, v] != 0)
                Console.Write("Edge already present");
            else
            {
                Adj[u, v] = wt;
                Adj[v, u] = wt;
                E++;
            }
        }

        //Get vertex names
        public List<string> GetVertexNames()
        {
            List<string> names = new List<string>();

            foreach (var v in VertexList)
                names.Add(v.name);

            return names;
        }

        //Get edge by one vertex
        public Dictionary<string, int> GetAdjacentVertices(String s1, UndirectedWeightedGraph g)
        {
            //List<string> indices = new List<string>();
            Dictionary<string, int> edges = new Dictionary<string, int>();
            int u = GetIndex(s1);

            for (int i = 0; i < g.N; i++)
                if (Adj[u, i] != 0)
                {
                    edges.Add(VertexList[i].name, Adj[u, i]);
                    //Rhino.RhinoApp.WriteLine("Hi" + edges.Count.ToString());
                }

            return edges;
        }

        //Delete the edge (s1,s2)
        public void DeleteEdge(String s1, String s2)
        {
            int u = GetIndex(s1);
            int v = GetIndex(s2);

            if (Adj[u, v] == 0)
                Rhino.RhinoApp.WriteLine("Edge not present in the graph");
            else
            {
                Adj[u, v] = 0;
                Adj[v, u] = 0;
                E--;
            }
        }

        //Returns number of edges going out/in from a vertex
        public int Degree(String s)
        {
            int u = GetIndex(s);
            int deg = 0;
            for (int v = 0; v < N; v++)
            {
                if (Adj[u, v] != 0)
                    deg++;
            }
            return deg;
        }

        public bool VertexExists(string s1, UndirectedWeightedGraph g)
        {
            foreach (Vertex vertex in g.VertexList)
                if (vertex.name.Equals(s1))
                    return true;
            return false;
        }

        //Graph to components using:
        //1. Depth First Search then
        //2. Get edges from base graph, comparing edges and indices names

        public List<UndirectedWeightedGraph> ConnectedComponents(UndirectedWeightedGraph graph)
        {
            List<UndirectedWeightedGraph> subGraphs = new List<UndirectedWeightedGraph>();

            List<List<int>> graphIndicesList = graph.DfsTraversal_All(); //Get all vertices

            for (int i = 0; i < graphIndicesList.Count; i++)
            {
                UndirectedWeightedGraph tempGraph = new UndirectedWeightedGraph(graphIndicesList[i].Count); ;

                //Insert vertices into graph
                for (int j = 0; j < graphIndicesList[i].Count; j++)
                {
                    tempGraph.InsertVertex(VertexList[graphIndicesList[i][j]].name);
                    //Rhino.RhinoApp.WriteLine(VertexList[graphIndicesList[i][j]].name);
                }


                //Insert edges into graph
                for (int j = 0; j < graphIndicesList[i].Count; j++)
                {

                    Dictionary<string, int> neighbours = graph.GetAdjacentVertices(VertexList[graphIndicesList[i][j]].name, graph);
                    //Rhino.RhinoApp.WriteLine(j.ToString() + " " + neighbours.Count);

                    foreach (var item in neighbours)
                        if (tempGraph.VertexExists(item.Key, tempGraph))
                            tempGraph.InsertEdge(VertexList[graphIndicesList[i][j]].name, item.Key, item.Value);
                }
                //Rhino.RhinoApp.WriteLine();

                subGraphs.Add(tempGraph);
            }
            return subGraphs;
        }

        private List<List<int>> DfsTraversal_All()
        {
            List<List<int>> vertexIndicesAll = new List<List<int>>();

            int v;
            for (v = 0; v < N; v++)
                VertexList[v].state = Initial;

            _vertexIndex = new List<int>();
            Dfs(0);
            vertexIndicesAll.Add(_vertexIndex);

            for (v = 0; v < N; v++)
                if (VertexList[v].state == Initial)
                {
                    _vertexIndex = new List<int>();
                    Dfs(v);
                    vertexIndicesAll.Add(_vertexIndex);
                }

            return vertexIndicesAll;
        }

        private void Dfs(int v)
        {
            _vertexIndex.Add(v);

            VertexList[v].state = Visited;

            for (int i = 0; i < N; i++)
                if (IsAdjacent(v, i) && VertexList[i].state == Initial)
                    Dfs(i);

            VertexList[v].state = Finished;
        }

        /// <summary>
        /// Dijkstra
        /// </summary>
        /// <param name="s"></param>

        private void Dijkstra(int s)
        {
            int v, c;

            for (v = 0; v < N; v++)
            {
                VertexList[v].state = Temporary;
                VertexList[v].pathLength = Infinity;
                VertexList[v].predecessor = Nil;
            }

            VertexList[s].pathLength = 0;

            while (true)
            {
                c = TempVertexMinPL();

                if (c == Nil)
                    return;

                VertexList[c].state = Permanent;

                for (v = 0; v < N; v++)
                {
                    if (IsAdjacent(c, v) && VertexList[v].state == Temporary)
                        if (VertexList[c].pathLength + Adj[c, v] < VertexList[v].pathLength)
                        {
                            VertexList[v].predecessor = c;
                            VertexList[v].pathLength = VertexList[c].pathLength + Adj[c, v];
                        }
                }
            }
        }

        private int TempVertexMinPL()
        {
            int min = Infinity;
            int x = Nil;
            for (int v = 0; v < N; v++)
            {
                if (VertexList[v].state == Temporary && VertexList[v].pathLength < min)
                {
                    min = VertexList[v].pathLength;
                    x = v;
                }
            }
            return x;
        }

        public List<List<string>> FindPaths(String source, List<String> v)
        {
            List<List<string>> indices = new List<List<string>>();

            int s = GetIndex(source);

            Dijkstra(s);

            // Console.WriteLine("Source Vertex : " + source + "\n");

            // for (int v = 0; v < n; v++)
            for (int i = 0; i < v.Count; i++)
            {
                int _v = GetIndex(v[i]);
                //Console.WriteLine("Destination Vertex : " + vertexList[v].name);
                if (VertexList[_v].pathLength == Infinity)
                    // Rhino.RhinoApp.WriteLine("There is no path from " + source + " to vertex " + vertexList[v[i]].name);
                    return null;
                else
                    indices.Add(FindPath(s, _v));
            }

            return indices;
        }

        public List<string> FindPath(int s, int v)
        {
            List<string> indices = new List<string>();

            int i, u;
            int[] path = new int[N];
            int sd = 0;
            int count = 0;

            //Start is not equal to start
            while (v != s )
            {
                
                count++;
                path[count] = v;
                u = VertexList[v].predecessor;
                sd += Adj[u, v];
                v = u;
                if(path.Length == count)
                    break;
            }

            count++;
            if (count != N)
                path[count] = s;
            else
                count -= 1;
                
            
           

            //if(count == N)
            //    for (i = count-1; i >= 1; i--)
            //        indices.Add(VertexList[path[i]].name);
            //else
                for (i = count ; i >= 1; i--)
                    indices.Add(VertexList[path[i]].name);



            return indices;
        }



    }
}
