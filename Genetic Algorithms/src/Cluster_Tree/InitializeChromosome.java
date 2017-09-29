package Cluster_Tree;

import java.util.ArrayList;
import java.util.Random;

public class InitializeChromosome {
	private Random r = new Random();

	public double[][] buildClusterWeightMatrix(double[][] weightMatrix, ArrayList<Integer> cluster) {
		int numberClusterVertices = cluster.size();
		double[][] clusterWeightMatrix = new double[numberClusterVertices][numberClusterVertices];
		for (int i = 0; i < numberClusterVertices; i++) {
			for (int j = 0; j < numberClusterVertices; j++) {
				clusterWeightMatrix[i][j] = weightMatrix[cluster.get(i)][cluster.get(j)];
			}
		}
		return clusterWeightMatrix;
	}

	public double[][] clusterPrimRST(double[][] weightMatrix, ArrayList<Cluster> clusters, int num_vertices) {
		int randomVertice = -1, indexRandomVertice = -1;
		double[][] tree = new double[num_vertices][num_vertices];
		double[][] clusterWeightMatrix;
		double[][] clusterSpanningTree;
		// build the spanning for each cluster
		for (int i = 0; i < ReadFiles.numberOfCluster; i++) {

			int numberClusterVertices = clusters.get(i).getCluster().size();
			clusterWeightMatrix = buildClusterWeightMatrix(weightMatrix, clusters.get(i).getCluster());
			clusterSpanningTree = primRST(clusterWeightMatrix, numberClusterVertices);

			for (int j = 0; j < numberClusterVertices; j++) {
				for (int k = 0; k < numberClusterVertices; k++) {
					if (clusterSpanningTree[j][k] > 0) {
						tree[clusters.get(i).getCluster().get(j)][clusters.get(i).getCluster()
								.get(k)] = clusterSpanningTree[j][k];
					}
				}
			}
		}

		// Each cluster is represented by one vertex
		Cluster clusterRepresentation = new Cluster();
		for (int i = 0; i < ReadFiles.numberOfCluster; i++) {
			int vertex = r.nextInt(clusters.get(i).getCluster().size());
			clusterRepresentation.getCluster().add(clusters.get(i).getCluster().get(vertex));
		}
		clusterWeightMatrix = buildClusterWeightMatrix(weightMatrix, clusterRepresentation.getCluster());
		clusterSpanningTree = primRST(clusterWeightMatrix, ReadFiles.numberOfCluster);

		// Transformation to spanning tree of G Graph
		for (int i = 0; i < ReadFiles.numberOfCluster; i++) {
			for (int j = 0; j < ReadFiles.numberOfCluster; j++) {
				if (clusterSpanningTree[i][j] > 0) {
					tree[clusterRepresentation.getCluster().get(i)][clusterRepresentation.getCluster()
							.get(j)] = clusterSpanningTree[i][j];
				}
			}
		}
		return tree;
	}

	// primRST
	public double[][] primRST(double[][] weightMatrix, int num_vertices) {
		int randomVertice = -1, indexRandomEdge = -1;
		ArrayList<Integer> contain = new ArrayList<Integer>();
		ArrayList<Edge> edgeListAdjacence = new ArrayList<Edge>();
		double[][] tempTable = new double[num_vertices][num_vertices];

		randomVertice = r.nextInt(num_vertices);
		contain.add(randomVertice);
		edgeListAdjacence.addAll(findEdge(num_vertices, randomVertice, contain, weightMatrix));

		while (contain.size() < num_vertices) {

			if (contain.size() == 0) {
				return null;
			}

			indexRandomEdge = r.nextInt(contain.size());
			Edge tempEdge = new Edge(edgeListAdjacence.get(indexRandomEdge).startVertice,
					edgeListAdjacence.get(indexRandomEdge).endVertice);
			edgeListAdjacence.remove(indexRandomEdge);

			if (!contain.contains(tempEdge.endVertice)) {

				tempTable[tempEdge.startVertice][tempEdge.endVertice] = 1.0;
				tempTable[tempEdge.endVertice][tempEdge.startVertice] = 1.0;

				contain.add(tempEdge.endVertice);
				edgeListAdjacence.addAll(findEdge(num_vertices, tempEdge.endVertice, contain, weightMatrix));
			}
		}
		return tempTable;
	}

	public ArrayList<Edge> findEdge(int num_vertices, int vertice, ArrayList<Integer> contain,
			double[][] weightMatrix) {
		ArrayList<Edge> listEdge = new ArrayList<Edge>();
		for (int i = 0; i < num_vertices; i++) {
			if (i != vertice) {
				if (!contain.contains(i)) {
					if (weightMatrix[vertice][i] > 0) {
						Edge tempEdge = new Edge(vertice, i);
						listEdge.add(tempEdge);

					}
				}
			}
		}
		return listEdge;
	}

	// get a Edge from the list of edges
	public Edge getEdge(int index, ArrayList<Edge> edgeList) {
		return edgeList.get(index);
	}

	public void findCycle(int start_vertice, int end_vertice, int number_vertices, double[][] weight_matrix,
			boolean[] visited, int[] pre) {
		visited[start_vertice] = true;
		for (int i = 0; i < number_vertices; i++) {
			if ((weight_matrix[start_vertice][i] > 0) && (i == end_vertice)) {
				pre[i] = start_vertice;
				return;
			}
			if ((weight_matrix[start_vertice][i] > 0) && (!visited[i])) {
				pre[i] = start_vertice;
				findCycle(i, end_vertice, number_vertices, weight_matrix, visited, pre);
			}
		}
	}

}
