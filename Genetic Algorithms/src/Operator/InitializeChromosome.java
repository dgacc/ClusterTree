package Cluster_Tree;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Queue;
import java.util.Random;

import javax.swing.JFrame;

public class InitializeChromosome {
	private Random a = new Random();
/**
 * copy the weightMatrix for each Cluster 
 * @param weightMatrix
 * @param cluster
 * @return cluster weightMatrix 
 */
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

/**
 * generate the  randomly tree based on PRimRST
 * first: generate the tree for each  cluster  then convert to the TRee
 * Second: generate the presentation tree: meaning each cluster is presented by a vertex then covert to the Tree  
 * @param weightMatrix
 * @param clusters
 * @param num_vertices
 * @return  Tree after do that stuff
 */
	public double[][] clusterPrimRST(double[][] weightMatrix, ArrayList<Cluster> clusters, int num_vertices) {
		int randomVertice = -1, indexRandomVertice = -1;
		double[][] tree = new double[num_vertices][num_vertices];
		double[][] clusterWeightMatrix;
		double[][] clusterSpanningTree;
		int[] indexPresentationVertex = new int[clusters.size()];
		// build the spanning for each cluster
		for (int i = 0; i < ReadFiles.numberOfCluster; i++) {

			int numberClusterVertices = clusters.get(i).getCluster().size();
//			if(i == 2){
//				indexPresentationVertex[i] = 14;
//			}else{
//				indexPresentationVertex[i] = a.nextInt(numberClusterVertices);
//		}
			clusterWeightMatrix = buildClusterWeightMatrix(weightMatrix, clusters.get(i).getCluster());
			clusterSpanningTree = primRST(clusterWeightMatrix, numberClusterVertices);
//	        clusterSpanningTree = BFSTree(clusterWeightMatrix, numberClusterVertices, indexPresentationVertex[i]);

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
			int vertex = a.nextInt(clusters.get(i).getCluster().size());
			clusterRepresentation.getCluster().add(clusters.get(i).getCluster().get(vertex));
		}
		clusterWeightMatrix = buildClusterWeightMatrix(weightMatrix, clusterRepresentation.getCluster());
     	clusterSpanningTree = primRST(clusterWeightMatrix, ReadFiles.numberOfCluster);
	//	clusterSpanningTree = BFSTree(clusterWeightMatrix, ReadFiles.numberOfCluster,indexPresentationVertex[a.nextInt(ReadFiles.numberOfCluster)] );

//		clusterSpanningTree = bigPrimRST(clusterWeightMatrix, ReadFiles.numberOfCluster);
		
		
		// Transform to spanning tree of G Graph
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

	
	/**
	 * this is algorithms to generate the  tree randomly
	 * first: apply "findEdge" method to find out the list of  addjence
	 * @param weightMatrix
	 * @param num_vertices
	 * @return
	 */
	public double[][] primRST(double[][] weightMatrix, int num_vertices) {
		int randomVertice = -1, indexRandomEdge = -1;
		ArrayList<Integer> contain = new ArrayList<Integer>();
		ArrayList<Edge> edgeListAdjacence = new ArrayList<Edge>();
		double[][] tempTable = new double[num_vertices][num_vertices];

		randomVertice = a.nextInt(num_vertices);
		contain.add(randomVertice);
		edgeListAdjacence.addAll(findEdge(num_vertices, randomVertice, contain, weightMatrix));

		while (contain.size() < num_vertices) {

			if (contain.size() == 0) {
				return null;
			}

			//indexRandomEdge = r.nextInt(contain.size());
			indexRandomEdge = a.nextInt(edgeListAdjacence.size());
			Edge tempEdge = new Edge(edgeListAdjacence.get(indexRandomEdge).startVertice,
					edgeListAdjacence.get(indexRandomEdge).endVertice);
			edgeListAdjacence.remove(indexRandomEdge);

			if (!contain.contains(tempEdge.endVertice)) {
//
//				tempTable[tempEdge.startVertice][tempEdge.endVertice] = 1.0f;
//				tempTable[tempEdge.endVertice][tempEdge.startVertice] = 1.0f;
				
				// initialize the  files instances
				Random r = new Random();
				double value = (double) (5 + r.nextInt(6)); 
				tempTable[tempEdge.startVertice][tempEdge.endVertice] = value;
				tempTable[tempEdge.endVertice][tempEdge.startVertice] = value;
				//end

				contain.add(tempEdge.endVertice);
				edgeListAdjacence.addAll(findEdge(num_vertices, tempEdge.endVertice, contain, weightMatrix));
			}
		}
		return tempTable;
	}
	/**
	 *  this is the algorithm to build the  weight matrix if have the edge from  cluster to cluster 
	 * @param weightMatrix
	 * @param num_vertices
	 * @return the tree with weight matrix;
	 */ 
	public double[][] bigPrimRST(double[][] weightMatrix, int num_vertices) {
		int randomVertice = -1, indexRandomEdge = -1;
		ArrayList<Integer> contain = new ArrayList<Integer>();
		ArrayList<Edge> edgeListAdjacence = new ArrayList<Edge>();
		double[][] tempTable = new double[num_vertices][num_vertices];

		randomVertice = a.nextInt(num_vertices);
		contain.add(randomVertice);
		edgeListAdjacence.addAll(findEdge(num_vertices, randomVertice, contain, weightMatrix));

		while (contain.size() < num_vertices) {

			if (contain.size() == 0) {
				return null;
			}

			//indexRandomEdge = r.nextInt(contain.size());
			indexRandomEdge = a.nextInt(edgeListAdjacence.size());
			Edge tempEdge = new Edge(edgeListAdjacence.get(indexRandomEdge).startVertice,
					edgeListAdjacence.get(indexRandomEdge).endVertice);
			edgeListAdjacence.remove(indexRandomEdge);

			if (!contain.contains(tempEdge.endVertice)) {
//
//				tempTable[tempEdge.startVertice][tempEdge.endVertice] = 1.0f;
//				tempTable[tempEdge.endVertice][tempEdge.startVertice] = 1.0f;
				// Initialize the weightMatrix 
				Random r = new Random();
				double value = (double) (20 + r.nextInt(6));
				tempTable[tempEdge.startVertice][tempEdge.endVertice] = value;
				tempTable[tempEdge.endVertice][tempEdge.startVertice] = value;
				// End 

				contain.add(tempEdge.endVertice);
				edgeListAdjacence.addAll(findEdge(num_vertices, tempEdge.endVertice, contain, weightMatrix));
			}
		}
		return tempTable;
	}

	
	/**
	 * find that addjence  with that edge
	 * @param num_vertices
	 * @param vertice
	 * @param contain
	 * @param weightMatrix
	 * @return the  list of  the addjence edge
	 */
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
/**
 * 
 * @param start_vertice
 * @param end_vertice
 * @param number_vertices
 * @param weight_matrix
 * @param visited
 * @param pre
 */
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
	
	
	/** 
	 *  find the tree based on BFS, make a tree  
	 * @param weightMatrix
	 * @param num_vertex
	 * @param startVertex
	 * @return BFS Tree
	 */
    public double[][] BFSTree(  double[][] weightMatrix, int num_vertex, int startVertex){
		  
		double[][] BFStree = new double[num_vertex][num_vertex];
		boolean[] mark = new boolean[num_vertex];
		Queue<Integer> queue = new LinkedList<>();
	    // initialize label for each vertex,  
		for( int i = 0; i < num_vertex; i++){
			mark[i] = true;
		}
		mark[startVertex] = false;
		queue.add(startVertex);
		
		 while(!queue.isEmpty()){
			 int considerVertex = queue.poll();
			 for( int i = 0; i < num_vertex; i++){
				 if(weightMatrix[considerVertex][i] > 0  && mark[i]){
					 queue.add(i);
					 mark[i] = false;
					 BFStree[considerVertex][i] = 1.0f;
					 BFStree[i][considerVertex] = 1.0f;
					 
				 }				 
			 }
			
		 }
	return BFStree;
	}
    /**
     * can sua dau vao cho bai toan
     * @param numCluster
     * @param weightMatrix
     * @param vertex_In_Cluster
     * @param idx
     * @param edMat
     * @return
     */

	public double[][] localSearch(int numCluster, double[][] weightMatrix, int[][] vertex_In_Cluster, int[] idx,
			double[][] edMat) {
		double temp[][] = new double[weightMatrix.length][weightMatrix.length];
		for (int i = 0; i < weightMatrix.length; i++) {
			for (int j = 0; j < weightMatrix.length; j++) {
				temp[i][j] = edMat[i][j];
			}
		}
		for (int i = 0; i < numCluster; i++) {
			for (int j = i + 1; j < numCluster; j++) {
				if (edMat[idx[i]][idx[j]] > 0) {
					int[] arr = findMinDistance(weightMatrix, vertex_In_Cluster[i], vertex_In_Cluster[j]);
					temp[idx[i]][idx[j]] = 0;
					idx[i] = arr[0];
					idx[j] = arr[1];
					temp[idx[i]][idx[j]] = 1;
				}
			}
		}
		return temp;
	}

	public int[] findMinDistance(double[][] weightMatrix, int[] a1, int a2[]) {
		int result[] = new int[2];
		double min = weightMatrix[a1[0]][a2[0]];
		for (int i = 0; i < a1.length; i++) {
			for (int j = 0; j < a2.length; j++) {
				if (weightMatrix[a1[i]][a2[j]] < min) {
					min = weightMatrix[a1[i]][a2[j]];
					result[0] = a1[i];
					result[1] = a2[j];
				}
			}
		}
		return result;
	}

}

