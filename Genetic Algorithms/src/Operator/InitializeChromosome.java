package Operator;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Queue;
import java.util.Random;

import javax.swing.JFrame;

import Dislay.Paint;
import Files_InOut.ReadFiles;
import Structures.Cluster;
import Structures.Edge;
import Structures.Individual;

public class InitializeChromosome {
	public static ArrayList<Cluster> maxClusters = new ArrayList<Cluster>();
	public static boolean isInstance1[];  // is number of  max vertices in cluster equals number of vertex in cluster instance1 
	public static boolean isMaxClusterInstance1;  //Return true if cluster in instance 1  have the number of vertex max. 
	private Random a = new Random();
    static	int minNumberOfCluster = 0;
	static int maxNumberOfCluster = 0;
	public static int maxNumberVertices ;
    public static int[] maxClusterVertices;
    public static int[] minClusterVertices;

	
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
			int vertex = a.nextInt(clusters.get(i).getCluster().size());
			clusterRepresentation.getCluster().add(clusters.get(i).getCluster().get(vertex));
		}
		clusterWeightMatrix = buildClusterWeightMatrix(weightMatrix, clusterRepresentation.getCluster());
     	clusterSpanningTree = primRST(clusterWeightMatrix, ReadFiles.numberOfCluster);
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

			indexRandomEdge = a.nextInt(edgeListAdjacence.size());
			Edge tempEdge = new Edge(edgeListAdjacence.get(indexRandomEdge).startVertice,
					edgeListAdjacence.get(indexRandomEdge).endVertice);
			edgeListAdjacence.remove(indexRandomEdge);

			if (!contain.contains(tempEdge.endVertice)) {
//
				tempTable[tempEdge.startVertice][tempEdge.endVertice] = 1.0f;
				tempTable[tempEdge.endVertice][tempEdge.startVertice] = 1.0f;
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
	/**
	 * After call this method all of the Information in the common gene going to be update
	 * Calculate  the  number of  vertices in cluster min for each cluster 
	 * Calculate  the  number of  vertices in cluster max for each cluster 
	 * And calculate the maxNumberVertices, then calculate  the number of cluster for each instance 
	 * Update the value of each the instance variables on top
	 * @param clusters1 : The list of vertices in instance 1
	 * @param clusters2 : The list of vertices in instance 2
	 */
	public void clutersVerticesInf(ArrayList<Cluster> clusters1, ArrayList<Cluster> clusters2){
		maxNumberVertices = 0;
	 if(clusters1.size() < clusters2.size()){
		 isMaxClusterInstance1 = false; 
		 minNumberOfCluster = clusters1.size();
		 maxNumberOfCluster = clusters2.size();
		 isInstance1 = new boolean[maxNumberOfCluster];
		 minClusterVertices = new int[maxNumberOfCluster];
		 maxClusterVertices = new int[maxNumberOfCluster];
		 
		 
		 for( int i = 0; i < minNumberOfCluster  ; i++){
				if(clusters1.get(i).getCluster().size() < clusters2.get(i).getCluster().size()){
					isInstance1[i] = false;
					minClusterVertices[i] = clusters1.get(i).getCluster().size();
					maxClusterVertices[i] = clusters2.get(i).getCluster().size();
					maxNumberVertices += maxClusterVertices[i];
				}else{
					isInstance1[i] = true;
					minClusterVertices[i] = clusters2.get(i).getCluster().size();
					maxClusterVertices[i] = clusters1.get(i).getCluster().size();
					maxNumberVertices += maxClusterVertices[i];
				}
			 }
		 
		 for(int j = minNumberOfCluster; j < clusters2.size(); j++){
			 isInstance1[j] = false;
			 minClusterVertices[j] = maxClusterVertices[j] = clusters2.get(j).getCluster().size();
			 maxNumberVertices += maxClusterVertices[j];
		 }
	 }else{
		 isMaxClusterInstance1 = true;
		 minNumberOfCluster = clusters2.size();
		 maxNumberOfCluster = clusters1.size();
		 
		 
		 for( int i = 0; i < minNumberOfCluster  ; i++){
				if(clusters1.get(i).getCluster().size() < clusters2.get(i).getCluster().size()){
					minClusterVertices[i] = clusters1.get(i).getCluster().size();
					maxClusterVertices[i] = clusters2.get(i).getCluster().size();
					maxNumberVertices += maxClusterVertices[i];
				}else{
					minClusterVertices[i] = clusters2.get(i).getCluster().size();
					maxClusterVertices[i] = clusters1.get(i).getCluster().size();
					maxNumberVertices += maxClusterVertices[i];
				}
			 }
		 
		 for(int j = minNumberOfCluster; j < clusters1.size(); j++){
			 minClusterVertices[j] = maxClusterVertices[j] = clusters1.get(j).getCluster().size();
			 maxNumberVertices += maxClusterVertices[j];
		 }
	 }
	 
	}
	/**
	 * 
	 * @param clusters : the cluster is mapped  to use for two instances
	 * @param r
	 * @return
	 */
	
	
	public double[][] decodeGroupGene(ArrayList<Cluster> clusters, Random r){
	
		
		
		double[][] treeMax  = new double[maxNumberVertices][maxNumberVertices];
		double[][] weightMatrix = new double[maxNumberVertices][maxNumberVertices];
		for( int i = 0; i < maxNumberVertices; i++ ){
			for( int j = 0; j < maxNumberVertices; j++){
				weightMatrix[i][j] = 1.0f;
			}
		}
		
		
		// initialize gene for each cluster 
		 for(int i = 0; i < maxNumberOfCluster; i++){
			 //small 
			double[][] bigClusterTree =  new double[ maxClusterVertices[i]][ maxClusterVertices[i]];
			double[][] smallClusterTree =  new double[ minClusterVertices[i]][ minClusterVertices[i]];
			for(int j = 0; j < minClusterVertices[i]; j++){
				for(int k = 0; k < minClusterVertices[i]; k++){
				smallClusterTree[j][k] = 1.0f;
				}
			}
		    
			smallClusterTree = primRST(smallClusterTree,minClusterVertices[i]);
			
			for(int j = 0; j < minClusterVertices[i]; j++){
				for(int k = 0; k < minClusterVertices[i]; k++){
				bigClusterTree[j][k] = smallClusterTree[j][k]; 
				}
			}
			
			
		    int tempNum =  minClusterVertices[i];
			while(tempNum < maxClusterVertices[i]){
				int randomIndex = r.nextInt(tempNum);
				bigClusterTree[randomIndex][tempNum] = 1.0f;
				bigClusterTree[tempNum][randomIndex] = 1.0f;
				tempNum++;
		    }
			// transform to the max tree
			for(int j = 0; j < maxClusterVertices[i]; j++){
				for(int k = 0; k < maxClusterVertices[i]; k++){
					if(bigClusterTree[j][k] > 0){
					treeMax[clusters.get(i).getCluster().get(j)][clusters.get(i).getCluster().get(k)]
							= bigClusterTree[j][k];
							}
				}
			}
			
		 }
		 
		 
		 //---------------------------End initialize gene for each cluster---------------------------//
		 
		 // Start to initialize for presentation cluster.
		 
		 // build up the Tree for  presentation Cluster in  instance which has the number of clusters less than. 
		 double[][] clusterWeightMatrix;
		 double[][] smallClusterSpanningTree;
	     // choose randomly the presentation  for  small presentation cluster
		 Cluster clusterRepresentation =  new Cluster();
		 for( int i = 0; i < minNumberOfCluster; i++ ){
			 int vertex = r.nextInt(minClusterVertices[i]);
			 clusterRepresentation.getCluster().add(clusters.get(i).getCluster().get(vertex));
	
		 }
		 
		 
		 clusterWeightMatrix = buildClusterWeightMatrix(weightMatrix, clusterRepresentation.getCluster());
	     smallClusterSpanningTree = primRST(clusterWeightMatrix,minNumberOfCluster);
	     
	     
	     // Choose randomly the presentation for big presentation cluster
	     for(int i = minNumberOfCluster; i < maxNumberOfCluster; i++ ){
	    	 int vertex = r.nextInt(minClusterVertices[i]);
	    	 clusterRepresentation.getCluster().add(clusters.get(i).getCluster().get(vertex));
	    	 
	     }

	     // Copy the matrix of the small presentation cluster to big presentation cluster. 
	     double[][] bigClusterSpanningTree = new double[maxNumberOfCluster][maxNumberOfCluster];
	     for(int j = 0; j < minNumberOfCluster; j++){
				for(int k = 0; k < minNumberOfCluster; k++){
					bigClusterSpanningTree[j][k] = 	 smallClusterSpanningTree[j][k]; 
				}
			}
	     
	     
	     //  build up the Tree for presentation Cluster in instance which has the number of clusters greater.
	     int tempNum = minNumberOfCluster;
			while(tempNum < maxNumberOfCluster){
				int randomIndex = r.nextInt(tempNum);
				bigClusterSpanningTree[randomIndex][tempNum] = 1.0f;
				bigClusterSpanningTree[tempNum][randomIndex] = 1.0f;
				tempNum++;
		    }
			
	     
	    // Transform presentation cluster to TreeMax
	     for (int i = 0; i < maxNumberOfCluster; i++) {
				for (int j = 0; j < maxNumberOfCluster; j++) {
					if (bigClusterSpanningTree[i][j] > 0) {
						treeMax[clusterRepresentation.getCluster().get(i)][clusterRepresentation.getCluster()
								.get(j)] = bigClusterSpanningTree[i][j];
					}
				}
			}
	     return treeMax;
	}
	/**
	 * we build cluster max following by:
	 *  the first cluster = (0, 1, 2 ..... (number of vertices in  1st cluster max))
	 *  the second cluster = ((number of vertices in 1st cluster max), (number of vertices in 1st cluster max + 1).....) 
	 *  ..etc...
	 *  the last cluster = (..........(number of vertices in last cluster max - 1)+(number of vertices in last  cluster max ))
	 */
	public static  void buildCluster(){

		int count = 0;
		
		for(int i = 0 ; i < maxNumberOfCluster; i++){
			Cluster cluster = new Cluster();
			for( int j = 0 ; j < maxClusterVertices[i]; j++){
				cluster.addElement(count, j);
			count++;
			}
			maxClusters.add(cluster);
		}
	}
	
/**
 * Tách từ cây biểu diễn chung ra các lời giải để dánh giá. 
 * @param maxTree
 * @param clusters1
 * @param clusters2
 * @param num_vertices1
 * @param num_vertices2
 * @return
 */
	public ArrayList<double[][]> devideTwoTree(double[][] maxTree, ArrayList<Cluster> clusters1,
			ArrayList<Cluster> clusters2,int num_vertices1,int num_vertices2 , Random r) {
		
		GraphMethods graphMethods = new GraphMethods();
		Crossover crossover = new Crossover();
		ArrayList<double[][]> list = new ArrayList<double[][]>(); 
		double[][] weightTree1 = new double[num_vertices1][num_vertices1];
		double[][] weightTree2 = new double[num_vertices2][num_vertices2];
		double[][] TreeInst1 = new double[num_vertices1][num_vertices1];
		double[][] TreeInst2 = new double[num_vertices2][num_vertices2];
		int numberOfcluster1 = clusters1.size();
		int numberOfcluster2 = clusters2.size();
	
			for( int i = 0; i < numberOfcluster1; i ++){
				  double[][] TreeOfEachCluster =  buildClusterWeightMatrix(maxTree, maxClusters.get(i).getCluster());
				  
			        // copy the tree in Each cluster for Tree in Instance1.
			        for( int j = 0 ; j < clusters1.get(i).getCluster().size(); j++){
			        	for(int k = 0; k < clusters1.get(i).getCluster().size(); k++){
			        		if(TreeOfEachCluster[j][k]  > 0){
			        			TreeInst1[clusters1.get(i).getCluster().get(j)][clusters1.get(i).getCluster().get(k)]=
			        					TreeInst1[clusters1.get(i).getCluster().get(k)][clusters1.get(i).getCluster().get(j)]= 1.0f;
			        		}
			        	}
			        }
			 
			}
			weightTree1 = findSpanningTreeBetweenClusters(maxTree,num_vertices1, maxClusters, clusters1);
			for( int j = 0 ; j < num_vertices1; j++){
	        	for(int k = 0; k < num_vertices1; k++){	
	        		if(weightTree1[j][k] > 0){
	        		TreeInst1[j][k] = weightTree1[j][k];
	        		}
	        		}
	        	}
			list.add(TreeInst1);
			
			
			
			for( int i = 0; i < numberOfcluster2; i ++){
				  double[][] TreeOfEachCluster =  buildClusterWeightMatrix(maxTree, maxClusters.get(i).getCluster());
			        // copy the tree in Each cluster for Tree in Instance2.
			        for( int j = 0 ; j < clusters2.get(i).getCluster().size(); j++){
			        	for(int k = 0; k <clusters2.get(i).getCluster().size(); k++){
			        		if(TreeOfEachCluster[j][k]  > 0){
			        			TreeInst2[clusters2.get(i).getCluster().get(j)][clusters2.get(i).getCluster().get(k)]=
			        					TreeInst2[clusters2.get(i).getCluster().get(k)][clusters2.get(i).getCluster().get(j)]= 1.0f;
			        		}
			        	}
			        }
			}
		
			
			weightTree2 = findSpanningTreeBetweenClusters(maxTree,num_vertices2, maxClusters, clusters2);
			for( int j = 0 ; j < num_vertices2; j++){
	        	for(int k = 0; k < num_vertices2; k++){	
	        		if(weightTree2[j][k] > 0){
	        		TreeInst2[j][k] = weightTree2[j][k];
	        		}
	        		}
	        	}
			list.add(TreeInst2);
			
//			  // ve
//					JFrame gf = new JFrame();
//					gf.setVisible(true);
//					gf.setSize(800, 800);
//					gf.setTitle(" cay 1");
//					gf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
//					gf.setVisible(true);
//					
//					Paint  p = new Paint();
//					p.weightMatrix = TreeInst2;
//					p.num_vertex = 76;
//				    gf.add(p);
//				    gf.setVisible(true);
					
	        return list;
			}
		
	
	/**
	 * 
	 * @param parent
	 * @param num_vertex
	 * @param clusters
	 * @return the spanning tree between clusters 
	 */
	public double[][] findSpanningTreeBetweenClusters(double[][] parent, int num_vertex, ArrayList<Cluster> clusters,  ArrayList<Cluster> clusters1 ){
		int numberOfCluster = clusters1.size();
		double[][]  tree = new double[num_vertex][num_vertex];
		int numberClusterVertex1 = 0;
		int numberClusterVertex2 = 0; 
		for( int i = 0; i < numberOfCluster; i++ ){
			for( int j = 0; j < numberOfCluster; j++){
				tree[i][j] = 0; 
				tree[j][i] = 0;
			}
		}
		 
		    for( int i = 0; i < numberOfCluster; i++){
			numberClusterVertex1 = clusters1.get(i).getCluster().size();
			for( int j = 0; j < numberClusterVertex1; j++){
				for( int k = i + 1; k < numberOfCluster; k++ ){
					numberClusterVertex2 = clusters1.get(k).getCluster().size();
					
					for( int t = 0; t < numberClusterVertex2; t++){
					if(parent[clusters.get(i).getCluster().get(j)][clusters.get(k).getCluster().get(t)] > 0){
						
						tree[clusters1.get(i).getCluster().get(j)][clusters1.get(k).getCluster().get(t)]  = 1.0f;
						tree[clusters1.get(k).getCluster().get(t)][clusters1.get(i).getCluster().get(j)]  = 1.0f;
	     }
					}
				}
			}
		}
		return tree;
	}
	

	

}

