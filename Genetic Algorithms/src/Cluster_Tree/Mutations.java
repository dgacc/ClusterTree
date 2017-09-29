package Cluster_Tree;

import java.util.ArrayList;
import java.util.Random;

public class Mutations {
	private static Random r = new Random();
	InitializeChromosome chromosome = new InitializeChromosome();
	GraphMethods graphMethods = new GraphMethods();
	
	/**
	 * 
	 * @param parents
	 * @param num_vertex
	 * @return
	 */
	public double[][] mutationForEachClusters(double[][] parents,int num_vertex){
		double[][] offspring = new double[num_vertex][num_vertex];
		for(  int i = 0; i < num_vertex; i++){
			for( int j = 0; j < num_vertex; j++){
				offspring[i][j] = parents[i][j]; 
			}
		}
		int startVertex = r.nextInt(num_vertex);
		int endVertex = r.nextInt(num_vertex);
		
		while((startVertex == endVertex) || (offspring[startVertex][endVertex] > 0 ))
		// two vertices are equal or adjacence > 0 then choose two vertices again
		{
			startVertex = r.nextInt(num_vertex);
			endVertex = r.nextInt(num_vertex);
		}
		// find the path from start_Vertice to end_vertice
		 // initialize  two matrix visited  matrix = false and pre Matrix = -1;
		boolean[] visited = new boolean[num_vertex];
		int[] pre = new int[num_vertex];
		for( int i = 0; i < num_vertex; i ++ ){
			visited[i] = false;
			pre[i] = -1;
		}
		
	   chromosome.findCycle(startVertex, endVertex, num_vertex, offspring, visited, pre);
		ArrayList<Integer> path = graphMethods.printPath(startVertex, endVertex, pre);
		// delete a edge from cycle
		
		int index1 = r.nextInt(path.size() - 1);
		int index2 =  index1 + 1;
		
		offspring[path.get(index1)][path.get(index2)] = 0f;
		offspring[path.get(index2)][path.get(index1)] = 0f;
		
		offspring[path.get(startVertex)][path.get(endVertex)] = 1f;
		offspring[path.get(endVertex)][path.get(startVertex)] = 1f;
		
		return offspring;
		
	}
	
	/**
	 * 
	 * @param parent
	 * @param num_vertex
	 * @param clusters
	 * @return
	 */
	public double[][] mutationClusterTree(double[][] parent, int num_vertex, ArrayList<Cluster> clusters){
		
		InitializeChromosome initializeChromosome = new InitializeChromosome();
		double[][] offspring = new double[num_vertex][num_vertex];
		int numberOfCluster  = clusters.size();
		int indexCluster = r.nextInt(numberOfCluster);
		
		for( int i = 0 ; i < num_vertex ; i ++){
			for( int j = 0; j < num_vertex; j++){
				offspring[i][j] = parent[i][j];
			}
		}
		//choose the cluster which have more then two vertex
		while(clusters.get(indexCluster).getCluster().size() < 3){
			indexCluster = r.nextInt(numberOfCluster);
		}
		
		int numberClusterVertex = clusters.get(indexCluster).getCluster().size();
		// create the weight matrix, and spanning tree for that  cluster then  do mutation
		double[][] clusterWeightMatrix = new double[numberClusterVertex][numberClusterVertex];
		double[][] clusterSpanningTree = new double[numberClusterVertex][numberClusterVertex];
		
		clusterWeightMatrix = initializeChromosome.buildClusterWeightMatrix(parent,
				clusters.get(indexCluster).getCluster());
		clusterSpanningTree = mutationForEachClusters(clusterWeightMatrix, numberClusterVertex);
		
		// 
		for( int i = 0; i < numberClusterVertex; i++ ){
			for( int j = 0; j < numberClusterVertex; j++){
				offspring[clusters.get(indexCluster).getCluster().get(i)]
						[clusters.get(indexCluster).getCluster().get(j)] = clusterSpanningTree[i][j];
			}
		}
		return offspring;
	}

 
}
