package Cluster_Tree;

import java.util.ArrayList;
import java.util.Random;

public class Mutations {
	private static Random r = new Random();
	InitializeChromosome chromosome = new InitializeChromosome();
	GraphMethods graphMethods = new GraphMethods();
	
	
	public int[][] Mutation_Clusters(int[][] parents,int num_vertex){
		int[][] offspring = new int[num_vertex][num_vertex];
		for(  int i = 0; i < num_vertex; i++){
			for( int j = 0; j < num_vertex; j++){
				offspring[i][j] = parents[i][j]; 
			}
		}
		int startVertex = r.nextInt(num_vertex);
		int endVertex = r.nextInt(num_vertex);
		
		while((startVertex == endVertex) || (offspring[num_vertex][num_vertex] > 0 ))
		// two vertices are equal or adjacence  then choose two vertices again
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
		
		offspring[path.get(index1)][path.get(index2)] = 0;
		offspring[path.get(index2)][path.get(index1)] = 0;
		
		offspring[path.get(startVertex)][path.get(endVertex)] = 1;
		offspring[path.get(endVertex)][path.get(startVertex)] = 1;
		
		return offspring;
		
	}

}
