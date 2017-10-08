package Cluster_Tree;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Queue;


public class Evaluation {
	GraphMethods graphMethodsClass  = new GraphMethods();
	/**
	 * evaluate the fitness for cluster tree
	 * @param weightMatrix
	 * @param edgeMatrix
	 * @param num_vertex
	 * @param startVertex
	 * @return  pathLength: meaning the distance from root to the others vertex; 
	 */
    public double clusterEvaluate(double[][] weightMatrix,double[][] edgeMatrix, int num_vertex, int startVertex){
     double pathLength = 0;
     int[] pre = new int[num_vertex];
     ArrayList<Integer> path = new ArrayList<Integer>();
     double[][] tempMatrix = new double[num_vertex][num_vertex];
     
     // initialize the value to  new weight matrix
     for(int i = 0 ; i < num_vertex; i++){
	   for( int j = 0; j < num_vertex; j++) { 	 	
		   if( edgeMatrix[i][j] > 0){
			   tempMatrix[i][j] = weightMatrix[i][j];
		   }else{
			   tempMatrix[i][j] = 0; 
		   }
	   }
     }
    // start to  the distances; 
     pre = graphMethodsClass.dijkstra(weightMatrix, num_vertex, startVertex);
     for( int i = 1 ; i < num_vertex; i++){
    	 path = graphMethodsClass.printPath(startVertex, i, pre);
    	 for( int j = path.size() -1; j > 0 ; j --){
    		 pathLength = pathLength  + tempMatrix[path.get(j)][path.get(j - 1)];
    		 
    	 }
     }
     return  pathLength;
     }
    
    public double evaluation(double[][] weightMatrix, double[][] tree, int num_vertex, int startVertex){
    	 double[]  distances  =  new double[num_vertex];// distance between root and the others
    	 double sum  = 0;
    	 distances[startVertex] = 0;
    	 boolean[] mark = new boolean[num_vertex];
    	 Queue<Integer> queue  = new LinkedList<>();
    	 for( int i = 0; i < num_vertex; i ++){
    		 mark[i] = true;
    	 }
    	 queue.add(startVertex);
    	 while(!queue.isEmpty()){
    		 int u = queue.poll();
    		 mark[u] = false;
    		 for( int i = 0; i < num_vertex; i++){
    			 if(tree[u][i] > 0 && mark[i] ){
    				queue.add(i);
    				mark[i] = false;
    				distances[i] =  distances[u] + weightMatrix[u][i];
    				sum += distances[i]; 
    				
    			 }
    		 }
    	 }
    	 return sum;
    }
    
}

