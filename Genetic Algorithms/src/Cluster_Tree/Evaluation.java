package Cluster_Tree;

import java.util.ArrayList;


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
     double pathLenght = 0;
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
    		 pathLenght = pathLenght  + tempMatrix[path.get(j)][path.get(j - 1)];
    		 
    	 }
     }
     return  pathLenght;
     }
    
    
}

