package Cluster_Tree;

import java.util.ArrayList;
import java.util.Random;


public class Crossover {
	InitializeChromosome initilizeChromosome = new InitializeChromosome();
	
	/**
	 * 
	 * @param parent
	 * @param num_vertex
	 * @param clusters
	 * @return the spanning tree between clusters 
	 */
	private double[][] findSpanningTreeBetweenClusters(double[][] parent, int num_vertex, ArrayList<Cluster> clusters ){
		int numberOfCluster = clusters.size();
		double[][]  tree = new double[numberOfCluster][numberOfCluster];
		int numberClusterVertex1 = 0;
		int numberClusterVertex2 = 0; 
		for( int i = 0; i < num_vertex; i++ ){
			for( int j = 0; j < num_vertex; j++){
				tree[i][j] = 0; 
				tree[j][i] = 0;
			}
		}
		
		for( int i = 0; i < numberOfCluster; i++){
			numberClusterVertex1 = clusters.get(i).getCluster().size();
			for( int j = 0; j < numberClusterVertex1; j++){
				for( int k = i + 1; k < numberOfCluster; k++ ){
					numberClusterVertex2 = clusters.get(k).getCluster().size();
					for( int t = 0; t < numberClusterVertex2; t++){
					if(parent[clusters.get(i).getCluster().get(j)][clusters.get(k).getCluster().get(t)] > 0){
						tree[j][t] = parent[clusters.get(i).getCluster().get(j)][clusters.get(k).getCluster().get(t)];
						tree[t][j] = tree[j][t];
					}
					}
				}
			}
		}
		return tree;
	}
	// crossover for big spanning  tree
	private double[][] primRSTcrossover(double[][] father, double[][] mother, int num_vertex){
		 double[][] G_cr = new double[num_vertex][num_vertex];
		 double[][] spanningTree = new double[num_vertex][num_vertex];
		 // initialize the G_cr Matrix = fatherMatrix + motherMatrix
		for(int i = 0; i < num_vertex; i++){
			for(int j = 0; j< num_vertex; j++){
				G_cr[i][j] = father[i][j] + mother[i][j];
			}
		}
		
		spanningTree = initilizeChromosome.primRST(G_cr, num_vertex);
	    return spanningTree;
	}
	// use crossover for each cluster
	public double[][] clusterCrossover(double[][] father,double[][] mother, int num_vertex, ArrayList<Cluster> clusters ){
		double[][] G_cr = new double[num_vertex][num_vertex];
		double[][] Tree = new double[num_vertex][num_vertex];
		double[][] clusterWeightMatrix1, clusterWeightMatrix2;
		double[][] spanningTreeOfCluster;
		int numberClusterVertex = 0;
		int numberOfCluster = clusters.size();
		for( int i = 0; i < numberOfCluster; i++){
			numberClusterVertex = clusters.get(i).getCluster().size();
			clusterWeightMatrix1 = initilizeChromosome.buildClusterWeightMatrix(father, clusters.get(i).getCluster());
			clusterWeightMatrix2 = initilizeChromosome.buildClusterWeightMatrix(mother, clusters.get(i).getCluster());
			
		    spanningTreeOfCluster = primRSTcrossover(clusterWeightMatrix1, clusterWeightMatrix2, numberClusterVertex);
		    // convert to the graph tree
		    for(int j = 0; j < numberClusterVertex; j++){
		    	for(int k = 0; k < numberClusterVertex; k++){
		    		Tree[clusters.get(i).getCluster().get(j)][clusters.get(i).getCluster().get(k)] = spanningTreeOfCluster[j][k];		
		    	}
		    }
		}
		    // build the vertex representation for each cluster
		    clusterWeightMatrix1 = findSpanningTreeBetweenClusters(father, num_vertex, clusters);
		    clusterWeightMatrix2 = findSpanningTreeBetweenClusters(mother, num_vertex, clusters);
		    
		    // generate the new offspring
		    spanningTreeOfCluster = primRSTcrossover(father, mother, num_vertex);
		    // convert to spanning tree of Graph
		    int[] indexCluster = new int[numberOfCluster];
		    Random r = new Random();
		    for(int i = 0; i < numberOfCluster; i++ ){
		    	int position = r.nextInt(clusters.get(i).getCluster().size());
		    	indexCluster[i] = clusters.get(i).getCluster().get(position);
		    	for( int j = 0; j < numberOfCluster; j++){
		    		for( int k = 0; k < numberOfCluster; k++){
		    			if( spanningTreeOfCluster[j][k] > 0){
		    				Tree[indexCluster[j]][indexCluster[k]] = spanningTreeOfCluster[j][k];
		    			}
		    		}
		    	}
		    }
		    return Tree;
		    
	}
	

}
