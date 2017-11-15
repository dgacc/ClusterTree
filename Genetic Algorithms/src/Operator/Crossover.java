package Cluster_Tree;

import java.util.ArrayList;
import java.util.Random;

import javax.swing.JFrame;


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
		for( int i = 0; i < numberOfCluster; i++ ){
			for( int j = 0; j < numberOfCluster; j++){
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
						tree[i][k] = parent[clusters.get(i).getCluster().get(j)][clusters.get(k).getCluster().get(t)];
						tree[k][i] = tree[i][k];
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
		double[][] Tree = new double[num_vertex][num_vertex];
		double[][] clusterWeightMatrix1, clusterWeightMatrix2;
		double[][] spanningTreeOfCluster;
		int numberClusterVertex = 0;
		int numberOfCluster = clusters.size();
		
		// 
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
		    spanningTreeOfCluster = primRSTcrossover(clusterWeightMatrix1, clusterWeightMatrix1, numberOfCluster);
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
	/**
	 * 
	 * @param father
	 * @param mother
	 * @param num_vertex
	 * @param startVertex
	 * @return
	 */
	public  double[][] BFSCrossover(double[][] father, double[][]  mother, int num_vertex, int startVertex){
		
		double[][] parentGroup  = new double[num_vertex][num_vertex];
		double[][] spanningTree = new double[num_vertex][num_vertex];
		// addition two matrix of father and mother
		
		for( int i = 0; i < num_vertex; i++){
			for( int  j = 0; j < num_vertex; j++){
				parentGroup[i][j] = father[i][j] + mother[i][j];
			}
		}
		spanningTree  = initilizeChromosome.BFSTree(parentGroup, num_vertex, startVertex);
		return spanningTree;
		
	}
	/**
	 *  apply BFS  for Crossover operation
	 *  The first:  choose the representation vertex for each cluster
	 *  The second: apply BFS to generate  the tree for each cluster
	 *  The third: apply BFS for big Tree then  copy the tree of cluster  and big tree to Tree.  
	 * @param father
	 * @param mother
	 * @param num_vertex
	 * @param clusters
	 * @return Tree
	 */
	public double[][] ClusterBFSCrossover(double[][] father, double[][] mother, int num_vertex, ArrayList<Cluster> clusters ){
		
		
		
		Random r = new Random();
		double[][] Tree = new double[num_vertex][num_vertex];
		double[][] clusterWeightMatrix1, clusterWeightMatrix2;
		double[][] spanningTreeOfCluster;
		int numberClusterVertex = 0;
		int numberOfCluster = clusters.size();
		// each Cluster is presented by one vertex in that cluster
		int[] presentationVertex = new int[numberOfCluster]; 
		int[] indexPresentationVertex = new int[numberOfCluster];
		// choose randomly a vertex in each cluster
		
		for( int i = 0; i < numberOfCluster; i ++){
			numberClusterVertex = clusters.get(i).getCluster().size();
			// we choose root is presentation vertex  for cluster which contains root vertex 
			if( i == this.Id_Cluster()){
				for( int j = 0; j < clusters.get(i).getCluster().size(); j++){
				      if(clusters.get(i).getCluster().get(j)  == ReadFiles.root){
				    	  indexPresentationVertex[i] = j;
				      }
				}
				}
		
			else{
			indexPresentationVertex[i] = r.nextInt(numberClusterVertex);
			}
			presentationVertex[i]  =  clusters.get(i).getCluster().get(indexPresentationVertex[i]);
			

			clusterWeightMatrix1 = initilizeChromosome.buildClusterWeightMatrix(father, clusters.get(i).getCluster());
			clusterWeightMatrix2 = initilizeChromosome.buildClusterWeightMatrix(mother, clusters.get(i).getCluster());
			
			spanningTreeOfCluster = BFSCrossover(clusterWeightMatrix1, clusterWeightMatrix2,
					numberClusterVertex, indexPresentationVertex[i]);
			
		    
			// copy to graph 
			for( int j = 0 ; j < numberClusterVertex ; j++){
				for( int k = 0; k < numberClusterVertex; k++){
					Tree[clusters.get(i).getCluster().get(j)][clusters.get(i).getCluster().get(k)] = spanningTreeOfCluster[j][k];
						}
			}	
			
		}
		// build the vertex representation for each cluster
	    clusterWeightMatrix1 = findSpanningTreeBetweenClusters(father, num_vertex, clusters);
	    clusterWeightMatrix2 = findSpanningTreeBetweenClusters(mother, num_vertex, clusters);
	    
	    // generate the new offspring
//	    spanningTreeOfCluster = BFSCrossover(clusterWeightMatrix1, clusterWeightMatrix1,
//	    		numberOfCluster, indexPresentationVertex[9]);
//		
		
	    spanningTreeOfCluster = primRSTcrossover(clusterWeightMatrix1, clusterWeightMatrix1,
	    		numberOfCluster);
	    
	    // convert to spanning tree of Graph
	    int[] indexCluster = new int[numberOfCluster];
	    for(int i = 0; i < numberOfCluster; i++ ){
	    	int position = indexPresentationVertex[i];
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
	/**
	 *  find out the index of cluster which contains root vertex.
	 * @return id cluster
	 */
	public int Id_Cluster(){
		 int id = 0;
		 boolean b = true;
		 for( int i = 0; i <ReadFiles.clusters.size(); i++ ){
			int clusterVertex = ReadFiles.clusters.get(i).getCluster().size();
			 for( int j = 0; j< clusterVertex; j ++ ){
				 if( ReadFiles.root == ReadFiles.clusters.get(i).getCluster().get(j)){
					 id = i;
					 b = false;
					break;
				 }
			 }
			 if( b == false ) break;
		 }
		 return id;
	 }
	/**\
	 * 
	 * @param crss_Para
	 * @param par_1
	 * @param par_2
	 * @param num_Genens
	 * @param rnd
	 * @return
	 */
	
	
	public Individual MFOCrossover(double crss_Para, Individual par_1, Individual par_2, int num_Genens, Random rnd) {
		Individual child = new Individual();
		child.setGene(primRSTcrossover(par_1.getGene(), par_2.getGene(), num_Genens));
		return child;
	}
	
	

   }