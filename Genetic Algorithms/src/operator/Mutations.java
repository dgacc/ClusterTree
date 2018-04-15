package operator;

import java.util.ArrayList;
import java.util.Random;

import javax.swing.JFrame;

import dislay.Paint;
import dislay.Windows;
import filesinout.ReadFiles;
import random.MyRandom;
import structures.Cluster;

public class Mutations {
	private InitializeChromosome chromosome = new InitializeChromosome();
	private GraphMethods graphMethods = new GraphMethods();
    
	
	/**
	 * 
	 * @param parents
	 * @param num_vertex
	 * @return
	 */
	public double[][] mutationForEachClusters(double[][] parents,int num_vertex , Random r){
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
		
		offspring[path.get(index1)][path.get(index2)] = 0.0f;
		offspring[path.get(index2)][path.get(index1)] = 0.0f;
		
		offspring[startVertex][endVertex] = 1.0f;
		offspring[endVertex][startVertex] = 1.0f;
//		
		return offspring;
		
	}
	
	
	public double[][] mutationForEachClusters1(double[][] parents, int minNum_vertex, int maxNum_vertex, Random r){
		double[][] offspring = new double[maxNum_vertex][maxNum_vertex];
		for(  int i = 0; i < maxNum_vertex; i++){
			for( int j = 0; j < maxNum_vertex; j++){
				offspring[i][j] = parents[i][j]; 
			}
		}  
		
		
		int startVertex = r.nextInt(maxNum_vertex);
		int endVertex = r.nextInt(maxNum_vertex);
		
		while((startVertex == endVertex) || (offspring[startVertex][endVertex] > 0 ))
		// two vertices are equal or adjacence > 0 then choose two vertices again
		{
			startVertex = r.nextInt(maxNum_vertex);
			endVertex = r.nextInt(maxNum_vertex);
		}
		
		
		// find the path from start_Vertice to end_vertice
		 // initialize  two matrix visited  matrix = false and pre Matrix = -1;
		boolean[] visited = new boolean[maxNum_vertex];
		int[] pre = new int[maxNum_vertex];
		
	
		for( int i = 0; i < maxNum_vertex; i ++ ){
			visited[i] = false;
			pre[i] = -1;
		}
				
	   chromosome.findCycle(startVertex, endVertex, maxNum_vertex, offspring, visited, pre);
		ArrayList<Integer> path = graphMethods.printPath(startVertex, endVertex, pre);
		// delete a edge from cycle
		int index1 = 0;
		int index2 = 0;
		if(startVertex < minNum_vertex && endVertex <  minNum_vertex){
			 index1 = r.nextInt(path.size() - 1);
			 index2 =  index1 + 1;
			
		}else{
			
			index1 = r.nextInt(path.size());
			while((path.get(index1)  < minNum_vertex) ){
				index1 = r.nextInt(path.size());	
				}
			if(index1 == 0){
			index2 = index1 + 1;
			}else if(index1 == path.size() - 1){
				index2 = index1 - 1;
			}else{
				index2 = index1 + 1; 
				
			}
		}
			
		
		offspring[path.get(index1)][path.get(index2)] = 0.0f;
		offspring[path.get(index2)][path.get(index1)] = 0.0f;
		
		offspring[startVertex][endVertex] = 1.0f;
		offspring[endVertex][startVertex] = 1.0f;
		
		return offspring;
		
	}
	
	
	/**
	 * 
	 * @param parent
	 * @param num_vertex
	 * @param clusters
	 * @return
	 */
	public double[][] mutationClusterTree(double[][] parent, int num_vertex, ArrayList<Cluster> clusters, 
			int[] minClusterVertices, int[] maxClusterVertices, Random r, double muRate ){
		Crossover crossover = new Crossover();
		InitializeChromosome initializeChromosome = new InitializeChromosome();
		double[][] offspring = new double[num_vertex][num_vertex];
		int numberOfCluster  = clusters.size();
		int indexCluster = r.nextInt(numberOfCluster);
		double[][] presentationTree = new double[numberOfCluster][numberOfCluster];
		
		for( int i = 0 ; i < num_vertex ; i ++){
			for( int j = 0; j < num_vertex; j++){
				offspring[i][j] = parent[i][j];
			}
		}
		int numGene = (int)(1.0/muRate);
		for(int t= 0; t < numGene; t++){
		if(MyRandom.r.nextDouble() < muRate){
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

		clusterSpanningTree = mutationForEachClusters1(clusterWeightMatrix, 
				minClusterVertices[indexCluster], maxClusterVertices[indexCluster], r);
	
		
		for( int i = 0; i < numberClusterVertex; i++ ){
			for( int j = 0; j < numberClusterVertex; j++){
				offspring[clusters.get(indexCluster).getCluster().get(i)]
						[clusters.get(indexCluster).getCluster().get(j)] = clusterSpanningTree[i][j];
			}
		}
		
		presentationTree =  crossover.findSpanningTreeBetweenClusters(parent, num_vertex, clusters);
		
		
//		Windows windows = new Windows();
//		windows.runWindow("cay dai dien");
//		Paint  p = new Paint();
//		 p.setPaint(presentationTree, ReadFiles.vertices1, ReadFiles.clusters1, numberOfCluster, 
//		  0, 0, ReadFiles.root1);
//		 windows.addPaint(p);
//		
		
		
		presentationTree = mutationForEachClusters1(presentationTree, 
				ReadFiles.numberOfCluster, ReadFiles.numberOfCluster1, r);
	
		int[] indexCluster1 = new int[numberOfCluster];
		indexCluster1 = findPresentationVertex(parent, num_vertex, clusters);
     	for( int j = 0; j < numberOfCluster; j++){
	    		for( int k = 0; k < numberOfCluster; k++){
//	    			if(presentationTree[j][k] > 0){
	    				offspring[indexCluster1[j]][indexCluster1[k]] = presentationTree[j][k];
	    		
	    			}
	    	}
		}
		}
		return offspring;
	}
	
	public double[][] mutationClusterTreeGA(double[][] parent, int num_vertex, ArrayList<Cluster> clusters, Random r, double muRate){
		Crossover crossover = new Crossover();
		InitializeChromosome initializeChromosome = new InitializeChromosome();
		double[][] offspring = new double[num_vertex][num_vertex];
		int numberOfCluster  = clusters.size();
		
		double[][] presentationTree = new double[numberOfCluster][numberOfCluster];
		
		for( int i = 0 ; i < num_vertex ; i ++){
			for( int j = 0; j < num_vertex; j++){
				offspring[i][j] = parent[i][j];
			}
		}
		int temp= ReadFiles.num_vertex;
		for(int t = 0; t < temp; t++){
		   if(r.nextDouble() < muRate){
		int indexCluster = r.nextInt(numberOfCluster);
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
		clusterSpanningTree = mutationForEachClusters(clusterWeightMatrix, numberClusterVertex, r);
//		
		
		for( int i = 0; i < numberClusterVertex; i++ ){
			for( int j = 0; j < numberClusterVertex; j++){
				offspring[clusters.get(indexCluster).getCluster().get(i)]
						[clusters.get(indexCluster).getCluster().get(j)] = clusterSpanningTree[i][j];
			}
		}
		
		presentationTree =  crossover.findSpanningTreeBetweenClusters(parent, num_vertex, clusters);
		
		presentationTree = mutationForEachClusters(presentationTree, ReadFiles.numberOfCluster, r);
		int[] indexCluster1 = new int[numberOfCluster];
		indexCluster1 = findPresentationVertex(parent, num_vertex, clusters);
     	for( int j = 0; j < numberOfCluster; j++){
	    		for( int k = 0; k < numberOfCluster; k++){
	    				offspring[indexCluster1[j]][indexCluster1[k]] = presentationTree[j][k];
	    		
	    			}
	    	}
			}
		}
			
		return offspring;
	}
	
	/**
	 * 
	 * @param parent
	 * @param num_vertex
	 * @param clusters
	 * @return
	 */
	
	public int[] findPresentationVertex(double[][] parent, int num_vertex, ArrayList<Cluster> clusters ){
		int numberOfCluster = clusters.size();
		int[] presentationvertex = new int[numberOfCluster];
		int numberClusterVertex1 = 0;
		int numberClusterVertex2 = 0; 
		    for( int i = 0; i < numberOfCluster; i++){
			numberClusterVertex1 = clusters.get(i).getCluster().size();
			for( int j = 0; j < numberClusterVertex1; j++){
				for( int k = i + 1; k < numberOfCluster; k++ ){
					numberClusterVertex2 = clusters.get(k).getCluster().size();
					
					for( int t = 0; t < numberClusterVertex2; t++){
					if(parent[clusters.get(i).getCluster().get(j)][clusters.get(k).getCluster().get(t)] > 0){
						presentationvertex[i] = clusters.get(i).getCluster().get(j);
						presentationvertex[k] = clusters.get(k).getCluster().get(t);
					}
					}
				}
			}
		}
		return presentationvertex;
	}
	
	
	public int[] mutationPrufer(ArrayList<Cluster> clusters, Random rnd, int[] pruferNumberCode, int[] startpoint, double muRate){
		int numberOfClusters = clusters.size();
		int geneNum = (int)(1.0/muRate);
		int[] temp = pruferNumberCode;
		for(int i = 0; i< geneNum; i++){
		int indexCluster  = rnd.nextInt(numberOfClusters);
		int numberClusterVertices = clusters.get(indexCluster).getCluster().size();
	
	
		if(rnd.nextDouble() < muRate){
		if(indexCluster == 0){
		 temp =  swap( rnd.nextInt(numberClusterVertices -2), 
			    	rnd.nextInt(numberClusterVertices -2), pruferNumberCode);
		
		}else{
	     temp =  swap( (startpoint[indexCluster - 1] + rnd.nextInt(numberClusterVertices -2)), 
	    		(startpoint[indexCluster - 1] + rnd.nextInt(numberClusterVertices -2)), pruferNumberCode);
		}
		 temp[startpoint[numberOfClusters - 1] + numberOfClusters - 2 + indexCluster ] = 
				 rnd.nextInt(InitializeChromosome.minClusterVertices[indexCluster]);
	 
		}
		}
		return temp;
	}
	/**
	 * 
	 * @param clusters
	 * @param rnd
	 * @param pruferNumberCode
	 * @param startpoint
	 * @return
	 */

	public int[] mutationPrufer1(ArrayList<Cluster> clusters, Random rnd, int[] pruferNumberCode, int[] startpoint, double muRate){
		int numberOfClusters = clusters.size();
		int geneNum = (int)(1.0/muRate);
		int[] temp = pruferNumberCode;
		for(int i = 0; i< geneNum; i++){
		int indexCluster  = rnd.nextInt(numberOfClusters);
		int numberClusterVertices = clusters.get(indexCluster).getCluster().size();
		if(indexCluster == 0){
		 temp =  swap( rnd.nextInt(numberClusterVertices -2), 
			    	rnd.nextInt(numberClusterVertices -2), pruferNumberCode);
		
		}else{
	     temp =  swap( (startpoint[indexCluster - 1] + rnd.nextInt(numberClusterVertices -2)), 
	    		(startpoint[indexCluster - 1] + rnd.nextInt(numberClusterVertices -2)), pruferNumberCode);
		}
		 temp[startpoint[numberOfClusters - 1] + numberOfClusters - 2 + indexCluster ] = 
				 rnd.nextInt(numberClusterVertices);
		}
		return temp;	
	}
	public int[] swap(int a, int b, int[] pruferNumberCode) {
		int temp = pruferNumberCode[a];
		 pruferNumberCode[a] =  pruferNumberCode[b];
		 pruferNumberCode[b] =  temp;
		 return pruferNumberCode;
	}
	

}

