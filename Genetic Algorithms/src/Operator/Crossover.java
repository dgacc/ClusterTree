package Operator;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Queue;
import java.util.Random;
import java.util.Stack;

import javax.swing.JFrame;

import Dislay.Paint;
import Files_InOut.ReadFiles;
import Structures.Cluster;
import Structures.Individual;


public class Crossover {
	InitializeChromosome initilizeChromosome = new InitializeChromosome();
	
	
	/**
	 * 
	 * @param parent
	 * @param num_vertex
	 * @param clusters
	 * @return the spanning tree between clusters 
	 */
	public double[][] findSpanningTreeBetweenClusters(double[][] parent, int num_vertex, ArrayList<Cluster> clusters ){
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
	public double[][] clusterCrossover(double[][] father,double[][] mother, int num_vertex, ArrayList<Cluster> clusters , Random r){
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
		    for(int i = 0; i < numberOfCluster; i++ ){
		    	int position = r.nextInt(clusters.get(i).getCluster().size());
		    	indexCluster[i] = clusters.get(i).getCluster().get(position);
		    }
		    	for( int j = 0; j < numberOfCluster; j++){
		    		for( int k = 0; k < numberOfCluster; k++){
		    			if( spanningTreeOfCluster[j][k] > 0){
		    				Tree[indexCluster[j]][indexCluster[k]] = spanningTreeOfCluster[j][k];
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
	 * @param r TODO
	 * @return Tree
	 */
	public double[][] ClusterBFSCrossover(double[][] father, double[][] mother, int num_vertex, ArrayList<Cluster> clusters, Random r ){
		
		
		
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
	    }
	    	for( int j = 0; j < numberOfCluster; j++){
	    		for( int k = 0; k < numberOfCluster; k++){
	    			if( spanningTreeOfCluster[j][k] > 0){
	    				Tree[indexCluster[j]][indexCluster[k]] = spanningTreeOfCluster[j][k];
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
	/**
	 * 
	 * @param father
	 * @param mother
	 * @param num_vertex:  this is the number of vertex max.
	 * @param clusters
	 * @param r
	 * @return
	 */
	public double[][] clusterCrossover1(double[][] father,double[][] mother, int num_vertex, ArrayList<Cluster> clusters ,
		 int[] minClusterVertices ,int[] maxClusterVertices, Random r){
		double[][] Tree = new double[num_vertex][num_vertex];
		double[][] clusterWeightMatrix1, clusterWeightMatrix2;
		double[][] spanningTreeOfCluster;
		int numberClusterVertex = 0;
		int numberOfCluster = clusters.size();
		
		// 
		for( int i = 0; i < numberOfCluster; i++){
			numberClusterVertex = maxClusterVertices[i];
			clusterWeightMatrix1 = initilizeChromosome.buildClusterWeightMatrix(father, clusters.get(i).getCluster());
			clusterWeightMatrix2 = initilizeChromosome.buildClusterWeightMatrix(mother, clusters.get(i).getCluster());
			
		    spanningTreeOfCluster = primRSTcrossover(clusterWeightMatrix1, clusterWeightMatrix2,  minClusterVertices[i],maxClusterVertices[i], r);
		
    
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
		    spanningTreeOfCluster = primRSTcrossover(clusterWeightMatrix1, clusterWeightMatrix1,numberOfCluster);
		    // convert to spanning tree of Graph
		    int[] indexCluster = new int[numberOfCluster];
		    for(int i = 0; i < numberOfCluster; i++ ){
		    	int position = r.nextInt(clusters.get(i).getCluster().size());
		    	indexCluster[i] = clusters.get(i).getCluster().get(position);
		    }
		    	for( int j = 0; j < numberOfCluster; j++){
		    		for( int k = 0; k < numberOfCluster; k++){
		    			if( spanningTreeOfCluster[j][k] > 0){
		    				Tree[indexCluster[j]][indexCluster[k]] = spanningTreeOfCluster[j][k];
		    			}
		    		}
		    	}
		    
		    return Tree;
	}
	private double[][] primRSTcrossover(double[][] father, double[][] mother, int minNum_vertex, int maxNum_vertex , Random r){
		 double[][] G_cr = new double[minNum_vertex][minNum_vertex];
		 double[][] spanningTree = new double[minNum_vertex][minNum_vertex];
		 Queue<Integer> queue = new LinkedList();
		 boolean[] mask = new boolean[maxNum_vertex];
		 GraphMethods graphMethods = new GraphMethods();
		 
		 
		 
		 // initialize the G_cr Matrix = fatherMatrix + motherMatrix
		for(int i = 0; i < minNum_vertex; i++){
			for(int j = 0; j< minNum_vertex; j++){
				G_cr[i][j] = father[i][j] + mother[i][j];
				
			}
		}
		
		spanningTree = initilizeChromosome.primRST(G_cr, minNum_vertex);
	
		
		 double[][] G_cr1 = new double[maxNum_vertex][maxNum_vertex];
		 double[][] spanningTree1 = new double[maxNum_vertex][maxNum_vertex];
		
		for(int i = 0; i < maxNum_vertex; i++){
			mask[i] = true;
			for( int j = 0; j < maxNum_vertex; j++){

				if((i < minNum_vertex )&&(j < minNum_vertex)){
					G_cr1[i][j] = spanningTree[i][j];
				}
				else{
				G_cr1[i][j] = father[i][j] + mother[i][j];
			}
		}
	   }
		 // ve
//		JFrame gf = new JFrame();
//		gf.setVisible(true);
//		gf.setSize(800, 800);
//		gf.setTitle(" cay 1");
//		gf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
//		gf.setVisible(true);
//		
//		Paint  p = new Paint();
//		p.weightMatrix = G_cr1;
//		p.num_vertex = 14;
//	    gf.add(p);
//	    gf.setVisible(true);
		
		for(int i = 0; i < minNum_vertex; i++){
			for( int j = minNum_vertex;  j < maxNum_vertex; j++){
				if(G_cr1[i][j] > 0  && mask[i] ){
					queue.add(i);
					mask[i] = false;
				}
		}
		}
		for( int j = minNum_vertex;  j < maxNum_vertex; j++){
				queue.add(j);
			}
		
	
	
		while(!queue.isEmpty()){
			int startVertex = queue.poll();
//			System.out.println("queue:");
		
			
			ArrayList<Integer> path = new ArrayList<Integer>();
			path  = dfs(G_cr1,maxNum_vertex, startVertex);
			
			// delete from cycle
			if( (path.size() > 1)){	
				int index1 = r.nextInt(path.size());
				while((path.get(index1)  < minNum_vertex) ){
					index1 = r.nextInt(path.size());	
					}
				int index2;
				if(index1 == 0){
				index2 = index1 + 1;
				}else if(index1 == path.size() - 1){
					index2 = index1 - 1;
				}else{
					index2 = index1 + 1; 
				}
				
				
				G_cr1[path.get(index1)][path.get(index2)] = 0.0f;
				G_cr1[path.get(index2)][path.get(index1)] = 0.0f;
				
				}
			System.out.println();
		}
		

	
		
//		}
		JFrame gf1 = new JFrame();
		gf1.setVisible(true);
		gf1.setSize(800, 800);
		gf1.setTitle(" cay sau xoa");
		gf1.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		gf1.setVisible(true);
		
		Paint  p1 = new Paint();
		p1.weightMatrix = G_cr1;
		p1.num_vertex = 14;
	    gf1.add(p1);
	    gf1.setVisible(true);
    return G_cr1;
    }
	/**
	 * this method gonna find the cycles based on the BFS algorithms. 
	 * @param weightMatrix
	 * @param num_vertex
	 * @param startVertex
	 * @return the list of path we  use that for delete the edges.
	 */
	 public LinkedList<Integer> findTheCycle(  double[][] weightMatrix, int num_vertex, int startVertex){
		  
			
			boolean[] mark = new boolean[num_vertex];
			LinkedList<Integer> preVertex = new LinkedList<Integer>();
			Stack<Integer> stack = new Stack<Integer>();
		    // initialize label for each vertex,  
			for( int i = 0; i < num_vertex; i++){
				mark[i] = true;
			}
			mark[startVertex] = false;
			stack.push(startVertex);
			
			 while(!stack.isEmpty()){
				 int considerVertex = stack.peek();
				 int i = 0;
				 while (i < num_vertex){
					 if(weightMatrix[considerVertex][i] > 0  && !preVertex.contains(i)){
						 
						 stack.push(i);
						 preVertex.add(i);
						 mark[i] = false;
						 considerVertex = stack.pop();
						 i = 0;
					 }	else{
						 i++;
					 }
				 }
				 stack.pop();
			 }
		return preVertex;
		}
	  public  ArrayList<Integer> dfs(double[][] g_cr1,   int number_of_nodes ,int source)
	    {
		    Stack<Integer> stack = new Stack<Integer>();
		    int[] preVertex = new int[number_of_nodes];
			ArrayList<Integer> path = new ArrayList<Integer>();
			ArrayList<Integer> TreeDFS = new ArrayList<Integer>();
	        int visited[] = new int[number_of_nodes];		
	        int element = source;		
	        int i = source;	
	        boolean flag = true;
	        System.out.print(element + "\t");		
	        visited[source] = 1;		
	        stack.push(source);
	      
	        TreeDFS.add(source);
	        for( int k = 0; k < number_of_nodes; k++){
	        	preVertex[k] = -1;
	        }
	        preVertex[source] = source;
	 
	        while (!stack.isEmpty() && flag)
	        {
	            element = stack.peek();
	            i = 0;	
		    while (i < number_of_nodes)
		    {
	     	        if (g_cr1[element][i] == 1 && visited[i] == 0)
		        {       preVertex[i] = element;
	                    stack.push(i);
	                    visited[i] = 1;
	                    element = i;
	                    TreeDFS.add(element);
	                    for( int j = 0; j< number_of_nodes; j++){
	                    	if((preVertex[j] != -1 )&&(g_cr1[preVertex[j]][element] == 1) && (preVertex[j] != preVertex[element])){
	                    		int temp = element;
	                    		path.add(element);
	                    		while( preVertex[temp] != preVertex[j]){
	                    			path.add(preVertex[temp]);
	                    			temp = preVertex[temp];
	                    		}
	                    		path.add(preVertex[j]);
	                    		flag = false;
	                    	    break;
	                    	}
	                    	if(!flag){
	                    		break;
	                    	}
	                    }
	                    i = 0;
	                    System.out.print(element + "\t");;
        
	                }else{
	                i++;
	                }
		    }
	            stack.pop();         
	        }
	    	System.out.println();
	        for( int k = 0; k < number_of_nodes; k++){
	        
	        	System.out.print(preVertex[k] + "\t");
	        }
	        return path;
	    }
	 
}