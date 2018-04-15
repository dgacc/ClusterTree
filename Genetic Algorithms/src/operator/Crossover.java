package operator;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Queue;
import java.util.Random;
import java.util.Stack;

import javax.swing.JFrame;

import dislay.Paint;
import dislay.Windows;
import filesinout.ReadFiles;
import random.MyRandom;
import structures.Cluster;
import structures.Individual;


public class Crossover {
	private InitializeChromosome initilizeChromosome = new InitializeChromosome();
    private Mutations mutations = new Mutations();
	
	
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
 * @param father : tree
 * @param mother : tree
 * @param num_vertex : number of vertices of tree
 * @param maxClusters : clusters 
 * @param numberOfCluster1 : numberOfCluster1  < numberOfCluster2
 * @param numberOfCluster2 : numberOfCluster1  < numberOfCluster2
 * @param minClusterVertices  : //array of number of vertex in cluster have less than  the other.
 * @param maxClusterVertices   //array of number of vertex in cluster have greater than  the other.

 * @param r
 * @return
 */
	public double[][] clusterCrossover1(double[][] father,double[][] mother, int num_vertex,
			ArrayList<Cluster> maxClusters ,int numberOfCluster1, int numberOfCluster2,
		int[] minClusterVertices ,int[] maxClusterVertices, Random r , Windows windows1){
		double[][] Tree = new double[num_vertex][num_vertex];
		double[][] clusterWeightMatrix1, clusterWeightMatrix2;
		double[][] spanningTreeOfCluster;
		int numberClusterVertex = 0;
		int numberOfCluster = maxClusters.size();
		

		
		// 
		for( int i = 0; i < numberOfCluster; i++){
			numberClusterVertex = maxClusterVertices[i];
			clusterWeightMatrix1 = initilizeChromosome.buildClusterWeightMatrix(father, maxClusters.get(i).getCluster());
			clusterWeightMatrix2 = initilizeChromosome.buildClusterWeightMatrix(mother, maxClusters.get(i).getCluster());
			
		    spanningTreeOfCluster = primRSTcrossover(clusterWeightMatrix1, clusterWeightMatrix2,  minClusterVertices[i],maxClusterVertices[i], r);
		   
		    //debug
//		    windows1.runWindow(" cluster thu " + i);
//			Paint  p = new Paint();
//			p.setPaint(spanningTreeOfCluster, ReadFiles.vertices, ReadFiles.clusters, minClusterVertices[i],0, 0, ReadFiles.root);
//			windows1.addPaint(p);
			
			
			
           
		    // convert to the graph tree
		    for(int j = 0; j < numberClusterVertex; j++){
		    	for(int k = 0; k < numberClusterVertex; k++){
		    		Tree[maxClusters.get(i).getCluster().get(j)][maxClusters.get(i).getCluster().get(k)] = spanningTreeOfCluster[j][k];		
		    	}
		    }
		 
		}
		    // build the vertex representation for each cluster
		    clusterWeightMatrix1 = findSpanningTreeBetweenClusters(father, num_vertex, maxClusters);
		    clusterWeightMatrix2 = findSpanningTreeBetweenClusters(mother, num_vertex, maxClusters);
		  
		    spanningTreeOfCluster = primRSTcrossover(clusterWeightMatrix1, clusterWeightMatrix1,numberOfCluster1, numberOfCluster2, r);
		    // convert to spanning tree of Graph
		    int[] indexCluster = new int[numberOfCluster];
		    for(int i = 0; i < numberOfCluster; i++ ){
		    	int position = r.nextInt(minClusterVertices[i]);
		    	indexCluster[i] = maxClusters.get(i).getCluster().get(position);
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
	 * apply primRSTcrossOver 
	 * @param father 
	 * @param mother
	 * @param minNum_vertex : number of vertices in cluster have number of vertices less than the other. 
	 * @param maxNum_vertex : number of vertices in cluster have number of vertices greater than the other.
	 * @param r
	 * @return tree after apply  primRSTcrossover
	 */
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
		// generate  smaller tree from G-cr matrix which
		spanningTree = initilizeChromosome.primRST(G_cr, minNum_vertex);
	
		
		 double[][] G_cr1 = new double[maxNum_vertex][maxNum_vertex];
		 
		 //copy the  small  tree to big graph. 
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
		
		
		// Initialize label  and add  vertex  in smaller tree
		// which is connected to the others vertices aren't in  smaller tree  to queue
		for(int i = 0; i < minNum_vertex; i++){
			for( int j = minNum_vertex;  j < maxNum_vertex; j++){
				if(G_cr1[i][j] > 0  && mask[i] ){
					queue.add(i);
					mask[i] = false;
				}
		}
		}
		
		
		// add all of vertices which aren't in smaller tree to queue
 		for( int j = minNum_vertex;  j < maxNum_vertex; j++){
				queue.add(j);
				mask[j] = false;
			}
		
	
       // consider vertex in queue to find the cycle and delete a edge 
 	   //in cycles and generate bigger tree	
		while(!queue.isEmpty()){
			int startVertex = queue.poll();		
			
			ArrayList<Integer> path = new ArrayList<Integer>();
			path  = findCycleDFS(G_cr1,maxNum_vertex, startVertex);
			
			// choose edge is not in tree  then delete that.
			
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
				
				// delete the edge in that cycle
				
				G_cr1[path.get(index1)][path.get(index2)] = 0.0f;
				G_cr1[path.get(index2)][path.get(index1)] = 0.0f;
				
				}
		}
    return G_cr1;
    }
	
	
      /**
       *  find  cycle in graph using DFS algorithm,
       *  we traversal tree by Using the BFS,
       *  if we have the reserving edges in graph then we have the cycle in that graph, then stop Traversal Tree
       *  after find the cycle already yet, we  use the previous vertex of the end vertex of that  cycle
       *  to find  the others vertex based on the previous vertex( recur), over and over again until 
       *  the start vertex of that cycle, then return the path.
       * 
       * @param Graph : that graph  which have cycle inside
       * @param num_vertices : the number of vertices in that graph
       * @param startVertex : the start vertex to travel tree
       * @return The path  in cycle 
       */
	
	  public  ArrayList<Integer> findCycleDFS(double[][] Graph,   int num_vertices ,int startVertex)
	    {
		    Stack<Integer> stack = new Stack<Integer>();
		    int[] preVertex = new int[num_vertices];
			ArrayList<Integer> path = new ArrayList<Integer>();
			ArrayList<Integer> TreeDFS = new ArrayList<Integer>();
	        int visited[] = new int[num_vertices];		
	        int element = startVertex;		
	        int i = startVertex;	
	        boolean flag = true;
//	        System.out.print(element + "\t");		
	        visited[startVertex] = 1;		
	        stack.push(startVertex);
	     // initialize the previous vertex for each vertex in graph
	        
	        TreeDFS.add(startVertex);
	        for( int k = 0; k < num_vertices; k++){
	        	preVertex[k] = -1;
	        }
	        preVertex[startVertex] = startVertex;
	 
	        // start to   DFS  tree traversal use stack structure. 
	        
	        while (!stack.isEmpty() && flag)
	        {
	            element = stack.peek();
	            i = 0;	
		    while (i < num_vertices)
		    {
	     	        if (Graph[element][i] > 0 && visited[i] == 0)
		        {       preVertex[i] = element;
	                    stack.push(i);
	                    visited[i] = 1;
	                    element = i;
	                    TreeDFS.add(element);
	                    
	                    // traversal tree DFS  until detect the reserving edges then stop. return the path.
	                    
	                    for( int j = 0; j< num_vertices; j++){
	                    	if((preVertex[j] != -1 )&&(Graph[preVertex[j]][element] > 0) && (preVertex[j] != preVertex[element])){
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
//	                    System.out.print(element + "\t");;
        
	                }else{
	                i++;
	                }
		    }
	            stack.pop();         
	        }
	        return path;
	    }
	 /**
	  *  
	  * @param father
	  * @param mother
	  * @param genLength
	  * @param rnd
	  * @return
	  */
	  public ArrayList<int[]>  pruferNumberCrossover(int[] father, int[] mother,  int genLength, Random rnd){
	  int point = rnd.nextInt(genLength);
	  int[] offspring1 = new int[genLength];
	  int[] offspring2 = new int[genLength]; 
	  ArrayList<int[]> offsprings = new ArrayList<int[]>();
	  for( int i = 0; i < genLength; i ++){
		  if(i < point){
		  offspring1[i] = father[i];
		  offspring2[i] = mother[i];
		  }else{
			  offspring1[i] = mother[i];
			  offspring2[i] = father[i]; 
		  }
	  }
	  
	  offsprings.add(offspring1);
	  offsprings.add(offspring2);
	  return  offsprings;
	  }
	  
	  
	  
	  }
	  
		  
	 
