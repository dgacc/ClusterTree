package filesinout;

//import java.util.ArrayList;
//import java.util.Random;
//
//import structures.Edge;

public class bug {
	///*-----------------------------------------Build the new instance--------------------------------------------------------------*/
//	//optimal tree
//	 double[][] optimalTree = new double[num_vertex][num_vertex];
//	 // build the weight matrix from optimalTree
//	 double[][] weightMatrix1 = new double[num_vertex][num_vertex]; 
//	 Individual individual1 = new Individual();
//	 individual1.inintilizeIndividual();
//	 optimalTree = individual1.getGene();
//	 ArrayList<Cluster> clusters = ReadFiles.clusters;
//	 
//	 for(int i = 0; i < clusters.size(); i ++){
//			int numberClusterVertex = clusters.get(i).getCluster().size();
//			for( int j = 0; j < numberClusterVertex; j++ ){
//				for( int k = 0; k < numberClusterVertex; k++){
//					double value = 11 + r.nextInt(10);
//					weightMatrix1[clusters.get(i).getCluster().get(j)][clusters.get(i).getCluster().get(k)] = value;
//					weightMatrix1[clusters.get(i).getCluster().get(k)][clusters.get(i).getCluster().get(j)] = value;
//				}
//			}
//		}
//	 for( int i = 0; i < num_vertex; i++){
//		 for(int j = 0; j < num_vertex; j++){
//			if(optimalTree[i][j] > 0){
//			weightMatrix1[i][j] = optimalTree[i][j];
//			
//			}
//		 }
//	 }
//	 for( int i = 0; i < num_vertex; i++){
//		 for(int j = 0; j < num_vertex; j++){
//			if(weightMatrix1[i][j] == 0){
//				double value = 26 + r.nextInt(75);
//			    weightMatrix1[i][j] = value; 
//			}
//		 }
//	 }
//	 
//	 for( int i = 0; i < num_vertex; i ++){
//		 for( int j = 0; j < num_vertex; j++){
//			System.out.print(" " +weightMatrix1[i][j]);
//		 }
//		 System.out.println();
//	 }
//	 
//	 Evaluation evaluation = new Evaluation();
////	 Paint  p = new Paint();
////		p.weightMatrix = optimalTree;
////		p.fitness = evaluation.distanceEvaluate(weightMatrix1, optimalTree, num_vertex, ReadFiles.root);
////		p.num_vertex = num_vertex;
////	    gf.add(p);
////	    gf.setVisible(true);
//	
//	   // print files:
//	    CreateFile g = new CreateFile();
//	    g.cost =  evaluation.distanceEvaluate(weightMatrix1, optimalTree, num_vertex, ReadFiles.root);
//	    g.weightMatrix = weightMatrix1;
//     	g.writeFile();
//     	
//	    
///*--------------------------------------------------------END-------------------------------------------------------------*/

//	/**
//	 *  this is the algorithm to build the  weight matrix if have the edge from  cluster to cluster 
//	 * @param weightMatrix
//	 * @param num_vertices
//	 * @return the tree with weight matrix;
//	 */ 
//	public double[][] bigPrimRST(double[][] weightMatrix, int num_vertices) {
//		int randomVertice = -1, indexRandomEdge = -1;
//		ArrayList<Integer> contain = new ArrayList<Integer>();
//		ArrayList<Edge> edgeListAdjacence = new ArrayList<Edge>();
//		double[][] tempTable = new double[num_vertices][num_vertices];
//
//		randomVertice = a.nextInt(num_vertices);
//		contain.add(randomVertice);
//		edgeListAdjacence.addAll(findEdge(num_vertices, randomVertice, contain, weightMatrix));
//
//		while (contain.size() < num_vertices) {
//
//			if (contain.size() == 0) {
//				return null;
//			}
//
//			//indexRandomEdge = r.nextInt(contain.size());
//			indexRandomEdge = a.nextInt(edgeListAdjacence.size());
//			Edge tempEdge = new Edge(edgeListAdjacence.get(indexRandomEdge).startVertice,
//					edgeListAdjacence.get(indexRandomEdge).endVertice);
//			edgeListAdjacence.remove(indexRandomEdge);
//
//			if (!contain.contains(tempEdge.endVertice)) {
////
////				tempTable[tempEdge.startVertice][tempEdge.endVertice] = 1.0f;
////				tempTable[tempEdge.endVertice][tempEdge.startVertice] = 1.0f;
//				// Initialize the weightMatrix 
//				Random r = new Random();
//				double value = (double) (20 + r.nextInt(6));
//				tempTable[tempEdge.startVertice][tempEdge.endVertice] = value;
//				tempTable[tempEdge.endVertice][tempEdge.startVertice] = value;
//				// End 
//
//				contain.add(tempEdge.endVertice);
//				edgeListAdjacence.addAll(findEdge(num_vertices, tempEdge.endVertice, contain, weightMatrix));
//			}
//		}
//		return tempTable;
//	}

	
}
