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
    
    
    /**
     * 
     * @param weightMatrix
     * @param tree
     * @param num_vertex
     * @param startVertex
     * @return
     */
    public double  distanceEvaluate(double[][] weightMatrix, double[][] tree, int num_vertex, int startVertex){
     int considerVertex;
     double sum = 0;
     boolean[] mark = new boolean[ num_vertex];
     Queue<Integer> queue =  new LinkedList<>();
     // initialize for : distances root,  mark = true, and etc.....
     for( int  i = 0; i < num_vertex; i++){
    	 mark[i] = true;
     }
     queue.add(startVertex);
     mark[startVertex] = false;
     while(!queue.isEmpty()){
    	 considerVertex = queue.poll(); 
    	for( int i= 0; i< num_vertex; i++){
    		if( tree[considerVertex][i] > 0 && mark[i]){
    			sum = sum + weightMatrix[considerVertex][i];
    			queue.add(i);
    			mark[i] = false;
    		}
    	}
    	
     }
    return sum;
    }
	/*********************************************************************************************************************************************
	 * Decoding cho bai toan cay khung Giam tu n dinh ve thanh m dinh 01. Xóa
	 * các đỉnh và cạnh liên thuộc với nó > m 02. Tìm các thành phần liên thông
	 * 03. Nối các thành phần liên thông thứ i -> i + 1 + Chọn ngẫu nhiên ở mỗi
	 * thành phần liên thông 1 đỉnh + Chọn đỉnh đầu tiên
	 ********************************************************************************************************************************************/
	public  double[][] decodingMFOVertexInSubGraph(double[][] ind_Matrix, int max_Genes, int num_Gen_of_Task_j) {
		double[][] tmp_Matrix = new double[max_Genes][max_Genes];
		// 01. Tìm các đỉnh lớn hơn num_Gen_of_Task_j và xóa nó cùng cạnh liên
		// thuộc
		for (int i = 0; i < num_Gen_of_Task_j; i++) {
			for (int j = 0; j < num_Gen_of_Task_j; j++) {
				if ((ind_Matrix[i][j] > 0) && (ind_Matrix[i][j] < Double.MAX_VALUE)) {
					tmp_Matrix[i][j] = 1.0f;
				} else {
					tmp_Matrix[i][j] = 0.0f;
				}
			}
		}
		// 02. Tìm các thành phần liên thông va lay moi tp lien thong 1 dinh
		int[] tp_LT = graphMethodsClass.get_Vertex_In_Each_SubGraph(tmp_Matrix, num_Gen_of_Task_j);

		// 03. Nối các thành phần liên thông thứ i -> i + 1
		for (int i = 0; i < tp_LT.length - 1; i++) {
			tmp_Matrix[tp_LT[i]][tp_LT[i + 1]] = 1.0f;
			tmp_Matrix[tp_LT[i + 1]][tp_LT[i]] = 1.0f;
		}

		// Tao ra ma tran ket qua
		double[][] final_Matrix = new double[num_Gen_of_Task_j][num_Gen_of_Task_j];
		for (int i = 0; i < num_Gen_of_Task_j; i++) {
			for (int j = 0; j < num_Gen_of_Task_j; j++) {
				final_Matrix[i][j] = (int) tmp_Matrix[i][j];
			}
		}
		return final_Matrix;
	}
}

