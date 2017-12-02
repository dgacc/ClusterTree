package Operator;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Queue;

public class GraphMethods {
	public ArrayList<Integer> printPath( int start_vertice, int end_vertice, int[] pre){
		ArrayList<Integer> tmp_path = new ArrayList<Integer>();
		
		if(pre[end_vertice] == -1){
		}
		else{
			int j = end_vertice;
			while(j != start_vertice){
				tmp_path.add(j);
				j = pre[j];
			}
			tmp_path.add(start_vertice);
		}
		return tmp_path;
	}
	// nguyen duc nghia
	public int[]  dijkstra(double[][] weightMatrix, int num_vertex, int startVertex){
		
		int[] truoc =  new int[num_vertex]; // vertex before
		double[] d =  new double[num_vertex]; // the list of vertex consider
		ArrayList<Integer> Tree = new ArrayList<Integer>(); // the weight Matrix = distances  between startVertex  and the others vertex in consider vertex list
		int u = 0;
		double maxp = 10000000000000000f;
		boolean[] cuoi = new boolean[num_vertex];
		// initialize label
		for( int v = 0; v < num_vertex; v++){
			if( (weightMatrix[startVertex -1][v] <= 0) ||(weightMatrix[startVertex -1][v]) == maxp)
			{
				d[v] = maxp;
			}
			else{
 				d[v] = weightMatrix[startVertex -1][v];
			}
			truoc[v] = startVertex - 1;
			Tree.add(v);
				
		}
		// xong buoc khoi tao 
		// khoi tao cac gia tri cho nut goc
		
		truoc[startVertex -1] = -1; 
		d[startVertex -1] = 0;
		Tree.remove(startVertex - 1);
// buocs lap 
        while (Tree.size() > 0)
        {
            //Tim u la dinh co nhan tam thoi nho nhat
        	double minp = 10000000000000000f;
        	u = -1;
          
            for(int v : Tree){
         
            	if(minp > d[v] ){
            		u = v;
            		minp =d[v]; 
            	}
            }
            Tree.remove(u);
            
            for(int v: Tree){

                if ((d[u] == maxp) || (weightMatrix[u][v] == maxp))
                {
                    break;
                }else{
            	if(((d[u] + weightMatrix[u][v]) < d[v]) && (weightMatrix[u][v] > 0)){
            		d[v] = d[u] + weightMatrix[u][v];
            		truoc[v] = u;
            	}
                }
            }
        }
	return truoc;
	}
	public int[] get_Vertex_In_Each_SubGraph(double[][] weigh_Matrix, int num_Vertex) {
		int[] chua_Xet = new int[num_Vertex];
		int solt = 0;
		for (int i = 0; i < num_Vertex; i++) {
			chua_Xet[i] = -1;
		}
		for (int i = 0; i < num_Vertex; i++) {
			if (chua_Xet[i] == -1) {

				BFS(weigh_Matrix, num_Vertex, chua_Xet, solt, i);
				solt++;
			}
		}

		int[] vertex_In_SubGraph = new int[solt];
		for (int i = 0; i < solt; i++) {
			for (int j = 0; j < num_Vertex; j++) {
				if (chua_Xet[j] == i) {
					vertex_In_SubGraph[i] = j;
				}
			}
		}

		return vertex_In_SubGraph;
	}
	
	
	public int[] BFS(double[][] weigh_Matrix, int num_Vertex, int[] chua_Xet, int solt, int start_Vertex) {
		int[] QUEUE = new int[num_Vertex];
		// bool[] chua_Xet = new bool[num_Vertex];
		int[] truoc = new int[num_Vertex];// ko can vi ko can tim duoc di
		int u, dauQ, cuoiQ;

		for (int i = 0; i < num_Vertex; i++) {
			// chua_Xet[i] = true;
			truoc[i] = start_Vertex;
			QUEUE[i] = -1;
		}
		truoc[start_Vertex] = -1;
		dauQ = 0;
		cuoiQ = 0;
		QUEUE[cuoiQ] = start_Vertex;
		chua_Xet[start_Vertex] = solt;

		while (dauQ <= cuoiQ) {
			u = QUEUE[dauQ];
			dauQ++;
			for (int i = 0; i < num_Vertex; i++) {

				if (weigh_Matrix[u][i] == Double.MAX_VALUE) {
					continue;
				} else {
					if ((weigh_Matrix[u][i] > 0) && (chua_Xet[i] == -1)) {
						cuoiQ++;
						QUEUE[cuoiQ] = i;
						chua_Xet[i] = solt;
						truoc[i] = u;
					}
				}
			} // for
		} // while
		return truoc;
	}

	//
	

}
