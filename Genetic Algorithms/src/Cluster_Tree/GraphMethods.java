package Cluster_Tree;

import java.util.ArrayList;

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
		int[] truoc =  new int[num_vertex];
		double[] d =  new double[num_vertex];
		ArrayList<Integer> Tree = new ArrayList<Integer>();
		int u = 0;
		double minp = 10000000000000000f;
		
		// initialize lable
		for( int v = 0; v < num_vertex; v++){
			if( weightMatrix[startVertex][v] > 0){
				d[v] = weightMatrix[startVertex][v];
			}
			truoc[v] = startVertex;
			Tree.add(v);
				
		}
		
		truoc[startVertex] = -1; 
		d[startVertex] = 0;
		Tree.remove(startVertex);

        while (Tree.size() > 0)
        {
            //Tim u la dinh co nhan tam thoi nho nhat
            u = -1;
          
            for(int v : Tree){
         
            	if(minp > d[v] ){
            		u = v;
            		minp =d[v]; 
            	}
            }
            Tree.remove(u);
            
            for(int v: Tree){
            	if((d[u] + weightMatrix[u][v] < d[v]) && (weightMatrix[u][v] > 0)){
            		d[v] = d[u] + weightMatrix[u][v];
            		truoc[v] = u;
            	}
            }
        }
	return truoc;
	}
	
	
}
