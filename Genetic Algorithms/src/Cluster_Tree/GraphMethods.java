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
	
}
