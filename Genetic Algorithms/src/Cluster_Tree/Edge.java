package Cluster_Tree;

public class Edge {
	int startVertice;
	int endVertice;
	public Edge(int s_vertice, int e_vertice){
		startVertice = s_vertice;
		endVertice = e_vertice;
	}
	//getter 
	public int getStartVertice(){
		return startVertice;
	}
	//setter 
	public void setStartVertice(int value){
		this.startVertice = value;
	}
	//getter
	public int getEndVertice(){
		return endVertice;
	}
	// setter 
	public void setEndVertice(int value){
		this.endVertice = value;
	}

}
