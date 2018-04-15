package clustertree;


import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Queue;
import java.util.Scanner;
import javax.swing.JFrame;

import dislay.Paint;
import filesinout.ReadFiles;
import structures.Cluster;



public class Heuristic_CTMS {
 int root = ReadFiles.root;
 ArrayList<Cluster> clusters = ReadFiles.clusters;
 double[][] weightMatrix = ReadFiles.weightMatrix;
 int[][] tree;
  int numberVertices = ReadFiles.num_vertex;


 public int Id_Cluster(){
	 int id = 0;
	 boolean b = true;
	 for( int i = 0; i <clusters.size(); i++ ){
		int clusterVertex = clusters.get(i).getCluster().size();
		 for( int j = 0; j< clusterVertex; j ++ ){
			 if( root == clusters.get(i).getCluster().get(j)){
				 id = i;
				 b = false;
				break;
			 }
		 }
		 if( b == false ) break;
	 }
	 return id;
 }

 public int getId(ArrayList<Integer> list){
     double min = 1000000000000000f;
     int id = 0;
	 for( int i = 0; i < list.size(); i++){
		 double d = 0;
		 for( int j = 0; j < list.size(); j ++){
			 d+= weightMatrix[list.get(i)][list.get(j)];
		 }
		 d = d + list.size()*weightMatrix[root][list.get(i)];
		 if( d < min){
			 min = d;
			 id = i;
		 }
	 }
	 return list.get(id);
 }


public void makeTree(){
 System.out.println(numberVertices);
 tree = new int[numberVertices][numberVertices];
 int id = Id_Cluster();

 for( int i = 0 ; i <clusters.get(id).getCluster().size(); i ++ ){
	 tree[root][clusters.get(id).getCluster().get(i)]  = 1;
 }
 for( int i = 0; i < clusters.size(); i ++){
	 if(i != id){
		 int pointMin = getId(clusters.get(i).getCluster());
		 tree[root][pointMin] = 1;
		 for( int j = 0; j < clusters.get(i).getCluster().size(); j++){
			 int u = clusters.get(i).getCluster().get(j);
			 tree[pointMin][u] = 1;
		 }
	 }
 }
}
 public double solve(){
	 this.makeTree();
//	 Evaluation evaluation = new Evaluation();
	 double cost = this.evaluation(weightMatrix, this.readFilesTree(), numberVertices, root);
	 System.out.println(" cost  = " + cost);
	 return cost;


//	    JFrame gf = new JFrame();
//		gf.setVisible(true);
//		gf.setSize(800, 800);
//		gf.setTitle("the best Individual");
//		gf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
//		gf.setVisible(true);
//
//		Paint  p = new Paint();
//		//	p.weightMatrix = population.getIndividual(bestFinessIndex).getGene();
//			p.weightMatrix = this.readFilesTree();
////			p.weightMatrix = tree;
////			p.fitness = evaluation(weightMatrix, this.readFilesTree(), numberVertices, root);
//			p.fitness = cost;
//			//System.out.println("cost" + p.fitness);
//			//p.fitness1 = popFitness1[bestFinessIndex];
//			p.num_vertex = numberVertices;
//		    gf.add(p);
//		    gf.setVisible(true);
 }

 public double evaluation(double[][] weightMatrix, int[][] tree, int num_vertex, int startVertex){
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
 public int[][]  readFilesTree(){
	 int[][] tree = null;
	 int num_vertex ;
	 String fileName = "C:/Users/TrungTB/Desktop/tex.txt";
		Scanner sc;
		try {
			sc = new Scanner(new File(fileName));
			num_vertex = sc.nextInt();
			tree = new int[num_vertex][num_vertex];
			for (int i = 0; i < num_vertex; i++)
				for (int j = 0; j < num_vertex; j++)
					tree[i][j] = (int) sc.nextDouble();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} // string to store data from file
		return tree;
 }


 public static void main(String[] args){
 // read file

//	 ReadFiles.clusterReadFiles() ;
//	 System.out.println(args[0]);
//	 ReadFiles.clusterReadFiles("C:/Users/TrungTB/Desktop/test/" +args[0]+ ".clt");
	 Heuristic_CTMS heuristic_CTMS = new Heuristic_CTMS();
	// int[][] weight = heuristic_CTMS.readFilesTree();
	 double  cost = heuristic_CTMS.solve();

			    int [][] check = new int[ heuristic_CTMS.numberVertices][ heuristic_CTMS.numberVertices];
			    check = heuristic_CTMS.readFilesTree();
			    for(int i = 0; i <  heuristic_CTMS.numberVertices; i++){
			    	for( int j = 0 ; j <  heuristic_CTMS.numberVertices; j++){
			    		if( check[i][j] != check[j][i]){
			    			System.out.println(i+" "+j);
			    		}
			    	}
			    }

//	 PrintWriter pw = null;
//	try {
//		pw = new PrintWriter(new FileWriter(new File(args[1]), true));
//	} catch (IOException e) {
//		e.printStackTrace();
//	}
//	 pw.println(cost);
//	 pw.close();




 }


 }
