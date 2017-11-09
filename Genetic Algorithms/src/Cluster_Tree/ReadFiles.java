package Cluster_Tree;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Scanner;

public class ReadFiles {
	public static ArrayList<Cluster> clusters = new ArrayList<Cluster>();
	public static Vertex[] vertices = new Vertex[5000];
	//public static int num_vertex = clusterReadFiles();
	public static int num_vertex = clusterReadFiles();
//	public static double weightMatrix[][] = readWeightMatrix();
	public static double weightMatrix[][];
	public static int numberOfCluster; 
    public static int root;
	
    
	
	public static int clusterReadFiles() {
//    public static int clusterReadFiles(String fileName) {
		String fileName = "C:/Users/TrungTB/Desktop/test/5eil51.clt";
		BufferedReader br = null; // string to store data from file
		try {
			String sCurrentLine = null;
			br = new BufferedReader(new FileReader(fileName));
			// read lines 1..4
			for (int j = 0; j < 3; j++) {
				sCurrentLine = br.readLine();
			}
//			System.out.println(" str =  "+str[1]);
			String[] str = sCurrentLine.split(": ");
			System.out.println(" str =  "+str[1]);
			num_vertex = Integer.parseInt(str[1]); 
			
//			weightMatrix = new double[num_vertex][num_vertex];
			weightMatrix = new double[num_vertex][num_vertex];
			sCurrentLine = br.readLine();
			str = sCurrentLine.split(": ");
			numberOfCluster = Integer.parseInt(str[1]);
			sCurrentLine = br.readLine(); 
			sCurrentLine = br.readLine();
			
			// read the detail of the vertex
			for (int j = 0; j < num_vertex; j++) {
				sCurrentLine = br.readLine();
				str = sCurrentLine.split("\\s+");
				vertices[j] = new Vertex();
				// set coordinates to city
				vertices[j].setX(Double.parseDouble(str[1]));
				vertices[j].setY(Double.parseDouble(str[2]));
				// calculate distances
				for (int i = 0; i <= j; i++) {
					if (i == j) {
						weightMatrix[j][i] = 0;
					} else {
						weightMatrix[j][i] = weightMatrix[i][j] = Math.sqrt(Math.pow((vertices[j].getX()
								- vertices[i].getX()), 2)
								+ Math.pow((vertices[j].getY() - vertices[i].getY()), 2));
					}
				}
			}
		
			sCurrentLine = br.readLine();
			sCurrentLine = br.readLine();
			str = sCurrentLine.split(": ");
			root = Integer.parseInt(str[1]);
			for(int i = 0; i < numberOfCluster; i++ ){
				int  numberClusterVertex;
				int arrayCluster;
				sCurrentLine = br.readLine();
				str = sCurrentLine.split(" ");
				numberClusterVertex = str.length;
				Cluster cluster = new Cluster();
				for(int j = 0; j < numberClusterVertex -2; ++j ){
					arrayCluster = Integer.parseInt(str[j+1]);
				     cluster.addElement(arrayCluster, j);
				}
//				System.out.println();
				clusters.add(cluster);
			}
			
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				if (br != null)
					br.close();
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		}
		return num_vertex; // return number of vertex
	}
	public static double[][] readWeightMatrix(){
		String filename= "C:/Users/TrungTB/git/ClusterTree/Genetic Algorithms/4eil76NON.txt";
		Scanner sc;
		double optimalCost;
		double[][] weghtMatrix2 = new double [num_vertex][num_vertex];
		try {
			sc = new Scanner(new File(filename));
			sc.nextLine();
			sc.nextLine();
			sc.nextLine();
//	     	optimalCost = sc.nextDouble();
			for( int i = 0; i < num_vertex; i++ ){
				for( int j = 0; j < num_vertex; j++){
					weghtMatrix2[i][j] = sc.nextDouble();
					System.out.print(" " + weghtMatrix2[i][j]);
				}
				System.out.println();
			}
	
		} catch (FileNotFoundException e) {
		
			e.printStackTrace();
		}
		
		return weghtMatrix2;
		
	}
}
