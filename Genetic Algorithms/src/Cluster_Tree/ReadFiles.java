package Cluster_Tree;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;

public class ReadFiles {
	public static ArrayList<Cluster> clusters = new ArrayList<Cluster>();
	private  Cluster cluster = new Cluster();
	public static int num_vertex = clusterReadFiles();
	public static double weightMatrix[][];
	public static int numberOfCluster; 
    public static int root;
    // files path: C:\Users\TrungTB\Desktop\16pr76.tsp
	
	public static int clusterReadFiles() {
		String fileName = null;
		// read filename from keyboard
		try {
			BufferedReader buf = new BufferedReader(new InputStreamReader(System.in));
			System.out.println("Input file name:");
			fileName = buf.readLine();
		} catch (IOException ex) {
		}
		BufferedReader br = null; // string to store data from file
		int num_vertex = 0;
		try {
			String sCurrentLine = null;
			br = new BufferedReader(new FileReader(fileName));
			// read lines 1..4
			for (int j = 0; j < 4; j++) {
				sCurrentLine = br.readLine();
			}
			String[] str = sCurrentLine.split(": ");
			num_vertex = Integer.parseInt(str[1]); 
			weightMatrix = new double[num_vertex][num_vertex];
			sCurrentLine = br.readLine();
			str = sCurrentLine.split(": ");
			numberOfCluster = Integer.parseInt(str[1]);
			sCurrentLine = br.readLine(); 
			sCurrentLine = br.readLine();
			Vertex[] vertices = new Vertex[num_vertex];
			
			// read the detail of the vertex
			for (int j = 0; j < num_vertex; j++) {
				sCurrentLine = br.readLine();
				str = sCurrentLine.split(" ");
				vertices[j] = new Vertex();
				// set coordinates to city
				vertices[j].setX(Double.parseDouble(str[1]));
				vertices[j].setY(Double.parseDouble(str[2]));
				// calculate distances
				for (int i = 0; i <= j; i++) {
					if (i == j) {
						weightMatrix[j][i] = 0;
					} else {
						weightMatrix[j][i] = weightMatrix[i][j] = Math.sqrt(Math.pow((vertices[j].getX() - vertices[i].getX()), 2)
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

}
