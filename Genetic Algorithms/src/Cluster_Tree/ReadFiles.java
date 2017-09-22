package Cluster_Tree;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;

public class ReadFiles {
	
	public static int num_vertex = clusterReadFiles();
	public static double distances[][];
	private static  Cluster cluster = new Cluster();
	public static ArrayList<Cluster> clusters = new ArrayList<Cluster>();
	public static int numberOfCluster = clusters.size(); 
	public static int root = 0;


	public static int clusterReadFiles() {
		String fileName = null;
		try {
			BufferedReader buf = new BufferedReader(new InputStreamReader(System.in));
			System.out.println("Input file: ");
			fileName = buf.readLine();
		} catch (IOException ex) {
		}
		BufferedReader br = null; 
		int num_vertex = 0;
		int numberOfCluster = 0;
		try {
			String sCurrentLine = null;
			br = new BufferedReader(new FileReader(fileName));
			for (int j = 0; j < 4; j++) {
				sCurrentLine = br.readLine();
			}
			String[] str = sCurrentLine.split(": ");
			num_vertex = Integer.parseInt(str[1]); 
			distances = new double[num_vertex][num_vertex];
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
						distances[j][i] = 0;
					} else {
						distances[j][i] = distances[i][j] = Math.sqrt(Math.pow((vertices[j].getX() - vertices[i].getX()), 2)
								+ Math.pow((vertices[j].getY() - vertices[i].getY()), 2));
					}
				}
			}
			sCurrentLine = br.readLine();
			str = sCurrentLine.split(": ");
			root = Integer.parseInt(str[1]);
			sCurrentLine = br.readLine();
			for(int i = 0; i < numberOfCluster; i++ ){
				int  numberClusterVertex;
				int arrayCluster;
				sCurrentLine = br.readLine();
				str = sCurrentLine.split(" ");
				numberClusterVertex = str.length;
				Cluster cluster = new Cluster();
				for(int j = 0; j < numberClusterVertex -1; j++ ){
					arrayCluster = Integer.parseInt(str[j]);
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
